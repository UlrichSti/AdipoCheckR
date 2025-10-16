#!/usr/bin/env python3
"""
Deconvolution backend (CIBERSORT-like) using NNLS.

Inputs:
  - Signature matrix (CSV/TSV): rows = genes, columns = cell types
  - Bulk RNA-seq expression (CSV/TSV): rows = genes, columns = samples

Outputs:
  - <out_prefix>_fractions.csv         : cell-type fractions per sample
  - <out_prefix>_diagnostics.csv       : per-sample fit metrics (RMSE, residual norm, genes used)
  - <out_prefix>_gene_match_report.txt : summary of gene matching

Usage example (from project root):
  venv\\Scripts\\activate
  python backend\\deconvolution.py ^
      --signature results\\tables\\signature_matrix_clean.csv ^
      --bulk path\\to\\bulk_matrix.csv ^
      --out-prefix results\\tables\\deconvolution

Notes:
  - Gene matching is case-insensitive and strips 'hg19-' prefixes automatically.
  - Bulk file can be CSV or TSV; gene symbols must be in the first column.
  - Fractions are non-negative and normalized to sum to 1 per sample.
"""

import argparse
import os
import sys
import re
from typing import Tuple, List

import numpy as np
import pandas as pd
from scipy.optimize import nnls


# --------------------------- Helpers ---------------------------

def _read_table_any(path: str) -> pd.DataFrame:
    """
    Read CSV or TSV into a DataFrame.
    Assumes first column contains gene symbols. Supports comma or tab separators.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")

    # Try comma, then tab
    try:
        df = pd.read_csv(path, sep=",", header=0)
        if df.shape[1] == 1:
            # probably a TSV misread as CSV; retry with tab
            df = pd.read_csv(path, sep="\t", header=0)
    except Exception:
        df = pd.read_csv(path, sep="\t", header=0)

    if df.shape[1] < 2:
        raise ValueError(f"{path}: expected at least 2 columns (gene + one sample/celltype).")

    # First column = gene symbols
    df = df.rename(columns={df.columns[0]: "gene"})
    df["gene"] = df["gene"].astype(str)
    df = df.dropna(subset=["gene"])
    df = df.set_index("gene")

    # Coerce numeric for the rest
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.fillna(0.0)

    return df


def _clean_gene_symbols(idx: pd.Index) -> pd.Index:
    """
    Clean gene symbols:
      - strip 'hg19-' prefix if present
      - uppercase for case-insensitive matching
      - strip whitespace
    """
    cleaned = []
    for g in idx.astype(str):
        g = re.sub(r"^hg19-", "", g, flags=re.IGNORECASE)
        g = g.strip().upper()
        cleaned.append(g)
    return pd.Index(cleaned, name="gene")


def load_signature(path: str) -> pd.DataFrame:
    """
    Load signature matrix: genes x cell_types.
    """
    sig = _read_table_any(path)
    sig.index = _clean_gene_symbols(sig.index)
    # Drop duplicate genes by keeping the mean (rare, but safer)
    sig = sig.groupby(sig.index).mean()
    # Remove all-zero genes
    sig = sig.loc[(sig.sum(axis=1) != 0)]
    return sig


def load_bulk(path: str) -> pd.DataFrame:
    """
    Load bulk matrix: genes x samples.
    """
    bulk = _read_table_any(path)
    bulk.index = _clean_gene_symbols(bulk.index)
    bulk = bulk.groupby(bulk.index).sum()
    # Remove all-zero genes
    bulk = bulk.loc[(bulk.sum(axis=1) != 0)]
    return bulk


def intersect_signature_bulk(sig: pd.DataFrame, bulk: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    """
    Intersect genes between signature and bulk.
    Returns: (sig_int, bulk_int, matched_genes)
    """
    genes = sig.index.intersection(bulk.index)
    sig_int = sig.loc[genes].copy()
    bulk_int = bulk.loc[genes].copy()
    return sig_int, bulk_int, list(genes)


def nnls_deconvolution(signature: pd.DataFrame, bulk: pd.Series) -> Tuple[np.ndarray, float]:
    """
    Solve non-negative least squares for one bulk sample.
    signature: genes x cell_types
    bulk:      genes (vector)
    Returns: (proportions, residual_norm)
    """
    S = signature.values  # shape G x C
    y = bulk.values       # shape G
    # NNLS on each sample
    coeffs, residual_norm = nnls(S, y)
    # Normalize to sum to 1 (avoid divide-by-zero)
    total = coeffs.sum()
    if total > 0:
        coeffs = coeffs / total
    return coeffs, residual_norm


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    err = y_true - y_pred
    return float(np.sqrt(np.mean(err ** 2)))


# --------------------------- Main ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Bulk RNA-seq deconvolution using NNLS.")
    ap.add_argument("--signature", required=True, help="Path to signature matrix CSV/TSV (genes x cell types).")
    ap.add_argument("--bulk", required=True, help="Path to bulk expression CSV/TSV (genes x samples).")
    ap.add_argument("--out-prefix", required=True, help="Output prefix, e.g., results/tables/deconvolution")
    ap.add_argument("--min-genes", type=int, default=50, help="Minimum overlapping genes required (default: 50).")
    ap.add_argument("--libsize-normalize", action="store_true",
                    help="Library-size normalize bulk to counts-per-million (CPM) before deconvolution.")
    args = ap.parse_args()

    # Load inputs
    signature = load_signature(args.signature)
    bulk = load_bulk(args.bulk)

    # Optional library-size normalization (CPM) on bulk
    if args.libsize_normalize:
        libsize = bulk.sum(axis=0).replace(0, np.nan)
        bulk = bulk.divide(libsize, axis=1) * 1e6
        bulk = bulk.fillna(0.0)

    # Intersect genes
    sig_int, bulk_int, matched_genes = intersect_signature_bulk(signature, bulk)

    if len(matched_genes) < args.min_genes:
        print(f"[ERROR] Only {len(matched_genes)} overlapping genes < --min-genes ({args.min_genes}).", file=sys.stderr)
        print("        Check that your bulk gene symbols are standard HGNC (e.g., ADIPOQ) and first column is 'Gene'.",
              file=sys.stderr)
        sys.exit(1)

    # Deconvolve each sample
    cell_types = list(sig_int.columns)
    fractions = []
    diag_rows = []

    S = sig_int.values  # reuse for RMSE later

    for sample in bulk_int.columns:
        y = bulk_int[sample]
        coeffs, residual_norm = nnls_deconvolution(sig_int, y)
        # Predicted expression for RMSE diagnostic
        y_pred = S @ coeffs
        sample_rmse = rmse(y.values, y_pred)

        # Save
        fractions.append(pd.Series(coeffs, index=cell_types, name=sample))
        diag_rows.append({
            "sample": sample,
            "overlap_genes": len(matched_genes),
            "residual_norm": residual_norm,
            "rmse": sample_rmse
        })

    frac_df = pd.DataFrame(fractions)  # samples x cell_types
    diag_df = pd.DataFrame(diag_rows).set_index("sample")

    # Ensure rows are samples (common convention)
    frac_df.index.name = "sample"

    # Save outputs
    out_frac = f"{args.out_prefix}_fractions.csv"
    out_diag = f"{args.out_prefix}_diagnostics.csv"
    out_match = f"{args.out_prefix}_gene_match_report.txt"

    os.makedirs(os.path.dirname(out_frac), exist_ok=True)
    frac_df.to_csv(out_frac, index=True)
    diag_df.to_csv(out_diag, index=True)

    with open(out_match, "w", encoding="utf-8") as fh:
        fh.write("Deconvolution gene matching report\n")
        fh.write("----------------------------------\n")
        fh.write(f"Signature file : {args.signature}\n")
        fh.write(f"Bulk file      : {args.bulk}\n")
        fh.write(f"Matched genes  : {len(matched_genes)}\n\n")
        fh.write("Examples of matched genes (first 50):\n")
        for g in matched_genes[:50]:
            fh.write(f"  {g}\n")

    print(f"[OK] Fractions saved to   : {out_frac}")
    print(f"[OK] Diagnostics saved to : {out_diag}")
    print(f"[OK] Match report saved to: {out_match}")


if __name__ == "__main__":
    main()
