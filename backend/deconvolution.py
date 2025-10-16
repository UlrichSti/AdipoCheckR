#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
from scipy.optimize import nnls

"""
AdipoCheckR backend: NNLS deconvolution with hybrid markers + adaptive adipocyte boosting.

Features
- ALWAYS include curated markers from results/tables/top_markers.csv (columns: cluster, gene)
- Add auto-selected markers per cell type to reach a target size
- Curated markers get higher base weights than auto markers
- Z-score alignment so bulk & signature are comparable
- ADAPTIVE: if a sample shows adipocyte marker signal (in z-score space), boost adipocyte marker weights
  -> makes d8 more likely to show adipocytes while keeping d0 clean (no boost if no signal)
"""

# =========================
# IO helpers
# =========================

def read_table_any(path: str) -> pd.DataFrame:
    """Read CSV or TSV; do not assume index column in marker files."""
    if path.lower().endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")

def read_matrix(path: str) -> pd.DataFrame:
    """Read signature/bulk matrix: rows=genes, cols=samples/ct; first column is index."""
    sep = "," if path.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, index_col=0)

# =========================
# Gene processing & scaling
# =========================

def clean_gene_index(idx: pd.Index) -> pd.Index:
    out = []
    for g in idx.astype(str):
        g = g.strip().strip('"').strip("'")
        if g.lower().startswith("hg19-"):
            g = g[5:]
        # strip Ensembl version suffix like ".12"
        if "." in g:
            parts = g.split(".")
            if len(parts) == 2 and parts[1].isdigit():
                g = parts[0]
        out.append(g.upper())
    return pd.Index(out, name="gene")

def zscore_to_signature(sig: pd.DataFrame, bul: pd.DataFrame):
    """
    Z-score each gene using signature stats so bulk and signature are comparable.
    """
    mu = sig.mean(axis=1)
    sd = sig.std(axis=1).replace(0, np.nan)
    sig_z = (sig.sub(mu, axis=0)).div(sd, axis=0).fillna(0.0)
    bul_z = (bul.sub(mu, axis=0)).div(sd, axis=0).fillna(0.0)
    return sig_z, bul_z

# =========================
# Marker selection (hybrid)
# =========================

def load_curated_markers(markers_path: str, signature: pd.DataFrame) -> dict:
    """
    Load curated markers. Expected columns (case-insensitive):
      - 'cluster' (cell type; ASPC/Preadipocyte/Adipocyte)
      - 'gene'
    Extra columns (p_val, avg_log2FC, etc.) are ignored.
    Returns dict: { cell_type: set([genes...]) } filtered to signature genes.
    """
    df = read_table_any(markers_path)
    cols = {c.lower(): c for c in df.columns}
    ct_col = cols.get("cluster") or cols.get("cell_type") or cols.get("label")
    gene_col = cols.get("gene") or cols.get("symbol") or cols.get("gene_symbol")
    if ct_col is None or gene_col is None:
        raise ValueError("Curated markers must contain columns 'cluster' and 'gene' (case-insensitive).")

    genes_clean = clean_gene_index(df[gene_col].astype(str))
    df = df.copy()
    df[gene_col] = genes_clean

    sig_genes = set(signature.index)
    df = df[df[gene_col].isin(sig_genes)]

    curated = {}
    for ct, sub in df.groupby(ct_col):
        curated[str(ct)] = set(sub[gene_col].tolist())
    return curated

def select_markers_auto(signature: pd.DataFrame, per_class: int) -> dict:
    """
    Auto markers: for each cell type, pick top 'per_class' genes by
    absolute deviation from the across-CT mean.
    Returns dict {ct: [genes...]} (not union).
    """
    auto = {}
    mean_all = signature.mean(axis=1)
    for ct in signature.columns:
        delta = (signature[ct] - mean_all).abs()
        auto[ct] = list(delta.sort_values(ascending=False).head(per_class).index)
    return auto

def build_hybrid_marker_set(signature: pd.DataFrame,
                            markers_path: str | None,
                            per_class_total: int = 150,
                            curated_cap_per_ct: int = 120,
                            curated_weight: float = 3.0,
                            auto_weight: float = 1.0) -> tuple[pd.Index, dict, dict]:
    """
    Build hybrid markers:
      - ALWAYS include curated markers up to curated_cap_per_ct per CT (if provided)
      - THEN add auto-selected genes until per_class_total per CT is reached
      - Assign weights: curated_weight > auto_weight
    Returns:
      - selected gene Index (union across CTs)
      - weights dict {gene: weight}
      - curated_map dict {ct: set(genes)} for later adaptive boosting
    """
    weights = {}
    selected_union = set()

    curated_map = {}
    if markers_path is not None and os.path.isfile(markers_path):
        curated_map = load_curated_markers(markers_path, signature)

    auto = select_markers_auto(signature, per_class=per_class_total)

    for ct in signature.columns:
        curated_ct = list(curated_map.get(ct, []))
        curated_ct = curated_ct[:curated_cap_per_ct] if curated_ct else []

        need_auto = max(0, per_class_total - len(curated_ct))
        auto_ct = [g for g in auto[ct] if g not in curated_ct][:need_auto]

        selected_union.update(curated_ct)
        selected_union.update(auto_ct)

        for g in curated_ct:
            weights[g] = max(weights.get(g, 0.0), curated_weight)
        for g in auto_ct:
            weights.setdefault(g, auto_weight)

        # keep curated set per ct for adaptive boosting use
        curated_map[ct] = set(curated_ct)

    return pd.Index(sorted(selected_union)), weights, curated_map

# =========================
# Base NNLS
# =========================

def nnls_deconvolution(signature_df: pd.DataFrame, bulk_vector: pd.Series):
    S = signature_df.values
    y = bulk_vector.values
    coeffs, residual = nnls(S, y)
    return coeffs, residual

# =========================
# Main pipeline with adaptive adipocyte boost
# =========================

def prepare_and_deconvolve(signature: pd.DataFrame,
                           bulk: pd.DataFrame,
                           markers_path: str | None = None,
                           per_class_total: int = 150,
                           curated_cap_per_ct: int = 120,
                           curated_weight: float = 3.0,
                           auto_weight: float = 1.0,
                           adipocyte_adaptive_boost: bool = True,
                           boost_low_threshold: float = 0.20,
                           boost_high_threshold: float = 0.50,
                           boost_low_factor: float = 1.5,
                           boost_high_factor: float = 2.0) -> pd.DataFrame:
    """
    1) Hybrid markers (curated ALWAYS included; auto fills to target)
    2) Intersect with bulk
    3) Z-score to signature stats
    4) Apply per-gene weights (curated > auto)
    5) OPTIONAL: per-sample adaptive boost of ADIPOCYTE curated markers if marker signal present
    6) NNLS per sample + renormalize
    """
    marker_genes, weights, curated_map = build_hybrid_marker_set(
        signature=signature,
        markers_path=markers_path,
        per_class_total=per_class_total,
        curated_cap_per_ct=curated_cap_per_ct,
        curated_weight=curated_weight,
        auto_weight=auto_weight,
    )

    genes = signature.index.intersection(bulk.index).intersection(marker_genes)
    if len(genes) == 0:
        raise ValueError("No overlap after hybrid marker selection. Check markers and gene IDs.")

    S = signature.loc[genes]
    Y = bulk.loc[genes]

    S_z, Y_z = zscore_to_signature(S, Y)

    # Base per-gene weights (row-wise)
    if weights:
        w_vec = np.array([weights.get(g, 1.0) for g in S_z.index])
        S_z = (S_z.T * w_vec).T
        Y_z = (Y_z.T * w_vec).T

    # Adaptive adipocyte boosting (per sample, based on curated adipocyte markers signal)
    adipocyte_genes_curated = curated_map.get("Adipocyte", set())
    adipocyte_genes_curated = [g for g in adipocyte_genes_curated if g in S_z.index]

    if adipocyte_adaptive_boost and len(adipocyte_genes_curated) > 0:
        adipo_idx = S_z.index.get_indexer(adipocyte_genes_curated)
        adipo_idx = adipo_idx[adipo_idx >= 0]
        if len(adipo_idx) > 0:
            # For each sample, compute mean z-score of adipocyte curated markers
            adipo_scores = Y_z.iloc[adipo_idx, :].mean(axis=0)  # one score per sample
            # Build a per-sample multiplier for those genes
            multipliers = np.ones_like(Y_z.values)
            for j, sample in enumerate(Y_z.columns):
                score = float(adipo_scores.iloc[j])
                if score >= boost_high_threshold:
                    factor = boost_high_factor
                elif score >= boost_low_threshold:
                    factor = boost_low_factor
                else:
                    factor = 1.0
                if factor != 1.0:
                    # multiply ONLY adipocyte curated marker rows for this sample
                    multipliers[adipo_idx, j] *= factor
            # Apply sample-specific boosts to Y_z only (keeps signature fixed; we weight the evidence)
            Y_z = Y_z * multipliers

    rows = []
    for sample in Y_z.columns:
        coeffs, _ = nnls_deconvolution(S_z, Y_z[sample])
        s = coeffs.sum()
        if s > 0:
            coeffs = coeffs / s
        rows.append(pd.Series(coeffs, index=S_z.columns, name=sample))
    return pd.DataFrame(rows)

# =========================
# CLI (optional)
# =========================

def cpm(df: pd.DataFrame) -> pd.DataFrame:
    lib = df.sum(axis=0).replace(0, np.nan)
    return (df.div(lib, axis=1) * 1e6).fillna(0.0)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", required=True)
    parser.add_argument("--bulk", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--markers-path", default=None,
                        help="Optional curated markers CSV/TSV with columns: cluster,gene")
    parser.add_argument("--per-class-total", type=int, default=150)
    parser.add_argument("--curated-cap-per-ct", type=int, default=120)
    parser.add_argument("--curated-weight", type=float, default=3.0)
    parser.add_argument("--auto-weight", type=float, default=1.0)
    parser.add_argument("--no-adaptive", action="store_true",
                        help="Disable adaptive adipocyte boosting")
    args = parser.parse_args()

    sig = read_matrix(args.signature)
    bul = read_matrix(args.bulk)
    sig.index = clean_gene_index(sig.index)
    bul.index = clean_gene_index(bul.index)
    sig = sig.groupby(sig.index).mean()
    bul = bul.groupby(bul.index).sum()
    bul = cpm(bul)  # safe default

    frac = prepare_and_deconvolve(
        signature=sig,
        bulk=bul,
        markers_path=args.markers_path,
        per_class_total=args.per_class_total,
        curated_cap_per_ct=args.curated_cap_per_ct,
        curated_weight=args.curated_weight,
        auto_weight=args.auto_weight,
        adipocyte_adaptive_boost=not args.no_adaptive
    )
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    frac.to_csv(f"{args.out_prefix}_fractions.csv")

if __name__ == "__main__":
    main()
