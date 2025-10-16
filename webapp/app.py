import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sys
import os
import re
import csv

# ---- import backend ----
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "backend"))
from deconvolution import nnls_deconvolution

st.set_page_config(page_title="Adipocyte Deconvolution Tool", layout="wide")

# ---------- Helpers ----------
def read_any_table(uploaded_file):
    try:
        header = uploaded_file.readline().decode("utf-8")
    except:
        header = uploaded_file.readline()
    uploaded_file.seek(0)
    dialect = csv.Sniffer().sniff(header, delimiters=[",", ";", "\t"])
    df = pd.read_csv(uploaded_file, delimiter=dialect.delimiter)
    df = df.rename(columns={df.columns[0]: "gene"})
    return df.set_index("gene")

def clean_gene_index(idx):
    cleaned = []
    for g in idx.astype(str):
        g = re.sub(r"^hg19-", "", g, flags=re.IGNORECASE)
        g = re.sub(r"\.\d+$", "", g)  # remove version suffixes
        cleaned.append(g.strip().upper())
    return pd.Index(cleaned, name="gene")

def cpm(df):
    lib = df.sum(axis=0).replace(0, np.nan)
    res = df.divide(lib, axis=1) * 1e6
    return res.fillna(0.0)

# ---------- Load fixed signature ----------
SIGN_PATH = os.path.join(os.path.dirname(__file__), "..", "reference", "signature_matrix.csv")
signature = pd.read_csv(SIGN_PATH, index_col=0)
signature.index = clean_gene_index(signature.index)
signature = signature.groupby(signature.index).mean()

# ---------- UI ----------
st.title("AdipoCheckR")

bulk_file = st.file_uploader("Upload bulk RNA-seq file (CSV/TSV):", type=["csv", "tsv"])

if st.button("Run Deconvolution"):
    if bulk_file is None:
        st.error("Please upload a bulk RNA-seq file.")
    else:
        bulk = read_any_table(bulk_file)
        bulk.index = clean_gene_index(bulk.index)
        bulk = bulk.groupby(bulk.index).sum()

        # auto-scale raw counts -> CPM
        bulk = cpm(bulk)

        # intersect
        genes = signature.index.intersection(bulk.index)
        if len(genes) == 0:
            st.error("No overlapping genes found. Check gene symbols.")
        else:
            sig = signature.loc[genes]
            bul = bulk.loc[genes]

            # deconvolution
            results = []
            for sample in bul.columns:
                coeffs, _ = nnls_deconvolution(sig, bul[sample])
                coeffs = coeffs / coeffs.sum()
                results.append(pd.Series(coeffs, index=sig.columns, name=sample))

            frac_df = pd.DataFrame(results)

            # ---- OUTPUT TABLE ----
            st.subheader("Cell Type Fractions")
            st.dataframe(frac_df.style.format("{:.3f}"))

            # ---- STACKED BAR ----
            st.subheader("Composition Plot")
            color_map = {
                "Adipocyte": "#FFA500",
                "Preadipocyte": "#FFD700",
                "ASPC": "#808080"
            }
            fig = px.bar(
                frac_df,
                x=frac_df.index,
                y=frac_df.columns,
                barmode="stack",
                color_discrete_map=color_map
            )
            st.plotly_chart(fig, use_container_width=True)

            # ---- DOWNLOAD ----
            st.download_button("Download Results (CSV)", frac_df.to_csv().encode(), "fractions.csv")
