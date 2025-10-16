import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import sys, os, re, csv

# ---- import improved backend ----
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "backend"))
from deconvolution import prepare_and_deconvolve

st.set_page_config(page_title="AdipoCheckR", layout="wide")

# ---- light layout CSS ----
st.markdown("""
<style>
.block-container {padding-top: 2rem; padding-bottom: 3rem; max-width: 1200px;}
.card {border:1px solid #e6e6e6; border-radius:12px; padding:18px 20px; margin:16px 0; background:white;}
h1, h2, h3, p {text-align:left;}
footer {visibility:hidden;}
.small-dim {color:#777; font-size:0.85rem; margin-top:20px;}
</style>
""", unsafe_allow_html=True)

# ---- helpers ----
def read_any_table(uploaded_file):
    try:
        header = uploaded_file.readline().decode("utf-8")
    except Exception:
        header = uploaded_file.readline()
    uploaded_file.seek(0)
    dialect = csv.Sniffer().sniff(header, delimiters=[",",";","\\t"])
    df = pd.read_csv(uploaded_file, delimiter=dialect.delimiter)
    df = df.rename(columns={df.columns[0]: "gene"})
    return df.set_index("gene")

def clean_gene_index(idx):
    cleaned = []
    for g in idx.astype(str):
        g = re.sub(r"^hg19-", "", g, flags=re.IGNORECASE)
        g = re.sub(r"\\.[0-9]+$", "", g)
        cleaned.append(g.strip().upper())
    return pd.Index(cleaned, name="gene")

def cpm(df):
    lib = df.sum(axis=0).replace(0, np.nan)
    return (df.div(lib, axis=1) * 1e6).fillna(0.0)

# ---- load fixed signature ----
SIGN_PATH = os.path.join(os.path.dirname(__file__), "..", "reference", "signature_matrix.csv")
signature = pd.read_csv(SIGN_PATH, index_col=0)
signature.index = clean_gene_index(signature.index)
signature = signature.groupby(signature.index).mean()

# optional curated markers path
MARKERS_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "tables", "top_markers.csv")
if not os.path.isfile(MARKERS_PATH):
    MARKERS_PATH = None  # fallback to auto-only if file missing

# ---- header ----
st.title("AdipoCheckR")
st.write("Bulk RNA-seq adipocyte deconvolution into ASPC, Preadipocyte and Adipocyte fractions.")

# ---- upload panel ----
with st.container():
    st.markdown('<div class="card">', unsafe_allow_html=True)
    st.subheader("Upload bulk RNA-seq matrix")
    bulk_file = st.file_uploader("Upload CSV/TSV (rows = genes, columns = samples):", type=["csv","tsv"])
    run = st.button("Run Deconvolution")
    st.markdown("</div>", unsafe_allow_html=True)

# ---- main logic ----
if run:
    if bulk_file is None:
        st.error("Please upload a bulk RNA-seq file.")
    else:
        bulk = read_any_table(bulk_file)
        bulk.index = clean_gene_index(bulk.index)
        bulk = bulk.groupby(bulk.index).sum()
        bulk = cpm(bulk)

        genes = signature.index.intersection(bulk.index)
        if len(genes) == 0:
            st.error("No overlapping genes found between bulk data and signature.")
        else:
            frac_df = prepare_and_deconvolve(
                signature=signature.loc[genes],
                bulk=bulk.loc[genes],
                markers_path=MARKERS_PATH,
                per_class_total=150,
                curated_cap_per_ct=120,
                curated_weight=3.0,
                auto_weight=1.0,
            )

            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.subheader("Estimated Cell Fractions")
            st.dataframe(frac_df.style.format("{:.3f}"))

            color_map = {
                "Adipocyte": "#FFA500",
                "Preadipocyte": "#FFD700",
                "ASPC": "#808080"
            }
            st.subheader("Stacked Composition Plot")
            fig = px.bar(
                frac_df,
                x=frac_df.index,
                y=frac_df.columns,
                barmode="stack",
                color_discrete_map=color_map
            )
            st.plotly_chart(fig, use_container_width=True)

            st.download_button("Download Fractions (CSV)",
                               frac_df.to_csv().encode(), "fractions.csv")
            st.markdown("</div>", unsafe_allow_html=True)

# ---- footer ----
st.markdown('<p class="small-dim">Developed by Ulrich Sti · AdipoCheckR v0.1</p>', unsafe_allow_html=True)
