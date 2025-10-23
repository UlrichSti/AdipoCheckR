# AdipoCheckR
**Version 0.2**  
**Author:** Ulrich Stifel  
**License:** MIT

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://adipocheckr.streamlit.app/)


AdipoCheckR estimates adipogenic composition in **in vitro bulk RNA-seq** using a **hybrid signature-based NNLS deconvolution approach**. It quantifies the proportions of  
**ASPC (adipose stromal progenitor cells), Preadipocytes, and mature Adipocytes** in bulk samples derived from adipogenic differentiation experiments.

---

## Features
- Reference-based deconvolution using a fixed adipogenic signature matrix
- Hybrid marker strategy:
  - Curated adipogenesis marker genes from single-cell data
  - Supplementary automatically selected markers per cell type
- Adaptive boosting of adipocyte signal to improve mid-stage resolution (e.g., day 8 differentiation)
- Clean web interface (`Streamlit`) and CLI backend (`Python`)
- Designed for reproducible differentiation quantification in adipogenesis studies

---

## Installation

### Requirements
- Python 3.10+
- Git

### Setup
```bash
git clone https://github.com/UlrichSti/AdipoCheckR.git
cd AdipoCheckR
python -m venv venv
.\venv\Scripts\activate    # Windows
pip install -r requirements.txt

Usage
Option A: Web Interface

cd webapp
..\venv\Scripts\activate
streamlit run app.py

Upload your bulk RNA-seq matrix and view deconvolution results interactively.
Option B: Command Line (CLI)

python backend/deconvolution.py \
  --signature reference/signature_matrix.csv \
  --bulk path/to/bulk_matrix.csv \
  --out-prefix results/deconv_run \
  --markers-path results/tables/top_markers.csv

Input Format
Bulk RNA-seq matrix

    File format: .csv or .tsv

    Rows = genes, Columns = samples

    First column = gene names

    Example:

gene,Sample1,Sample2,Sample3
PDGFRA,120,98,102
ADIPOQ,0,15,2300
FABP4,5,90,1800

Signature matrix

Provided in repository:

reference/signature_matrix.csv

Curated markers

Provided in repository:

results/tables/top_markers.csv

Output
File	Description
fractions.csv	Estimated ASPC, Preadipocyte, Adipocyte fractions per sample

Example output:

Sample,ASPC,Preadipocyte,Adipocyte
d0_rep1,0.92,0.08,0.00
d8_rep2,0.43,0.37,0.20
d14_rep1,0.05,0.10,0.85

Method

AdipoCheckR uses non-negative least squares (NNLS) for reference-based cell proportion estimation.
The method combines curated biology and data-driven refinement:

    Hybrid marker strategy

        Curated markers from single-cell adipogenesis data (ASPC, Preadipocyte, Adipocyte)

        Supplementary auto-detected markers selected per cell type

    Gene expression harmonization

        Expression standardization via signature-based z-scoring

    Gene weighting

        Curated markers are weighted higher than auto markers

    Adaptive adipocyte boosting

        When adipocyte marker signal is detected in a sample, adipocyte genes are boosted

        Prevents underestimation of adipocyte fractions at intermediate timepoints (e.g. day 8)

    Final estimation

        NNLS solves:
        bulk ≈ signature × proportions,
        constrained to non-negative and normalized to sum to 1

Citation

If you use AdipoCheckR in your work, please cite this repository:

    Stifel U. AdipoCheckR: hybrid NNLS-based deconvolution for in vitro adipogenesis. 2025.
    GitHub: https://github.com/UlrichSti/AdipoCheckR

License

This project is released under the MIT License (see LICENSE file).
