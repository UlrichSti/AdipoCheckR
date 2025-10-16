\# Adipocyte Deconvolution Reference



This repository generates a cell-type signature matrix from single-cell RNA-seq data of human adipocyte differentiation (SGBS cell system). The signature matrix can be used for bulk RNA-seq deconvolution to estimate the proportion of:



\- Adipocyte Stem/Progenitor Cells (ASPC)

\- Preadipocytes

\- Mature Adipocytes



---



\## 🔬 Methods



The pipeline uses:



| Step | Tool |

|------|------|

| QC + Clustering | Seurat (SCTransform) |

| Trajectory inference | Monocle3 |

| Marker detection | Seurat FindAllMarkers |

| Signature Matrix | Seurat AverageExpression |



---



\## Repository Structure

adipocyte-deconvolution/

├── data/ # raw scRNA-seq input (not included)

├── scripts/ # R analysis scripts

├── results/ # output UMAPs, signature matrix

├── README.md # project documentation

└── LICENSE # MIT license for open use





---



\## Usage



To reproduce the signature matrix:



```bash

Rscript scripts/01\_SGBS\_processing.R



Output file:

results/tables/signature\_matrix\_clean.csv



Coming Next



Build deconvolution backend in Python

Deploy web app using Streamlit

Upload example bulk RNA-seq test files



Author



Created by Ulrich Stifel



License

This project is licensed under the MIT License – see LICENSE file for details.







