# Set your desired base path (adjust for your system)
setwd("D:/2024/3. Projekte/SGBS single cell")

# Create folders
dir.create("adipocyte-deconvolution")
setwd("adipocyte-deconvolution")

dir.create("data")
dir.create("scripts")
dir.create("results")
dir.create("results/figures")
dir.create("results/tables")

# Create empty placeholder files
file.create("README.md")
file.create(".gitignore")

# Copy your R script into scripts/
file.copy("path/to/your/current_script.R", "scripts/01_SGBS_processing.R")
