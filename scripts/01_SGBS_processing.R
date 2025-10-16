###############################################################################
# Title: SGBS Adipocyte Differentiation Analysis from GSE226365
# Author: Ulrich Stifel
# Description: Create Seurat object from D0/D8, perform QC, clustering,
#              pseudotime ordering, rename clusters, and export signature matrix
###############################################################################

# --- Load libraries ----------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(monocle3)
library(SeuratWrappers)

# --- Load 10X data -----------------------------------------------------------
SGBS_D0 <- Read10X("data/SGBS_D0")   # processed 10X folder
SGBS_D8 <- Read10X("data/SGBS_D8")

D0 <- CreateSeuratObject(counts = SGBS_D0, project = "D0")
D8 <- CreateSeuratObject(counts = SGBS_D8, project = "D8")

# Merge two Seurat objects (keep sample ID)
merged_seurat <- merge(D0, y = D8, add.cell.ids = c("D0", "D8"), project = "SGBS")
merged_seurat$sample <- sapply(strsplit(colnames(merged_seurat), "_"), `[`, 1)

# --- QC filtering ------------------------------------------------------------
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4500)

# --- Normalization, Dimensionality reduction, Clustering ---------------------
merged_seurat <- SCTransform(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat)
ElbowPlot(merged_seurat)

merged_seurat <- RunUMAP(merged_seurat, dims = 1:15)
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:15)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.6)

# --- Quick biological sanity checks -----------------------------------------
FeaturePlot(merged_seurat, features = c("hg19-PDGFRA", "hg19-ADIPOQ"), cols = c("lightgrey", "red"), ncol = 2)
VlnPlot(merged_seurat, features = c("hg19-PDGFRA", "hg19-ADIPOQ"), pt.size = 0)

# --- Pseudotime analysis using Monocle3 -------------------------------------
cds <- as.cell_data_set(merged_seurat)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")

# sync UMAP coordinates with Seurat
cds@int_colData@listData$reducedDims$UMAP <- merged_seurat@reductions$umap@cell.embeddings

cds <- learn_graph(cds, use_partition = FALSE)

# choose root cluster manually (adjust number if needed)

DefaultAssay(merged_seurat) <- "SCT"

# Plot UMAP with cluster numbers
DimPlot(
  merged_seurat,
  reduction = "umap",
  label = TRUE,        # adds big labels for each cluster
  label.size = 6,
  repel = TRUE         # avoids overlapping labels
) + ggtitle("Cluster UMAP - Choose root cluster for pseudotime")
root_cluster <- 2
cds <- order_cells(cds, root_cells = colnames(cds[, clusters(cds) == root_cluster]))

FeaturePlot(
  merged_seurat,
  features = c("hg19-PDGFRA", "hg19-PPARG", "hg19-ADIPOQ"),
  reduction = "umap",
  cols = c("lightgrey", "red"),
  ncol = 3
)

merged_seurat$pseudotime <- pseudotime(cds)

# visualize pseudotime
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE)

# --- Rename clusters into 3 biological states --------------------------------
# Rename individual clusters
merged_seurat <- RenameIdents(merged_seurat,
                              `2` = "ASPC",
                              `0` = "ASPC",
                              `5` = "ASPC",
                              `3` = "ASPC",
                              `7` = "Preadipocyte",        
                              `6` = "Preadipocyte",
                              `4` = "Preadipocyte",
                              `1` = "Adipocyte"
)

merged_seurat$CellType <- Idents(merged_seurat)
DefaultAssay(merged_seurat) <- "SCT"

DimPlot(merged_seurat, reduction = "umap", group.by = "CellType", label = TRUE) +
  ggtitle("Final Cell Type Annotation")

# --- Generate signature matrix -----------------------------------------------
# average expression per cell type
# Generate signature matrix
sig_matrix <- AverageExpression(merged_seurat, group.by = "CellType", slot = "data")$SCT
rownames(sig_matrix) <- gsub("^hg19-", "", rownames(sig_matrix))
# clean gene names (remove "hg19-" prefix)
rownames(sig_matrix) <- gsub("^hg19-", "", rownames(sig_matrix))
write.csv(sig_matrix, "results/tables/signature_matrix_clean.csv")

# marker genes for each cluster
merged_seurat <- PrepSCTFindMarkers(merged_seurat)
markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, logfc.threshold = 0.25)
top_markers <- markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

library(patchwork)

# Make sure assay is correct
DefaultAssay(merged_seurat) <- "SCT"

# 1. UMAP by cell type
p1 <- DimPlot(
  merged_seurat,
  reduction = "umap",
  group.by = "CellType",
  label = TRUE,
  repel = TRUE
) + ggtitle("Cell Type Annotation")

# 2. Pseudotime WITH trajectory
p2 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_branch_points = TRUE,
  label_roots = TRUE,
  label_leaves = TRUE,
  show_trajectory_graph = TRUE
) + ggtitle("Pseudotime + Trajectory")

# 3. Early marker (PDGFRA)
p3 <- FeaturePlot(
  merged_seurat,
  features = "hg19-PDGFRA",
  reduction = "umap",
  cols = c("lightgrey", "darkgreen")
) + ggtitle("Early Marker: PDGFRA")

# 4. Mid marker (PPARG)
p4 <- FeaturePlot(
  merged_seurat,
  features = "hg19-PPARG",
  reduction = "umap",
  cols = c("lightgrey", "purple")
) + ggtitle("Mid Marker: PPARG")

# 5. Late marker (ADIPOQ)
p5 <- FeaturePlot(
  merged_seurat,
  features = "hg19-ADIPOQ",
  reduction = "umap",
  cols = c("lightgrey", "red")
) + ggtitle("Late Marker: ADIPOQ")

# Combine as panel
combined_plot <- (p1 | p2) / (p3 | p4 | p5)

combined_plot

ggsave(
  filename = "results/figures/UMAP_combined_panel.png",
  plot = combined_plot,
  width = 13,
  height = 9,
  dpi = 300
)

# --- Save outputs ------------------------------------------------------------
if(!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)

write.csv(sig_matrix, "results/tables/sig_matrix.csv")
write.csv(top_markers, "results/tables/top_markers.csv")
saveRDS(merged_seurat, file = "results/SGBS_processed.rds")

message("✅ Analysis complete. Signature matrix and markers saved in results/tables/")
