# Load required libraries
library(Seurat)
library(Matrix)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(viridis)
library(ggplot2)
set.seed(123)

# Set working directory and load Seurat object
setwd("C:/Users/keiki/OneDrive/デスクトップ/Sara_project")
FL <- readRDS("FL_CITE.rds")

# Set assay and cluster identity
DefaultAssay(FL) <- "RNA"
Idents(FL) <- "clust.names"

# Convert to CellDataSet for Monocle3
cds <- as.cell_data_set(FL)
cds <- cluster_cells(cds)
fData(cds)$gene_short_name <- rownames(fData(cds))

# Force single partition
partition <- factor(rep(1, length(cds@colData@rownames)))
names(partition) <- cds@colData@rownames
cds@clusters$UMAP$partitions <- partition

# Set clusters and UMAP
cds@clusters$UMAP$clusters <- FL@active.ident
cds@int_colData@listData$reducedDims$UMAP <- FL@reductions$umap@cell.embeddings

# Learn graph and plot
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds,
           color_cells_by = "clust.names",
           label_cell_groups = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 1.5)

# Function to find root node
get_earliest_principal_node <- function(cds, time_bin = "HSC") {
  cell_ids <- which(colData(cds)[, "clust.names"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))
  ]
  return(root_pr_nodes)
}

# Order cells
root_pr_nodes <- get_earliest_principal_node(cds)
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

# Plot pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 1.5,
           label_principal_points = TRUE)

# Add pseudotime to Seurat object
pseudotime_values <- pseudotime(cds)
FL$monocle3 <- pseudotime_values

# Identify nodes for specific clusters
cluster_names <- c("Ly-II", "HSC", "GMP")
for (cluster_name in cluster_names) {
  cells <- which(colData(cds)$"clust.names" == cluster_name)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  node_counts <- table(closest_vertex[cells, ])
  selected_node <- names(which.max(node_counts))
  print(paste("Selected node for", cluster_name, ":", selected_node))
}

# Choose graph segments
l2 <- choose_graph_segments(cds, starting_pr_node = 'Y_243', ending_pr_nodes = 'Y_172', return_list = TRUE)
gmp <- choose_graph_segments(cds, return_list = TRUE)

# Extract pseudotime per trajectory
l2_ps <- pseudotime_values[l2$cells]
gmp_ps <- pseudotime_values[gmp$cells]

# Add to Seurat object
FL <- FL %>%
  AddMetaData(metadata = l2_ps, col.name = 'l2_ps') %>%
  AddMetaData(metadata = gmp_ps, col.name = 'gmp_ps')

# Plot trajectories
pdf("trajectories_FL_CITE.pdf")
FeaturePlot(FL, features = c('l2_ps'), reduction = 'umap', pt.size = 1) +
  scale_color_viridis_c() + coord_fixed() |
  FeaturePlot(FL, features = c('gmp_ps'), reduction = 'umap', pt.size = 1) +
  scale_color_viridis_c() + coord_fixed()
dev.off()

# Save Seurat object with pseudotime
# saveRDS(FL, file = "FL_CITE_trajectory.rds")
# identical(FL, FL_CITE_trajectory)
