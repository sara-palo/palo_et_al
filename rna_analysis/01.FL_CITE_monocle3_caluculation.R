# Load required libraries
library(Seurat)
library(Matrix)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(viridis)
library(ggplot2)
set.seed(123)

PRJ <- "C:/Users/keiki/OneDrive/デスクトップ/Sara_project"
setwd(PRJ)
out_dir <- sprintf("output_%s", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir

FL <- readRDS("rds/FL_CITE.rds")

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
#gmp <- choose_graph_segments(cds, return_list = TRUE)
gmp <- choose_graph_segments(cds, starting_pr_node = 'Y_243', ending_pr_nodes = 'Y_256', return_list = TRUE)
      gmp2 <- choose_graph_segments(cds, starting_pr_node = 'Y_243', ending_pr_nodes = 'Y_143', return_list = TRUE)
# identical(gmp,gmp2)
# same_cells <- setequal(gmp$cells,  gmp2$cells)
# same_nodes <- setequal(gmp$pr_nodes, gmp2$pr_nodes) 
# c(same_cells = same_cells, same_nodes = same_nodes)

# Extract pseudotime per trajectory
l2_ps <- pseudotime_values[l2$cells]
gmp_ps <- pseudotime_values[gmp$cells]
      gmp2_ps <- pseudotime_values[gmp2$cells]

# Add to Seurat object
FL <- FL %>%
  AddMetaData(metadata = l2_ps, col.name = 'l2_ps') %>%
  AddMetaData(metadata = gmp_ps, col.name = 'gmp_ps')%>%
      AddMetaData(metadata = gmp2_ps, col.name = 'gmp2_ps')


# Plot trajectories
pdf(paste0(out_dir,"/trajectories_FL_CITE.pdf"),width=12,height = 4)
FeaturePlot(FL, features = c('l2_ps'), reduction = 'umap', pt.size = 0.5) +
  scale_color_viridis_c() + coord_fixed() |
  FeaturePlot(FL, features = c('gmp_ps'), reduction = 'umap', pt.size = 0.5) +
  scale_color_viridis_c() + coord_fixed() |
  FeaturePlot(FL, features = c('gmp2_ps'), reduction = 'umap', pt.size = 0.5) +
  scale_color_viridis_c() + coord_fixed()
dev.off()

# Save Seurat object with pseudotime
# saveRDS(FL, file = "FL_CITE_trajectory_gmp2.rds")
# identical(FL, FL_CITE_trajectory)

          # library(dplyr)
          # library(ggrepel)
          # library(igraph)
          # 
          # cv <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
          # cv <- setNames(as.character(cv[colnames(cds), 1]), colnames(cds))
          # umap <- reducedDims(cds)$UMAP
          # colnames(umap) <- c("UMAP_1", "UMAP_2")
          # 
          # df_map <- data.frame(
          #   cell   = colnames(cds),
          #   node   = cv[colnames(cds)],
          #   UMAP_1 = umap[,1],
          #   UMAP_2 = umap[,2],
          #   stringsAsFactors = FALSE
          # )
          # node_pos <- df_map %>%
          #   group_by(node) %>%
          #   summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), n = n(), .groups="drop")
          # 
          # g <- monocle3::principal_graph(cds)[["UMAP"]]
          # deg <- igraph::degree(g)
          # node_pos$type <- case_when(
          #   node_pos$node %in% names(deg[deg == 1]) ~ "leaf",
          #   node_pos$node %in% names(deg[deg >= 3]) ~ "branch",
          #   TRUE                                     ~ "middle"
          # )
          # p <- plot_cells(cds,
          #                 label_principal_points = FALSE,
          #                 label_branch_points    = FALSE,
          #                 label_leaves           = FALSE)
          # p +
          #   geom_point(data = node_pos, aes(UMAP_1, UMAP_2), size = 1.5) +
          #   ggrepel::geom_text_repel(
          #     data = node_pos,
          #     aes(UMAP_1, UMAP_2, label = node),
          #     size = 3, max.overlaps = Inf, box.padding = 0.3,
          #     point.padding = 0.2, segment.size = 0.2
          #   )



# #============load
# library(tradeSeq)
# library(RColorBrewer)
# library(SingleCellExperiment)
# library(monocle3)
#       plot_cells(cds)
#       cells_lym <- l2$cells
#       cells_mye <- gmp$cells
#       cells_all <- union (cells_lym,cells_mye)
#       pt <- pseudotime(cds)
#       keep <- names(pt)[!is.na(pt)]
#       cells_keep <- intersect(cells_all, keep)
# 
# counts <- GetAssayData(FL, assay = "RNA", slot = "counts")[, cells_keep, drop = FALSE]
# co <- colnames(counts)
# 
# # 4) 所属フラグ
# in_lym <- co %in% cells_lym
# in_mye <- co %in% cells_mye
# both   <- which(in_lym & in_mye)
# 
# # 5) pseudotime 行列（NA禁止：非所属は 0 で埋める）
# PT <- cbind(
#   lym = ifelse(in_lym, pt[co], 0),
#   mye = ifelse(in_mye, pt[co], 0)
# )
# rownames(PT) <- co
# 
# # 6) cellWeights 行列
# W <- cbind(
#   lym = as.numeric(in_lym),
#   mye = as.numeric(in_mye)
# )
# if (length(both) > 0) {
#   W[both, "lym"] <- 0.5
#   W[both, "mye"] <- 0.5
# }
# rownames(W) <- co





# 
# BPPARAM <-BiocParallel::bpparam()
# BPPARAM$workers <- 2
# genes = VariableFeatures(FL)
# conditions<-factor(FL$clust.names)
# set.seed(123)
#       keep_genes <- Matrix::rowSums(counts >= 3) >= 10
#       counts_f <- counts[keep_genes, , drop = FALSE]
#       
#       set.seed(123)
#       pick <- sample(which(keep_genes), size = min(200, sum(keep_genes)))
# 
#       knot <- evaluateK(counts = counts_f[pick, , drop = FALSE],
#                         pseudotime = PT,
#                         cellWeights = W,
#                         k = 7:9,
#                         nGenes = length(pick),
#                         verbose = TRUE,
#                         BPPARAM = BPPARAM)
#       set.seed(123)
#       gam <- fitGAM(counts = counts,
#                     pseudotime = PT,
#                     cellWeights = W,
#                     genes=genes,
#                     nknots = 7, verbose = TRUE,
#                     BPPARAM = BPPARAM)
# 
#       #saveRDS(gam,file = "rds/FL_monocle_gam.rds")
# gam_fl <- gam
# gam_fl <- readRDS("FL_monocle_gam.rds")
#         colnames(PT)
#         metadata(gam_fl)$lineageNames <- colnames(PT)
#         assocRes <- associationTest(gam_fl, lineages = TRUE, l2fc = log2(2))
#         ln <- colnames(PT)
#         for (i in seq_along(ln)) {
#           old <- paste0("pvalue_", i)
#           new <- paste0("pvalue_", ln[i])
#           if (old %in% colnames(assocRes)) {
#             colnames(assocRes)[colnames(assocRes) == old] <- new
#           }
#         }
#         head(assocRes)
#         
#         
#         LYM_Genes_FL <-  rownames(assocRes)[
#           which(p.adjust(assocRes$pvalue_lym, "fdr") <= 0.05)
#         ]
#         GMP_Genes_FL <-  rownames(assocRes)[
#           which(p.adjust(assocRes$pvalue_mye, "fdr") <= 0.05)
#         ]
#         
#         library(UpSetR)
#         UpSetR::upset(fromList(list(FL_myeloid = GMP_Genes_FL , FL_lymphoid = LYM_Genes_FL)),mainbar.y.max = 400)
#         
#         
# # weights を Seurat に載せて可視化
# FL$w_lym <- 0; FL$w_mye <- 0
# FL$w_lym[colnames(counts)] <- W[, "lym"]
# FL$w_mye[colnames(counts)] <- W[, "mye"]
# 
# p_w1 <- FeaturePlot(FL, features = "w_lym", reduction = "umap") +
#   scale_color_viridis_c() + ggtitle("lineage 1 (= lym) weight")
# p_w2 <- FeaturePlot(FL, features = "w_mye", reduction = "umap") +
#   scale_color_viridis_c() + ggtitle("lineage 2 (= mye) weight")
# p_w1 | p_w2
# 
# 
# startRes <- startVsEndTest(gam_fl)
# oStart <- order(startRes$waldStat, decreasing = TRUE)
# sce <- as.SingleCellExperiment(FL)
# sigGeneStart <- names(sce)[oStart[3]]
# plotSmoothers(sce, counts, gene = sigGeneStart)
# 
# 
