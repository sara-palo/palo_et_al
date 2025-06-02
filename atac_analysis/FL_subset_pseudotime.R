source('scripts/atac_analysis/00.R')

library(harmony)
library(monocle3)
library(dendextend)
library(dendsort)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

Idents(FL) <- 'seurat_clusters'
DefaultAssay(FL) <- 'peaks'

# subset for primitive, ly and my clusters and make new UMAP
subs <- subset(FL, idents = c(0:8))
subs

# no features with 0 counts
table(rowSums(GetAssayData(object = subs, assay = 'peaks', slot = "counts")) > 0)

subs <- RunTFIDF(subs, assay = 'peaks')
subs <- FindTopFeatures(subs, assay = 'peaks')
subs <- RunSVD(subs, assay = 'peaks')

subs <-
  RunHarmony(
    subs,
    group.by.vars = "sample_number",
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F 
  )

subs <- RunUMAP(subs, reduction = 'harmony', dims = 2:20, return.model = TRUE, min.dist = 0.3, 
                spread = 0.9, n.neighbors = 30)

# invert UMAP coordinates
subs[['umap']]@cell.embeddings[,'UMAP_1'] <- - (subs[['umap']]@cell.embeddings[,'UMAP_1'])
subs[['umap']]@cell.embeddings[,'UMAP_2'] <- - (subs[['umap']]@cell.embeddings[,'UMAP_2'])

DimPlot(subs, label = TRUE, pt.size = 0.5, group.by = 'sample_number', label.size = 0, shuffle = T, cols = annot_colors$sample_numberFL) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/FLsubs_umap_harmony.pdf'))

DimPlot(subs, label = TRUE, pt.size = 0.6, group.by = 'seurat_clusters', cols = annot_colors$seurat_cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'sup3/FLsubs_umap.pdf'))

ncds <- new_cell_data_set(expression_data = subs[['peaks']]@counts, cell_metadata = subs@meta.data)

# Transfer UMAP embeddings
ncds@int_colData@listData$reducedDims$UMAP <- subs@reductions$umap@cell.embeddings
ncds@clusters$UMAP$clusters <- subs$seurat_clusters

# set partition to 1 for all cells
partition <- setNames(nm = ncds@colData@rownames, object = rep(1, length(ncds@colData@rownames))) %>% as.factor()
ncds@clusters$UMAP$partitions <- partition
ncds

ncds <- learn_graph(ncds, use_partition = T)

# plot to be able to manually select nodes
plot_cells(
  cds = ncds,
  show_trajectory_graph = TRUE, 
  cell_size = 0.7, 
  label_principal_points = T
)

# manually select root
ncds <- order_cells(ncds, reduction_method = 'UMAP')  
plot_cells(
  cds = ncds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, 
  cell_size = 0.7
) + coord_fixed()
ggsave(paste0(fig_path, 'sup3/FLsubs_monocle_ps.pdf'))

# manual node selection
l2 <- choose_graph_segments(ncds, return_list = T) 
m1 <- choose_graph_segments(ncds, return_list = T)

ps <- ncds@principal_graph_aux@listData$UMAP$pseudotime

# Add the pseudotimes for each trajectory to the cells in the seurat object
subs <- AddMetaData(object = subs,
                    metadata = ps[l2$cells],
                    col.name = 'l2_ps') %>% 
  AddMetaData(metadata = ps[m1$cells], col.name = 'm1_ps') 

(FeaturePlot(subs, features = c('l2_ps'), reduction = 'umap') + scale_color_viridis_c() + coord_fixed()) | 
  (FeaturePlot(subs, features = c('m1_ps'), reduction = 'umap') + scale_color_viridis_c() + coord_fixed()) 

FeaturePlot(subs, features = c('l2_ps'), reduction = 'umap') + scale_color_viridis_c() + coord_fixed()
ggsave(paste0(fig_path, 'fig3/FLsubs_l2_ps.pdf'))

FeaturePlot(subs, features = c('m1_ps'), reduction = 'umap') + scale_color_viridis_c() + coord_fixed()
ggsave(paste0(fig_path, 'fig3/FLsubs_m1_ps.pdf'))

saveRDS(subs, paste0(rds_path, 'FLsubs_monocle.rds'))

subs$clust.names <- subs$seurat_clusters
cols <- annot_colors$seurat_cluster_cols[c('4','8')]
names(cols) <- c('C4_UCell', 'C8_UCell')

tf_dotplot(subs, end_cluster = 'l2_ps', assay = 'metadata', features = c('C4_UCell', 'C8_UCell'), colors = cols, ymax = 1.1)
ggsave(paste0(fig_path, 'fig3/FLsubs_l2_ps_modscores.pdf'))

tf_dotplot(subs, end_cluster = 'm1_ps', assay = 'metadata', features = c('C4_UCell', 'C8_UCell'), colors = cols, ymax = 1.1)
ggsave(paste0(fig_path, 'fig3/FLsubs_m1_ps_modscores.pdf'))


subs$l2_ps_cut <- cut(subs$l2_ps, breaks = 10)
subs$m1_ps_cut <- cut(subs$m1_ps, breaks = 10)

l2_df <- subs@meta.data[!is.na(subs$l2_ps_cut), ]
l2_df <- l2_df[l2_df$l2_ps_cut == levels(subs$l2_ps_cut)[10], c('C4_UCell', 'C8_UCell', 'l2_ps')]
l2_df$diff <- l2_df$C4_UCell - l2_df$C8_UCell
hist(l2_df$diff)

m1_df <- subs@meta.data[!is.na(subs$m1_ps_cut), ]
m1_df <- m1_df[m1_df$m1_ps_cut == levels(subs$m1_ps_cut)[10], c('C4_UCell', 'C8_UCell', 'm1_ps')]
m1_df$diff <- m1_df$C8_UCell - m1_df$C4_UCell
hist(m1_df$diff)

wc <- wilcox.test(l2_df$diff, m1_df$diff)

mean(l2_df$diff)
mean(m1_df$diff)
median(l2_df$diff)
median(m1_df$diff)

stats_df <- data.frame('ly_ps_mean_diff' = mean(l2_df$diff), 'my_ps_mean_diff' = mean(m1_df$diff), 
                       'ly_ps_median_diff' = median(l2_df$diff), 'my_ps_median_diff' = median(m1_df$diff), 'wilcox_pval' = wc$p.value)
write.csv(stats_df, paste0(csv_path, 'sup/pseudotime_modscores_FL_Ly_my_stat.csv'))

levels(subs$l2_ps_cut)[10]
levels(subs$m1_ps_cut)[10]

tf_dotplot(subs, end_cluster = 'l2_ps', assay = 'metadata', features = c('C4_UCell', 'C8_UCell'), ymax = 1.1, mark_range = c(min(l2_df$l2_ps), max(l2_df$l2_ps)), colors = cols) 
ggsave(paste0(fig_path, 'fig3/FL_ly_ps_stat_ranges.pdf'))

tf_dotplot(subs, end_cluster = 'm1_ps', assay = 'metadata', features = c('C4_UCell', 'C8_UCell'), ymax = 1.1, mark_range = c(min(m1_df$m1_ps), max(m1_df$m1_ps)), colors = cols)
ggsave(paste0(fig_path, 'fig3/FL_my_ps_stat_ranges.pdf'))

plot_df <- rbind(l2_df[,c(1:2,4)], m1_df[,c(1:2,4)])
plot_df$trajectory <- c(rep('4', nrow(l2_df)), rep('8', nrow(m1_df)))
plot_df$trajectory <- as.factor(plot_df$trajectory)

ggplot(plot_df, aes(x = trajectory, y = diff, fill = trajectory)) + 
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9) + theme_classic() + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.1) +
  scale_fill_manual(values = annot_colors$seurat_cluster_cols)
ggsave(paste0(fig_path, 'sup3/FL_pseudotime_diffs_box.pdf'), width = 4, height = 5)


subs <- ScaleData(subs, assay = 'cisbpchromvar')

tf_dotplot(subs, end_cluster = 'l2_ps', assay = 'cisbpchromvar', features = c('HOXB7', 'TCF3', 'EBF1', 'SPIB'), 
           slot_to_use = 'scale.data', include.dots = F, ymin = -3.5, ymax = 7) 
ggsave(paste0(fig_path, 'fig3/FLsubs_l2_ps_TFs.pdf'))

tf_dotplot(subs, end_cluster = 'm1_ps', assay = 'cisbpchromvar', features = c('HOXB7', 'CEBPA', 'EBF1', 'SPIB'), 
           slot_to_use = 'scale.data', include.dots = F, ymin = -3.5, ymax = 7) 
ggsave(paste0(fig_path, 'fig3/FLsubs_m1_ps_TFs.pdf'))

up_FL <- read.csv(paste0(csv_path_ch, '/vs_prim/FL_perclust_up.csv'))
dn_FL <- read.csv(paste0(csv_path_ch, '/vs_prim/FL_perclust_dn.csv'))

subs2 <- subset(subs, cells = names(na.omit(subs$l2_ps)))
DimPlot(subs2, reduction = 'umap', cols = annot_colors$seurat_cluster_cols)
VlnPlot(subs2, 'l2_ps', group.by = 'seurat_clusters')
subs2 <- ScaleData(subs2, assay = 'cisbpchromvar')

annot_df <- data.frame('cluster' = subs2$seurat_clusters[order(subs2$l2_ps)])
column_ha <- HeatmapAnnotation(df = annot_df, col = list(
  'cluster' = annot_colors$seurat_cluster_cols[as.character(unique(subs2$seurat_clusters))]
))

mtx <- subs2[['cisbpchromvar']]@scale.data[na.omit(c(dn_FL$C4, up_FL$C4)), order(subs2$l2_ps)]
Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = T)

row_dend = dendsort(hclust(dist(mtx)), isReverse = T, type = 'average')
hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
              top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = F, cluster_rows = row_dend)

pdf(paste0(fig_path, 'sup3/FLsubs_l2_up_dn_C4_scaled.pdf'))
hm
dev.off()

mtx <- subs2[['cisbpchromvar']]@scale.data[na.omit(c(dn_FL$C8, up_FL$C8)), order(subs2$l2_ps)]
Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = T)

row_dend = dendsort(hclust(dist(mtx)), isReverse = T, type = 'average')
hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
              top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = F, cluster_rows = row_dend)

pdf(paste0(fig_path, 'additional_plots/FLsubs_l2_up_dn_C8_scaled.pdf'))
hm
dev.off()


subs2 <- subset(subs, cells = names(na.omit(subs$m1_ps)))
DimPlot(subs2, reduction = 'umap', cols = annot_colors$seurat_cluster_cols)
VlnPlot(subs2, 'm1_ps', group.by = 'seurat_clusters')
subs2 <- ScaleData(subs2, assay = 'cisbpchromvar')

annot_df <- data.frame('cluster' = subs2$seurat_clusters[order(subs2$m1_ps)])
column_ha <- HeatmapAnnotation(df = annot_df, col = list(
  'cluster' = annot_colors$seurat_cluster_cols[as.character(unique(subs2$seurat_clusters))]
))

mtx <- subs2[['cisbpchromvar']]@scale.data[na.omit(c(dn_FL$C8, up_FL$C8)), order(subs2$m1_ps)]
Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = T)

row_dend = dendsort(hclust(dist(mtx)), isReverse = T, type = 'average')
hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
              top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = F, cluster_rows = row_dend)

pdf(paste0(fig_path, 'additional_plots/FLsubs_m1_up_dn_C8_scaled.pdf'))
hm
dev.off()

mtx <- subs2[['cisbpchromvar']]@scale.data[na.omit(c(dn_FL$C4, up_FL$C4)), order(subs2$m1_ps)]
Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = T)

row_dend = dendsort(hclust(dist(mtx)), isReverse = F, type = 'average')
hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
              top_annotation = column_ha, row_names_gp = gpar(fontsize = 4), row_dend_reorder = F, cluster_rows = row_dend)

pdf(paste0(fig_path, 'additional_plots/FLsubs_m1_up_dn_C4_scaled.pdf'))
hm
dev.off()




