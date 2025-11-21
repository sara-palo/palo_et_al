source('scripts/atac_analysis/00.R')

library(monocle3)

merged <- readRDS(paste0(rds_path, 'merged.rds'))

Idents(merged) <- 'clust.names'

# create cell dataset
ncds <- new_cell_data_set(expression_data = merged[['peaks']]@counts, cell_metadata = merged@meta.data)

# Transfer Seurat UMAP embeddings
ncds@int_colData@listData$reducedDims$UMAP <- merged@reductions$mergedUMAP@cell.embeddings
ncds@clusters$UMAP$clusters <- merged$clust.names
plot_cells(ncds)

# set partition to 1 for all cells
partition <- setNames(nm = ncds@colData@rownames, object = rep(1, length(ncds@colData@rownames))) %>% as.factor()
ncds@clusters$UMAP$partitions <- partition
ncds

ncds <- learn_graph(ncds, use_partition = T)
plot_cells(ncds, color_cells_by = "partition", group_label_size = 3.5)

# plot to be able to manually select nodes
plot_cells(
  cds = ncds,
  show_trajectory_graph = TRUE, 
  cell_size = 0.7, 
  label_principal_points = T
)

ncds <- order_cells(ncds, reduction_method = 'UMAP', root_pr_nodes = c('Y_447') ) # 
plot_cells(
  cds = ncds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, 
  cell_size = 0.7
) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/merged_monocle_ps.pdf'))


l2 <- choose_graph_segments(ncds, starting_pr_node = 'Y_447', ending_pr_nodes = 'Y_460', return_list = T) 
m1 <- choose_graph_segments(ncds, starting_pr_node = 'Y_447', ending_pr_nodes = 'Y_423', return_list = T)
m2 <- choose_graph_segments(ncds, starting_pr_node = 'Y_447', ending_pr_nodes = 'Y_440', return_list = T)
e1 <- choose_graph_segments(ncds, starting_pr_node = 'Y_447', ending_pr_nodes = 'Y_266', return_list = T)
e2 <- choose_graph_segments(ncds, starting_pr_node = 'Y_447', ending_pr_nodes = 'Y_45', return_list = T)

ps <- ncds@principal_graph_aux@listData$UMAP$pseudotime

# Add the pseudotimes for each trajectory to the cells in the seurat object
merged <- AddMetaData(object = merged,
                      metadata = ps[l2$cells],
                      col.name = 'l2_ps') %>% 
  AddMetaData(metadata = ps[m1$cells], col.name = 'm1_ps') %>% 
  AddMetaData(metadata = ps[m2$cells], col.name = 'm2_ps') %>% 
  AddMetaData(metadata = ps[e1$cells], col.name = 'e1_ps') %>% 
  AddMetaData(metadata = ps[e2$cells], col.name = 'e2_ps')

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

BM <- AddMetaData(object = BM,
                  metadata = ps[l2$cells],
                  col.name = 'merged_l2_ps') %>% 
  AddMetaData(metadata = ps[m1$cells], col.name = 'merged_m1_ps') %>%
  AddMetaData(metadata = ps[m2$cells], col.name = 'merged_m2_ps') %>% 
  AddMetaData(metadata = ps[e1$cells], col.name = 'merged_e1_ps') %>% 
  AddMetaData(metadata = ps[e2$cells], col.name = 'merged_e2_ps')

FL <- AddMetaData(object = FL,
                  metadata = ps[l2$cells],
                  col.name = 'merged_l2_ps') %>% 
  AddMetaData(metadata = ps[m1$cells], col.name = 'merged_m1_ps') %>% 
  AddMetaData(metadata = ps[m2$cells], col.name = 'merged_m2_ps') %>% 
  AddMetaData(metadata = ps[e1$cells], col.name = 'merged_e1_ps') %>% 
  AddMetaData(metadata = ps[e2$cells], col.name = 'merged_e2_ps')

saveRDS(FL, paste0(rds_path, 'FL_only.rds'))
saveRDS(BM, paste0(rds_path, 'BM_only.rds'))
saveRDS(merged, paste0(rds_path, 'merged.rds'))

xmin <- min(c(BM[['umap']]@cell.embeddings[,1], FL[['BMumap']]@cell.embeddings[,1]))
xmax <- max(c(BM[['umap']]@cell.embeddings[,1], FL[['BMumap']]@cell.embeddings[,1]))

ymin <- min(c(BM[['umap']]@cell.embeddings[,2], FL[['BMumap']]@cell.embeddings[,2]))
ymax <- max(c(BM[['umap']]@cell.embeddings[,2], FL[['BMumap']]@cell.embeddings[,2]))


FeaturePlot(merged, features = c('l2_ps'), pt.size = 0.6) + scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'fig4/merged_l2_ps_umap.pdf'))

FeaturePlot(merged, features = c('m1_ps'), pt.size = 0.6) + scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/merged_m1_ps_umap.pdf'))

FeaturePlot(merged, features = c('m2_ps'), pt.size = 0.6) + scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/merged_m2_ps_umap.pdf'))

FeaturePlot(merged, features = c('e1_ps'), pt.size = 0.6) + scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/merged_e1_ps_umap.pdf'))

FeaturePlot(merged, features = c('e2_ps'), pt.size = 0.6) + scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/merged_e2_ps_umap.pdf'))


FeaturePlot(BM, features = c('merged_l2_ps'), pt.size = 0.6) + 
  scale_color_viridis_c(limits = c(0, max(merged$l2_ps, na.rm = T))) + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'fig4/BM_l2_ps.pdf'))

FeaturePlot(BM, features = c('merged_m1_ps'), pt.size = 0.6) + 
  scale_color_viridis_c(limits = c(0, max(merged$m1_ps, na.rm = T))) + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/BM_m1_ps.pdf'))

FeaturePlot(BM, features = c('merged_e1_ps'), pt.size = 0.6) + 
  scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/BM_e1_ps.pdf'))

FeaturePlot(BM, features = c('merged_e2_ps'), pt.size = 0.6) + 
  scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/BM_e2_ps.pdf'))

FeaturePlot(BM, features = c('merged_m2_ps'), pt.size = 0.6) + 
  scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/BM_m2_ps.pdf'))


FeaturePlot(FL, features = c('merged_l2_ps'), pt.size = 0.6, reduction = 'BMumap') +
  scale_color_viridis_c(limits = c(0, max(merged$l2_ps, na.rm = T))) + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'fig4/FL_mergedumap_l2ps.pdf'))

FeaturePlot(FL, features = c('merged_m1_ps'), pt.size = 0.6, reduction = 'BMumap') + 
  scale_color_viridis_c(limits = c(0, max(merged$m1_ps, na.rm = T))) + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_m1ps.pdf'))

FeaturePlot(FL, features = c('merged_e1_ps'), pt.size = 0.6, reduction = 'BMumap') + 
  scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_e1ps.pdf'))

FeaturePlot(FL, features = c('merged_e2_ps'), pt.size = 0.6, reduction = 'BMumap') + 
  scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_e2ps.pdf'))

FeaturePlot(FL, features = c('merged_m2_ps'), pt.size = 0.6, reduction = 'BMumap') + 
  scale_color_viridis_c() + coord_fixed() + xlim(xmin, xmax) + ylim(ymin, ymax)
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_m2ps.pdf'))


FL <- ScaleData(FL, assay = 'cisbpchromvar')
FL$clust.names <- FL$pred.id.BM

BM <- ScaleData(BM, assay = 'cisbpchromvar')
BM$clust.names <- BM$BM_annot


tf_dotplot(FL, end_cluster = 'merged_l2_ps', assay = 'cisbpchromvar', features = c('RUNX1', 'TCF3', 'EBF1', 'PAX5', 'CEBPA'), 
           slot_to_use = 'scale.data', include.dots = F, xmax = max(BM$merged_l2_ps, na.rm = T), ymin = -4.1, ymax = 10) 
ggsave(paste0(fig_path, 'fig4/FL_mergedumap_l2ps_motifs_B.pdf'))

tf_dotplot(BM, end_cluster = 'merged_l2_ps', assay = 'cisbpchromvar', features = c('RUNX1','TCF3', 'EBF1', 'PAX5', 'CEBPA'), 
           slot_to_use = 'scale.data', include.dots = F,  xmax = max(BM$merged_l2_ps, na.rm = T, ymin = -4.1, ymax = 10)) 
ggsave(paste0(fig_path, 'fig4/BM_mergedumap_l2ps_motifs_B.pdf'))


cols <- annot_colors$cluster_cols[c('Ly-II','GMP-II')]
names(cols) <- c('BM.Ly.II_UCell', 'BM.GMP.II_UCell')

tf_dotplot(FL, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('BM.Ly.II_UCell', 'BM.GMP.II_UCell'), include.dots = F, 
           xmax = max(BM$merged_l2_ps, na.rm = T), colors = cols, ymin = 0, ymax = 1.1) 
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_l2ps_BMmodscores.pdf'))

tf_dotplot(FL, end_cluster = 'merged_m1_ps', assay = 'metadata', features = c('BM.Ly.II_UCell', 'BM.GMP.II_UCell'), include.dots = F, 
           xmax = max(c(BM$merged_m1_ps, FL$merged_m1_ps), na.rm = T), colors = cols, ymin = 0, ymax = 1.1) 
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_m1ps_BMmodscores.pdf'))


names(cols) <- c('Ly.II_UCell', 'GMP.II_UCell')

tf_dotplot(BM, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('Ly.II_UCell', 'GMP.II_UCell'), include.dots = F, 
           xmax = max(BM$merged_l2_ps, na.rm = T), colors = cols, ymin = 0, ymax = 1.1) 
ggsave(paste0(fig_path, 'additional_plots/BM_mergedumap_l2ps_BMmodscores.pdf'))

tf_dotplot(BM, end_cluster = 'merged_m1_ps', assay = 'metadata', features = c('Ly.II_UCell', 'GMP.II_UCell'), include.dots = F, 
           xmax = max(c(BM$merged_m1_ps, FL$merged_m1_ps), na.rm = T), colors = cols, ymin = 0, ymax = 1.1) 
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_m1ps_BMmodscores.pdf'))

tf_dotplot(FL, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('Ly.II_UCell', 'GMP.II_UCell'), include.dots = F, 
           xmax = max(BM$merged_l2_ps, na.rm = T), colors = cols, ymin = 0, ymax = 1.1) 
ggsave(paste0(fig_path, 'additional_plots/FL_mergedumap_l2ps_FLmodscores.pdf'))

names(cols) <- c('FL.Ly.II_UCell', 'FL.GMP.II_UCell')

tf_dotplot(BM, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('FL.Ly.II_UCell', 'FL.GMP.II_UCell'), include.dots = F, 
           xmax = max(BM$merged_l2_ps, na.rm = T), colors = cols, ymin = 0, ymax = 1.1) 
ggsave(paste0(fig_path, 'additional_plots/BM_mergedumap_l2ps_FLmodscores.pdf'))


# compare end points (final 10% of the pseudotime)
FL$l2_ps_cut <- cut(FL$merged_l2_ps, breaks = 10)
# BM$l2_ps_cut <- cut(BM$merged_l2_ps, breaks = 10)

levels(FL$l2_ps_cut)[10]

# choose the last range for FL pseudotime 
fl_r10 <- WhichCells(FL, expression = l2_ps_cut == levels(FL$l2_ps_cut)[10])

FL_df <- FL@meta.data[!is.na(FL$l2_ps_cut), ]
FL_df <- FL_df[fl_r10, c('BM.Ly.II_UCell', 'BM.GMP.II_UCell', 'merged_l2_ps')]
FL_df$diff <- FL_df$BM.Ly.II_UCell - FL_df$BM.GMP.II_UCell
hist(FL_df$diff)

# select BM cells in the same range of pseudotime as FL
bm_r10 <- WhichCells(BM, expression = merged_l2_ps >= min(FL_df$merged_l2_ps) & merged_l2_ps <= max(FL_df$merged_l2_ps))
BM_df <- BM@meta.data[!is.na(BM$merged_l2_ps), ]
BM_df <- BM_df[bm_r10, c('Ly.II_UCell', 'GMP.II_UCell', 'merged_l2_ps')]
BM_df$diff <- BM_df$Ly.II_UCell - BM_df$GMP.II_UCell
hist(BM_df$diff)

wc <- wilcox.test(FL_df$diff, BM_df$diff)

mean(FL_df$diff)
mean(BM_df$diff)

stats_df <- data.frame('FL_mean_diff' = mean(FL_df$diff), 'BM_mean_diff' = mean(BM_df$diff),
                       'FL_median_diff' = median(FL_df$diff), 'BM_median_diff' = median(BM_df$diff), 
                       'wilcox_pval' = wc$p.value)
write.csv(stats_df, paste0(csv_path, 'sup/pseudotime_modscores_FLvsBM_stat.csv'))

cols <- annot_colors$cluster_cols[c('Ly-II','GMP-II')]
names(cols) <- c('BM.Ly.II_UCell', 'BM.GMP.II_UCell')

tf_dotplot(FL, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('BM.Ly.II_UCell', 'BM.GMP.II_UCell'), include.dots = T, 
           xmax = max(BM$merged_l2_ps, na.rm = T), ymin = 0, ymax = 1.1, mark_range = range(FL_df$merged_l2_ps), 
           colors = cols)  
ggsave(paste0(fig_path, 'fig4/FL_l2ps_stat_range.pdf'), width = 8, height = 4.3)

tf_dotplot(FL, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('BM.Ly.II_UCell', 'BM.GMP.II_UCell'), include.dots = F, 
           xmax = max(BM$merged_l2_ps, na.rm = T), ymin = 0, ymax = 1.1, mark_range = range(FL_df$merged_l2_ps), 
           colors = cols)  
ggsave(paste0(fig_path, 'additional_plots/FL_l2ps_stat_range_nodots.pdf'), width = 8, height = 4.3)

names(cols) <- c('Ly.II_UCell', 'GMP.II_UCell')

tf_dotplot(BM, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('Ly.II_UCell', 'GMP.II_UCell'), include.dots = T, 
           xmax = max(BM$merged_l2_ps, na.rm = T), ymin = 0, ymax = 1.1, mark_range = range(FL_df$merged_l2_ps), colors = cols)
ggsave(paste0(fig_path, 'fig4/BM_l2ps_stat_range.pdf'), width = 8, height = 4.3)

tf_dotplot(BM, end_cluster = 'merged_l2_ps', assay = 'metadata', features = c('Ly.II_UCell', 'GMP.II_UCell'), include.dots = F, 
           xmax = max(BM$merged_l2_ps, na.rm = T), ymin = 0, ymax = 1.1, mark_range = range(FL_df$merged_l2_ps), colors = cols)
ggsave(paste0(fig_path, 'additional_plots/BM_l2ps_stat_range_nodots.pdf'), width = 8, height = 4.3)


plot_df <- rbind(FL_df[,3:4], BM_df[,3:4])
plot_df$trajectory <- c(rep('FL', nrow(FL_df)), rep('BM', nrow(BM_df)))
plot_df$trajectory <- factor(plot_df$trajectory, levels = c('FL', 'BM'))

ggplot(plot_df, aes(x = trajectory, y = diff, fill = trajectory)) + 
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9) + theme_classic() + 
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.1) +
  scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'sup4/FL_BM_pseudotime_score_diffs_box.pdf'), width = 4, height = 5)

# stats for PAX5
goi <- c('TCF3', 'EBF1', 'PAX5', 'RUNX1', 'CEBPA')
FL_df <- data.frame(t(FL[['cisbpchromvar']]@scale.data[goi, fl_r10]))
FL_df$l2_ps <- FL@meta.data[fl_r10, 'merged_l2_ps']

BM_df <- data.frame(t(BM[['cisbpchromvar']]@scale.data[goi, bm_r10]))
BM_df$l2_ps <- BM@meta.data[bm_r10, 'merged_l2_ps']

res_df <- data.frame()
tfs <- c('TCF3', 'EBF1', 'PAX5', 'RUNX1')

for (tf in tfs){
  
  FLdiff <- FL_df[, tf] - FL_df[, 'CEBPA']
  BMdiff <- BM_df[, tf] - BM_df[, 'CEBPA']
  
  plot_df <- data.frame('diff'= c(FLdiff, BMdiff), 'trajectory' = c(rep('FL', nrow(FL_df)), rep('BM', nrow(BM_df))))
  plot_df$trajectory <- factor(plot_df$trajectory, levels = c('FL', 'BM'))
  
  wc <- wilcox.test(FLdiff, BMdiff)
  
  if(wc$p.value < 0.05){
    p <- ggplot(plot_df, aes(x = trajectory, y = diff, fill = trajectory)) + 
      geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9) + theme_classic() + 
      geom_jitter(width = 0.2, alpha = 0.5, size = 1.1) +
      scale_fill_manual(values = annot_colors$tissue_cols)
    print(p)
    ggsave(plot = p, filename = paste0(fig_path, 'sup4/FL_BM_pseudotime_', tf, '_diffs_box.pdf'), width = 4, height = 5)
  }
  
  res_df <- rbind(res_df, data.frame('FL_mean_diff' = mean(FLdiff), 'BM_mean_diff' = mean(BMdiff),
                                     'FL_median_diff' = median(FLdiff), 'BM_median_diff' = median(BMdiff), 
                                     'wilcox_pval' = wc$p.value))
}

rownames(res_df) <- tfs
res_df$p_adj <- p.adjust(res_df$wilcox_pval, method = 'bonf')

write.csv(res_df, paste0(csv_path, 'sup/pseudotime_motifs_FLvsBM_stat.csv'))


tf_dotplot(FL, end_cluster = 'merged_l2_ps', assay = 'cisbpchromvar', features = c('TCF3', 'EBF1', 'PAX5', 'CEBPA', 'RUNX1'), 
           slot_to_use = 'scale.data', include.dots = F, xmax = max(BM$merged_l2_ps, na.rm = T), 
           mark_range = range(FL_df$l2_ps), ymin = -4.1, ymax = 10) 
ggsave(paste0(fig_path, 'fig4/FL_l2ps_motifs_B_markrange.pdf'))

tf_dotplot(BM, end_cluster = 'merged_l2_ps', assay = 'cisbpchromvar', features = c('TCF3', 'EBF1', 'PAX5', 'CEBPA', 'RUNX1'), 
           slot_to_use = 'scale.data', include.dots = F,  xmax = max(BM$merged_l2_ps, na.rm = T), 
           mark_range = range(FL_df$l2_ps), ymin = -4.1, ymax = 10)
ggsave(paste0(fig_path, 'fig4/BM_l2ps_motifs_B_markrange.pdf'))

