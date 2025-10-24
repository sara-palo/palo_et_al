source('scripts/atac_analysis/00.R')

library(dendextend)
library(dendsort)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

FL_up <- read.csv(paste0(csv_path_ch, '/vs_prim/FL_perannot_up.csv'))
BM_up <- read.csv(paste0(csv_path_ch, '/vs_prim/BM_perclust_up.csv'))

FL_ps <- names(which(! is.na(FL$merged_l2_ps)))
DimPlot(FL, reduction = 'BMumap', cells.highlight = FL_ps)

table(FL$pred.id.BM[FL_ps]) / table(FL$pred.id.BM)

BM_ps <- names(which(! is.na(BM$merged_l2_ps)))
DimPlot(BM, cells.highlight = BM_ps)

table(BM$BM_annot[BM_ps]) / table(BM$BM_annot)

hm_motifs <- unique(na.omit(c(FL_up$Ly.I, FL_up$Ly.II, BM_up$Ly.I, BM_up$Ly.II)))

FL$clust.names <- FL$pred.id.BM
BM$clust.names <- BM$BM_annot

subsFL <- subset(FL, cells = names(na.omit(FL$merged_l2_ps)))
DimPlot(subsFL, reduction = 'BMumap', cols = annot_colors$cluster_cols)
subsFL <- ScaleData(subsFL, assay = 'cisbpchromvar')

subsBM <- subset(BM, cells = names(na.omit(BM$merged_l2_ps)))
DimPlot(subsBM, reduction = 'umap', cols = annot_colors$cluster_cols)
subsBM <- ScaleData(subsBM, assay = 'cisbpchromvar')

annot_df_BM <- data.frame('cluster' = subsBM$clust.names[order(subsBM$merged_l2_ps)])
column_ha_BM <- HeatmapAnnotation(df = annot_df_BM, col = list(
  'cluster' = annot_colors$cluster_cols[as.character(unique(subsBM$clust.names))]
))

annot_df_FL <- data.frame('cluster' = subsFL$clust.names[order(subsFL$merged_l2_ps)])
column_ha_FL <- HeatmapAnnotation(df = annot_df_FL, col = list(
  'cluster' = annot_colors$cluster_cols[as.character(unique(subsFL$clust.names))]
))



mtx <- subsBM[['cisbpchromvar']]@scale.data[hm_motifs, order(subsBM$merged_l2_ps)]

Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha_BM, row_names_gp = gpar(fontsize = 6), cluster_rows = F, use_raster = F)

hcl <- hclust(dist(mtx))
row_dend <- dendsort(hclust(dist(mtx)), isReverse = T, type = 'min')

hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
              top_annotation = column_ha_BM, row_names_gp = gpar(fontsize = 6), cluster_rows = row_dend, show_row_dend = T, use_raster = F)

hm <- draw(hm)
hm

pdf(paste0(fig_path, 'sup4/BMsubs_l2_lymotifs_scaled_withdend.pdf'))
hm
dev.off()

row_order <- row_order(hm)

mtxBM <- subsBM[['cisbpchromvar']]@scale.data[hm_motifs[row_order], order(subsBM$merged_l2_ps)]

pdf(paste0(fig_path, 'sup4/BMsubs_l2_lymotifs_scaled.pdf'))
Heatmap(mtxBM, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha_BM, row_names_gp = gpar(fontsize = 6), cluster_rows = F, show_row_dend = F, use_raster = F)
dev.off()

# plot FL using the same row order
mtxFL <- subsFL[['cisbpchromvar']]@scale.data[hm_motifs[row_order], order(subsFL$merged_l2_ps)]

pdf(paste0(fig_path, 'sup4/FLsubs_l2_lymotifs_scaled.pdf'))
Heatmap(mtxFL, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha_FL, row_names_gp = gpar(fontsize = 6), cluster_rows = F, show_row_dend = F, use_raster = F)
dev.off()



mtxBM <- subsBM[['cisbpchromvar']]@data[hm_motifs[row_order], order(subsBM$merged_l2_ps)]

pdf(paste0(fig_path, 'additional_plots/BMsubs_l2_lymotifs_unscaled.pdf'))
Heatmap(mtxBM, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha_BM, row_names_gp = gpar(fontsize = 6), cluster_rows = F, show_row_dend = F)
dev.off()

# plot FL using the same row order
mtxFL <- subsFL[['cisbpchromvar']]@data[hm_motifs[row_order], order(subsFL$merged_l2_ps)]

pdf(paste0(fig_path, 'additional_plots/FLsubs_l2_lymotifs_unscaled.pdf'))
Heatmap(mtxFL, cluster_columns = F, show_row_names = T, show_column_names = F, heatmap_legend_param = list(title = 'TFs'),
        top_annotation = column_ha_FL, row_names_gp = gpar(fontsize = 6), cluster_rows = F, show_row_dend = F)
dev.off()




