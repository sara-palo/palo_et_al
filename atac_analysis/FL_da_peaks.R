source('scripts/atac_analysis/00.R')

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

DefaultAssay(FL) <- 'peaks'
Idents(FL) <- 'seurat_clusters'

for (ident in c(1:13)) {
  
  da_peaks_FL <- FindMarkers(
    object = FL,
    ident.1 = ident,
    ident.2 = 0,
    test.use = 'LR',
    latent.vars = 'nCount_peaks', 
    assay = 'peaks'
  )
  
  da_peaks_FL <- cbind(da_peaks_FL, ClosestFeature(FL, rownames(da_peaks_FL)))
  write.csv(da_peaks_FL, paste0(csv_path_da, '/FL_vs0_', ident, '_clust_da_peaks.csv'))
  
}

for (ident in 0:13) {
  
  da_peaks_FL <- FindMarkers(
    object = FL,
    ident.1 = ident,
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    assay = 'peaks'
  )
  
  da_peaks_FL <- cbind(da_peaks_FL, ClosestFeature(FL, rownames(da_peaks_FL)))
  write.csv(da_peaks_FL, paste0(csv_path_da, '/FL_vsall_', ident, '_clust_da_peaks.csv'))
  
}

da_pos_FL <- list()
da_neg_FL <- list()

for (ident in as.character(1:13)){
  
  da_peaks <- read.csv(paste0(csv_path_da, '/FL_vs0_', ident, '_clust_da_peaks.csv'), row.names = 1)
  
  da_pos_FL[[ident]] <- rownames(subset(da_peaks, subset = avg_log2FC > 2 & p_val_adj < 0.05 & pct.1 > 0.1)) 
  da_neg_FL[[ident]] <- rownames(subset(da_peaks, subset = avg_log2FC < -2 & p_val_adj < 0.05 & pct.2 > 0.1)) 
}

# make data frame with number of DA peaks per cluster for plots
da_up <- lapply(da_pos_FL, length) %>% unlist()
da_dn <- lapply(da_neg_FL, length) %>% unlist()

df_da <- rbind(data.frame('cluster' = names(da_up), 'count' = unname(da_up), 'up_or_dn' = 'up'), 
               data.frame('cluster' = names(da_dn), 'count' = unname(da_dn), 'up_or_dn' = 'down'))

df_da$cluster <- factor(df_da$cluster, levels = 1:13)
df_da$up_or_dn <- factor(df_da$up_or_dn, levels = c('up', 'down'))

ggplot(df_da, aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.7, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('Number of peaks up- or downregulated in FL compared cluster 0')
ggsave(paste0(fig_path, 'additional_plots/FLvsc0_da_peaks.pdf'), width = 9)

ggplot(subset(df_da, cluster %in% c('4', '8', '13')), aes(x = cluster, y = count, fill = up_or_dn)) + 
  geom_bar(stat = 'identity', width = 0.4, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))
ggsave(paste0(fig_path, 'sup3/FLvsc0_da_peaks_subs.pdf'), width = 6)

