source('scripts/atac_analysis/00.R')

library(ComplexHeatmap)
library(dendextend)
library(dendsort)
library(openxlsx)

merged <- readRDS(paste0(rds_path, 'merged.rds'))

Idents(merged) <- 'clust.names'
table(merged$clust.names)
DefaultAssay(merged) <- 'peaks'

da_pos_FL <- list()
da_neg_FL <- list()

for (ident in unique(merged$clust.names)){
  
  if(ident == 'Ba/M'){
    da_peaks <- read.csv(paste0(csv_path_da, '/FLvsBM/', 'BaM', '.csv'), row.names = 1)
  } else {
    da_peaks <- read.csv(paste0(csv_path_da, '/FLvsBM/', ident, '.csv'), row.names = 1)
  }  
  
  da_pos_FL[[ident]] <- rownames(subset(da_peaks, subset = avg_log2FC > 2 & p_val_adj < 0.05)) 
  da_neg_FL[[ident]] <- rownames(subset(da_peaks, subset = avg_log2FC < -2 & p_val_adj < 0.05)) 
}

# make cellradar plot for genes associated with DA peaks
up <- da_pos_FL[['HSC']]
dn <- da_neg_FL[['HSC']]

up_feat <- ClosestFeature(merged, regions = up)
dn_feat <- ClosestFeature(merged, regions = dn)

up_genes <- unique(up_feat[up_feat$distance < 2000, 'gene_name'])
dn_genes <- unique(dn_feat[dn_feat$distance < 2000, 'gene_name'])

print(table(up_genes %in% dn_genes))

res1 <- prep_data('human', up_genes)
res2 <- prep_data('human', dn_genes)

orca(p = plot_cellradar(res_list = list(res1, res2), colors = unlist(annot_colors$tissue_cols, use.names = F),
                        names = c('FL', 'BM'), dash = c('solid', '4px')), 
     file = paste0(fig_path, 'fig2/da_peaks_genes_HSC_cr.pdf'), more_args = '--disable-gpu')

table(up_genes %in% res1$nonexisting_genes)
table(dn_genes %in% res2$nonexisting_genes)


df_list <- list(up_feat, dn_feat)
names(df_list) = c('up in FL HSCs', 'up in BM HSCs')

wb <- createWorkbook('table S2')

for(name in names(df_list)){
  
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = df_list[[name]], keepNA = F)
}
saveWorkbook(wb, file = paste0(csv_path, 'sup/tableS2.xlsx'), overwrite = TRUE)

meta.feature <- GetAssayData(merged, assay = 'peaks', slot = 'meta.features')

enr.res.FL <- list()
enr.motifs.FL <- list()
enr.res.BM <- list()
enr.motifs.BM <- list()

for (name in names(da_pos_FL)){
  
  cells <- WhichCells(merged, expression = clust.names == name)
  open.peaks <- AccessiblePeaks(merged, cells = cells, min.cells = length(cells)*0.1)
  
  regs1 <- da_pos_FL[[name]]
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ], 
    query.feature = meta.feature[regs1, ],
    n = 20000 
  )
  
  enriched.motifs <- FindMotifs(
    object = merged,
    features = regs1, 
    background = peaks.matched
  )
  enr.res.FL[[name]] <- enriched.motifs
  enr.motifs.FL[[name]] <- rownames(subset(enriched.motifs, subset = fold.enrichment > 1.5 & p.adjust < 0.05)) 
  
  regs2 <- da_neg_FL[[name]]
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ], 
    query.feature = meta.feature[regs2, ],
    n = 20000 
  )
  
  enriched.motifs <- FindMotifs(
    object = merged,
    features = regs2, 
    background = peaks.matched
  )
  enr.res.BM[[name]] <- enriched.motifs
  enr.motifs.BM[[name]] <- rownames(subset(enriched.motifs, subset = fold.enrichment > 1.5 & p.adjust < 0.05)) 
}

hsc.motifs <- unique(c(enr.motifs.FL[['HSC']], enr.motifs.BM[['HSC']]))

HSC.df <- data.frame('FL' = enr.res.FL[['HSC']][hsc.motifs, 'fold.enrichment'], 
                     'BM' = enr.res.BM[['HSC']][hsc.motifs, 'fold.enrichment'], 
                     row.names = hsc.motifs)

mtx <- as.matrix(HSC.df)
Heatmap(mtx)

hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = T, heatmap_legend_param = list(title = 'Fold enrichment'), column_names_side = 'top', column_names_rot = 0,
              column_names_centered = T, row_names_gp = gpar(fontsize = 10), row_dend_reorder = F)

pdf(paste0(fig_path, 'fig2/da_peak_motifs_HSC_w4.pdf'), width = 4)
hm
dev.off()

ly.i.motifs <- unique(c(enr.motifs.FL[['Ly-I']], enr.motifs.BM[['Ly-I']]))

ly.i.df <- data.frame('FL' = enr.res.FL[['Ly-I']][ly.i.motifs, 'fold.enrichment'], 
                     'BM' = enr.res.BM[['Ly-I']][ly.i.motifs, 'fold.enrichment'], 
                     row.names = ly.i.motifs)

mtx <- as.matrix(ly.i.df)
Heatmap(mtx)

hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = T, heatmap_legend_param = list(title = 'Fold enrichment'), column_names_side = 'top', column_names_rot = 0,
              column_names_centered = T, row_names_gp = gpar(fontsize = 10), row_dend_reorder = F)


pdf(paste0(fig_path, 'additional_plots/da_peak_motifs_LyI.pdf'), width = 4)
hm
dev.off()

ly.ii.motifs <- unique(c(enr.motifs.FL[['Ly-II']], enr.motifs.BM[['Ly-II']]))

ly.ii.df <- data.frame('FL' = enr.res.FL[['Ly-II']][ly.ii.motifs, 'fold.enrichment'], 
                      'BM' = enr.res.BM[['Ly-II']][ly.ii.motifs, 'fold.enrichment'], 
                      row.names = ly.ii.motifs)

mtx <- as.matrix(ly.ii.df)
Heatmap(mtx)

hm <- Heatmap(mtx, cluster_columns = F, show_row_names = T, show_column_names = T, heatmap_legend_param = list(title = 'Fold enrichment'), column_names_side = 'top', column_names_rot = 0,
              column_names_centered = T, row_names_gp = gpar(fontsize = 10), row_dend_reorder = F)


pdf(paste0(fig_path, 'additional_plots/da_peak_motifs_LyII.pdf'), width = 4)
hm
dev.off()

# make data frame with number of DA peaks per cluster
da_FL <- lapply(da_pos_FL, length) %>% unlist()
da_BM <- lapply(da_neg_FL, length) %>% unlist()

df_da <- rbind(data.frame('cluster' = names(da_FL), 'count' = unname(da_FL), 'up_or_dn' = 'up'), 
               data.frame('cluster' = names(da_BM), 'count' = unname(da_BM), 'up_or_dn' = 'down'))

df_da$cluster <- factor(df_da$cluster, levels = clust_order)
df_da$up_or_dn <- factor(df_da$up_or_dn, levels = c('up', 'down'))

ggplot(df_da, aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.7, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('Number of peaks up- or downregulated in FL compared to BM')
ggsave(paste0(fig_path, 'additional_plots/FLvsBM_da_peaks.pdf'), width = 9)

ggplot(subset(df_da, cluster %in% c('HSC')), aes(x = cluster, y = count, fill = up_or_dn)) + 
  geom_bar(stat = 'identity', width = 0.3, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))
ggsave(paste0(fig_path, 'additional_plots/FLvsBM_da_peaks_subs.pdf'), width = 3, height = 2.5)

write.table(hsc.motifs, paste0(csv_path, 'FLvsBM_enr_motifs.txt'), row.names = F, col.names = F, quote = F)

write.csv(enr.res.FL[['HSC']][hsc.motifs, ], paste0(sup_table_path, '/FL_HSC_motif_enrichment.csv'), row.names = T)
write.csv(enr.res.BM[['HSC']][hsc.motifs, ], paste0(sup_table_path, '/BM_HSC_motif_enrichment.csv'), row.names = T)
