source('scripts/atac_analysis/00.R')

library(BSgenome.Hsapiens.UCSC.hg38)

merged <- readRDS(paste0(rds_path, 'merged_with_ranz.rds'))

table(merged$tissue, merged$dataset)

Idents(merged) <- 'tissue'

da_peaks <- read.csv(paste0(csv_path_da, '/FLvsBM/merged_ranz_da_peaks_HSCs_vs_adult.csv'))
da_peaks <- subset(da_peaks, p_val_adj < 0.05)
da_peaks <- subset(da_peaks, avg_log2FC > 2 |  avg_log2FC < -2)

table(da_peaks$tissue, da_peaks$avg_log2FC > 2)

da_peaks <- cbind(da_peaks, ClosestFeature(merged, regions = da_peaks$peak))

# subset for peaks within 2kb of a gene for cellradar
da_peaks2 <- subset(da_peaks, distance < 2000)

res1 <- prep_data('human', unique(da_peaks2[da_peaks2$tissue == 'CS22_FL' & da_peaks2$avg_log2FC > 2, 'gene_name']))
res2 <- prep_data('human', unique(da_peaks2[da_peaks2$tissue == 'trim2_FL' & da_peaks2$avg_log2FC > 2, 'gene_name']))
res3 <- prep_data('human', unique(da_peaks2[da_peaks2$tissue == 'trim2_BM' & da_peaks2$avg_log2FC > 2, 'gene_name']))

res4 <- prep_data('human', unique(da_peaks2[da_peaks2$tissue == 'CS22_FL' & da_peaks2$avg_log2FC < -2, 'gene_name']))
res5 <- prep_data('human', unique(da_peaks2[da_peaks2$tissue == 'trim2_FL' & da_peaks2$avg_log2FC < -2, 'gene_name']))
res6 <- prep_data('human', unique(da_peaks2[da_peaks2$tissue == 'trim2_BM' & da_peaks2$avg_log2FC < -2, 'gene_name']))

orca(plot_cellradar(res_list = list(res1, res4), colors = c(annot_colors$tissue_cols[['FL']], annot_colors$tissue_cols[['BM']]), names = c('CS22 FL', 'Adult BM'), 
                         dash = c('solid', '4px')), file = paste0(fig_path, 'fig3/da_peaks_genes_CS22_HSC_cr.pdf'), more_args = '--disable-gpu')

orca(plot_cellradar(res_list = list(res2, res5), colors = c('#5949EB', annot_colors$tissue_cols[['BM']]), names = c('2nd trim FL', 'Adult BM'),
                    dash = c('solid', '4px')), file = paste0(fig_path, 'sup3/da_peaks_genes_trim2_FL_HSC_cr.pdf'), more_args = '--disable-gpu')

orca(plot_cellradar(res_list = list(res3, res6), colors = c('#E8AB09', annot_colors$tissue_cols[['BM']]), names = c('2nd trim BM', 'Adult BM'),
                    dash = c('solid', '4px')), file = paste0(fig_path, 'sup3/da_peaks_genes_trim2_BM_HSC_cr.pdf'), more_args = '--disable-gpu')

missinggenes <- data.frame('CS22_up' = str_flatten(res1$nonexisting_genes, collapse = ','), 
                           'trim2_FL_up' = str_flatten(res2$nonexisting_genes, collapse = ','),
                           'trim2_BM_up' = str_flatten(res3$nonexisting_genes, collapse = ','),
                           'CS22_dn' = str_flatten(res4$nonexisting_genes, collapse = ','),
                           'trim2_FL_dn' = str_flatten(res5$nonexisting_genes, collapse = ','),
                           'trim2_BM_dn' = str_flatten(res6$nonexisting_genes, collapse = ','))

write.table(missinggenes, paste0(sup_table_path, '/HSC_CR_missinggenes.txt'))

table(toupper(unique(da_peaks2[da_peaks2$tissue == 'CS22_FL' & da_peaks2$avg_log2FC > 2, 'gene_name'])) %in% res1$nonexisting_genes)
table(toupper(unique(da_peaks2[da_peaks2$tissue == 'trim2_FL' & da_peaks2$avg_log2FC > 2, 'gene_name'])) %in% res2$nonexisting_genes)
table(toupper(unique(da_peaks2[da_peaks2$tissue == 'trim2_BM' & da_peaks2$avg_log2FC > 2, 'gene_name'])) %in% res3$nonexisting_genes)

table(toupper(unique(da_peaks2[da_peaks2$tissue == 'CS22_FL' & da_peaks2$avg_log2FC < -2, 'gene_name'])) %in% res4$nonexisting_genes)
table(toupper(unique(da_peaks2[da_peaks2$tissue == 'trim2_FL' & da_peaks2$avg_log2FC < -2, 'gene_name'])) %in% res5$nonexisting_genes)
table(toupper(unique(da_peaks2[da_peaks2$tissue == 'trim2_BM' & da_peaks2$avg_log2FC < -2, 'gene_name'])) %in% res6$nonexisting_genes)


merged <- RegionStats(merged, assay = 'peaks', genome = BSgenome.Hsapiens.UCSC.hg38)
meta.feature <- GetAssayData(merged, assay = 'peaks', slot = 'meta.features')

hm_list <- list()

for (name in c('CS22_FL', 'trim2_FL', 'trim2_BM')){
  
  # find peaks open in at least 5% of the relevant HSCs for background selection
  cells <- WhichCells(merged, expression = clust.names == 'HSC' & tissue %in% c('ABM', name))
  open.peaks <- AccessiblePeaks(merged, cells = cells, min.cells = round(length(cells) * 0.05))
  
  print(length(open.peaks))
  
  regs1 <- da_peaks[da_peaks$tissue == name & da_peaks$avg_log2FC > 2, 'peak']
  
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
  
  enr.res.pos <- enriched.motifs
  enr.motifs.pos <- rownames(subset(enriched.motifs, subset = fold.enrichment > 1.5 & p.adjust < 0.05)) 
  
  regs2 <- da_peaks[da_peaks$tissue == name & da_peaks$avg_log2FC < -2, 'peak']
  
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
  
  enr.res.neg <- enriched.motifs
  enr.motifs.neg <- rownames(subset(enriched.motifs, subset = fold.enrichment > 1.5 & p.adjust < 0.05)) 
  
  plot.motifs <- unique(c(enr.motifs.pos, enr.motifs.neg))
  
  plot.df <- data.frame('fetal' = enr.res.pos[plot.motifs, 'fold.enrichment'], 
                        'BM' = enr.res.neg[plot.motifs, 'fold.enrichment'], 
                        row.names = plot.motifs)
  
  colnames(plot.df)[1] <- name
  
  mtx <- as.matrix(plot.df)
  hm_list[[name]] <- mtx
  
  write.csv(enr.res.pos, paste0(sup_table_path, '/', name, '_HSC_motif_enrichment_pos.csv'), row.names = T)
  write.csv(enr.res.pos, paste0(sup_table_path, '/', name, '_HSC_motif_enrichment_neg.csv'), row.names = T)
}

pdf(paste0(fig_path, 'fig3/da_peak_motifs_', 'CS22_FL', '_HSC.pdf'), height = nrow(hm_list[['CS22_FL']])*0.25, width = 5)
Heatmap(hm_list[['CS22_FL']], cluster_columns = F, show_row_names = T, show_column_names = T, heatmap_legend_param = list(title = 'Fold enrichment'),
        column_names_side = 'top', column_names_rot = 0,
        column_names_centered = T, row_dend_reorder = F, row_names_gp = gpar(fontsize = 12))
dev.off()


pdf(paste0(fig_path, 'sup3/da_peak_motifs_', 'trim2_FL', '_HSC.pdf'), height = nrow(hm_list[['trim2_FL']])*0.25, width = 4)
Heatmap(hm_list[['trim2_FL']], cluster_columns = F, show_row_names = T, show_column_names = T, heatmap_legend_param = list(title = 'Fold enrichment'),
        column_names_side = 'top', column_names_rot = 0,
        column_names_centered = T, row_dend_reorder = F, row_names_gp = gpar(fontsize = 12))
dev.off()

pdf(paste0(fig_path, 'sup3/da_peak_motifs_', 'trim2_BM', '_HSC.pdf'), height = nrow(hm_list[['trim2_BM']])*0.25, width = 4)
Heatmap(hm_list[['trim2_BM']], cluster_columns = F, show_row_names = T, show_column_names = T, heatmap_legend_param = list(title = 'Fold enrichment'),
        column_names_side = 'top', column_names_rot = 0, 
        column_names_centered = T, row_dend_reorder = F, row_names_gp = gpar(fontsize = 12))
dev.off()



