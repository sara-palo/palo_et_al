source('scripts/atac_analysis/00.R')

merged <- readRDS(paste0(rds_path, 'merged.rds'))
Idents(merged) <- 'clust.names'

subs <- subset(merged, idents = 'Ly-II')
Idents(subs) <- 'sample'

CoveragePlot(subs, 'RAG2', group.by = 'sample') & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'sup4/LyII_RAG_coverageplot.pdf'))

CoveragePlot(subs, 'DNTT', group.by = 'sample') & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'sup4/LyII_DNTT_coverageplot.pdf'))

CoveragePlot(subs, 'IL7R', group.by = 'sample') & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/LyII_IL7R_coverageplot.pdf'))



subs <- subset(merged, idents = 'Ly-I')
Idents(subs) <- 'sample'

CoveragePlot(subs, 'RAG2', group.by = 'sample') & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/LyI_RAG_coverageplot.pdf'))

CoveragePlot(subs, 'DNTT', group.by = 'sample') & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/LyI_DNTT_coverageplot.pdf'))


Idents(merged) <- 'sample'

CoveragePlot(merged, 'LIN28B', group.by = 'sample') & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/all_LIN28B_coverageplot.pdf'))




