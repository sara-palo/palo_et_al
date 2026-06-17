source('scripts/atac_analysis/00.R')

merged <- readRDS(paste0(rds_path, 'merged.rds'))
Idents(merged) <- 'sample'

ly_subs <- subset(merged, clust.names == 'Ly-II')
hsc_subs <- subset(merged, clust.names == 'HSC')

CoveragePlot(ly_subs, 'LIN28B', extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/covplot_lyII_lin28b.pdf'), width = 7.2, height = 4)

CoveragePlot(ly_subs, 'RAG2', extend.upstream = 70000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'sup4/covplot_lyII_rag.pdf'), width = 7.2, height = 4)

CoveragePlot(ly_subs, 'DNTT', extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'sup4/covplot_lyII_dntt.pdf'), width = 7.2, height = 4)

CoveragePlot(merged, 'DNTT', extend.upstream = 5000, extend.downstream = 5000)
ggsave(paste0(fig_path, 'additional_plots/covplot_merged_dntt.pdf'), width = 8, height = 6)



CoveragePlot(hsc_subs, 'LIN28B', extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/covplot_hsc_lin28b.pdf'), width = 7.2, height = 4)

CoveragePlot(hsc_subs, 'SPI1', extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/covplot_hsc_spi1.pdf'), width = 7.2, height = 4)


CoveragePlot(hsc_subs, 'GATA1', extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/covplot_hsc_gata1.pdf'), width = 7.2, height = 4)

CoveragePlot(hsc_subs, 'GATA2', extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/covplot_hsc_gata2.pdf'), width = 7.2, height = 4)

CoveragePlot(hsc_subs, 'GATA3', extend.upstream = 5000, extend.downstream = 5000) & scale_fill_manual(values = annot_colors$tissue_cols)
ggsave(paste0(fig_path, 'additional_plots/covplot_hsc_gata3.pdf'), width = 7.2, height = 4)


merged <- readRDS(paste0(rds_path, 'merged_with_mard.rds'))
Idents(merged) <- 'tissue'

ly_subs <- subset(merged, clust.names == 'Ly-II')

CoveragePlot(ly_subs, 'RAG2', extend.upstream = 70000, extend.downstream = 5000)
ggsave(paste0(fig_path, 'additional_plots/covplot_withmard_lyII_rag.pdf'), width = 7.2, height = 5)

CoveragePlot(ly_subs, 'DNTT', extend.upstream = 5000, extend.downstream = 5000)
ggsave(paste0(fig_path, 'additional_plots/covplot_withmard_lyII_dntt.pdf'), width = 7.2, height = 5)


