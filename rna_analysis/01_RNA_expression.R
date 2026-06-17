source('scripts/atac_analysis/00.R')

library(UCell)

FL <- readRDS(paste0(rds_path, 'FL_CITE.rds'))
BM <- readRDS(paste0(rds_path, 'BM_CITE.rds'))

FL$tissue <- 'FL'
BM$tissue <- 'BM'

BM <- RenameIdents(BM, c('HSC-I' = 'HSC', 'MEP-I' = 'MEP', 'MEP-II' = 'MEP'))
BM$clust.names <- Idents(BM)

merge <- merge(FL, BM)
merge <- NormalizeData(merge, 'RNA')

merge$plotclusters <- paste0(merge$tissue, '_', merge$clust.names)
Idents(merge) <- 'plotclusters'

subs <- subset(merge, idents = c('FL_Ly-I', 'FL_Ly-II', 'BM_Ly-I', 'BM_Ly-II'))
levels(subs) <- c('FL_Ly-I', 'FL_Ly-II', 'BM_Ly-I', 'BM_Ly-II')

goi <- rev(c('TCF3', 'EBF1', 'PAX5', 'MYC'))

DotPlot(subs, features = goi) + coord_flip() + RotatedAxis()
ggsave(paste0(fig_path, 'fig4/CITEseq_dotplot.pdf'), width = 5, height = 4)

goi <- rev(c('TCF3', 'EBF1', 'PAX5', 'MYC', 'LIN28B', 'RUNX1'))

DotPlot(subs, features = goi) + coord_flip() + RotatedAxis()
ggsave(paste0(fig_path, 'additional_plots/CITEseq_dotplot_incLIN28B.pdf'), width = 5, height = 4)

goi <- rev(c('JCHAIN', 'DNTT', 'RAG1', 'RAG2'))
DotPlot(subs, features = goi) + coord_flip() + RotatedAxis()
ggsave(paste0(fig_path, 'additional_plots/CITEseq_dotplot_Bgenes.pdf'), width = 5, height = 4)


subs3 <- subset(merge, idents = c('FL_HSC', 'FL_MPP-I', 'FL_MPP-II', 'BM_HSC', 'BM_MPP-I', 'BM_MPP-II'))
levels(subs3) <- c('FL_HSC', 'FL_MPP-I', 'FL_MPP-II', 'BM_HSC', 'BM_MPP-I', 'BM_MPP-II')

goi <- rev(c('GATA1', 'SPI1', 'RUNX1', 'IRF8'))

DotPlot(subs3, features = goi) + coord_flip() + RotatedAxis()
ggsave(paste0(fig_path, 'additional_plots/CITEseq_dotplot_GATA_etc.pdf'), width = 5, height = 4)


subs4 <- subset(merge, idents = c('FL_HSC', 'FL_MPP-I', 'FL_MEP', 'BM_HSC', 'BM_MPP-II', 'BM_MEP'))
levels(subs4) <- c('FL_HSC', 'FL_MPP-I', 'FL_MEP', 'BM_HSC', 'BM_MPP-II', 'BM_MEP')

goi <- rev(c('GATA1', 'GATA2', 'GATA3', 'GATA5'))

DotPlot(subs4, features = goi) + coord_flip() + RotatedAxis()
ggsave(paste0(fig_path, 'additional_plots/CITEseq_dotplot_GATAs.pdf'), width = 5, height = 4)


aksoz_data <- read.csv(paste0(ext_data_path, 'aksoz_s3.csv'), na.strings = "")

aksoz_data <- list('Plt_old' = na.omit(aksoz_data$Platelet_bias_hHSC_oldDonor), 
                   'ML_old' = na.omit(aksoz_data$Multilineage_output_hHSC_oldDonor), 
                   'Plt_young' = na.omit(aksoz_data$Platelet_bias_hHSC_youngDonor),
                   'ML_young' = na.omit(aksoz_data$Multilineage_output_hHSC_youngDonor),
                   'Plt_ind' = na.omit(aksoz_data$Platelet_bias_hHSC_ageIndependent),
                   'ML_ind' = na.omit(aksoz_data$Multilineage_output_hHSC_ageIndependent))

merge <- AddModuleScore_UCell(merge, features = aksoz_data, name = '_score', assay = 'RNA', maxRank = 2200)

VlnPlot(merge, idents = c('FL_HSC', 'BM_HSC'), features = paste0(c('Plt_young', 'Plt_old', 'Plt_ind', 'ML_young', 'ML_old', 'ML_ind'), '_score'), 
        cols = unname(annot_colors$tissue_cols))
ggsave(paste0(fig_path, 'additional_plots/Plt_ML_scores_all.pdf'), width = 7.4, height = 6)

VlnPlot(merge, idents = c('FL_HSC', 'BM_HSC'), features = paste0(c('Plt_ind', 'ML_ind'), '_score'), cols = unname(annot_colors$tissue_cols))
ggsave(paste0(fig_path, 'fig1/Plt_ML_scores.pdf'), width = 6.5, height = 4.1)

VlnPlot(merge, idents = c('FL_HSC', 'BM_HSC'), features = 'Plt_ind_score', cols = unname(annot_colors$tissue_cols)) + ylim(0,1.1) + NoLegend() |
  VlnPlot(merge, idents = c('FL_HSC', 'BM_HSC'), features = 'ML_ind_score', cols = unname(annot_colors$tissue_cols)) + ylim(0,1.1) + NoLegend()
ggsave(paste0(fig_path, 'fig1/Plt_ML_scores_ylim.pdf'), width = 5, height = 4.1)



