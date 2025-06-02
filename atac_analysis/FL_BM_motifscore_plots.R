source('scripts/atac_analysis/00.R')

library(ggfortify)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

clusters_oi <- c('HSC', 'Ly-II', 'GMP-II')
Idents(FL) <- 'pred.id.BM'

VlnPlot(BM, 'Ly.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1, idents = clusters_oi) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'sup5/BM_LyIIscore_vln.pdf'), width = 2.5, height = 3.1)

VlnPlot(BM, 'GMP.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1, idents = clusters_oi) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'sup5/BM_GMPIIscore_vln.pdf'), width = 2.5, height = 3.1)

VlnPlot(BM, 'Ly.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/BM_full_LyIIscore_vln.pdf'), width = 8)

VlnPlot(BM, 'GMP.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/BM_full_GMPIIscore_vln.pdf'), width = 8)

VlnPlot(BM, 'GMP.I_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/BM_GMPIscore_vln.pdf'), width = 8)

VlnPlot(BM, 'Ery.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/BM_EryIIscore_vln.pdf'), width = 8)

VlnPlot(FL, 'BM.Ly.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1, idents = clusters_oi) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'sup5/FL_BMLyIIscore_vln.pdf'), width = 2.5, height = 3.1)

VlnPlot(FL, 'BM.GMP.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1, idents = clusters_oi) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'sup5/FL_BMGMPIIscore_vln.pdf'), width = 2.5, height = 3.1)

VlnPlot(FL, 'BM.Ly.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/FL_full_BMLyIIscore_vln.pdf'), width = 8)

VlnPlot(FL, 'BM.GMP.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/FL_full_BMGMPIIscore_vln.pdf'), width = 8)

VlnPlot(FL, 'BM.GMP.I_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/FL_BMGMPIscore_vln.pdf'), width = 8)

VlnPlot(BM, 'FL.Ly.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/BM_full_FLLyIIscore_vln.pdf'), width = 8)

VlnPlot(BM, 'FL.GMP.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/BM_full_FLGMPIIscore_vln.pdf'), width = 8)

VlnPlot(FL, 'GMP.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/FL_full_FLGMPIIscore_vln.pdf'), width = 8)

VlnPlot(FL, 'Ly.II_UCell', cols = annot_colors$cluster_cols, pt.size = 0.1) + ylim(c(-0.1,1)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/FL_full_FLLyIIscore_vln.pdf'), width = 8)


# PCA for only Ly and GMP-II
clusters_oi <- c('Ly-I', 'Ly-II', 'GMP-II')

FLsubs <- subset(FL, pred.id.BM %in% clusters_oi)
BMsubs <- subset(BM, BM_annot %in% clusters_oi)

FLvars <- apply(FLsubs[['cisbpchromvar']]@data, 1, var)
hist(FLvars, breaks = 500)
hist(FLvars, breaks = 500, xlim = c(0,3))
table(FLvars > 1.5)

BMvars <- apply(BMsubs[['cisbpchromvar']]@data, 1, var)
hist(BMvars, breaks = 500)
hist(BMvars, breaks = 500, xlim = c(0,3))
table(BMvars > 1.6)

FLvars_ord <- FLvars[order(FLvars, decreasing = T)]
BMvars_ord <- BMvars[order(BMvars, decreasing = T)]

FLsubs <- ScaleData(FLsubs, assay = 'cisbpchromvar', features = names(FLvars_ord[1:150]))
BMsubs <- ScaleData(BMsubs, assay = 'cisbpchromvar', features = names(BMvars_ord[1:150]))

FLsubsmotifs <- FLsubs[['cisbpchromvar']]@scale.data
dim(FLsubsmotifs)
FLsubsmotifs <- data.frame(t(FLsubsmotifs))
FLsubsmotifs$cluster <- FLsubs$pred.id.BM
FLsubsmotifs$cluster <- factor(FLsubsmotifs$cluster)

pcaFLsubs <- prcomp(FLsubsmotifs[, 1:(dim(FLsubsmotifs)[2]-1)])

autoplot(object = pcaFLsubs, data = FLsubsmotifs, color = 'cluster', size = 1.1) +
  scale_color_manual(values = annot_colors$cluster_cols[clusters_oi]) + 
  theme_bw()
ggsave(paste0(fig_path, 'fig4/FL_lymy_top150motifs_pca.pdf'))

loadings <- as.data.frame(pcaFLsubs$rotation)

PC1 <- loadings$PC1
names(PC1) <- rownames(loadings)
PC1 <- PC1[order(PC1, decreasing = T)]

PC1_pos <- head(PC1, 25)
PC1_neg <- tail(PC1, 25)

PC2 <- loadings$PC2
names(PC2) <- rownames(loadings)
PC2 <- PC2[order(PC2, decreasing = T)]

PC2_pos <- head(PC2, 25)
PC2_neg <- tail(PC2, 25)

write.csv(data.frame('PC1_pos' = names(PC1_pos), 'PC1_pos_vals' = unname(PC1_pos), 'PC1_neg' = names(PC1_neg), 'PC1_neg_vals' = unname(PC1_neg),
                     'PC2_pos' = names(PC2_pos), 'PC2_pos_vals' = unname(PC2_pos), 'PC2_neg' = names(PC2_neg), 'PC2_neg_vals' = unname(PC2_neg)), 
          paste0(sup_table_path, '/FLmotifs_PC_components.csv'))

BMsubsmotifs <- BMsubs[['cisbpchromvar']]@scale.data
dim(BMsubsmotifs)
BMsubsmotifs <- data.frame(t(BMsubsmotifs))
BMsubsmotifs$cluster <- BMsubs$BM_annot
BMsubsmotifs$cluster <- factor(BMsubsmotifs$cluster)

pcaBMsubs <- prcomp(BMsubsmotifs[, 1:(dim(BMsubsmotifs)[2]-1)])

trng <- rep('triangle', length(clusters_oi))
names(trng) <- clusters_oi

autoplot(object = pcaBMsubs, data = BMsubsmotifs, color = 'cluster', shape = 'cluster', size = 1.1) +
  scale_color_manual(values = annot_colors$cluster_cols[clusters_oi]) +
  scale_shape_manual(values = trng) +
  theme_bw()
ggsave(paste0(fig_path, 'fig4/BM_lymy_top150motifs_pca.pdf'))

loadings <- as.data.frame(pcaBMsubs$rotation)

PC1 <- loadings$PC1
names(PC1) <- rownames(loadings)
PC1 <- PC1[order(PC1, decreasing = T)]

PC1_pos <- head(PC1, 25)
PC1_neg <- tail(PC1, 25)

PC2 <- loadings$PC2
names(PC2) <- rownames(loadings)
PC2 <- PC2[order(PC2, decreasing = T)]

PC2_pos <- head(PC2, 25)
PC2_neg <- tail(PC2, 25)

write.csv(data.frame('PC1_pos' = names(PC1_pos), 'PC1_pos_vals' = unname(PC1_pos), 'PC1_neg' = names(PC1_neg), 'PC1_neg_vals' = unname(PC1_neg),
                     'PC2_pos' = names(PC2_pos), 'PC2_pos_vals' = unname(PC2_pos), 'PC2_neg' = names(PC2_neg), 'PC2_neg_vals' = unname(PC2_neg)), 
          paste0(sup_table_path, '/BMmotifs_PC_components.csv'))

# PCA for selected clusters
clusters_oi <- c('HSC', 'Ly-II', 'GMP-II', 'MEP-II', 'Ery-II')

FLsubs <- subset(FL, pred.id.BM %in% clusters_oi)
BMsubs <- subset(BM, BM_annot %in% clusters_oi)

FLvars <- apply(FLsubs[['cisbpchromvar']]@data, 1, var)
hist(FLvars, breaks = 500)
hist(FLvars, breaks = 1000, xlim = c(0,3))
table(FLvars > 1.5)

BMvars <- apply(BMsubs[['cisbpchromvar']]@data, 1, var)
hist(BMvars, breaks = 500)
hist(BMvars, breaks = 1000, xlim = c(0,3))
table(BMvars > 1.6)

FLvars_ord <- FLvars[order(FLvars, decreasing = T)]
BMvars_ord <- BMvars[order(BMvars, decreasing = T)]

FLsubs <- ScaleData(FLsubs, assay = 'cisbpchromvar', features = names(FLvars_ord[1:150]))
BMsubs <- ScaleData(BMsubs, assay = 'cisbpchromvar', features = names(BMvars_ord[1:150]))

FLsubsmotifs <- FLsubs[['cisbpchromvar']]@scale.data
dim(FLsubsmotifs)
FLsubsmotifs <- data.frame(t(FLsubsmotifs))
FLsubsmotifs$cluster <- FLsubs$pred.id.BM
FLsubsmotifs$cluster <- factor(FLsubsmotifs$cluster)

pcaFLsubs <- prcomp(FLsubsmotifs[, 1:(dim(FLsubsmotifs)[2]-1)])

autoplot(object = pcaFLsubs, data = FLsubsmotifs, color = 'cluster', size = 1.1) +
  scale_color_manual(values = annot_colors$cluster_cols[clusters_oi]) + 
  theme_bw()
ggsave(paste0(fig_path, 'fig2/FL_top150motifs_pca.pdf'))

BMsubsmotifs <- BMsubs[['cisbpchromvar']]@scale.data
dim(BMsubsmotifs)
BMsubsmotifs <- data.frame(t(BMsubsmotifs))
BMsubsmotifs$cluster <- BMsubs$BM_annot
BMsubsmotifs$cluster <- factor(BMsubsmotifs$cluster)

pcaBMsubs <- prcomp(BMsubsmotifs[, 1:(dim(BMsubsmotifs)[2]-1)])

trng <- rep('triangle', length(clusters_oi))
names(trng) <- clusters_oi

autoplot(object = pcaBMsubs, data = BMsubsmotifs, color = 'cluster', shape = 'cluster', size = 1.1) +
  scale_color_manual(values = annot_colors$cluster_cols[clusters_oi]) +
  scale_shape_manual(values = trng) +
  theme_bw()
ggsave(paste0(fig_path, 'fig2/BM_top150motifs_pca.pdf'))


# PCA on all cells
FLvars <- apply(FL[['cisbpchromvar']]@data, 1, var)
hist(FLvars, breaks = 500)
hist(FLvars, breaks = 1000, xlim = c(0,3))
table(FLvars > 1.4)

BMvars <- apply(BM[['cisbpchromvar']]@data, 1, var)
hist(BMvars, breaks = 500)
hist(BMvars, breaks = 1000, xlim = c(0,3))
table(BMvars > 1.6)

FLvars_ord <- FLvars[order(FLvars, decreasing = T)]
BMvars_ord <- BMvars[order(BMvars, decreasing = T)]

FL <- ScaleData(FL, assay = 'cisbpchromvar', features = names(FLvars_ord[1:150]))
BM <- ScaleData(BM, assay = 'cisbpchromvar', features = names(BMvars_ord[1:150]))

FLmotifs <- FL[['cisbpchromvar']]@scale.data
dim(FLmotifs)
FLmotifs <- data.frame(t(FLmotifs))
FLmotifs$cluster <- FL$pred.id.BM
FLmotifs$cluster <- factor(FLmotifs$cluster)

pcaFL <- prcomp(FLmotifs[, 1:(dim(FLmotifs)[2]-1)])

autoplot(object = pcaFL, data = FLmotifs, color = 'cluster', size = 1.1) +
  scale_color_manual(values = annot_colors$cluster_cols) + 
  theme_bw()
ggsave(paste0(fig_path, 'additional_plots/FL_top150motifs_pca.pdf'))

BMmotifs <- BM[['cisbpchromvar']]@scale.data
dim(BMmotifs)
BMmotifs <- data.frame(t(BMmotifs))
BMmotifs$cluster <- BM$BM_annot
BMmotifs$cluster <- factor(BMmotifs$cluster)

pcaBM <- prcomp(BMmotifs[, 1:(dim(BMmotifs)[2]-1)])

trng <- rep('triangle', length(clust_order))
names(trng) <- clust_order

autoplot(object = pcaBM, data = BMmotifs, color = 'cluster', shape = 'cluster', size = 1.1) +
  scale_color_manual(values = annot_colors$cluster_cols) +
  scale_shape_manual(values = trng) +
  theme_bw()
ggsave(paste0(fig_path, 'additional_plots/BM_top150motifs_pca.pdf'))


# scatter plots
FL <- ScaleData(FL, assay = 'cisbpchromvar')
BM <- ScaleData(BM, assay = 'cisbpchromvar')

fl_chrmv <- data.frame(t(FL[['cisbpchromvar']]@scale.data[, WhichCells(FL, expression = pred.id.BM == 'HSC')]))
fl_chrmv$sample <- 'FL'
bm_chrmv <- data.frame(t(BM[['cisbpchromvar']]@scale.data[, WhichCells(BM, expression = BM_annot == 'HSC')]))
bm_chrmv$sample <- 'BM'

ggplot(fl_chrmv, aes(y = GATA1, x = SPI1, color = sample)) +
  scale_color_manual(values = annot_colors$tissue_cols) +
  ylim(min(c(fl_chrmv$GATA1, bm_chrmv$GATA1)), max(c(fl_chrmv$GATA1, bm_chrmv$GATA1))) +
  xlim(min(c(fl_chrmv$SPI1, bm_chrmv$SPI1)), max(c(fl_chrmv$SPI1, bm_chrmv$SPI1))) +
  geom_point() + theme_bw() + coord_fixed()
ggsave(paste0(fig_path, 'fig2/FL_GATA1_SPI1_scatter.pdf'))

ggplot(bm_chrmv, aes(y = GATA1, x = SPI1, color = sample)) +
  scale_color_manual(values = annot_colors$tissue_cols) +
  ylim(min(c(fl_chrmv$GATA1, bm_chrmv$GATA1)), max(c(fl_chrmv$GATA1, bm_chrmv$GATA1))) +
  xlim(min(c(fl_chrmv$SPI1, bm_chrmv$SPI1)), max(c(fl_chrmv$SPI1, bm_chrmv$SPI1))) +
  geom_point() + theme_bw() + coord_fixed()
ggsave(paste0(fig_path, 'fig2/BM_GATA1_SPI1_scatter.pdf'))

ggplot(fl_chrmv, aes(x = EBF1, y = SPI1, color = sample)) +
  scale_color_manual(values = annot_colors$tissue_cols) +
  xlim(min(c(fl_chrmv$EBF1, bm_chrmv$EBF1)), max(c(fl_chrmv$EBF1, bm_chrmv$EBF1))) +
  ylim(min(c(fl_chrmv$SPI1, bm_chrmv$SPI1)), max(c(fl_chrmv$SPI1, bm_chrmv$SPI1))) +
  geom_point() + theme_bw() + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/FL_EBF1_SPI1_scatter.pdf'))

ggplot(bm_chrmv, aes(x = EBF1, y = SPI1, color = sample)) +
  scale_color_manual(values = annot_colors$tissue_cols) +
  xlim(min(c(fl_chrmv$EBF1, bm_chrmv$EBF1)), max(c(fl_chrmv$EBF1, bm_chrmv$EBF1))) +
  ylim(min(c(fl_chrmv$SPI1, bm_chrmv$SPI1)), max(c(fl_chrmv$SPI1, bm_chrmv$SPI1))) +
  geom_point() + theme_bw() + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_EBF1_SPI1_scatter.pdf'))


