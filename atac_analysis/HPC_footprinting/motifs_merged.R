set.seed(12345)

packages <- c('Signac', 'Seurat', 'BSgenome.Hsapiens.UCSC.hg38', 'ggplot2', 'stringr')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/atac/rds_files/'
fig_path <- '/home/spalo/atac/figures/'

merged <- readRDS(paste0(rds_path, 'merged.rds'))

DefaultAssay(merged) <- 'peaks'
merged$plotclusters <- paste0(merged$sample, '_', merged$clust.names)
Idents(merged) <- 'plotclusters'

# update fragment object due to change in file path
Fragments(merged) <- NULL
Fragments(merged) <- CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = Cells(merged))

motifs_oi_1 <- c('RUNX1', 'IRF8', 'SPIB')
motifs_oi_2 <- c('EBF1', 'HOXB7', 'IRF8', 'PAX5')
motifs_oi_3 <- c('RUNX1', 'IRF8', 'FOS', 'JUN', 'SPIB', 'FOSL2')
motifs_oi_4 <- c('NFIB', 'NFE2')
motifs_oi_5 <- 'HOXB7'

motifs_oi <- unique(c(motifs_oi_1, motifs_oi_2, motifs_oi_3, motifs_oi_4, motifs_oi_5))

for (motif in motifs_oi){
  merged <- Footprint(object = merged, genome = BSgenome.Hsapiens.UCSC.hg38, motif.name = motif, compute.expected = T, in.peaks = T, assay = 'peaks')
 }

for (motif in motifs_oi_1){
  p <- PlotFootprint(merged, motif, idents = c('FL_HSC', 'BM_HSC', 'FL_Ly-II', 'BM_Ly-II', 'FL_GMP-II', 'BM_GMP-II'), group.by = 'plotclusters')
  ggsave(plot = p, filename = paste0(fig_path, 'merged_', motif, '_footprint_lymy.pdf'), width = 8, height = 6)   
}

for (motif in motifs_oi_2){
  p <- PlotFootprint(merged, motif, idents = c('FL_HSC', 'BM_HSC', 'FL_Ly-II', 'BM_Ly-II'), group.by = 'plotclusters')
  ggsave(plot = p, filename = paste0(fig_path, 'merged_', motif, '_footprint_ly.pdf'), width = 8, height = 6)   
}

for (motif in motifs_oi_3){
  p <- PlotFootprint(merged, motif, idents = c('FL_HSC', 'BM_HSC', 'FL_GMP-II', 'BM_GMP-II'), group.by = 'plotclusters')
  ggsave(plot = p, filename = paste0(fig_path, 'merged_', motif, '_footprint_gmp.pdf'), width = 8, height = 6)   
}

for (motif in motifs_oi_4){
  p <- PlotFootprint(merged, motif, idents = c('FL_HSC', 'BM_HSC', 'FL_Ery-II', 'BM_Ery-II'), group.by = 'plotclusters')
  ggsave(plot = p, filename = paste0(fig_path, 'merged_', motif, '_footprint_ery.pdf'), width = 8, height = 6)   
}

for (motif in motifs_oi_4){
  p <- PlotFootprint(merged, motif, idents = c('FL_HSC', 'BM_HSC', 'FL_MEP-II', 'BM_MEP-II'), group.by = 'plotclusters')
  ggsave(plot = p, filename = paste0(fig_path, 'merged_', motif, '_footprint_mep.pdf'), width = 8, height = 6)   
}

for (motif in motifs_oi_5){
  p <- PlotFootprint(merged, motif, idents = c('FL_HSC', 'BM_HSC'), group.by = 'plotclusters')
  ggsave(plot = p, filename = paste0(fig_path, 'merged_', motif, '_footprint_hsc.pdf'), width = 8, height = 6)   
}

saveRDS(merged[['peaks']], paste0(rds_path, 'merged_peaks_withfootprint.rds'))

