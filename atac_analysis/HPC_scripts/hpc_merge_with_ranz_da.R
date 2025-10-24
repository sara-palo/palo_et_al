
set.seed(12345)

packages <- c('Signac', 'Seurat', 'stringr', 'ggplot2', 'GenomicRanges')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/atac/rds_files/'
fig_path <- '/home/spalo/atac/figures/'

merged <- readRDS(paste0(rds_path, 'merged_with_ranz.rds'))

Idents(merged) <- 'clust.names'
subs <- subset(merged, idents = 'HSC')

Idents(subs) <- 'tissue'

comb_df <- data.frame()

for(ident in c('CS22_FL', 'trim2_FL', 'trim2_BM')){

da_peaks <- FindMarkers(
  object = subs,
  ident.1 = ident,
  ident.2 = 'ABM',
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  assay = 'peaks'
 )
 
 da_peaks$tissue <- ident
 da_peaks$peak <- rownames(da_peaks)
 comb_df <- rbind(comb_df, da_peaks)

}

write.csv(comb_df, paste0(fig_path, 'merged_ranz_da_peaks_HSCs_vs_adult.csv'))


