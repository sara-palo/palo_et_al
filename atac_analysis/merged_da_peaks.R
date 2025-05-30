source('scripts/atac_analysis/00.R')

merged <- readRDS(paste0(rds_path, 'merged.rds'))

Idents(merged) <- 'clust.names'
table(merged$clust.names)

for (ident in unique(merged$clust.names)){

  da_peaks <- FindMarkers(
    object = merged,
    ident.1 = 'FL',
    ident.2 = 'BM',
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    assay = 'peaks',
    subset.ident = ident,
    group.by = 'sample'
  )

  if(ident == 'Ba/M'){
    write.csv(da_peaks, paste0(csv_path_da, '/FLvsBM/', 'BaM', '.csv'))
  } else {
    write.csv(da_peaks, paste0(csv_path_da, '/FLvsBM/', ident, '.csv'))
  }
}

# Also try FL cluster 0 vs BM HSCs
subs <- subset(merged, cells = c(WhichCells(merged, expression = sample == 'BM' & clust.names == 'HSC'), 
                                 WhichCells(merged, expression = seurat_clusters == 0 & sample == 'FL')))
Idents(subs) <- 'sample'

da_peaks <- FindMarkers(
  object = subs,
  ident.1 = 'FL',
  ident.2 = 'BM',
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  assay = 'peaks'
)

write.csv(da_peaks, paste0(csv_path_da, '/FLvsBM/', 'FLc0_vs_BMHSC', '.csv'))






