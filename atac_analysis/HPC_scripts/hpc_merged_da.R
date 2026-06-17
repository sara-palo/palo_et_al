set.seed(12345)

packages <- c('Signac', 'Seurat', 'stringr', 'ggplot2', 'GenomicRanges')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/lu2025-12-50/spalo/atac/rds_files/'
csv_path <- '/home/spalo/lu2025-12-50/spalo/atac/output/da_peaks/'

merged <- readRDS(paste0(rds_path, 'merged.rds'))
Idents(merged) <- 'clust.names'

clusters <- as.character(unique(Idents(merged)))
print(clusters)

for(clust in clusters){
  
  print(clust)
  
  da_peaks <- FindMarkers(
    object = merged,
    ident.1 = 'FL',
    ident.2 = 'BM',
    group.by = 'sample',
    subset.ident = clust,
    only.pos = FALSE,
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    assay = 'peaks'
  )
  
  if(clust == 'Ba/M'){
    write.csv(da_peaks, paste0(csv_path, 'merged_da_BaM.csv'))
  } else {
    write.csv(da_peaks, paste0(csv_path, 'merged_da_', clust, '.csv'))
  }
}

quit(save = 'no')


