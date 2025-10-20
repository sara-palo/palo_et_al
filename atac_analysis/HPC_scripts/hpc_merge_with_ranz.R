
set.seed(12345)

packages <- c('Signac', 'Seurat', 'stringr', 'ggplot2', 'GenomicRanges')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/atac/rds_files/'
fig_path <- '/home/spalo/atac/figures/'

CS22 <-  readRDS(paste0(rds_path, 'FL_only.rds'))
ABM <-  readRDS(paste0(rds_path, 'BM_only.rds'))
trim2.1 <- readRDS(paste0(rds_path, 'ranzoni.rds'))

DefaultAssay(CS22) <- 'peaks'
DefaultAssay(trim2.1) <- 'peaks'
DefaultAssay(ABM) <- 'peaks'

CS22$clust.names <- as.character(CS22$pred.id.BM)
trim2.1$clust.names <- as.character(trim2.1$pred.id.BM)
ABM$clust.names <- as.character(ABM$BM_annot)

UMAP_CS22 <- CS22[['BMumap']]@cell.embeddings 
colnames(UMAP_CS22) <- c('UMAP_1', 'UMAP_2')

UMAP_trim2.1 <- trim2.1[['BMumap']]@cell.embeddings 
colnames(UMAP_trim2.1) <- c('UMAP_1', 'UMAP_2')

UMAP_ABM <- ABM[['umap']]@cell.embeddings
colnames(UMAP_ABM) <- c('UMAP_1', 'UMAP_2')

mergedUMAP <- rbind(UMAP_CS22, UMAP_trim2.1, UMAP_ABM) %>% as.matrix()

CS22$tissue <- 'CS22_FL'
ABM$tissue <- 'ABM'

trim2.1$dataset <- 'ranzoni'
CS22$dataset <- 'palo'
ABM$dataset <- 'sommarin'

merged_metadata <- rbind(CS22@meta.data[, c('tissue', 'sample_number', 'seurat_clusters', 'clust.names', 'dataset')],
                         trim2.1@meta.data[, c('tissue', 'sample_number', 'seurat_clusters', 'clust.names', 'dataset')],
                         ABM@meta.data[, c('tissue', 'sample_number', 'seurat_clusters', 'clust.names', 'dataset')])


gr_merge <- reduce(c(granges(CS22[['peaks']]), granges(trim2.1[['peaks']]), granges(ABM[['peaks']])))
gr_merge

# Filter out possible too short/long peaks
peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]
gr_merge

# update fragment object due to change in file path
Fragments(CS22) <- NULL
Fragments(CS22) <- CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = Cells(CS22))

Fragments(ABM) <- NULL
Fragments(ABM) <- CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = Cells(ABM))

Fragments(trim2.1) <- NULL
Fragments(trim2.1) <- CreateFragmentObject('/home/spalo/atac/ranzoni_frags/frags2.tsv.gz', cells = Cells(trim2.1))

merged.counts <- FeatureMatrix(
  fragments = c(Fragments(CS22), Fragments(trim2.1), Fragments(ABM)),
  features = gr_merge, 
  cells = c(Cells(CS22), Cells(trim2.1), Cells(ABM))
)

peaks <- CreateChromatinAssay(
  counts = merged.counts,
  fragments = c(Fragments(trim2.1), Fragments(CS22), Fragments(ABM)),
  min.cells = 10, 
  genome = 'hg38',
)

merged <- CreateSeuratObject(
  counts = peaks, 
  assay = 'peaks'
)

merged

Annotation(merged) <- Annotation(CS22)

merged <- AddMetaData(merged, metadata = merged_metadata)

merged[['mergedUMAP']] <-
  CreateDimReducObject(embeddings = mergedUMAP,
                       key = "mergedUMAP_",
                       assay = DefaultAssay(merged))

merged$clust.names <- factor(merged$clust.names)
Idents(merged) <- 'clust.names'

merged <- RunTFIDF(merged)

saveRDS(merged, paste0(rds_path, 'merged_with_ranz.rds'))

DimPlot(merged, reduction = 'mergedUMAP', pt.size = 0.6, group.by = 'tissue', label = F, shuffle = T) + coord_fixed()
ggsave(paste0(fig_path, 'merge_ranz_BMumap.pdf'), width = 8, height = 5)


