source('scripts/atac_analysis/00.R')

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

table(duplicated(c(Cells(FL), Cells(BM)))) # no barcodes present in both objects

# First make a merged umap between the BM umap and the FL projection onto it
FL$clust.names <- as.character(FL$pred.id.BM)
BM$clust.names <- as.character(BM$BM_annot)

# create matrix of umap coordinates 
UMAP_FL <- FL[['BMumap']]@cell.embeddings 
colnames(UMAP_FL) <- c('UMAP_1', 'UMAP_2')

UMAP_BM <- BM[['umap']]@cell.embeddings
colnames(UMAP_BM) <- c('UMAP_1', 'UMAP_2')

mergedUMAP <- rbind(UMAP_FL, UMAP_BM) %>% as.matrix()

merged_metadata <- rbind(FL@meta.data[, c('sample', 'sample_number', 'seurat_clusters', 'clust.names')], 
                         BM@meta.data[, c('sample', 'sample_number', 'seurat_clusters', 'clust.names')])

# create a shared peak set
gr_merge <- reduce(c(granges(FL[['peaks']]), granges(BM[['peaks']])))
gr_merge

# Filter out possible too short/long peaks
peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]
gr_merge

# make a new peaks assay and seurat object
fragments <- CreateFragmentObject(
  path = paste0(rawdata_path, 'fragments.tsv.gz'),
  cells = c(Cells(FL), Cells(BM))
)

merged.counts <- FeatureMatrix(
  fragments = fragments,
  features = gr_merge, 
  cells = c(Cells(FL), Cells(BM))
)

peaks <- CreateChromatinAssay(
  counts = merged.counts,
  fragments = fragments,
  min.cells = 10, 
  genome = 'hg38',
)

merged <- CreateSeuratObject(
  counts = peaks, 
  assay = 'peaks'
)

merged

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(merged) <- annotations

merged <- AddMetaData(merged, metadata = merged_metadata)

merged[['mergedUMAP']] <-
  CreateDimReducObject(embeddings = mergedUMAP,
                       key = "mergedUMAP_",
                       assay = DefaultAssay(merged))

merged$clust.names <- factor(merged$clust.names, levels = clust_order)
merged$sample <- factor(merged$sample, levels = c('FL', 'BM'))
Idents(merged) <- 'clust.names'

DimPlot(merged, reduction = 'mergedUMAP', pt.size = 0.6, cols = annot_colors$cluster_cols, label = T) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/merged_umap_clust.pdf'), width = 8, height = 5)

DimPlot(merged, reduction = 'mergedUMAP', pt.size = 0.5, group.by = 'sample', shuffle = T, cols = annot_colors$tissue_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/merged_umap.pdf'), width = 8, height = 5)

merged <- RunTFIDF(merged)

VlnPlot(merged, 'nCount_peaks', group.by = 'clust.names', split.by = 'sample', cols = annot_colors$tissue_cols, pt.size = 0.01)
ggsave(paste0(fig_path, 'additional_plots/merged_peak_counts.pdf'), width = 7)

saveRDS(merged, paste0(rds_path, 'merged.rds'))

gene.activities <- GeneActivity(merged, assay = 'peaks')
saveRDS(gene.activities, paste0(rds_path, 'merged_geneact_mtx.rds'))


