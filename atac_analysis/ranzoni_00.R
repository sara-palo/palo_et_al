source('scripts/atac_analysis/00.R')

library(harmony)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
annotations <- Annotation(FL[['peaks']])
rm(FL)

ranzoni.frags <- CreateFragmentObject(path = paste0(ext_data_path, 'ranzoni/frags2.tsv.gz')) #this one has the "chr" added to the chromosomes
ranzoni.md.c <- read.csv(paste0(ext_data_path, 'ranzoni/mergedATAC-metadata.csv'), row.names = 'X')
ranzoni.umap <- read.csv(paste0(ext_data_path, 'ranzoni/mergedATAC-UMAP.csv'), row.names = 1)

# extract cell IDs that can be matched to fastq files
sampleinfo <- read.table(paste0(ext_data_path, 'ranzoni/filereport_read_run_PRJEB38346_tsv_2.txt'), header = T)
sampleinfo$samplename <- str_split_i(sampleinfo$submitted_ftp, pattern = '/', i = 6)
sampleinfo$samplename <- str_split_i(sampleinfo$samplename, pattern = "[.]", i = 1)
sampleinfo$samplename <- str_replace(sampleinfo$samplename, pattern = '_', replacement = '-')

all_na <- sapply(ranzoni.md.c, \(x) all(is.na(x)))
ranzoni.md.c <- ranzoni.md.c[, !all_na]
colnames(ranzoni.md.c)

ranzoni.md.c[1:5, c(1:5, 7, 14, 36)]
ranzoni.md.c[1:5]
ranzoni.md <- ranzoni.md.c[, c(1:5, 7, 13, 14, 36, 45)] 
ranzoni.md$paper_cluster <- ranzoni.md$seurat_clusters
rownames(ranzoni.md) <- str_remove(rownames(ranzoni.md), pattern = '.clean.bam')

ranzoni.md$fastq <- rownames(sampleinfo)[match(rownames(ranzoni.md), sampleinfo$samplename)]
rownames(ranzoni.md) <- ranzoni.md$fastq

identical(rownames(ranzoni.umap), rownames(ranzoni.md.c))
rownames(ranzoni.umap) <- rownames(ranzoni.md)

ranzoni.md$tissue <- NA

ranzoni.md$tissue <- ifelse(str_detect(ranzoni.md$SUPPLIER.SAMPLE.NAME, pattern = 'femur'), 'trim2_BM', ranzoni.md$tissue)
ranzoni.md$tissue <- ifelse(str_detect(ranzoni.md$SUPPLIER.SAMPLE.NAME, pattern = 'liver'), 
                                          'trim2_FL', ranzoni.md$tissue)

ranzoni.md$sample_number <- ranzoni.md$SAMPLE

# 1st round of peak calling
peaks <- CallPeaks(ranzoni.frags, macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2')

peaks
hist(width(peaks))

peaks <- peaks[peaks@seqnames %in% paste0('chr', c(1:22))]
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]

counts <- FeatureMatrix(
  fragments = ranzoni.frags,
  features = peaks
)

peaks <- CreateChromatinAssay(
  counts = counts,
  fragments = ranzoni.frags,
  genome = 'hg38',
  min.cells = 10
)

ranz <- CreateSeuratObject(
  counts = peaks, 
  assay = 'peaks'
)


Annotation(ranz) <- annotations

ranz <- TSSEnrichment(ranz, assay = 'peaks')

ranz$blacklist_fraction <- FractionCountsInRegion(
  object = ranz, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

# keep only cells that are included in the metadata from the paper, for information about tissue type etc
ranz <- subset(ranz, cells = ranzoni.md$fastq)

VlnPlot(ranz, c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction'))

ranz <- subset(ranz, subset = nCount_peaks > 5000 & nCount_peaks < 55000 & TSS.enrichment > 2.5)

ranz

VlnPlot(ranz, c('nCount_peaks', 'TSS.enrichment', 'blacklist_fraction'))

ranz <- AddMetaData(object = ranz, metadata = ranzoni.md)

table(ranz$paper_cluster)

ranz <- FindTopFeatures(ranz) %>% RunTFIDF() %>% RunSVD()
ElbowPlot(ranz, reduction = 'lsi')

ranz <- FindNeighbors(object = ranz, reduction = 'lsi', dims = 2:10)
ranz <- FindClusters(object = ranz, verbose = FALSE, algorithm = 3, resolution = 0.5)
ranz <- RunUMAP(ranz, reduction = 'lsi', dims = 2:10, return.model = TRUE,
                spread = 1, n.neighbors = 30)

DimPlot(ranz) | DimPlot(ranz, group.by = 'sample_number')

ranz <-
  RunHarmony(
    ranz,
    group.by.vars = 'sample_number',
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F,
    reduction.save = 'harmony'
  )

ranz <- FindNeighbors(object = ranz, reduction = 'harmony', dims = 2:10)
ranz <- FindClusters(object = ranz, verbose = FALSE, algorithm = 3, resolution = 0.5)
ranz <- RunUMAP(ranz, reduction = 'harmony', dims = 2:10, return.model = TRUE,
                spread = 1, n.neighbors = 30)

DimPlot(ranz) | DimPlot(ranz, group.by = 'sample_number') | DimPlot(ranz, group.by = 'paper_cluster')
ggsave(paste0(fig_path, 'additional_plots/ranzoni_firstpeakcall.pdf'), width = 12)

ranzoni.umap <- ranzoni.umap[Cells(ranz), ]
ranz[['paper_UMAP']] <- CreateDimReducObject(embeddings = as.matrix(ranzoni.umap), assay = 'peaks')

DimPlot(ranz, group.by = 'paper_cluster', reduction = 'paper_UMAP')

Idents(ranz) <- 'paper_cluster'
ranz <- RenameIdents(ranz, '0' = '7')
levels(ranz) <- 1:7
ranz$paper_cluster <- Idents(ranz)
DimPlot(ranz, group.by = 'paper_cluster', reduction = 'paper_UMAP')

table(ranz$seurat_clusters)

macs2peaks <- CallPeaks(
  ranz,
  macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2',
  group.by = 'seurat_clusters',
  combine.peaks = TRUE, 
  idents = 0:5
)

macs2peaks
hist(width(macs2peaks))

macs2peaks <- macs2peaks[macs2peaks@seqnames %in% paste0('chr', c(1:22))]
peakwidths <- width(macs2peaks)
macs2peaks <- macs2peaks[peakwidths < 10000 & peakwidths > 20]

counts <- FeatureMatrix(
  fragments = ranzoni.frags,
  features = macs2peaks, 
  cells = Cells(ranz)
)

peaks <- CreateChromatinAssay(
  counts = counts,
  fragments = ranzoni.frags,
  min.cells = 10, 
  genome = 'hg38',
)

ranz[['peaks']] <- peaks

Annotation(ranz) <- annotations
# ranz <- BinarizeCounts(ranz, assay = 'peaks')

ranz <- FindTopFeatures(ranz) %>% RunTFIDF() %>% RunSVD()
ElbowPlot(ranz, reduction = 'lsi')

ranz <- FindNeighbors(object = ranz, reduction = 'lsi', dims = 2:10)
ranz <- FindClusters(object = ranz, verbose = FALSE, algorithm = 3, resolution = 0.5)
ranz <- RunUMAP(ranz, reduction = 'lsi', dims = 2:10, return.model = TRUE,
                spread = 1, n.neighbors = 30)

DimPlot(ranz) | DimPlot(ranz, group.by = 'sample_number')

ranz <-
  RunHarmony(
    ranz,
    group.by.vars = 'sample_number',
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F,
    reduction.save = 'harmony'
  )

ElbowPlot(ranz, reduction = 'harmony')

ranz <- FindNeighbors(object = ranz, reduction = 'harmony', dims = 2:10)
ranz <- FindClusters(object = ranz, verbose = FALSE, algorithm = 3, resolution = 0.5)
ranz <- RunUMAP(ranz, reduction = 'harmony', dims = 2:10, return.model = TRUE,
                spread = 1, n.neighbors = 30)

DimPlot(ranz, reduction = 'umap', group.by = 'seurat_clusters') | 
  DimPlot(ranz, group.by = 'sample_number', reduction = 'umap', shuffle = T) | 
  DimPlot(ranz, group.by = 'paper_cluster', reduction = 'umap')
ggsave(paste0(fig_path, 'additional_plots/ranzoni_umap.pdf'), width = 14, height = 4.5)

VlnPlot(ranz, 'nCount_peaks', group.by = 'seurat_clusters')

VlnPlot(ranz, 'nCount_peaks', group.by = 'paper_cluster') + 
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/ranzoni_ncount_peaks_paperclust.pdf'))

VlnPlot(ranz, 'nCount_peaks', group.by = 'seurat_clusters') + 
  theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/ranzoni_ncount_peaks.pdf'))


DimPlot(ranz, reduction = 'umap', group.by = 'tissue') | DimPlot(ranz, reduction = 'paper_UMAP', group.by = 'tissue')
ggsave(paste0(fig_path, 'additional_plots/ranzoni_umap_tissue.pdf'))

saveRDS(ranz, paste0(rds_path, 'ranzoni.rds'))

VlnPlot(ranz, 'nCount_peaks', group.by = 'paper_cluster', split.by = 'tissue')
ggsave(paste0(fig_path, 'additional_plots/ranzoni_ncount_peaks_split.pdf'))
