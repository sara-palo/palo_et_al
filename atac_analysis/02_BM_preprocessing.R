source('scripts/atac_analysis/00.R')

library(harmony)

BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

# run a first round of harmony to get clusters for MACS2 peak calling
BM <-
  RunHarmony(
    BM,
    group.by.vars = "sample_number",
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F 
  )

BM <- FindNeighbors(object = BM, reduction = 'harmony', dims = 2:20)
BM <- RunUMAP(BM, reduction = "harmony", dims = 2:20, return.model = TRUE, min.dist = 0.3, 
              spread = 1, n.neighbors = 30)
BM <- FindClusters(object = BM, verbose = FALSE, algorithm = 3, resolution = 0.7)
DimPlot(BM, label = TRUE, pt.size = 0.85) + coord_fixed()

# smallest cluster is 45 cells - keep for peak calling
table(BM$seurat_clusters)

DimPlot(BM, label = TRUE, pt.size = 0.7, group.by = 'sample_number', label.size = 0) + coord_fixed() | 
  DimPlot(BM, label = TRUE, pt.size = 0.7, group.by = 'seurat_clusters') + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_firstharmony.pdf'))

Idents(BM) <- 'seurat_clusters'

clust_peaks <- CallPeaks(
  object = BM,
  macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2',
  group.by = 'seurat_clusters',
  idents = as.character(0:12)
)

clust_peaks

# remove X and Y chromosomes since there is likely one male and one female sample
clust_peaks <- clust_peaks[! clust_peaks@seqnames %in% c('chrY', 'chrX')]

# check peak widths
peakwidths <- width(clust_peaks)
clust_peaks <- clust_peaks[peakwidths < 10000 & peakwidths > 20]


mtx <- FeatureMatrix(fragments = Fragments(BM), features = clust_peaks, cells = Cells(BM))

BM[['peaks']] <- CreateChromatinAssay(counts = mtx, annotation = Annotation(BM), fragments = Fragments(BM))

saveRDS(BM, paste0(rds_path, 'BM_only.rds'))
