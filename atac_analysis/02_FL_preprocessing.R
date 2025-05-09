source('scripts/atac_analysis/00.R')

library(harmony)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

# run a first round of harmony to get clusters for MACS2 peak calling
FL <-
  RunHarmony(
    FL,
    group.by.vars = "sample_number",
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F 
  )

FL <- FindNeighbors(object = FL, reduction = 'harmony', dims = 2:20)
FL <- RunUMAP(FL, reduction = "harmony", dims = 2:20, return.model = TRUE, min.dist = 0.3, 
              spread = 1, n.neighbors = 30)

FL <- FindClusters(object = FL, verbose = FALSE, algorithm = 3, resolution = 0.7)
DimPlot(FL, label = TRUE, pt.size = 0.85) + coord_fixed()

table(FL$seurat_clusters)

DimPlot(FL, label = TRUE, pt.size = 0.7, group.by = 'sample_number', label.size = 0) + coord_fixed() | 
  DimPlot(FL, label = TRUE, pt.size = 0.7, group.by = 'seurat_clusters') + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/FL_firstharmony.pdf'), width = 7, height = 4)

table(FL$seurat_clusters)
Idents(FL) <- 'seurat_clusters'

clust_peaks <- CallPeaks(
  object = FL,
  macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2',
  group.by = 'seurat_clusters',
  idents = as.character(0:13)
)

clust_peaks

# remove X and Y chromosomes since there is likely one male and one female sample
clust_peaks <- clust_peaks[! clust_peaks@seqnames %in% c('chrY', 'chrX')]

# check peak widths
peakwidths <- width(clust_peaks)
clust_peaks <- clust_peaks[peakwidths < 10000 & peakwidths > 20]


mtx <- FeatureMatrix(fragments = Fragments(FL), features = clust_peaks, cells = Cells(FL))

FL[['peaks']] <- CreateChromatinAssay(counts = mtx, annotation = Annotation(FL), fragments = Fragments(FL))

saveRDS(FL, paste0(rds_path, 'FL_only.rds'))
