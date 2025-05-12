source('scripts/atac_analysis/00.R')

library(harmony)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
annotations <- Annotation(FL[['peaks']])
rm(FL)

md <- read.table(file = 'external_data/marderstein/E-MTAB-13070.sdrf.txt', sep='\t', header = T, fill = T, check.names = F, quote = "")
table(md$`Characteristics[disease]`)
samples_oi <- md[md$`Characteristics[disease]` == 'normal', ]
table(samples_oi$`Characteristics[age]`)
samples_oi <- md[md$`Characteristics[age]` == '12', ]
table(samples_oi$`Characteristics[individual]`) #H5 and H6

# call peaks separately on each fragment file
frags_h5a <- CreateFragmentObject(path = paste0(ext_data_path,'marderstein/multi/H_5_A-atac_fragments.tsv.gz'))
frags_h5b <- CreateFragmentObject(path = paste0(ext_data_path,'marderstein/multi/H_5_B-atac_fragments.tsv.gz'))
frags_h6a <- CreateFragmentObject(path = paste0(ext_data_path,'marderstein/multi/H_6_A-atac_fragments.tsv.gz'))

h5a_peaks <- CallPeaks(object = frags_h5a, macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2')
h5b_peaks <- CallPeaks(object = frags_h5b, macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2')
h6a_peaks <- CallPeaks(object = frags_h6a, macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2')

# remove non-chromosomes and the X and Y chromosomes
h5a_peaks <- h5a_peaks[h5a_peaks@seqnames %in% paste0('chr', 1:22)]
h5b_peaks <- h5b_peaks[h5b_peaks@seqnames %in% paste0('chr', 1:22)]
h6a_peaks <- h6a_peaks[h6a_peaks@seqnames %in% paste0('chr', 1:22)]

# Filter out possible too short/long peaks
peakwidths <- width(h5a_peaks)
h5a_peaks <- h5a_peaks[peakwidths < 10000 & peakwidths > 20]

peakwidths <- width(h5b_peaks)
h5b_peaks <- h5b_peaks[peakwidths < 10000 & peakwidths > 20]

peakwidths <- width(h6a_peaks)
h6a_peaks <- h6a_peaks[peakwidths < 10000 & peakwidths > 20]

h5a_mtx <- FeatureMatrix(fragments = frags_h5a, features = h5a_peaks)
h5b_mtx <- FeatureMatrix(fragments = frags_h5b, features = h5b_peaks)
h6a_mtx <- FeatureMatrix(fragments = frags_h6a, features = h6a_peaks)

h5a <- CreateSeuratObject(counts = CreateChromatinAssay(counts = h5a_mtx, fragments = frags_h5a, min.cells = 10, min.features = 100), assay = 'peaks')
VlnPlot(h5a, 'nCount_peaks')
h5a <- subset(h5a, nCount_peaks > 3000 & nCount_peaks < 50000)

h5b <- CreateSeuratObject(counts = CreateChromatinAssay(counts = h5b_mtx, fragments = frags_h5b, min.cells = 10, min.features = 100), assay = 'peaks')
VlnPlot(h5b, 'nCount_peaks')
h5b <- subset(h5b, nCount_peaks > 3000 & nCount_peaks < 53000)

h6a <- CreateSeuratObject(counts = CreateChromatinAssay(counts = h6a_mtx, fragments = frags_h6a, min.cells = 10, min.features = 100), assay = 'peaks')
VlnPlot(h6a, 'nCount_peaks')
h6a <- subset(h6a, nCount_peaks > 3000 & nCount_peaks < 40000)

Annotation(h5a) <- annotations
Annotation(h5b) <- annotations
Annotation(h6a) <- annotations

h5a <- NucleosomeSignal(object = h5a) %>% TSSEnrichment(fast = TRUE)
h5b <- NucleosomeSignal(object = h5b) %>% TSSEnrichment(fast = TRUE)
h6a <- NucleosomeSignal(object = h6a) %>% TSSEnrichment(fast = TRUE)

VlnPlot(h5a, c('TSS.enrichment', 'nucleosome_signal')) /
  VlnPlot(h5b, c('TSS.enrichment', 'nucleosome_signal')) /
  VlnPlot(h6a, c('TSS.enrichment', 'nucleosome_signal'))

h5a <- subset(h5a, TSS.enrichment > 2 & nucleosome_signal < 3)
h5b <- subset(h5b, TSS.enrichment > 2 & nucleosome_signal < 3)
h6a <- subset(h6a, TSS.enrichment > 2 & nucleosome_signal < 3)

h5a <- RenameCells(h5a, add.cell.id = 'h5a')
h5b <- RenameCells(h5b, add.cell.id = 'h5b')
h6a <- RenameCells(h6a, add.cell.id = 'h6a')

h5a$sample <- 'h5a'
h5b$sample <- 'h5b'
h6a$sample <- 'h6a'

saveRDS(h5a, paste0(rds_path, 'cvejic_h5a.rds'))
saveRDS(h5b, paste0(rds_path, 'cvejic_h5b.rds'))
saveRDS(h6a, paste0(rds_path, 'cvejic_h6a.rds'))

# make a merged object
mg_peaks <- reduce(c(h5a[['peaks']]@ranges, h5b[['peaks']]@ranges, h6a[['peaks']]@ranges))
mg_peaks <- mg_peaks[mg_peaks@seqnames %in% paste0('chr', 1:22)]

# Filter out possible too short/long peaks
peakwidths <- width(mg_peaks)
mg_peaks <- mg_peaks[peakwidths < 10000 & peakwidths > 20]

mtx <- FeatureMatrix(fragments = c(Fragments(h5a), Fragments(h5b), Fragments(h6a)), features = mg_peaks, 
                     cells = c(Cells(h5a), Cells(h5b), Cells(h6a))) 

cvejic <- CreateSeuratObject(counts = CreateChromatinAssay(counts = mtx, 
                                                           fragments = c(Fragments(h5a), Fragments(h5b), Fragments(h6a)), 
                                                           min.cells = 10), 
                             meta.data = rbind(h5a@meta.data, h5b@meta.data, h6a@meta.data), 
                             assay = 'peaks')

Annotation(cvejic) <- annotations

cvejic <- RunTFIDF(cvejic, assay = 'peaks') %>% FindTopFeatures(min.cutoff = 'q75') %>% RunSVD()

cvejic <-
  RunHarmony(
    cvejic,
    group.by.vars = 'sample',
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F,
    reduction.save = 'harmony'
  )

cvejic <- FindNeighbors(object = cvejic, dims = 2:20, reduction = 'harmony')
cvejic <- FindClusters(object = cvejic, verbose = FALSE, algorithm = 3, resolution = 0.6)
cvejic <- RunUMAP(cvejic, reduction = 'harmony', dims = 2:20, return.model = TRUE, min.dist = 0.3,
                  spread = 1, n.neighbors = 30)

DimPlot(cvejic, label = TRUE, pt.size = 0.6) + coord_fixed() + NoLegend() 
DimPlot(cvejic, label = F, pt.size = 0.4, group.by = 'sample', shuffle = T, cols = annot_colors$sample_number2nd) + coord_fixed() 

table(cvejic$seurat_clusters)

# now re-run peak calling on a per-cluster basis
cvejic_peaks <- CallPeaks(cvejic, group.by = 'seurat_clusters',
                          macs2.path = '/Users/sa5488pa/Library/Python/3.9/bin/macs2', idents = 0:12)

peakwidths <- width(cvejic_peaks)
cvejic_peaks <- cvejic_peaks[peakwidths < 10000 & peakwidths > 20]
cvejic_peaks <- cvejic_peaks[cvejic_peaks@seqnames %in% paste0('chr', 1:22)]

mtx <- FeatureMatrix(fragments = Fragments(cvejic), features = cvejic_peaks, cells = Cells(cvejic)) 
cvejic[['peaks']] <- CreateChromatinAssay(counts = mtx, fragments = Fragments(cvejic), min.cells = 10)


cvejic <- RunTFIDF(cvejic, assay = 'peaks') %>% FindTopFeatures(min.cutoff = 'q75') %>% RunSVD()

cvejic <-
  RunHarmony(
    cvejic,
    group.by.vars = 'sample',
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F,
    reduction.save = 'harmony'
  )

cvejic <- FindNeighbors(object = cvejic, dims = 2:20, reduction = 'harmony')
cvejic <- FindClusters(object = cvejic, verbose = FALSE, algorithm = 3, resolution = 0.6)
cvejic <- RunUMAP(cvejic, reduction = 'harmony', dims = 2:20, return.model = TRUE, min.dist = 0.3,
                  spread = 1, n.neighbors = 30)

saveRDS(cvejic, paste0(rds_path, 'cvejic_merge.rds'))

DimPlot(cvejic, label = TRUE, pt.size = 0.6) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/cvejic_umap_presubset.pdf'))

DimPlot(cvejic, label = F, pt.size = 0.4, group.by = 'sample', shuffle = T, cols = annot_colors$sample_number2nd) + coord_fixed() 
ggsave(paste0(fig_path, 'additional_plots/cvejic_umap_presubset_samples.pdf'))
