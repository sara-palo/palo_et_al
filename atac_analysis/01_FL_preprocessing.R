source('scripts/atac_analysis/00.R')

counts <- Read10X_h5(filename = paste0(rawdata_path, 'filtered_peak_bc_matrix.h5'))

metadata <- read.csv(
  file = paste0(rawdata_path, 'singlecell.csv'),
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  genome = 'hg38',
  sep = c(":", "-"),
  fragments = paste0(rawdata_path, 'fragments.tsv.gz'),
  min.cells = 10,
  min.features = 200
)

sobj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata 
)

#add sample information to metadata
sobj$sample_number <- str_sub((Cells(sobj)), start=-1)

table(sobj$sample_number)

sobj$sample <- ifelse(sobj$sample_number <= 2, "FL", "BM") #based on cellranger report
table(sobj$sample_number, sobj$sample)

# subset to only contain FL
FL <- subset(sobj, sample == 'FL')
FL$sample_number <- paste0('FL', FL$sample_number)
table(FL$sample_number, FL$sample)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style 
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(FL) <- annotations

saveRDS(FL, paste0(rds_path, 'FL_only_preQC.rds'))

# remove features with 0 counts
FL <- FL[names(which(rowSums(GetAssayData(object = FL, slot = "counts")) > 0))]

# compute nucleosome signal
FL <- NucleosomeSignal(object = FL)

# compute TSS enrichment score
FL <- TSSEnrichment(object = FL, fast = TRUE)

# add blacklist ratio and fraction reads in peaks
FL$pct_reads_in_peaks <- FL$peak_region_fragments / FL$passed_filters * 100
FL$blacklist_ratio <- FL$blacklist_region_fragments / FL$peak_region_fragments

FL$nucleosome_group <- ifelse(FL$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = FL, group.by = 'nucleosome_group')

VlnPlot(
  object = FL,
  group.by = 'sample',
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5, 
  split.by = 'sample_number'
) & theme(text = element_text(size = 8))
ggsave(paste0(fig_path, 'additional_plots/FL_preQC_violin.pdf'), width = 12, height = 4.8)


FL <- subset(
  x = FL,
  subset = peak_region_fragments > 5000 &
    peak_region_fragments < 55000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.005 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2.5
)

FL

VlnPlot(
  object = FL,
  group.by = 'sample',
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5, 
  split.by = 'sample_number'
) & theme(text = element_text(size = 8))
ggsave(paste0(fig_path, 'additional_plots/FL_postQC_violin.pdf'), width = 12, height = 4.8)

FL <- RunTFIDF(FL)
FL <- FindTopFeatures(FL, min.cutoff = 'q75')
FL <- RunSVD(FL) 

DepthCor(FL)

#non-linear dimension reduction and clustering
FL <- RunUMAP(object = FL, reduction = 'lsi', dims = 2:30)
FL <- FindNeighbors(object = FL, reduction = 'lsi', dims = 2:30)
FL <- FindClusters(object = FL, verbose = FALSE, algorithm = 3, resolution = 1.3)

DimPlot(FL, group.by = 'sample_number') 

saveRDS(FL, paste0(rds_path, 'FL_only.rds'))
