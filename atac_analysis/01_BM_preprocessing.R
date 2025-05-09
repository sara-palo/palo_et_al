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

# subset to only contain BM
BM <- subset(sobj, sample == 'BM')
table(BM$sample_number, BM$sample)
BM$sample_number <- ifelse(BM$sample_number == 3, 'BM1', 'BM2')
table(BM$sample_number, BM$sample)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style 
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(BM) <- annotations

#read bone marrow metadata from previous project
BM_meta <- read.csv(file = paste0(rawdata_path, 'BM_metadata.csv'), row.names = 1)
BM_meta$orig_cell_id <- rownames(BM_meta)
rownames(BM_meta) <- str_remove(rownames(BM_meta), pattern = 'atac_')

head(BM_meta)

#convert sample naming system to the one of the current object
rownames(BM_meta) <- gsub(rownames(BM_meta), pattern = '-1_BM10_34', replacement = '-3')
rownames(BM_meta) <- gsub(rownames(BM_meta), pattern = '-1_BM12_34', replacement = '-4')
colnames(BM_meta) <- c('BM_annot', 'old_UMAP1', 'old_UMAP2', 'orig_cell_id')

head(BM_meta)

# remove cells not in the BM metadata
BM <- subset(BM, cells = rownames(BM_meta))
BM <- AddMetaData(BM, metadata = BM_meta)

BM

head(BM@meta.data)

saveRDS(BM, paste0(rds_path, 'BM_only_preQC.rds'))

# remove features with 0 counts
BM <- BM[names(which(rowSums(GetAssayData(object = BM, slot = "counts")) > 0))]

# compute nucleosome signal
BM <- NucleosomeSignal(object = BM)

# compute TSS enrichment score
BM <- TSSEnrichment(object = BM, fast = TRUE)

# add blacklist ratio and fraction reads in peaks
BM$pct_reads_in_peaks <- BM$peak_region_fragments / BM$passed_filters * 100
BM$blacklist_ratio <- BM$blacklist_region_fragments / BM$peak_region_fragments

BM$nucleosome_group <- ifelse(BM$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = BM, group.by = 'nucleosome_group')

VlnPlot(
  object = BM,
  group.by = 'sample',
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5, 
  split.by = 'sample_number'
) & theme(text = element_text(size = 8))
ggsave(paste0(fig_path, 'additional_plots/BM_preQC_violin.pdf'), width = 12, height = 4.8)


BM <- subset(
  x = BM,
  subset = peak_region_fragments > 5000 &
    peak_region_fragments < 60000 &
    pct_reads_in_peaks > 45 &
    blacklist_ratio < 0.005 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2.5
)

BM

VlnPlot(
  object = BM,
  group.by = 'sample',
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5, 
  split.by = 'sample_number'
) & theme(text = element_text(size = 8))
ggsave(paste0(fig_path, 'additional_plots/BM_postQC_violin.pdf'), width = 12, height = 4.8)

BM <- BM[names(which(rowSums(GetAssayData(object = BM, slot = "counts")) > 0))]

BM

BM <- RunTFIDF(BM)
BM <- FindTopFeatures(BM, min.cutoff = 'q75')
BM <- RunSVD(BM) 

DepthCor(BM)

#non-linear dimension reduction and clustering
BM <- RunUMAP(object = BM, reduction = 'lsi', dims = 2:30)
BM <- FindNeighbors(object = BM, reduction = 'lsi', dims = 2:30)
BM <- FindClusters(object = BM, verbose = FALSE, algorithm = 3, resolution = 1.3)

DimPlot(BM, group.by = 'sample_number') 

saveRDS(BM, paste0(rds_path, 'BM_only.rds'))

