# create separate fragment files etc for the two FL donors

source('scripts/atac_analysis/00.R')

counts <- Read10X_h5(filename = paste0(rawdata_path, 'filtered_peak_bc_matrix.h5'))

# make a seurat object containing all cells and then extract metadata, cellranger peaks and fragments for each FL donor
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
  min.cells = 0,
  min.features = 0
)

sobj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata 
)

sobj$sample_number <- str_sub((Cells(sobj)), start=-1)
table(sobj$sample_number)

sobj$sample <- ifelse(sobj$sample_number <= 2, "FL", "BM") 
table(sobj$sample_number, sobj$sample)

FL <- subset(sobj, sample == 'FL')
FL$sample_number <- paste0('FL', FL$sample_number)

Idents(FL) <- 'sample_number'
table(Idents(FL))

# remove peaks with 0 counts
FL <- FL[names(which(rowSums(GetAssayData(object = FL, slot = "counts")) > 0))]

dim(FL[['peaks']]@counts)
dim(counts)

counts[1:5,1:5]
counts1 <- FL[['peaks']]@counts[, WhichCells(FL, idents = 'FL1')]
counts2 <- FL[['peaks']]@counts[, WhichCells(FL, idents = 'FL2')]

dim(counts1)
dim(counts2)

md1 <- metadata[WhichCells(FL, idents = 'FL1'), ]
md2 <- metadata[WhichCells(FL, idents = 'FL2'), ]

write.table(x = counts1, file = 'FL_data_perdonor/feature_bc_matrix_FL1.tsv', col.names = colnames(counts1))
write.table(x = counts2, file = 'FL_data_perdonor/feature_bc_matrix_FL2.tsv', col.names = colnames(counts2))

write.csv(md1, 'FL_data_perdonor/metadata_FL1.csv')
write.csv(md2, 'FL_data_perdonor/metadata_FL2.csv')

SplitFragments(FL, assay = 'peaks', group.by = 'sample_number', outdir = 'FL_data_perdonor')
