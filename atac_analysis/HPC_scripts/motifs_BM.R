set.seed(12345)

packages <- c('Signac', 'Seurat', 'BSgenome.Hsapiens.UCSC.hg38', 'ggplot2', 'stringr')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/atac/rds_files/'
fig_path <- '/home/spalo/atac/figures/'

cluster_cols <- c(
  'HSC' = '#E9CB53',
  'MPP-I' = '#52BE80',
  'MPP-III' = '#16A085',
  'Ly-I' = '#A569BD',
  'Ly-II' = '#793885',
  'DC' = '#4e3666',
  'GMP-II' = '#2874A6',
  'GMP-I' = '#56ace5',
  'MPP-II' = '#E7A235',
  'MEP-I' = '#DC7633',
  'MEP-II' = '#E46457',
  'Ery-I' = '#B03A2E',
  'Ery-II' = '#78281F',
  'Ba/M' = '#F06292'
)

BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

DefaultAssay(BM) <- 'peaks'
Idents(BM) <- 'BM_annot'

# update fragment object due to change in file path
Fragments(BM) <- NULL
Fragments(BM) <- CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = Cells(BM))

motifs_oi <- c('RUNX1', 'EBF1', 'HOXB7', 'IRF8', 'PAX5', 'FOS', 'JUN', 'SPIB')

for (motif in motifs_oi){
  BM <- Footprint(object = BM, genome = BSgenome.Hsapiens.UCSC.hg38, motif.name = motif, compute.expected = T, in.peaks = T, assay = 'peaks')

  p <- PlotFootprint(BM, motif, idents = c('HSC', 'Ly-II', 'GMP-II'), group.by = 'BM_annot') & scale_color_manual(values = cluster_cols)
  ggsave(plot = p, filename = paste0(fig_path, 'BM_', motif, '_footprint.pdf'), width = 8, height = 6)   
 }



