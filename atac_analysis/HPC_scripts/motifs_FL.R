set.seed(12345)

packages <- c('Signac', 'Seurat', 'BSgenome.Hsapiens.UCSC.hg38', 'ggplot2', 'stringr')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/atac/rds_files/'
fig_path <- '/home/spalo/atac/figures/'

seurat_cluster_cols <- c(
  '0' = '#E9CB53',
  '1' = '#52BE80',
  '2' = '#16A085',
  '3' = '#A569BD',
  '4' = '#793885',
  '5' = '#E7A235',
  '6' = '#3498DB',
  '7' = '#58558f',
  '8' = '#2874A6',
  '9' = '#F06292',
  '10' = '#DC7633',
  '11' = '#E46457',
  '12' = '#78281F',
  '13' = '#4e342e'
)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

DefaultAssay(FL) <- 'peaks'
Idents(FL) <- 'seurat_clusters'

# update fragment object due to change in file path
Fragments(FL) <- NULL
Fragments(FL) <- CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = Cells(FL))

motifs_oi <- c('RUNX1', 'EBF1', 'HOXB7', 'IRF8', 'PAX5', 'FOS', 'JUN', 'SPIB')

for (motif in motifs_oi){
  FL <- Footprint(object = FL, genome = BSgenome.Hsapiens.UCSC.hg38, motif.name = motif, compute.expected = T, in.peaks = T, assay = 'peaks')

  p <- PlotFootprint(FL, motif, idents = c(0,4,8), group.by = 'seurat_clusters') & scale_color_manual(values = seurat_cluster_cols)
  ggsave(plot = p, filename = paste0(fig_path, 'FL_', motif, '_footprint.pdf'), width = 8, height = 6)   
 }
#
saveRDS(FL[['peaks']], paste0(rds_path, 'FL_peaks_withfootprint.rds'))

sessionInfo()

