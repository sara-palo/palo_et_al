set.seed(12345)

packages <- c('Signac', 'Seurat', 'stringr', 'ggplot2', 'GenomicRanges', 'BSgenome.Hsapiens.UCSC.hg38')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/atac/rds_files/'
fig_path <- '/home/spalo/atac/figures/'

merged <- readRDS(paste0(rds_path, 'merged_with_ranz.rds'))

frags <- list(CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = WhichCells(merged, expression = tissue == 'CS22_FL')),
              CreateFragmentObject(path = '/home/spalo/atac/ranzoni_frags/frags2.tsv.gz', cells = WhichCells(merged, expression = dataset == 'ranzoni')),
              CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = WhichCells(merged, expression = tissue == 'ABM'))
)

Fragments(merged) <- NULL
Fragments(merged) <- frags

merged$plotclusters <- paste0(merged$tissue, '_', merged$clust.names)
Idents(merged) <- 'plotclusters'

motifs_oi <- c('HOXB7', 'HOXB8', 'PBX1', 'NFIB', 'NFIA', 'NFIX', 'NFE2', 'JUN', 'FOS', 'BACH1', 'BACH2', 'FOSL2', 'SPIB')

merged <- Footprint(object = merged, genome = BSgenome.Hsapiens.UCSC.hg38, motif.name = motifs_oi, compute.expected = T, in.peaks = T, assay = 'peaks')

for (motif in motifs_oi){
  
  p <- PlotFootprint(merged, motif, group.by = 'plotclusters', idents = c('CS22_FL_HSC', 'trim2_FL_HSC', 'ABM_HSC')) & 
    scale_color_manual(values = c('CS22_FL_HSC' = '#1c6494', 'trim2_FL_HSC' = '#E8AB09', 'ABM_HSC' = '#1b4b1f')) 
  ggsave(plot = p, filename = paste0(fig_path, 'merged_ranz_', motif, '_footprint_HSC_inpeaksT.pdf'), width = 8, height = 6)   
  
    p <- PlotFootprint(merged, motif, group.by = 'plotclusters', idents = c('trim2_FL_HSC', 'trim2_BM_HSC')) & 
    scale_color_manual(values = c('trim2_FL_HSC' = '#E8AB09', 'trim2_BM_HSC' = '#5949EB')) 
  ggsave(plot = p, filename = paste0(fig_path, 'merged_ranz_', motif, '_footprint_HSC_inpeaksT_trim2.pdf'), width = 8, height = 6)   
}

warnings()

quit(save = 'no')
