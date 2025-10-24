set.seed(12345)

packages <- c('Signac', 'Seurat', 'stringr', 'ggplot2', 'GenomicRanges', 'BSgenome.Hsapiens.UCSC.hg38')
invisible(lapply(packages, library, character.only = TRUE))

rds_path <- '/home/spalo/atac/rds_files/'
fig_path <- '/home/spalo/atac/figures/'

merged <- readRDS(paste0(rds_path, 'merged_with_mard.rds'))

h5acells <- WhichCells(merged, expression = sample_number == 'h5a')
names(h5acells) <- str_remove(h5acells, pattern = 'h5a_')
h5bcells <- WhichCells(merged, expression = sample_number == 'h5b')
names(h5bcells) <- str_remove(h5bcells, pattern = 'h5b_')
h6acells <- WhichCells(merged, expression = sample_number == 'h6a')
names(h6acells) <- str_remove(h6acells, pattern = 'h6a_')

frags <- list(CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = WhichCells(merged, expression = tissue == 'CS22_FL')),
              CreateFragmentObject(path = '~/atac/cvejic_frags/H_5_A-atac_fragments.tsv.gz', cells = names(h5acells)),
              CreateFragmentObject(path = '~/atac/cvejic_frags/H_5_B-atac_fragments.tsv.gz', cells = names(h5bcells)),
              CreateFragmentObject(path = '~/atac/cvejic_frags/H_6_A-atac_fragments.tsv.gz', cells = names(h6acells)),
              CreateFragmentObject(path = '~/atac/merged_FL12_BM12/outs/fragments.tsv.gz', cells = WhichCells(merged, expression = tissue == 'ABM'))
)

frags[[2]] <- RenameCells(frags[[2]], h5acells)
frags[[3]] <- RenameCells(frags[[3]], h5bcells)
frags[[4]] <- RenameCells(frags[[4]], h6acells)

Fragments(merged) <- NULL
Fragments(merged) <- frags

merged$plotclusters <- paste0(merged$tissue, '_', merged$clust.names)
Idents(merged) <- 'plotclusters'

motifs_oi <- c('EBF1', 'IRF8', 'PAX5', 'TCF3', 'IKZF1', 'RUNX1', 'SPI1', 'JUN', 'FOS', 'SPIB', 'IRF4')

for (motif in motifs_oi){
  
  merged <- Footprint(object = merged, genome = BSgenome.Hsapiens.UCSC.hg38, motif.name = motif, compute.expected = T, in.peaks = T, assay = 'peaks')
  
  p <- PlotFootprint(merged, motif, group.by = 'plotclusters', idents = c('CS22_FL_Ly-II', 'trim2_FL_Ly-II', 'ABM_Ly-II')) & 
    scale_color_manual(values = c('CS22_FL_Ly-II' = '#1c6494', 'trim2_FL_Ly-II' = '#E8AB09', 'ABM_Ly-II' = '#1b4b1f')) 
  ggsave(plot = p, filename = paste0(fig_path, 'merged_mard_', motif, '_footprint_LyII_inpeaksT.pdf'), width = 8, height = 6)   
}

