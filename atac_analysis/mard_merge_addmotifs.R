source('scripts/atac_analysis/00.R')

library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(motifmatchr)
library(GenomeInfoDb)

merged <- readRDS(paste0(rds_path, 'merged_with_mard.rds'))

data('human_pwms_v2')
cisbp <- human_pwms_v2

motifnames <- lapply(cisbp, function(x) x@name[[1]]) %>% unlist(use.names = F)
names(cisbp) <- motifnames

cisbp

# fix the misnaming of some motifs
errors <- grep('LINE', motifnames)

for(err in errors){
  print(cisbp[[err]]@ID)
  correctname <- str_split_i(cisbp[[err]]@ID, '_', i = 1)
  names(cisbp)[err] <- correctname
}

merged <- AddMotifs(object = merged,
                    genome = BSgenome.Hsapiens.UCSC.hg38,
                    pfm = cisbp,
                    assay = 'peaks')

h5acells <- WhichCells(merged, expression = sample_number == 'h5a')
names(h5acells) <- str_remove(h5acells, pattern = 'h5a_')
h5bcells <- WhichCells(merged, expression = sample_number == 'h5b')
names(h5bcells) <- str_remove(h5bcells, pattern = 'h5b_')
h6acells <- WhichCells(merged, expression = sample_number == 'h6a')
names(h6acells) <- str_remove(h6acells, pattern = 'h6a_')

frags <- list(CreateFragmentObject(path = paste0(rawdata_path, 'fragments.tsv.gz'), cells = WhichCells(merged, expression = tissue == 'CS22_FL')),
              CreateFragmentObject(path = paste0(rawdata_path, 'fragments.tsv.gz'), cells = WhichCells(merged, expression = tissue == 'ABM')),
              CreateFragmentObject(path = paste0(ext_data_path, 'marderstein/multi/H_5_A-atac_fragments.tsv.gz'), cells = names(h5acells)),
              CreateFragmentObject(path = paste0(ext_data_path, 'marderstein/multi/H_5_B-atac_fragments.tsv.gz'), cells = names(h5bcells)),
              CreateFragmentObject(path = paste0(ext_data_path, 'marderstein/multi/H_6_A-atac_fragments.tsv.gz'), cells = names(h6acells))
)

frags[[3]] <- RenameCells(frags[[3]], h5acells)
frags[[4]] <- RenameCells(frags[[4]], h5bcells)
frags[[5]] <- RenameCells(frags[[5]], h6acells)

Fragments(merged) <- NULL
Fragments(merged) <- frags

saveRDS(merged, paste0(rds_path, 'merged_with_mard.rds'))
