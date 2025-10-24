source('scripts/atac_analysis/00.R')

library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(motifmatchr)
library(GenomeInfoDb)

merged <- readRDS(paste0(rds_path, 'merged_with_ranz.rds'))

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

frags <- list(CreateFragmentObject(path = paste0(rawdata_path, 'fragments.tsv.gz'), cells = WhichCells(merged, expression = tissue == 'CS22_FL')),
              CreateFragmentObject(path = paste0(rawdata_path, 'fragments.tsv.gz'), cells = WhichCells(merged, expression = tissue == 'ABM')),
              CreateFragmentObject(path = paste0(ext_data_path, 'ranzoni/frags2.tsv.gz'), cells = WhichCells(merged, expression = dataset == 'ranzoni'))
)

Fragments(merged) <- NULL
Fragments(merged) <- frags

saveRDS(merged, paste0(rds_path, 'merged_with_ranz.rds'))
