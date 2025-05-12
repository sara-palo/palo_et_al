source('scripts/atac_analysis/00.R')

library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(GenomeInfoDb)
library(chromVARmotifs)

merged <- readRDS(paste0(rds_path, 'merged.rds'))

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

merged <- RegionStats(merged, assay = 'peaks', genome = BSgenome.Hsapiens.UCSC.hg38)

saveRDS(merged, paste0(rds_path, 'merged.rds'))
