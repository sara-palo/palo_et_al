source('scripts/atac_analysis/00.R')

FL <- readRDS(paste0(rds_path, 'FL_CITE.rds'))
BM <- readRDS(paste0(rds_path, 'BM_CITE.rds'))

FL$tissue <- 'FL'
BM$tissue <- 'BM'

# merge BM HSC clusters into one
BM <- RenameIdents(BM, c('HSC-I' = 'HSC', 'HSC-II' = 'HSC', 'HSC-III' = 'HSC'))
BM$clust.names <- Idents(BM)

merge <- merge(FL, BM)
merge <- NormalizeData(merge, 'RNA')

merge$plotclusters <- paste0(merge$tissue, '_', merge$clust.names)
Idents(merge) <- 'plotclusters'

subs <- subset(merge, idents = c('FL_Ly-I', 'FL_Ly-II', 'BM_Ly-I', 'BM_Ly-II'))
levels(subs) <- c('FL_Ly-I', 'FL_Ly-II', 'BM_Ly-I', 'BM_Ly-II')

goi <- rev(c('TCF3', 'EBF1', 'PAX5', 'MYC'))

DotPlot(subs, features = goi) + coord_flip() + RotatedAxis()
ggsave(paste0(fig_path, 'fig4/CITEseq_dotplot.pdf'), width = 5, height = 4)
