source('scripts/atac_analysis/00.R')
library(SeuratDisk)

merged_atac <- readRDS(paste0(rds_path, 'merged.rds'))
FL_CITE <- readRDS(paste0(rds_path, 'FL_CITE.rds'))
BM_CITE <- readRDS(paste0(rds_path, 'BM_CITE.rds'))

# Fix cell type annotation and Label FL cluster 0 as HSCs
FL_HSC <- WhichCells(merged_atac, expression = seurat_clusters == 0 & sample == 'FL')

merged_atac$scenic_ident <- paste0(merged_atac$sample, '_atac_', merged_atac$clust.names)

Idents(merged_atac) <- 'scenic_ident'
merged_atac <- RenameIdents(merged_atac, c('BM_atac_HSC' = 'BMHSC', 'BM_atac_Ly-II' = 'BMLy', 'FL_atac_Ly-II' = 'FLLy', 'BM_atac_Ba/M' = 'BM_atac_BaM',
                                           'FL_atac_Ba/M' = 'FL_atac_BaM', 'BM_atac_Ery-I' = 'BMEry', 'BM_atac_Ery-II' = 'BMEry', 
                                           'FL_atac_Ery-I' = 'FLEry', 'FL_atac_Ery-II' = 'FLEry', 'FL_atac_MEP-II' = 'FLEry',
                                           'FL_atac_GMP-II' = 'FLMy', 'BM_atac_GMP-II' = 'BMMy'))

merged_atac$scenic_ident <- as.character(Idents(merged_atac))

merged_atac@meta.data[FL_HSC, 'scenic_ident'] <- 'FLHSC'
table(merged_atac$scenic_ident)
Idents(merged_atac) <- 'scenic_ident'

DimPlot(merged_atac, split.by = 'sample', label = T, pt.size = 0.8, shuffle = T)

FL_CITE$tissue <- 'FL'
BM_CITE$tissue <- 'BM'

# rename relevant clusters to have same name in CITE and ATAC
BM_newID <-  c('HSC-I' = 'BMHSC', 'Ly-II' = 'BMLy', 'GMP' = 'BMMy', 'MEP-I' = 'BMEry', 'MEP-II' = 'BMEry')
BM_notused <- unique(BM_CITE$paper_clusters)[! unique(BM_CITE$paper_clusters) %in% names(BM_newID)]
BM_notused_nm <- paste0('BM_cite_', BM_notused)
names(BM_notused_nm) <- BM_notused

BM_CITE <- RenameIdents(BM_CITE, c(BM_newID, BM_notused_nm))
BM_CITE$scenic_ident <- Idents(BM_CITE)


FL_newID <-  c('Ly-II' = 'FLLy', 'DC-Mono' = 'FLMy', 'MEP' = 'FLEry', 'HSC' = 'FLHSC')
FL_notused <- unique(FL_CITE$clust.names)[! unique(FL_CITE$clust.names) %in% names(FL_newID)]
FL_notused_nm <- paste0('FL_cite_', FL_notused)
names(FL_notused_nm) <- FL_notused

FL_CITE <- RenameIdents(FL_CITE, c(FL_newID, FL_notused_nm))
FL_CITE$scenic_ident <- Idents(FL_CITE)

plot_colors <- c('FLHSC' = '#EDCB5A', 'BMHSC' = '#E0A65A', 'FLLy' = '#674791', 'BMLy' = '#853870',
                 'FLMy' = '#3E8DD6', 'BMMy' = '#126D82', 'FLEry' = '#BF4236', 'BMEry' = '#964536')

DimPlot(FL_CITE, reduction = 'umap', label = T, cols = plot_colors, pt.size = 0.8) + NoLegend() + coord_fixed()
ggsave(paste0(fig_path, 'sup3/FLCITE_scenic_clusts.pdf'), width = 8.5, height = 5)

DimPlot(BM_CITE, group.by = 'scenic_ident', reduction = 'umap', label = T, cols = plot_colors, pt.size = 0.8) + NoLegend() + coord_fixed()
ggsave(paste0(fig_path, 'sup3/BMCITE_scenic_clusts.pdf'), width = 8.5, height = 5)

DimPlot(merged_atac, label = T, cols = plot_colors, pt.size = 0.8, cells = WhichCells(merged_atac, expression = sample == 'FL')) + NoLegend() + coord_fixed()
ggsave(paste0(fig_path, 'sup3/FLATAC_scenic_clusts.pdf'), width = 8.5, height = 5)

DimPlot(merged_atac, label = T, cols = plot_colors, pt.size = 0.8, cells = WhichCells(merged_atac, expression = sample == 'BM')) + NoLegend() + coord_fixed()
ggsave(paste0(fig_path, 'sup3/BMATAC_scenic_clusts.pdf'), width = 8.5, height = 5)

merged_CITE <- merge(FL_CITE, BM_CITE)
Idents(merged_CITE) <- 'scenic_ident'

table(merged_CITE$scenic_ident)

FLc <- data.frame(table(merged_CITE$scenic_ident)[c('FLHSC', 'FLLy', 'FLMy', 'FLEry')])
FLc$Freq <- FLc$Freq/length(Cells(FL_CITE))*100

FLa <- data.frame(table(merged_atac$scenic_ident)[c('FLHSC', 'FLLy', 'FLMy', 'FLEry')])
FLa$Freq <- FLa$Freq/length(WhichCells(merged_atac, expression = sample == 'FL'))*100

fldf <- rbind(FLa, FLc)
fldf$mod <- c(rep('ATAC', 4), rep('CITE', 4))

ggplot(fldf, aes(x = Var1, y = Freq, fill = mod)) + geom_bar(stat = 'identity', position = 'dodge2', width = 0.5) + 
  theme_bw() + scale_fill_manual(values = c('ATAC' = '#363636', 'CITE' = '#BFBFBF'))
ggsave(paste0(fig_path, 'sup3/FL_scenic_clusts.pdf'))

BMc <- data.frame(table(merged_CITE$scenic_ident)[c('BMHSC', 'BMLy', 'BMMy', 'BMEry')])
BMc$Freq <- BMc$Freq/length(Cells(BM_CITE))*100

BMa <- data.frame(table(merged_atac$scenic_ident)[c('BMHSC', 'BMLy', 'BMMy', 'BMEry')])
BMa$Freq <- BMa$Freq/length(WhichCells(merged_atac, expression = sample == 'BM'))*100

bmdf <- rbind(BMa, BMc)
bmdf$mod <- c(rep('ATAC', 4), rep('CITE', 4))

ggplot(bmdf, aes(x = Var1, y = Freq, fill = mod)) + geom_bar(stat = 'identity', position = 'dodge2', width = 0.5) + 
  theme_bw() + scale_fill_manual(values = c('ATAC' = '#363636', 'CITE' = '#BFBFBF'))
ggsave(paste0(fig_path, 'sup3/BM_scenic_clusts.pdf'))

# export metadata and counts matrix for SCENIC+
write.table(merged_atac@meta.data, file = 'rds_files/scenicplus/merged_ATAC_md.tsv', quote = F, row.names = T, col.names = T, sep = '\t')
write.table(merged_CITE@meta.data, file = 'rds_files/scenicplus/merged_CITE_md.tsv', quote = F, row.names = T, col.names = T, sep = '\t')

# save peak matrix in two halves due to size
dim(merged_atac[['peaks']]@counts)
write.table(merged_atac[['peaks']]@counts[,1:10000], file = 'rds_files/scenicplus/merged_ATAC_rawc_1.tsv', quote = F, row.names = T, col.names = T, sep = '\t')
write.table(merged_atac[['peaks']]@counts[,10001:20366], file = 'rds_files/scenicplus/merged_ATAC_rawc_2.tsv', quote = F, row.names = T, col.names = T, sep = '\t')

# overwrite normalized and scaled data with raw counts
merged_CITE[['RNA']]@data <- merged_CITE[['RNA']]@counts

merged_CITE[['RNA']]@counts[1:5,1:3]
merged_CITE[['RNA']]@data[1:5,1:3]

SaveH5Seurat(merged_CITE, filename = 'rds_files/scenicplus/merged_CITE')
Convert('rds_files/scenicplus/merged_CITE.h5Seurat', dest = 'h5ad')

merged_peaks <- merged_atac[['peaks']]@ranges
merged_peaks <- as.data.frame(merged_peaks)

merged_bed <- merged_peaks[,c('seqnames', 'start', 'end', 'strand')]
colnames(merged_bed) <- c('chrom', 'chromStart', 'chromEnd', 'strand')
write.table(merged_bed, 'rds_files/scenicplus/mergedpeaks.bed', quote = F, row.names = F, sep = '\t', col.names = F)

CITE_clustorder <- c('HSC-I', 'HSC-II', 'HSC-III', 'MPP-I', 'MPP-II', 'MPP-III', 'Ly-I', 'Ly-II', 'DC-I', 'DC-II', 'GMP', 'MEP-I', 'MEP-II', 'MB', 'Cyc-I', 'Cyc-II')
CITE_clustcolors <- c(
  'HSC-I' = '#E9CB53',
  'HSC-II' = '#B5A253',
  'HSC-III' = '#C1C750',
  'MPP-I' = '#52BE80',
  'MPP-III' = '#16A085',
  'Ly-I' = '#A569BD',
  'Ly-II' = '#793885',
  'DC-I' = '#5C4078',
  'DC-II' = '#392257',
  'GMP' = '#2874A6',
  'MPP-II' = '#E7A235',
  'MEP-I' = '#E46457',
  'MEP-II' = '#B03A2E',
  'MB' = '#F06292',
  'Cyc-I' = '#3A6929',
  'Cyc-II' = '#1A4A08'
)

BM_CITE$paper_clusters <- factor(BM_CITE$paper_clusters, levels = CITE_clustorder)

DimPlot(BM_CITE, pt.size = 0.7, reduction = 'umap', group.by = 'paper_clusters', cols = CITE_clustcolors) + coord_fixed()
ggsave(paste0(fig_path, 'sup3/BM_CITE_newumap.pdf'), width = 8, height = 5)


