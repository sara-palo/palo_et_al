source('scripts/atac_analysis/00.R')

library(monocle3)
library(cicero)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(readxl)

merged <- readRDS(paste0(rds_path, 'merged.rds'))

Idents(merged) <- 'sample'

FLs <- subset(merged, idents = 'FL')
BMs <- subset(merged, idents = 'BM')

# create cell dataset
FLcds <- new_cell_data_set(expression_data = FLs[['peaks']]@counts, cell_metadata = FLs@meta.data)
BMcds <- new_cell_data_set(expression_data = BMs[['peaks']]@counts, cell_metadata = BMs@meta.data)

# Transfer Seurat UMAP embeddings
FLcds@int_colData@listData$reducedDims$UMAP <- FLs@reductions$mergedUMAP@cell.embeddings
FLcds@clusters$UMAP$clusters <- FLs$clust.names
plot_cells(FLcds)

BMcds@int_colData@listData$reducedDims$UMAP <- BMs@reductions$mergedUMAP@cell.embeddings
BMcds@clusters$UMAP$clusters <- BMs$clust.names
plot_cells(BMcds)

# run with lower k since there are quite few cells in e.g. HSC cluster
FLccds <- make_cicero_cds(FLcds, reduced_coordinates = reducedDims(FLcds)$UMAP, k = 30)
BMccds <- make_cicero_cds(BMcds, reduced_coordinates = reducedDims(BMcds)$UMAP, k = 30)

hg38 <- data.frame('chr' = seqnames(BSgenome.Hsapiens.UCSC.hg38), 'lengths' = seqlengths(BSgenome.Hsapiens.UCSC.hg38))

FLconns <- run_cicero(FLccds, hg38, sample_num = 100)
head(FLconns)

BMconns <- run_cicero(BMccds, hg38, sample_num = 100)
head(BMconns)

ccanFL <- generate_ccans(FLconns)
ccanBM <- generate_ccans(BMconns)

saveRDS(FLconns, paste0(rds_path, 'FLconns.rds'))
saveRDS(BMconns, paste0(rds_path, 'BMconns.rds'))
saveRDS(ccanFL, paste0(rds_path, 'ccanFL.rds'))
saveRDS(ccanBM, paste0(rds_path, 'ccanBM.rds'))

# look for distal peaks associated with e.g. lymphoid gene promoters:

# define the CollapseToLongestTranscript function from the Signac package to make it accessible
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

# Define promoter regions as the 2kb just upstream of the full, protein coding genes. Define distal peaks as peaks >2500 bp from any annotated region, protein coding or not. 
genes <- Annotation(merged)

full_genes <- CollapseToLongestTranscript(ranges = genes)
full_genes

# Filter for protein coding genes
full_genes <- full_genes[full_genes$gene_biotype == 'protein_coding']

table(duplicated(full_genes$gene_name))
range(width(full_genes))

promoters <- promoters(full_genes, upstream = 2000, downstream = 0) 
head(promoters)

feats <- ClosestFeature(merged, regions = rownames(merged[['peaks']]))
dist <- subset(feats, distance > 2500)

# Filter the co-accessible peaks using the same cutoff as what was automatically selected when finding CCANs
# default co-accessibility cutoffs were 0.22 for FL and 0.3 for BM
hist(FLconns$coaccess, breaks = 100, xlim = c(0.01, 1))
abline(v = 0.22)

hist(BMconns$coaccess, breaks = 100, xlim = c(0.01, 1))
abline(v = 0.3)

FLconns <- subset(FLconns, coaccess > 0.22) 
BMconns <- subset(BMconns, coaccess > 0.3)


# Find promoter and distal peaks
FLprom <- ClosestFeature(merged, regions = StringToGRanges(unique(FLconns$Peak1)), annotation = promoters) %>% subset(distance == 0)
FLconns$peak1_promoter <- ifelse(FLconns$Peak1 %in% FLprom$query_region, 'promoter', 'no')
FLconns$peak2_distal <- ifelse(FLconns$Peak2 %in% dist$query_region, 'distal', 'no')

BMprom <- ClosestFeature(merged, regions = StringToGRanges(unique(BMconns$Peak1)), annotation = promoters) %>% subset(distance == 0)
BMconns$peak1_promoter <- ifelse(BMconns$Peak1 %in% BMprom$query_region, 'promoter', 'no')
BMconns$peak2_distal <- ifelse(BMconns$Peak2 %in% dist$query_region, 'distal', 'no')

# filter for promoter-distal connections
prom_dist_FL <- FLconns[FLconns$peak1_promoter == 'promoter' & FLconns$peak2_distal == 'distal', ]
prom_dist_BM <- BMconns[BMconns$peak1_promoter == 'promoter' & BMconns$peak2_distal == 'distal', ]

prom_dist_FL <- cbind(prom_dist_FL, ClosestFeature(merged, regions = prom_dist_FL$Peak1, annotation = promoters))
prom_dist_BM <- cbind(prom_dist_BM, ClosestFeature(merged, regions = prom_dist_BM$Peak1, annotation = promoters))

table(prom_dist_FL$gene_name %in% prom_dist_BM$gene_name)


# read DA peaks between FL and BM for each cluster, filter for log2FC and p-adj
files <- list.files(paste0(csv_path_da, '/FLvsBM'), pattern = 'merged_da_', full.names = T)
da_peaks <- lapply(files, read.csv)

names <- lapply(files, function(x) str_split_i(x, '_', i = 5))
names(da_peaks) <- lapply(names, function(x) str_remove(x, pattern = '.csv'))

da_peaks <- lapply(da_peaks, function (x) subset(x, p_val_adj < 0.05 & abs(avg_log2FC) > 2))

# make data frame and annotate by which promoter(s) the distal peak is connected to
FL_df <- data.frame(matrix(0, nrow = nrow(prom_dist_FL), ncol = 6 + length(da_peaks)))
colnames(FL_df) <- c(colnames(prom_dist_FL[i, c(1:3, 6:8)]), names(da_peaks))

for (i in 1:nrow(prom_dist_FL)){
  
  data <- prom_dist_FL[i, c(1:3, 6:8)]
  
  da_list <- rep('no', length(names(da_peaks)))
  names(da_list) <- names(da_peaks)
  
  for (clust in names(da_peaks)){
    
    if(as.character(prom_dist_FL[i, 'Peak2']) %in% da_peaks[[clust]]$X){
      da_list[clust] <- ifelse(da_peaks[[clust]][da_peaks[[clust]]$X == as.character(prom_dist_FL[i, 'Peak2']), 'avg_log2FC'] > 0, 'up_FL', 'up_BM')
    }
  }
  
  # fix weird problem with peak column
  cols <- c(data, da_list)
  cols[['Peak2']] <- as.character(cols[['Peak2']])
  
  FL_df[i,] <- unlist(cols, use.names = F)
}


BM_df <- data.frame(matrix(0, nrow = nrow(prom_dist_BM), ncol = 6 + length(da_peaks)))
colnames(BM_df) <- c(colnames(prom_dist_BM[i, c(1:3, 6:8)]), names(da_peaks))

for (i in 1:nrow(prom_dist_BM)){
  
  data <- prom_dist_BM[i, c(1:3, 6:8)]
  
  da_list <- rep('no', length(names(da_peaks)))
  names(da_list) <- names(da_peaks)
  
  for (clust in names(da_peaks)){
    
    if(as.character(prom_dist_BM[i, 'Peak2']) %in% da_peaks[[clust]]$X){
      da_list[clust] <- ifelse(da_peaks[[clust]][da_peaks[[clust]]$X == as.character(prom_dist_BM[i, 'Peak2']), 'avg_log2FC'] > 0, 'up_FL', 'up_BM')
    }
  }
  
  cols <- c(data, da_list)
  cols[['Peak2']] <- as.character(cols[['Peak2']])
  
  BM_df[i,] <- unlist(cols, use.names = F)
}

write.csv(FL_df, paste0(sup_table_path, '/cicero_FL_da.csv'))
write.csv(BM_df, paste0(sup_table_path, '/cicero_BM_da.csv'))

# subset to keep only peaks which are DA in any cluster
FL_df_da <- FL_df[apply(FL_df[, 7:20], 1, function(x) any(! x == 'no')), ]
BM_df_da <- BM_df[apply(BM_df[, 7:20], 1, function(x) any(! x == 'no')), ]

FL_df[FL_df$`Ly-II` == 'up_FL' & FL_df$gene_name == 'PAX5', ]
BM_df[BM_df$`Ly-II` == 'up_BM' & BM_df$gene_name == 'PAX5', ]

da_pax5 <- BM_df_da[BM_df_da$gene_name == 'PAX5', 'Peak2']

subs <- subset(merged, clust.names == 'Ly-II')
Idents(subs) <- 'sample'

# ChIP-seq tracks downloaded from ENCODE
bw_me1 <- list('B cell H3K4me1-1' = paste0(ext_data_path, 'chip_seq/ENCFF798CGE_27F_H3K4me1_pval.bigWig'),
               'B cell H3K4me1-2' = paste0(ext_data_path, 'chip_seq/ENCFF775QXE_21M_H3K4me1_pval.bigWig'))

bw_ac <- list('B cell H3K27ac-1' = paste0(ext_data_path, 'chip_seq/ENCFF576QGU_27F_H3K27ac_pval.bigWig'),
              'B cell H3K27ac-2' = paste0(ext_data_path, 'chip_seq/ENCFF435QXY_37M_H3K27ac_pval.bigWig'))

bw_me3 <- list('B cell H3K27me3-1' = paste0(ext_data_path, 'chip_seq/ENCFF960VJJ_27F_H3K27me3_pval.bigWig'),
               'B cell H3K27me3-2' = paste0(ext_data_path, 'chip_seq/ENCFF480VDY_21M_H3K27me3_pval.bigWig'),
               'B cell H3K27me3-3' = paste0(ext_data_path, 'chip_seq/ENCFF951MIP_37M_H3K27me3_pval.bigWig'))

bw_pax5 <- list('PAX5-1' = paste0(ext_data_path, 'chip_seq/ENCFF886UGF_PAX5_GM12878_pval.bigWig'),
                'PAX5-2' = paste0(ext_data_path, 'chip_seq/ENCFF914QGY_PAX5_GM12878_pval.bigWig'),
                'PAX5-3' = paste0(ext_data_path, 'chip_seq/ENCFF809RDU_PAX5_GM12892_pval.bigWig'),
                'PAX5-4' = paste0(ext_data_path, 'chip_seq/ENCFF125FTI_PAX5_GM12891_pval.bigWig'))

bw_ebf1 <- list('EBF1-1' = paste0(ext_data_path, 'chip_seq/ENCFF389SRN_GM12878_EBF1_pval.bigWig'),
                'EBF1-2' = paste0(ext_data_path, 'chip_seq/ENCFF810XRY_GM12878_EBF1_pval.bigWig'))


CoveragePlot(subs, full_genes[full_genes$gene_name == 'PAX5'])  & scale_fill_manual(values = annot_colors$tissue_cols) 
ggsave(paste0(fig_path, 'fig4/FL_BM_LyII_PAX5.pdf'), width = 6, height = 4)

a <- CoveragePlot(subs, full_genes[full_genes$gene_name == 'PAX5'])  & scale_fill_manual(values = annot_colors$tissue_cols) 
b <- BigwigTrack(bigwig = bw_ac, full_genes[full_genes$gene_name == 'PAX5'], bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K27ac-1' = '#242323', 'B cell H3K27ac-2' = '#9E9E9E'))
c <- BigwigTrack(bigwig = bw_me1, full_genes[full_genes$gene_name == 'PAX5'], bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K4me1-1' = '#242323', 'B cell H3K4me1-2' = '#9E9E9E'))
d <- BigwigTrack(bigwig = bw_pax5, full_genes[full_genes$gene_name == 'PAX5'], bigwig.scale = 'separate')
e <- BigwigTrack(bigwig = bw_ebf1, full_genes[full_genes$gene_name == 'PAX5'], bigwig.scale = 'separate')

CombineTracks(plotlist = list(a, b, c, d, e), heights = c(6,1,1,2,1))
ggsave(paste0(fig_path, 'additional_plots/FL_BM_LyII_PAX5_chip.pdf'), width = 10, height = 7)

f <- BigwigTrack(bigwig = bw_me3, region = full_genes[full_genes$gene_name == 'PAX5'], extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')

CombineTracks(plotlist = list(a, b, c, f, d, e), heights = c(6,1,1,2,2,1))
ggsave(paste0(fig_path, 'additional_plots/FL_BM_LyII_PAX5_chip_me3.pdf'), width = 10, height = 8)


CoveragePlot(subs, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000,  
             region.highlight = StringToGRanges(da_pax5)) & scale_fill_manual(values = annot_colors$tissue_cols) 
ggsave(paste0(fig_path, 'fig4/FL_BM_LyII_PAX5_dist.pdf'), width = 6, height = 4)

a <- CoveragePlot(subs, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000,
                  region.highlight = StringToGRanges(da_pax5)) & scale_fill_manual(values = annot_colors$tissue_cols) 
b <- BigwigTrack(bigwig = bw_ac, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K27ac-1' = '#242323', 'B cell H3K27ac-2' = '#9E9E9E'))
c <- BigwigTrack(bigwig = bw_me1, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K4me1-1' = '#242323', 'B cell H3K4me1-2' = '#9E9E9E'))
d <- BigwigTrack(bigwig = bw_pax5, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')
e <- BigwigTrack(bigwig = bw_ebf1, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')

CombineTracks(plotlist = list(a, b, c, d, e), heights = c(6,1,1,2,1))
ggsave(paste0(fig_path, 'fig4/FL_BM_LyII_PAX5_dist_chip.pdf'), width = 10, height = 7)

f <- BigwigTrack(bigwig = bw_me3, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')

CombineTracks(plotlist = list(a, b, c, f, d, e), heights = c(6,1,1,2,2,1))
ggsave(paste0(fig_path, 'additional_plots/FL_BM_LyII_PAX5_dist_chip_me3.pdf'), width = 10, height = 8)

subs2 <- subset(merged, clust.names == 'DC')
Idents(subs2) <- 'sample'

a <- CoveragePlot(subs2, region = FL_df_da[FL_df_da$gene_name == 'PAX5', 'Peak2'], extend.upstream = 10000, extend.downstream = 10000,
                  region.highlight = StringToGRanges(FL_df_da[FL_df_da$gene_name == 'PAX5', 'Peak2'])) & scale_fill_manual(values = annot_colors$tissue_cols) 
b <- BigwigTrack(bigwig = bw_ac, region = FL_df_da[FL_df_da$gene_name == 'PAX5', 'Peak2'], extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K27ac-1' = '#242323', 'B cell H3K27ac-2' = '#9E9E9E'))
c <- BigwigTrack(bigwig = bw_me1, region = FL_df_da[FL_df_da$gene_name == 'PAX5', 'Peak2'], extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K4me1-1' = '#242323', 'B cell H3K4me1-2' = '#9E9E9E'))
d <- BigwigTrack(bigwig = bw_pax5, region = FL_df_da[FL_df_da$gene_name == 'PAX5', 'Peak2'], extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')
e <- BigwigTrack(bigwig = bw_ebf1, region = FL_df_da[FL_df_da$gene_name == 'PAX5', 'Peak2'], extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')
f <- BigwigTrack(bigwig = bw_me3, region =FL_df_da[FL_df_da$gene_name == 'PAX5', 'Peak2'], extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')

CombineTracks(plotlist = list(a, b, c, f, d, e), heights = c(6,1,1,2,2,1))
ggsave(paste0(fig_path, 'additional_plots/FL_BM_DC_PAX5_dist_chip.pdf'), width = 10, height = 8)

# include 2nd trimester Ly cells in coverage plot
merged_mard <- readRDS(paste0(rds_path, 'merged_with_mard.rds'))

frags <- Fragments(merged_mard)

h5acells <- WhichCells(merged_mard, expression = sample_number == 'h5a')
names(h5acells) <- str_remove(h5acells, pattern = 'h5a_')

h5bcells <- WhichCells(merged_mard, expression = sample_number == 'h5b')
names(h5bcells) <- str_remove(h5bcells, pattern = 'h5b_')

h6acells <- WhichCells(merged_mard, expression = sample_number == 'h6a')
names(h6acells) <- str_remove(h6acells, pattern = 'h6a_')

frags[[1]] <- CreateFragmentObject(path = paste0(rawdata_path, 'fragments.tsv.gz'), WhichCells(merged_mard, expression = tissue == 'CS22_FL'))
frags[[2]] <- CreateFragmentObject(path = paste0(ext_data_path, 'marderstein/multi/H_5_A-atac_fragments.tsv.gz'), cells = names(h5acells))
frags[[3]] <- CreateFragmentObject(path =  paste0(ext_data_path, 'marderstein/multi/H_5_B-atac_fragments.tsv.gz'), cells = names(h5bcells))
frags[[4]] <- CreateFragmentObject(path =  paste0(ext_data_path, 'marderstein/multi/H_6_A-atac_fragments.tsv.gz'), cells = names(h6acells))
frags[[5]] <- CreateFragmentObject(path = paste0(rawdata_path, 'fragments.tsv.gz'), WhichCells(merged_mard, expression = tissue == 'ABM'))

frags[[2]] <- RenameCells(frags[[2]], h5acells)
frags[[3]] <- RenameCells(frags[[3]], h5bcells)
frags[[4]] <- RenameCells(frags[[4]], h6acells)

Fragments(merged_mard) <- NULL
Fragments(merged_mard) <- frags

subs <- subset(merged_mard, clust.names == 'Ly-II')
Idents(subs) <- 'tissue'

plot_cols <- c('CS22_FL' = '#1c6494', 'trim2_FL' = '#DE9E50', 'ABM' = '#1b4b1f')

a <- CoveragePlot(subs, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000,
                  region.highlight = StringToGRanges(da_pax5)) & scale_fill_manual(values = plot_cols) 
b <- BigwigTrack(bigwig = bw_ac, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K27ac-1' = '#242323', 'B cell H3K27ac-2' = '#9E9E9E'))
c <- BigwigTrack(bigwig = bw_me1, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate') & 
  scale_fill_manual(values = c('B cell H3K4me1-1' = '#242323', 'B cell H3K4me1-2' = '#9E9E9E'))
d <- BigwigTrack(bigwig = bw_pax5, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')
e <- BigwigTrack(bigwig = bw_ebf1, region = 'chr9-37369534-37375442', extend.upstream = 10000, extend.downstream = 10000, bigwig.scale = 'separate')

CombineTracks(plotlist = list(a, b, c, d, e), heights = c(6,1,1,2,1))
ggsave(paste0(fig_path, 'fig4/FL_BM_mard_LyII_PAX5_dist_chip.pdf'), width = 10, height = 7)

CoveragePlot(subs, full_genes[full_genes$gene_name == 'PAX5'])  & scale_fill_manual(values = plot_cols) 
ggsave(paste0(fig_path, 'fig4/FL_BM_mard_LyII_PAX5.pdf'), width = 7, height = 5)

CoveragePlot(subs, full_genes[full_genes$gene_name == 'DNTT'])  & scale_fill_manual(values = plot_cols) 
ggsave(paste0(fig_path, 'fig4/FL_BM_mard_LyII_DNTT.pdf'), width = 7, height = 5)
