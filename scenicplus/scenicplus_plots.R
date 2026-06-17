source('scripts/atac_analysis/00.R')

library(hdf5r)
library(MuData)
library(MultiAssayExperiment)
library(dendextend)
library(ggrepel)

mdata <- readH5MU(paste0(fig_path, 'scenicplus_outs/scplusmdata.h5mu'))

gene_AUC_dir <- assay(mdata[['direct_gene_based_AUC']])
reg_AUC_dir <- assay(mdata[['direct_region_based_AUC']])

reg <- read.table(paste0(fig_path, 'scenicplus_outs/eRegulon_direct.tsv'), header = T)

annot_df <- data.frame('cluster' = stringr::str_split_i(colnames(gene_AUC_dir), pattern = '_', i = 1))
column_ha <- HeatmapAnnotation(df = annot_df, col = list('cluster' = c('FLHSC' = '#F2C25A', 'BMHSC' = '#D99466', 'FLLy' = '#6A4791', 'BMLy' = '#853870',
                                                                       'FLMy' = '#3E8DD6', 'BMMy' = '#127182', 'FLEry' = '#BF4236', 'BMEry' = '#964536')))

# filter for +/+  and -/- interactions
act <- reg[grep(pattern = '_\\+\\/\\+', x = reg$eRegulon_name), ]
rep <- reg[grep(pattern = '_\\-\\/\\-', x = reg$eRegulon_name), ]

act_genes <- unique(act$Gene_signature_name)
rep_genes <- unique(rep$Gene_signature_name)
act_regs <- unique(act$Region_signature_name)
rep_regs <- unique(rep$Region_signature_name)

act_genemat <- gene_AUC_dir[act_genes,]
act_regmat <- reg_AUC_dir[act_regs,]
rep_genemat <- gene_AUC_dir[rep_genes,]
rep_regmat <- reg_AUC_dir[rep_regs,]

# remove rows with mostly (> 95%) 0 AUC scores and reorder columns
act_genemat <- act_genemat[rowSums(act_genemat > 0) > ncol(act_genemat)*0.05, ]
act_regmat <- act_regmat[rowSums(act_regmat > 0) > ncol(act_regmat)*0.05, ]

# trim regulon names
rownames(act_genemat) <- sub("_.*", "", rownames(act_genemat))
rownames(act_regmat) <- sub("_.*", "", rownames(act_regmat))

act_genemat <- act_genemat[intersect(rownames(act_genemat), rownames(act_regmat)), ]
act_regmat <- act_regmat[intersect(rownames(act_genemat), rownames(act_regmat)), ]

act_genemat_sc <- scale(t(act_genemat))
act_genemat_sc <- t(act_genemat_sc)
act_regmat_sc <- scale(t(act_regmat))
act_regmat_sc <- t(act_regmat_sc)

# plot all scaled and unscaled heatmaps, row order based on the region score clustering
hm <- Heatmap(act_regmat_sc, top_annotation = column_ha, show_column_names = F)
hm <- draw(hm)
hm

# row_order <- row_order(hm)
row_dend = as.dendrogram(hclust(dist(act_regmat_sc)))

pdf(paste0(fig_path, 'sup3/scenic_act_reg_hm_dendr.pdf'), width = 7.5, height = 6.5)
Heatmap(act_regmat_sc, top_annotation = column_ha, show_column_names = F, cluster_rows = row_dend)
dev.off()

row_dend = as.dendrogram(hclust(dist(act_genemat_sc)))

pdf(paste0(fig_path, 'sup3/scenic_act_gene_hm_dendr.pdf'), width = 7.5, height = 6.5)
Heatmap(act_genemat_sc, top_annotation = column_ha, show_column_names = F, cluster_rows = row_dend)
dev.off()



pdf(paste0(fig_path, 'additional_plots/scenic_unsc_act_reg_hm.pdf'), width = 7.5, height = 6.5)
Heatmap(act_regmat, top_annotation = column_ha, show_column_names = F)
dev.off()

pdf(paste0(fig_path, 'additional_plots/scenic_unsc_act_gene_hm.pdf'), width = 7.5, height = 6.5)
Heatmap(act_genemat, top_annotation = column_ha, show_column_names = F)
dev.off()



rep_genemat <- rep_genemat[rowSums(rep_genemat > 0) > ncol(rep_genemat)*0.05, ]
rep_regmat <- rep_regmat[rowSums(rep_regmat > 0) > ncol(rep_regmat)*0.05, ]

rownames(rep_genemat) <- sub("_.*", "", rownames(rep_genemat))
rownames(rep_regmat) <- sub("_.*", "", rownames(rep_regmat))

rep_genemat <- rep_genemat[intersect(rownames(rep_genemat), rownames(rep_regmat)), ]
rep_regmat <- rep_regmat[intersect(rownames(rep_genemat), rownames(rep_regmat)), ]

rep_genemat_sc <- scale(t(rep_genemat))
rep_genemat_sc <- t(rep_genemat_sc)
rep_regmat_sc <- scale(t(rep_regmat))
rep_regmat_sc <- t(rep_regmat_sc)

pdf(paste0(fig_path, 'additional_plots/scenic_rep_reg_hm.pdf'), width = 6.5, height = 4)
Heatmap(rep_regmat_sc, top_annotation = column_ha, show_column_names = F)
dev.off()

pdf(paste0(fig_path, 'additional_plots/scenic_rep_gene_hm.pdf'), width = 6.5, height = 4)
Heatmap(rep_genemat_sc, top_annotation = column_ha, show_column_names = F)
dev.off()

pdf(paste0(fig_path, 'additional_plots/scenic_unsc_rep_reg_hm.pdf'), width = 6.5, height = 4)
Heatmap(rep_regmat, top_annotation = column_ha, show_column_names = F)
dev.off()

pdf(paste0(fig_path, 'additional_plots/scenic_unsc_rep_gene_hm.pdf'), width = 6.5, height = 4)
Heatmap(rep_genemat, top_annotation = column_ha, show_column_names = F)
dev.off()

# HSC heatmaps
hsc_genes <- act_genemat[, grep('HSC', colnames(act_genemat))]
hsc_regs <- act_regmat[, grep('HSC', colnames(act_regmat))]

gene_diff <- c()
pval_gene <- c()

for(i in 1:nrow(hsc_genes)){
  fl <- hsc_genes[i, grep('FL', colnames(hsc_genes))]
  bm <- hsc_genes[i, grep('BM', colnames(hsc_genes))]
  mean_diff <- mean(fl) - mean(bm)
  gene_diff[i] <- mean_diff
  pval_gene[i] <- wilcox.test(fl, bm)$p.val
}

gene_stats <- data.frame('mean_diff' = gene_diff, 'p_val' = pval_gene, 'p_adj' = p.adjust(pval_gene, method = 'bonf'), row.names = rownames(hsc_genes))
write.csv(gene_stats, paste0(sup_table_path, '/scenic_hsc_gene_stats.csv'))


reg_diff <- c()
pval_reg <- c()

for(i in 1:nrow(hsc_regs)){
  fl <- hsc_regs[i, grep('FL', colnames(hsc_regs))]
  bm <- hsc_regs[i, grep('BM', colnames(hsc_regs))]
  mean_diff <- mean(fl) - mean(bm)
  reg_diff[i] <- mean_diff
  pval_reg[i] <- wilcox.test(fl, bm)$p.val
}

reg_stats <- data.frame('mean_diff' = reg_diff, 'p_val' = pval_reg, 'p_adj' = p.adjust(pval_reg, method = 'bonf'), row.names = rownames(hsc_regs))
write.csv(reg_stats, paste0(sup_table_path, '/scenic_hsc_reg_stats.csv'))

# subset for heatmaps
gene_stats <- gene_stats[gene_stats$p_adj < 0.05, ]
reg_stats <- reg_stats[reg_stats$p_adj < 0.05, ]

hm_genes <- rownames(gene_stats)[c(head(order(gene_stats$mean_diff, decreasing = T), 5), head(order(gene_stats$mean_diff, decreasing = F), 5))]
hm_regs <- rownames(reg_stats)[c(head(order(reg_stats$mean_diff, decreasing = T), 5), head(order(reg_stats$mean_diff, decreasing = F), 5))]

hsc_genemat <- hsc_genes[hm_genes, ]
hsc_regmat <- hsc_regs[hm_regs, ]

annot_df <- data.frame('cluster' = stringr::str_split_i(colnames(hsc_genemat), pattern = '_', i = 1))
column_ha <- HeatmapAnnotation(df = annot_df, col = list('cluster' = c('FLHSC' = '#F2C25A', 'BMHSC' = '#D99466')))

pdf(paste0(fig_path, 'additional_plots/scenic_unsc_hsc_gene_hm.pdf'), width = 4.5, height = 3)
Heatmap(hsc_genemat, top_annotation = column_ha, show_column_names = F, 
        column_order = c(grep('FL', colnames(hsc_genes)), grep('BM', colnames(hsc_genes))))
dev.off()

pdf(paste0(fig_path, 'additional_plots/scenic_unsc_hsc_reg_hm.pdf'), width = 4.5, height = 3)
Heatmap(hsc_regmat, top_annotation = column_ha, show_column_names = F, 
        column_order = c(grep('FL', colnames(hsc_genes)), grep('BM', colnames(hsc_genes))))
dev.off()

hsc_genemat <- scale(t(hsc_genemat))
hsc_genemat <- t(hsc_genemat)
hsc_regmat <- scale(t(hsc_regmat))
hsc_regmat <- t(hsc_regmat)


row_dend = as.dendrogram(hclust(dist(hsc_genemat)))

pdf(paste0(fig_path, 'fig3/scenic_hsc_gene_hm.pdf'), width = 4.5, height = 3)
Heatmap(hsc_genemat, top_annotation = column_ha, show_column_names = F, 
        column_order = c(grep('FL', colnames(hsc_genes)), grep('BM', colnames(hsc_genes))), cluster_rows = rev(row_dend))
dev.off()

row_dend = as.dendrogram(hclust(dist(hsc_regmat)))

pdf(paste0(fig_path, 'fig3/scenic_hsc_reg_hm.pdf'), width = 4.5, height = 3)
Heatmap(hsc_regmat, top_annotation = column_ha, show_column_names = F, 
        column_order = c(grep('FL', colnames(hsc_regs)), grep('BM', colnames(hsc_regs))), cluster_rows = rev(row_dend))
dev.off()


rss <- read.csv(paste0(fig_path, 'scenicplus_outs/rss_gene_region.csv'), row.names = 1, check.names = F)
rss <- as.matrix(rss)
rss <- rss[, grep(pattern = '_\\+\\/\\+', x = colnames(rss))]

FL_df <- data.frame('FLHSC_score' = rss['FLHSC', order(rss['FLHSC',], decreasing = T)], 'FLHSC_rank' = 1:ncol(rss), 'regulon' = names(rss['FLHSC', order(rss['FLHSC',], decreasing = T)]))
FL_df[! FL_df$regulon %in% c(grep('SPI', rownames(FL_df), value = T), grep('GATA', rownames(FL_df), value = T)), 'regulon'] <- ''

BM_df <- data.frame( 'BMHSC_score' = rss['BMHSC', order(rss['BMHSC',], decreasing = T)], 'BMHSC_rank' = 1:ncol(rss), 'regulon' = names(rss['BMHSC', order(rss['BMHSC',], decreasing = T)]))
BM_df[! BM_df$regulon %in% c(grep('SPI', rownames(BM_df), value = T), grep('GATA', rownames(BM_df), value = T)), 'regulon'] <- ''

ggplot(data = FL_df, mapping = aes(x = FLHSC_rank, y = FLHSC_score, label = regulon)) + geom_point() + theme_bw() +
  geom_label_repel(max.overlaps = 100)
ggsave(paste0(fig_path, 'additional_plots/SCENIC_FLHSC_rss_ranks.pdf'), width = 5.5, height = 4.5)

ggplot(data = BM_df, mapping = aes(x = BMHSC_rank, y = BMHSC_score, label = regulon)) + geom_point() + theme_bw() +
  geom_label_repel(max.overlaps = 100)
ggsave(paste0(fig_path, 'additional_plots/SCENIC_BMHSC_rss_ranks.pdf'), width = 5.5, height = 4.5)



FL_df <- data.frame('FLHSC_score' = rss['FLHSC', order(rss['FLHSC',], decreasing = T)], 'FLHSC_rank' = 1:ncol(rss), 'regulon' = names(rss['FLHSC', order(rss['FLHSC',], decreasing = T)]))
FL_df[11:nrow(FL_df), 'regulon'] <- ''

BM_df <- data.frame( 'BMHSC_score' = rss['BMHSC', order(rss['BMHSC',], decreasing = T)], 'BMHSC_rank' = 1:ncol(rss), 'regulon' = names(rss['BMHSC', order(rss['BMHSC',], decreasing = T)]))
BM_df[11:nrow(BM_df), 'regulon'] <- ''

ggplot(data = FL_df, mapping = aes(x = FLHSC_rank, y = FLHSC_score, label = regulon)) + geom_point() + theme_bw() +
  geom_label_repel(max.overlaps = 100) + ylim(c(0, max(c(FL_df$FLHSC_score, BM_df$BMHSC_score), na.rm = T) + 0.03)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/SCENIC_FLHSC_rss_ranks_top10.pdf'), width = 4.5, height = 3.8)

ggplot(data = BM_df, mapping = aes(x = BMHSC_rank, y = BMHSC_score, label = regulon)) + geom_point() + theme_bw() +
  geom_label_repel(max.overlaps = 100) + ylim(c(0, max(c(FL_df$FLHSC_score, BM_df$BMHSC_score), na.rm = T) + 0.03)) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/SCENIC_BMHSC_rss_ranks_top10.pdf'), width = 4.5, height = 3.8)


ggplot(data = FL_df, mapping = aes(x = FLHSC_rank, y = FLHSC_score, color = (FLHSC_rank > 10))) + geom_point() + theme_bw() +
  ylim(c(0, max(c(FL_df$FLHSC_score, BM_df$BMHSC_score), na.rm = T) + 0.03)) + scale_color_manual(values = c('red', 'black')) + NoLegend()
ggsave(paste0(fig_path, 'fig3/SCENIC_FLHSC_rss_ranks_coltop10.pdf'), width = 4.5, height = 3.8)

ggplot(data = BM_df, mapping = aes(x = BMHSC_rank, y = BMHSC_score, color = (BMHSC_rank > 10))) + geom_point() + theme_bw() +
  ylim(c(0, max(c(FL_df$FLHSC_score, BM_df$BMHSC_score), na.rm = T) + 0.03)) + scale_color_manual(values = c('red', 'black')) + NoLegend()
ggsave(paste0(fig_path, 'fig3/SCENIC_BMHSC_rss_ranks_coltop10.pdf'), width = 4.5, height = 3.8)

