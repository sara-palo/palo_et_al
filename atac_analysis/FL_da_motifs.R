source('scripts/atac_analysis/00.R')

library(ggvenn)
library(openxlsx)
library(UCell)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

DefaultAssay(FL) <- 'cisbpchromvar'
Idents(FL) <- 'seurat_clusters'

clusters_oi <- c(1:13)

da_pos_FL <- list()
da_neg_FL <- list()

FL_df_up <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'up')
FL_df_dn <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'down')

for(i in 1:length(clusters_oi)) {
  
  ident <- as.character(clusters_oi[i])
  
  da_motifs <- FindMarkers(
    object = FL,
    assay = 'cisbpchromvar',
    ident.1 = ident,
    ident.2 = 0,
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff'
  )
  
  write.csv(da_motifs, paste0(csv_path_ch, '/vs_prim/FL_', i, '_vsc0.csv'))
  
  da_pos_FL[[ident]] <- subset(da_motifs, subset = avg_diff > 1.5 & p_val_adj < 0.05)
  da_neg_FL[[ident]] <- subset(da_motifs, subset = avg_diff < -1.5 & p_val_adj < 0.05)
  
  FL_df_up[i,'count'] <- dim(da_pos_FL[[ident]])[1]
  FL_df_up[i,'cluster'] <- ident
  
  FL_df_dn[i,'count'] <- dim(da_neg_FL[[ident]])[1]
  FL_df_dn[i,'cluster'] <- ident
}

df_da <- rbind(FL_df_up, FL_df_dn)
df_da$cluster <- factor(df_da$cluster, levels = clusters_oi)
df_da$up_or_dn <- factor(df_da$up_or_dn, levels = c('up', 'down'))

ggplot(df_da, aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.7, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('Number of TF motifs up- or downregulated compared to HSC-like cells')
ggsave(paste0(fig_path, 'sup3/FL_tf_da_vs_c0.pdf'), width = 9)

ggplot(subset(df_da,cluster %in% c('4', '8', '13')), aes(x = cluster, y = count, fill = up_or_dn)) + 
  geom_bar(stat = 'identity', width = 0.4, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('TF motifs up- or downregulated compared to HSC-like cells')
ggsave(paste0(fig_path, 'fig3/FL_tf_da_vs_c0_subs.pdf'), width = 6)


pos_list <- lapply(da_pos_FL, function(x) rownames(x)) 
pos_df <- data.frame(sapply(pos_list, function(x) x[1:max(lengths(pos_list))]))
colnames(pos_df) <- paste0('C', as.character(clusters_oi))
write.csv(pos_df, paste0(csv_path_ch, '/vs_prim/FL_perclust_up.csv'), row.names = F)

neg_list <- lapply(da_neg_FL, function(x) rownames(x)) 
neg_df <- data.frame(sapply(neg_list, function(x) x[1:max(lengths(neg_list))]))
colnames(neg_df) <- paste0('C', as.character(clusters_oi))
write.csv(neg_df, paste0(csv_path_ch, '/vs_prim/FL_perclust_dn.csv'), row.names = F)

# save the enriched motifs to excel sheets
df_list <- list(pos_df, neg_df)
names(df_list) = c('up in FL clusters vs C0', 'down in FL clusters vs C0')

wb <- createWorkbook('table S3')

for(name in names(df_list)){
  
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = df_list[[name]], keepNA = F)
}
saveWorkbook(wb, file = paste0(csv_path, 'sup/tableS3.xlsx'), overwrite = TRUE)

# module scores for upregulated motifs
DefaultAssay(FL) <- 'cisbpchromvar'
FL <- AddModuleScore_UCell(FL, features = list('C4' = na.omit(pos_df$C4), 'C8' = na.omit(pos_df$C8), 'C13' = na.omit(pos_df$C13), 
                                               'C2' = na.omit(pos_df$C2), 'C6' = na.omit(pos_df$C6)), 
                           assay = 'cisbpchromvar', slot = 'data', maxRank = 200)

saveRDS(FL, paste0(rds_path, 'FL_only.rds'))

VlnPlot(FL, 'C4_UCell', cols = annot_colors$seurat_cluster_cols, pt.size = 0.1, idents = c(0,4,8)) + ylim(-0.01, 1) +
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'fig3/FL_c4score_vln.pdf'), width = 4)

VlnPlot(FL, 'C8_UCell', cols = annot_colors$seurat_cluster_cols, pt.size = 0.1, idents = c(0,4,8)) + ylim(-0.01, 1) +
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'fig3/FL_c8score_vln.pdf'), width = 4)

VlnPlot(FL, 'C4_UCell', cols = annot_colors$seurat_cluster_cols, pt.size = 0.1) + ylim(-0.01, 1) +
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'sup3/FL_c4score_vln_full.pdf'), width = 6)

VlnPlot(FL, 'C8_UCell', cols = annot_colors$seurat_cluster_cols, pt.size = 0.1) + ylim(-0.01, 1) +
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'sup3/FL_c8score_vln_full.pdf'), width = 6)

VlnPlot(FL, 'C13_UCell', cols = annot_colors$seurat_cluster_cols, pt.size = 0.1) + ylim(-0.01, 1) +
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'additional_plots/FL_c13score_vln_full.pdf'), width = 6)

VlnPlot(FL, 'C6_UCell', cols = annot_colors$seurat_cluster_cols, pt.size = 0.1) + ylim(-0.01, 1) +
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'additional_plots/FL_c6score_vln_full.pdf'), width = 6)

VlnPlot(FL, 'C2_UCell', cols = annot_colors$seurat_cluster_cols, pt.size = 0.1) + ylim(-0.01, 1) +
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'additional_plots/FL_c2score_vln_full.pdf'), width = 6)

fl_df <- FL@meta.data[, c('seurat_clusters', 'C4_UCell', 'C8_UCell')]
fl_df <- subset(fl_df, seurat_clusters %in% c(4,8))
fl_df$seurat_clusters <- factor(fl_df$seurat_clusters, levels = c(4,8))

ggplot(fl_df, aes(x = C8_UCell, y = C4_UCell, color = seurat_clusters)) +
  scale_color_manual(values = annot_colors$seurat_cluster_cols[c('4','8')]) +
  geom_point() + xlim(0,1) + ylim(0,1) + theme_bw()
ggsave(paste0(fig_path, 'additional_plots/FL_C4_C8_score_scatter.pdf'))

ggvenn(list('cluster 4' = pos_df$C4, 'cluster 8' = pos_df$C8, 
            'cluster 13' = pos_df$C13), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$seurat_cluster_cols[c('4', '8', '13')])) +
  labs(title = 'TF motifs increased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'fig3/FL_tf_up_venn.pdf'), width = 7, height = 6)

ggvenn(list('cluster 4' = neg_df$C4, 'cluster 8' = neg_df$C8, 
            'cluster 13' = neg_df$C13), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$seurat_cluster_cols[c('4', '8', '13')])) +
  labs(title = 'TF motifs decreased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'sup3/FL_tf_dn_venn.pdf'), width = 7, height = 6)

res2 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C2))
res3 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C3))
res4 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C4))
res6 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C6))
res7 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C7))
res8 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C8))
res13 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C13))

orca(p = plot_cellradar(res_list = list(res4, res8), colors = unname(annot_colors$seurat_cluster_cols[c('4', '8')]), names = c('cluster 4', 'cluster 8')),
     paste0(fig_path, 'fig3/FL_CR_c4_8.pdf'))

orca(p = plot_cellradar(res_list = list(res6, res7, res8), colors = unname(annot_colors$seurat_cluster_cols[c('6', '7', '8')]), names = c('cluster 6', 'cluster 7', 'cluster 8')),
     paste0(fig_path, 'additional_plots/FL_CR_c6_7_8.pdf'))

orca(p = plot_cellradar(res_list = list(res2, res3, res4), colors = unname(annot_colors$seurat_cluster_cols[c('2', '3', '4')]), names = c('cluster 2', 'cluster 3', 'cluster 4')),
     paste0(fig_path, 'additional_plots/FL_CR_c2_3_4.pdf'), more_args = '--disable-gpu') # for some reason this specific plot gives some errors if not disabling gpu

orca(p = plot_cellradar(res_list = list(res4, res8, res13), colors = unname(annot_colors$seurat_cluster_cols[c('4', '8', '13')]), names = c('cluster 4', 'cluster 8', 'cluster 13')),
     paste0(fig_path, 'additional_plots/FL_CR_c4_8_13.pdf'))


missinggenes <- data.frame('C2' = str_flatten(res2$nonexisting_genes, collapse = ','), 
                           'C3' = str_flatten(res3$nonexisting_genes, collapse = ','),
                           'C4' = str_flatten(res4$nonexisting_genes, collapse = ','),
                           'C6' = str_flatten(res6$nonexisting_genes, collapse = ','),
                           'C7' = str_flatten(res7$nonexisting_genes, collapse = ','),
                           'C8' = str_flatten(res8$nonexisting_genes, collapse = ','),
                           'C13' = str_flatten(res13$nonexisting_genes, collapse = ','))

write.table(missinggenes, paste0(sup_table_path, '/FL_CR_missinggenes.txt'))

table(na.omit(pos_df$C4) %in% res4$nonexisting_genes)
table(na.omit(pos_df$C8) %in% res8$nonexisting_genes)


