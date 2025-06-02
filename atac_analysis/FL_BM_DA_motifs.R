source('scripts/atac_analysis/00.R')

library(ggvenn)
library(openxlsx)
library(UCell)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

DefaultAssay(FL) <- 'cisbpchromvar'

clusters_oi <- clust_order[2:14]

da_pos_FL <- list()
da_neg_FL <- list()

Idents(FL) <- 'pred.id.BM'

FL_df_up <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'up')
FL_df_dn <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'down')

for(i in 1:length(clusters_oi)) {
  
  ident <- as.character(clusters_oi[i])
  
  da_motifs <- FindMarkers(
    object = FL,
    assay = 'cisbpchromvar',
    ident.1 = ident,
    ident.2 = 'HSC',
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff'
  )
  
  if(ident == 'Ba/M'){
    write.csv(da_motifs, paste0(csv_path_ch, '/vs_prim/FL_', 'BaM', '_vsHSC.csv'))
  } else {
    write.csv(da_motifs, paste0(csv_path_ch, '/vs_prim/FL_', ident, '_vsHSC.csv'))
  }  
  
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
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.8, vjust = 0.8), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('TF motifs up- or downregulated compared to HSCs')
ggsave(paste0(fig_path, 'additional_plots/FL_tf_da_vs_hsc.pdf'), width = 9)


pos_list <- lapply(da_pos_FL, function(x) rownames(x)) 
pos_df <- data.frame(sapply(pos_list, function(x) x[1:max(lengths(pos_list))]))
write.csv(pos_df, paste0(csv_path_ch, '/vs_prim/FL_perannot_up.csv'), row.names = F)

neg_list <- lapply(da_neg_FL, function(x) rownames(x)) 
neg_df <- data.frame(sapply(neg_list, function(x) x[1:max(lengths(neg_list))]))
write.csv(neg_df, paste0(csv_path_ch, '/vs_prim/FL_perannot_dn.csv'), row.names = F)

ggvenn(list('Ly-II' = pos_df$Ly.II, 'GMP-II' = pos_df$GMP.II, 
            'Ery-II' = pos_df$Ery.II), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II', 'Ery-II')])) +
  labs(title = 'TF motifs increased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/FL_tf_up_venn_labeltr.pdf'), width = 7, height = 6)

ggvenn(list('Ly-II' = neg_df$Ly.II, 'GMP-II' = neg_df$GMP.II, 
            'Ery-II' = neg_df$Ery.II), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II', 'Ery-II')])) +
  labs(title = 'TF motifs decreased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/FL_tf_dn_venn_labeltr.pdf'), width = 7, height = 6)

resly1 <- prep_data(dataset = 'human', genes = na.omit(pos_df$Ly.I))
resly2 <- prep_data(dataset = 'human', genes = na.omit(pos_df$Ly.II))
resgm1 <- prep_data(dataset = 'human', genes = na.omit(pos_df$GMP.I))
resgm2 <- prep_data(dataset = 'human', genes = na.omit(pos_df$GMP.II))
reser2 <- prep_data(dataset = 'human', genes = na.omit(pos_df$Ery.II))

orca(p = plot_cellradar(res_list = list(resly1, resly2), colors = unname(annot_colors$cluster_cols[c('Ly-I', 'Ly-II')]), names = c('Ly-I', 'Ly-II')),
     paste0(fig_path, 'additional_plots/FL_LyI_LyII.pdf'), more_args = '--disable-gpu')

orca(p = plot_cellradar(res_list = list(resgm1, resgm2), colors = unname(annot_colors$cluster_cols[c('GMP-I', 'GMP-II')]), names = c('GMP-I', 'GMP-II')),
     paste0(fig_path, 'additional_plots/FL_GMPI_GMPII.pdf'), more_args = '--disable-gpu')

orca(p = plot_cellradar(res_list = list(resly2, resgm2), colors = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II')]), names = c('Ly-II', 'GMP-II')),
     paste0(fig_path, 'fig4/FL_LyII_GMPII_cellradar.pdf'), more_args = '--disable-gpu')

orca(p = plot_cellradar(res_list = list(resly2, resgm2, reser2), colors = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II', 'Ery-II')]), names = c('Ly-II', 'GMP-II', 'Ery-II')),
     paste0(fig_path, 'additional_plots/FL_LyII_GMPII_EryII.pdf'), more_args = '--disable-gpu')

missinggenes <- data.frame('ly1' = str_flatten(resly1$nonexisting_genes, collapse = ','), 
                           'ly2' = str_flatten(resly2$nonexisting_genes, collapse = ','),
                           'gmp1' = str_flatten(resgm1$nonexisting_genes, collapse = ','),
                           'gmp2' = str_flatten(resgm2$nonexisting_genes, collapse = ','),
                           'ery2' = str_flatten(reser2$nonexisting_genes, collapse = ','))

write.table(missinggenes, paste0(sup_table_path, '/FL_annot_CR_missinggenes.txt'))

table(na.omit(pos_df$Ly.II) %in% resly2$nonexisting_genes)
table(na.omit(pos_df$GMP.II) %in% resgm2$nonexisting_genes)

# do the same for adult BM
da_pos_BM <- list()
da_neg_BM <- list()

Idents(BM) <- 'BM_annot'

BM_df_up <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'up')
BM_df_dn <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'down')

for(i in 1:length(clusters_oi)) {
  
  ident <- as.character(clusters_oi[i])
  
  da_motifs <- FindMarkers(
    object = BM,
    assay = 'cisbpchromvar',
    ident.1 = ident,
    ident.2 = 'HSC',
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff'
  )
  
  if(ident == 'Ba/M'){
    write.csv(da_motifs, paste0(csv_path_ch, '/vs_prim/BM_', 'BaM', '_vsHSC.csv'))
  } else {
    write.csv(da_motifs, paste0(csv_path_ch, '/vs_prim/BM_', ident, '_vsHSC.csv'))
  }  
  
  da_pos_BM[[ident]] <- subset(da_motifs, subset = avg_diff > 1.5 & p_val_adj < 0.05)
  da_neg_BM[[ident]] <- subset(da_motifs, subset = avg_diff < -1.5 & p_val_adj < 0.05)
  
  BM_df_up[i,'count'] <- dim(da_pos_BM[[ident]])[1]
  BM_df_up[i,'cluster'] <- ident
  
  BM_df_dn[i,'count'] <- dim(da_neg_BM[[ident]])[1]
  BM_df_dn[i,'cluster'] <- ident
}

df_da_BM <- rbind(BM_df_up, BM_df_dn)
df_da_BM$cluster <- factor(df_da_BM$cluster, levels = clusters_oi)
df_da_BM$up_or_dn <- factor(df_da_BM$up_or_dn, levels = c('up', 'down'))

ggplot(df_da_BM, aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.7, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.8, vjust = 0.8), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('TF motifs up- or downregulated compared to HSCs')
ggsave(paste0(fig_path, 'additional_plots/BM_tf_da_vs_hsc.pdf'), width = 9)

df_da$sample <- 'FL'
df_da_BM$sample <- 'BM'
df_da_both <- rbind(df_da, df_da_BM)
df_da_both$x <- paste0(df_da_both$up_or_dn, ' in ', df_da_both$sample)
df_da_both$x <- factor(df_da_both$x, levels = c('up in FL', 'down in FL', 'up in BM', 'down in BM'))

ggplot(subset(df_da_both, cluster %in% c('Ly-II', 'GMP-II', 'Ery-II')), aes(x = cluster, y = count, fill = x, color = x)) + 
  geom_col(data = subset(df_da_both, cluster %in% c('Ly-II', 'GMP-II', 'Ery-II') & sample == 'FL'), width = 0.35, just = 1.05, position = position_dodge2(width = 0.7)) +
  geom_col(data = subset(df_da_both, cluster %in% c('Ly-II', 'GMP-II', 'Ery-II') & sample == 'BM'), width = 0.35, just = -0.05, position = position_dodge2(width = 0.7)) +
  scale_fill_manual('up or down', values = c("up in FL"= alpha('#840000', 1), "down in FL"= alpha('#081b80', 1), 
                                             "up in BM"= alpha('#840000', 0.4), "down in BM"= alpha('#081b80', 0.4))) + 
  scale_color_manual(values = annot_colors$da_cols2, guide = 'none') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('TF motifs up- or downregulated compared to HSCs')
ggsave(paste0(fig_path, 'fig4/both_tf_da_vs_hsc_subs.pdf'), width = 6)


pos_list <- lapply(da_pos_BM, function(x) rownames(x)) 
pos_df_BM <- data.frame(sapply(pos_list, function(x) x[1:max(lengths(pos_list))]))
write.csv(pos_df_BM, paste0(csv_path_ch, '/vs_prim/BM_perclust_up.csv'), row.names = F)

neg_list <- lapply(da_neg_BM, function(x) rownames(x)) 
neg_df_BM <- data.frame(sapply(neg_list, function(x) x[1:max(lengths(neg_list))]))
write.csv(neg_df_BM, paste0(csv_path_ch, '/vs_prim/BM_perclust_dn.csv'), row.names = F)

df_list <- list(pos_df, neg_df, pos_df_BM, neg_df_BM)
names(df_list) = c('up in FL clusters vs HSC', 'down in FL clusters vs HSC', 'up in BM clusters vs HSC', 'down in BM clusters vs HSC')

wb <- createWorkbook('table S6')

for(name in names(df_list)){
  
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = df_list[[name]], keepNA = F)
}
saveWorkbook(wb, file = paste0(sup_table_path, '/tableS6.xlsx'), overwrite = TRUE)


FL <- AddModuleScore_UCell(FL,
                           features = list('BM.Ly.II' = na.omit(pos_df_BM$Ly.II), 'BM.GMP.II' = na.omit(pos_df_BM$GMP.II), 
                                           'BM.Ery.II' = na.omit(pos_df_BM$Ery.II), 'BM.GMP.I' = na.omit(pos_df_BM$GMP.I), 'BM.MEP.II' = na.omit(pos_df_BM$MEP.II)), 
                           assay = 'cisbpchromvar', slot = 'data', maxRank = 200)

BM <- AddModuleScore_UCell(BM, 
                           features = list('Ly.II' = na.omit(pos_df_BM$Ly.II), 'GMP.II' = na.omit(pos_df_BM$GMP.II), 
                                           'Ery.II' = na.omit(pos_df_BM$Ery.II), 'GMP.I' = na.omit(pos_df_BM$GMP.I), 'MEP.II' = na.omit(pos_df_BM$MEP.II)), 
                           assay = 'cisbpchromvar', slot = 'data', maxRank = 200)

FL <- AddModuleScore_UCell(FL,
                           features = list('Ly.II' = na.omit(pos_df$Ly.II), 'GMP.II' = na.omit(pos_df$GMP.II), 
                                           'Ery.II' = na.omit(pos_df$Ery.II), 'GMP.I' = na.omit(pos_df$GMP.I), 'MEP.II' = na.omit(pos_df$MEP.II)), 
                           assay = 'cisbpchromvar', slot = 'data', maxRank = 200)

BM <- AddModuleScore_UCell(BM,
                           features = list('FL.Ly.II' = na.omit(pos_df$Ly.II), 'FL.GMP.II' = na.omit(pos_df$GMP.II), 
                                           'FL.Ery.II' = na.omit(pos_df$Ery.II), 'FL.GMP.I' = na.omit(pos_df$GMP.I), 'FL.MEP.II' = na.omit(pos_df$MEP.II)), 
                           assay = 'cisbpchromvar', slot = 'data', maxRank = 200)

saveRDS(FL, paste0(rds_path, 'FL_only.rds'))
saveRDS(BM, paste0(rds_path, 'BM_only.rds'))

ggvenn(list('Ly-II' = pos_df_BM$Ly.II, 'GMP-II' = pos_df_BM$GMP.II, 
            'Ery-II' = pos_df_BM$Ery.II), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II', 'Ery-II')])) +
  labs(title = 'TF motifs increased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/BM_tf_up_venn.pdf'), width = 7, height = 6)

ggvenn(list('Ly-II' = neg_df_BM$Ly.II, 'GMP-II' = neg_df_BM$GMP.II, 
            'Ery-II' = neg_df_BM$Ery.II), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II', 'Ery-II')])) +
  labs(title = 'TF motifs decreased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/BM_tf_dn_venn.pdf'), width = 7, height = 6)

ggvenn(list('FL Ly-II' = pos_df$Ly.II, 'BM Ly-II' = pos_df_BM$Ly.II), 
       show_percentage = FALSE, text_size = 14, fill_alpha = 0.7,  
       fill_color = unname(annot_colors$tissue_cols)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = 1))
ggsave(filename = paste0(fig_path, 'sup5/shared_LyII_up_venn.pdf'), width = 7, height = 6)

df <- data.frame(sapply(list('up FL only' = setdiff(pos_df$Ly.II, pos_df_BM$Ly.II), 
                             'up in both' = na.omit(intersect(pos_df$Ly.II, pos_df_BM$Ly.II)), 
                             'up BM only' = setdiff(pos_df_BM$Ly.II, pos_df$Ly.II)), 
                        function(x) x[1:max(lengths(neg_list))]))
write.csv(df, paste0(sup_table_path, '/up_LyII_FLBM.csv'), row.names = F)

ggvenn(list('FL Ly-II' = neg_df$Ly.II, 'BM Ly-II' = neg_df_BM$Ly.II), 
       show_percentage = FALSE, text_size = 14, fill_alpha = 0.7,  
       fill_color = unname(annot_colors$tissue_cols)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = 1))
ggsave(filename = paste0(fig_path, 'sup5/shared_LyII_dn_venn.pdf'), width = 7, height = 6)

ggvenn(list('FL GMP-II' = pos_df$GMP.II, 'BM GMP-II' = pos_df_BM$GMP.II), 
       show_percentage = FALSE, text_size = 14, fill_alpha = 0.7,  
       fill_color = unname(annot_colors$tissue_cols)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/shared_GMPII_up_venn.pdf'), width = 7, height = 6)

df <- data.frame(sapply(list('up FL only' = setdiff(pos_df$GMP.II, pos_df_BM$GMP.II), 
                             'up in both' = na.omit(intersect(pos_df$GMP.II, pos_df_BM$GMP.II)), 
                             'up BM only' = setdiff(pos_df_BM$GMP.II, pos_df$GMP.II)), 
                        function(x) x[1:max(lengths(neg_list))]))
write.csv(df, paste0(sup_table_path, '/up_GMPII_FLBM.csv'), row.names = F)

df <- data.frame(sapply(list('dn FL only' = setdiff(neg_df$Ly.II, neg_df_BM$Ly.II), 
                             'dn in both' = na.omit(intersect(neg_df$Ly.II, neg_df_BM$Ly.II)), 
                             'dn BM only' = setdiff(neg_df_BM$Ly.II, neg_df$Ly.II)), 
                        function(x) x[1:max(lengths(neg_list))]))
write.csv(df, paste0(sup_table_path, '/dn_LyII_FLBM.csv'), row.names = F)

df <- data.frame(sapply(list('dn FL only' = setdiff(neg_df$GMP.II, neg_df_BM$GMP.II), 
                             'dn in both' = na.omit(intersect(neg_df$GMP.II, neg_df_BM$GMP.II)), 
                             'dn BM only' = setdiff(neg_df_BM$GMP.II, neg_df$GMP.II)), 
                        function(x) x[1:max(lengths(neg_list))]))
write.csv(df, paste0(sup_table_path, '/dn_GMPII_FLBM.csv'), row.names = F)

ggvenn(list('FL GMP-II' = neg_df$GMP.II, 'BM GMP-II' = neg_df_BM$GMP.II), 
       show_percentage = FALSE, text_size = 14, fill_alpha = 0.7,  
       fill_color = unname(annot_colors$tissue_cols)) + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/shared_GMPII_dn_venn.pdf'), width = 7, height = 6)

resly1 <- prep_data(dataset = 'human', genes = na.omit(pos_df_BM$Ly.I))
resly2 <- prep_data(dataset = 'human', genes = na.omit(pos_df_BM$Ly.II))
resgm1 <- prep_data(dataset = 'human', genes = na.omit(pos_df_BM$GMP.I))
resgm2 <- prep_data(dataset = 'human', genes = na.omit(pos_df_BM$GMP.II))
reser2 <- prep_data(dataset = 'human', genes = na.omit(pos_df_BM$Ery.II))

orca(p = plot_cellradar(res_list = list(resly1, resly2), colors = unname(annot_colors$cluster_cols[c('Ly-I', 'Ly-II')]), names = c('Ly-I', 'Ly-II')),
     paste0(fig_path, 'additional_plots/BM_LyI_LyII.pdf'), more_args = '--disable-gpu')

orca(p = plot_cellradar(res_list = list(resgm1, resgm2), colors = unname(annot_colors$cluster_cols[c('GMP-I', 'GMP-II')]), names = c('GMP-I', 'GMP-II')),
     paste0(fig_path, 'additional_plots/BM_GMPI_GMPII.pdf'), more_args = '--disable-gpu')

orca(p = plot_cellradar(res_list = list(resly2, resgm2), colors = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II')]), 
                        names = c('Ly-II', 'GMP-II'), opac = c(0.8, 0.8), dash = c('4px', '4px')),
     paste0(fig_path, 'fig4/BM_LyII_GMPII.pdf'), more_args = '--disable-gpu')

orca(p = plot_cellradar(res_list = list(resly2, resgm2, reser2), colors = unname(annot_colors$cluster_cols[c('Ly-II', 'GMP-II', 'Ery-II')]), 
                        names = c('Ly-II', 'GMP-II', 'Ery-II')),
     paste0(fig_path, 'additional_plots/BM_LyII_GMPII_EryII.pdf'), more_args = '--disable-gpu')

missinggenes <- data.frame('ly1' = str_flatten(resly1$nonexisting_genes, collapse = ','), 
                           'ly2' = str_flatten(resly2$nonexisting_genes, collapse = ','),
                           'gmp1' = str_flatten(resgm1$nonexisting_genes, collapse = ','),
                           'gmp2' = str_flatten(resgm2$nonexisting_genes, collapse = ','),
                           'ery2' = str_flatten(reser2$nonexisting_genes, collapse = ','))

write.table(missinggenes, paste0(sup_table_path, '/BM_annot_CR_missinggenes.txt'))

table(na.omit(pos_df_BM$Ly.II) %in% resly2$nonexisting_genes)
table(na.omit(pos_df_BM$GMP.II) %in% resgm2$nonexisting_genes)
