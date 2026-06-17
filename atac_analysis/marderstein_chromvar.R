source('scripts/atac_analysis/00.R')

library(chromVARmotifs)
library(motifmatchr)
library(chromVAR)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggvenn)
library(ggfortify)

cvejic <- readRDS(paste0(rds_path, 'cvejic_subset.rds'))

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

cvejic <- AddMotifs(object = cvejic, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = cisbp, assay = 'peaks')
cvejic <- RunChromVAR(object = cvejic, genome = BSgenome.Hsapiens.UCSC.hg38, assay = 'peaks', new.assay.name = 'cisbpchromvar')

saveRDS(cvejic, paste0(rds_path, 'cvejic_subset.rds'))

DefaultAssay(cvejic) <- 'cisbpchromvar'
cvejic <- ScaleData(cvejic, assay = 'cisbpchromvar')
motifs_oi <- c('EBF1', 'HOXB7', 'CEBPA', 'GATA1', 'PAX5', 'IRF8', 'SPI1', 'NFE2', 'TCF3')

FeaturePlot(cvejic, features = motifs_oi, min.cutoff = 'q50', ncol = 3)
ggsave(paste0(fig_path, 'additional_plots/mard_motifs.pdf'), width = 12, height = 10)

# first run findmarker for each cluster using the seurat clusters to make a dot plot
Idents(cvejic) <- 'seurat_clusters'
levels(cvejic) <- rev(levels(cvejic))

motifs_chromv <- data.frame()
motifs_chromv_list <- list()

for(i in levels(cvejic)) {
  
  da_motifs <- FindMarkers(
    object = cvejic,
    assay = 'cisbpchromvar',
    ident.1 = i,
    only.pos = TRUE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff'
  )
  
  if (min(da_motifs$p_val_adj) < 0.05 &&
      max(da_motifs$avg_diff) > 1.5) {
    top_motifs <- da_motifs[da_motifs$p_val_adj < 0.05, ]
    top_motifs <- da_motifs[da_motifs$avg_diff > 1.5, ]
    
    motifs_chromv_list[[i]] <- top_motifs
    
    top_motifs_diff <-
      top_motifs[order(top_motifs$avg_diff, decreasing = TRUE), ]
    
    if (length(rownames(top_motifs_diff)) > 3) {
      df <- top_motifs_diff[1:3, ]
      
    } else {
      df <- top_motifs_diff
    }
    
    df$cluster_num <- i
    df$motif <- rownames(df)
    rownames(df) <- NULL
    
    motifs_chromv <- rbind(motifs_chromv, df)
  }
}

DotPlot(cvejic, features = rev(unique(motifs_chromv$motif)), cols = annot_colors$expr_cols_2) + RotatedAxis() + 
  theme(panel.background = element_rect(fill = 'white')) + 
  plot_annotation(title = 'DA TF motifs', theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(paste0(fig_path, 'additional_plots/mard_dotplot.pdf'), width = 15, height = 7.5)


clusters_oi <- c(1:13)

da_pos_mard <- list()
da_neg_mard <- list()

Idents(cvejic) <- 'pred.id.CS22'

mard_df_up <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'up')
mard_df_dn <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'down')

for(i in 1:length(clusters_oi)) {
  
  ident <- as.character(clusters_oi[i])
  
  da_motifs <- FindMarkers(
    object = cvejic,
    assay = 'cisbpchromvar',
    ident.1 = ident,
    ident.2 = 0,
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff'
  )
  
  write.csv(da_motifs, paste0(csv_path_ch, '/vs_prim/mard_', i, '_vsc0.csv'))
  
  da_pos_mard[[ident]] <- subset(da_motifs, subset = avg_diff > 1.5 & p_val_adj < 0.05)
  da_neg_mard[[ident]] <- subset(da_motifs, subset = avg_diff < -1.5 & p_val_adj < 0.05)
  
  mard_df_up[i,'count'] <- dim(da_pos_mard[[ident]])[1]
  mard_df_up[i,'cluster'] <- ident
  
  mard_df_dn[i,'count'] <- dim(da_neg_mard[[ident]])[1]
  mard_df_dn[i,'cluster'] <- ident
}

df_da <- rbind(mard_df_up, mard_df_dn)
df_da$cluster <- factor(df_da$cluster, levels = clusters_oi)
df_da$up_or_dn <- factor(df_da$up_or_dn, levels = c('up', 'down'))

ggplot(df_da, aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.7, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('Number of TF motifs up- or downregulated compared to HSC-like cells')
ggsave(paste0(fig_path, 'additional_plots/mard_tf_da_vs_c0.pdf'), width = 9)

ggplot(subset(df_da, cluster %in% c(4,8,13)), aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.4, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = annot_colors$da_cols) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, vjust = 1), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  ggtitle('Number of TF motifs up- or downregulated compared to HSC-like cells')
ggsave(paste0(fig_path, 'additional_plots/mard_tf_da_vs_c0_subs.pdf'), width = 6)

pos_list <- lapply(da_pos_mard, function(x) rownames(x)) 
pos_df <- data.frame(sapply(pos_list, function(x) x[1:max(lengths(pos_list))]))
colnames(pos_df) <- paste0('C', as.character(clusters_oi))
write.csv(pos_df, paste0(csv_path_ch, '/vs_prim/mard_perclust_up.csv'), row.names = F)

neg_list <- lapply(da_neg_mard, function(x) rownames(x)) 
neg_df <- data.frame(sapply(neg_list, function(x) x[1:max(lengths(neg_list))]))
colnames(neg_df) <- paste0('C', as.character(clusters_oi))
write.csv(neg_df, paste0(csv_path_ch, '/vs_prim/mard_perclust_dn.csv'), row.names = F)

# save the enriched motifs to excel sheets
library(openxlsx)

df_list <- list(pos_df, neg_df)
names(df_list) = c('up in 12pcw clusters vs C0', 'down in 12pcw clusters vs C0')

wb <- createWorkbook('marderstein_da_vs_C0')

for(name in names(df_list)){
  
  addWorksheet(wb, sheetName = name)
  writeData(wb, sheet = name, x = df_list[[name]], keepNA = F)
}
saveWorkbook(wb, file = paste0(csv_path, 'sup/marderstein_da_vs_C0.xlsx'), overwrite = TRUE)


ggvenn(list('cluster 4' = pos_df$C4, 'cluster 8' = pos_df$C8, 
            'cluster 13' = pos_df$C13), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$seurat_cluster_cols[c('4', '8', '13')])) +
  labs(title = 'TF motifs increased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/mard_tf_up_venn.pdf'), width = 7, height = 6)

ggvenn(list('cluster 4' = neg_df$C4, 'cluster 8' = neg_df$C8, 
            'cluster 13' = neg_df$C13), 
       show_percentage = FALSE, text_size = 8, fill_alpha = 0.65,  
       fill_color = unname(annot_colors$seurat_cluster_cols[c('4', '8', '13')])) +
  labs(title = 'TF motifs decreased in ly/my/ery clusters') +
  theme(plot.title = element_text(hjust = 0.5, size = 22, vjust = 1))
ggsave(filename = paste0(fig_path, 'additional_plots/mard_tf_dn_venn.pdf'), width = 7, height = 6)


res3 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C3))
res4 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C4))
res6 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C6))
res8 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C8))
res13 <- prep_data(dataset = 'human', genes = na.omit(pos_df$C13))

orca(p = plot_cellradar(res_list = list(res4, res8), colors = unname(annot_colors$seurat_cluster_cols[c('4', '8')]), names = c('cluster 4', 'cluster 8')),
     paste0(fig_path, 'additional_plots/mard_c4_8.pdf'), more_args = '--disable-gpu')

orca(p = plot_cellradar(res_list = list(res6, res3, res13), colors = unname(annot_colors$seurat_cluster_cols[c('6', '3', '13')]), names = c('cluster 6', 'cluster 3', 'cluster 13')),
     paste0(fig_path, 'additional_plots/mard_CR_c6_3_13.pdf'), more_args = '--disable-gpu')

missinggenes <- data.frame('C3' = str_flatten(res3$nonexisting_genes, collapse = ','),
                           'C4' = str_flatten(res4$nonexisting_genes, collapse = ','),
                           'C6' = str_flatten(res6$nonexisting_genes, collapse = ','),
                           'C8' = str_flatten(res8$nonexisting_genes, collapse = ','),
                           'C13' = str_flatten(res13$nonexisting_genes, collapse = ','))

write.table(missinggenes, paste0(sup_table_path, '/mard_CR_missinggenes.txt'))

table(na.omit(pos_df$C4) %in% res4$nonexisting_genes)
table(na.omit(pos_df$C8) %in% res8$nonexisting_genes)

# make PC plot of Ly/GMP 
clusters_oi <- c('Ly-I', 'Ly-II', 'GMP-II')

pca_cols <- c('Ly-I' = '#BC89CC', 'Ly-II' = '#63226E', 'GMP-II' = '#2874A6')

cvsubs <- subset(cvejic, pred.id.BM %in% clusters_oi)

cvvars <- apply(cvsubs[['cisbpchromvar']]@data, 1, var)
hist(cvvars, breaks = 500)
hist(cvvars, breaks = 500, xlim = c(0,3))
table(cvvars > 2.5)

cvvars_ord <- cvvars[order(cvvars, decreasing = T)]

cvsubs <- ScaleData(cvsubs, assay = 'cisbpchromvar', features = names(cvvars_ord[1:150]))

cvsubsmotifs <- cvsubs[['cisbpchromvar']]@scale.data
dim(cvsubsmotifs)
cvsubsmotifs <- data.frame(t(cvsubsmotifs))
cvsubsmotifs$cluster <- cvsubs$pred.id.BM
cvsubsmotifs$cluster <- factor(cvsubsmotifs$cluster)

pcacvsubs <- prcomp(cvsubsmotifs[, 1:(dim(cvsubsmotifs)[2]-1)])

autoplot(object = pcacvsubs, data = cvsubsmotifs, color = 'cluster', size = 1.1) +
  scale_color_manual(values = pca_cols) + 
  theme_bw()
ggsave(paste0(fig_path, 'fig4/mard_lymy_top150motifs_pca.pdf'), width = 7.61, height = 4.15)

# Flip x axis to match other plots
pcacvsubs$x[,1] <- - pcacvsubs$x[,1]
pcacvsubs$rotation[,1] <- -pcacvsubs$rotation[,1]

autoplot(object = pcacvsubs, data = cvsubsmotifs, color = 'cluster', size = 1.1) +
  scale_color_manual(values = pca_cols) + 
  theme_bw()
ggsave(paste0(fig_path, 'fig4/mard_lymy_top150motifs_pca_xflip.pdf'), width = 7.61, height = 4.15)


