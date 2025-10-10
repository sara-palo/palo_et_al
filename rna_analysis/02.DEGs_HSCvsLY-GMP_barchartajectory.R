library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)

set.seed(123)

PRJ <- "C:/Users/keiki/OneDrive/デスクトップ/Sara_project"
setwd(PRJ)
out_dir <- sprintf("output_%s", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir

FL <- readRDS(file="rds/FL_CITE_trajectory_gmp2.rds")
da_pos_FL <- list()
da_neg_FL <- list()

DefaultAssay(FL)<-"RNA"
Idents(FL) <- 'clust.names'
clusters_oi <- c("Ly-II","GMP")

cite_cols = c('HSC' = '#E9CB53','MPP-II' = '#52BE80', 'Ly-I' = '#A569BD', 'Ly-II' = '#793885', 'DC-I' = '#58558f',
              'Ly-III' = '#4e3666', 'DC-Mono' = '#2874A6','GMP' = '#56ace5','MPP-I' = '#E7A235','MEP' = '#78281F','Cyc' = '#533511')

FL <- SetIdent(FL, value = factor(Idents(FL), 
                                  levels = c("HSC", "MPP-I", "MPP-II", "Ly-I", "Ly-II", 
                                             "Ly-III","GMP","DC-I", "DC-Mono",  "MEP", "Cyc")))
pdf(file=paste0(out_dir,"/FL_UMAP_Dimplot.pdf"),width=8,height = 4)
DimPlot(FL, cols = cite_cols,pt.size=1,reduction="umap")+NoLegend()+coord_fixed(ratio=1)
dev.off()

FL_df_up <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'up')
FL_df_dn <- data.frame('count' = rep(NA, length(clusters_oi)), 'cluster' = clusters_oi, 'up_or_dn' = 'down')

for (i in 1:length(clusters_oi)){
  ident <- clusters_oi[i]
  de_genes <- FindMarkers(FL,ident.1 = ident,ident.2 = "HSC",only.pos = F)
  da_pos_FL[[ident]] <- subset(de_genes, subset = avg_log2FC > 1 & p_val_adj < 0.05)
  da_neg_FL[[ident]] <- subset(de_genes, subset = avg_log2FC < -1 & p_val_adj < 0.05)
  write.csv(da_pos_FL$`Ly-II`,file = paste0(out_dir,"/FLsub_LyII_df_up_log2_1.csv")) 
  write.csv(da_neg_FL$`Ly-II`,file = paste0(out_dir,"/FLsub_LyII_df_dn_log2_1.csv") )
  write.csv(da_pos_FL$`GMP`,file = paste0(out_dir,"/FLsub_GMP_df_up_log2_1.csv") )
  write.csv(da_neg_FL$`GMP`,file = paste0(out_dir,"/FLsub_GMP_df_dn_log2_1.csv"))
  FL_df_up[i,'count'] <- dim(da_pos_FL[[ident]])[1]
  FL_df_up[i,'cluster'] <- ident
  FL_df_dn[i,'count'] <- dim(da_neg_FL[[ident]])[1]
  FL_df_dn[i,'cluster'] <- ident
}

df_de <- rbind(FL_df_up, FL_df_dn)
df_de$cluster <- factor(df_de$cluster, levels = clusters_oi)
df_de$up_or_dn <- factor(df_de$up_or_dn, levels = c('up', 'down'))
pdf(paste0(out_dir,"/FL_subset_harmony_genes_l2fc_1.pdf"))
ggplot(df_de, aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.7, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = c("darkred","darkblue")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, vjust = 1), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.8, vjust = 0.8), 
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 85)) +
  ggtitle('FL_subset_harmony;up- or down to HSC, log2fc=1.0')
dev.off()


pdf(paste0(out_dir,"/FL_subset_harmony_genes_l2fc_1 with number_v2.pdf"))
ggplot(df_de, aes(x = cluster, y = count, fill = up_or_dn)) + geom_bar(stat = 'identity', width = 0.7, position = 'dodge') + 
  scale_fill_manual('up_or_dn', values = c("darkred","darkblue")) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7), vjust = -0.5, size = 4) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, vjust = 1), 
        axis.text.x = element_text(size = 12, angle = 45, hjust = 0.8, vjust = 0.8), 
        axis.text.y = element_text(size = 12)) +
  scale_y_continuous(limits = c(0, 85)) +
  ggtitle('FL_subset_harmony;up- or down to HSC, log2fc=1.0')
dev.off()





