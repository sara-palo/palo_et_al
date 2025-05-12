source('scripts/atac_analysis/00.R')

library(chromVARmotifs)
library(motifmatchr)
library(chromVAR)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)

BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

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

BM <- AddMotifs(object = BM,
                genome = BSgenome.Hsapiens.UCSC.hg38,
                pfm = cisbp,
                assay = 'peaks')


BM <- RunChromVAR(object = BM, genome = BSgenome.Hsapiens.UCSC.hg38, assay = 'peaks', new.assay.name = 'cisbpchromvar')

saveRDS(BM, paste0(rds_path, 'BM_only.rds'))

Idents(BM) <- 'seurat_clusters'
levels(BM) <- rev(c(0:11))

DefaultAssay(BM) <- 'cisbpchromvar'

motifsBM_chromv <- data.frame()
motifsBM_chromv_list <- list()

table(BM$seurat_clusters)

for(i in levels(BM)) {
  
  da_motifs <- FindMarkers(
    object = BM,
    assay = 'cisbpchromvar',
    ident.1 = i,
    only.pos = TRUE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff'
  )
  
  write.csv(da_motifs, paste0(csv_path_ch, '/vs_all/BM_', i, '_vsall.csv'))
  
  if (min(da_motifs$p_val_adj) < 0.05 &&
      max(da_motifs$avg_diff) > 1.5) {
    top_motifs <- da_motifs[da_motifs$p_val_adj < 0.05, ]
    top_motifs <- da_motifs[da_motifs$avg_diff > 1.5, ]
    
    motifsBM_chromv_list[[i]] <- top_motifs
    
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
    
    motifsBM_chromv <- rbind(motifsBM_chromv, df)
  }
}

DotPlot(BM, features = rev(unique(motifsBM_chromv$motif)), cols = annot_colors$expr_cols_2) + RotatedAxis() + 
  theme(panel.background = element_rect(fill = 'white')) + 
  plot_annotation(title = 'DA TF motifs', theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(paste0(fig_path, 'additional_plots/BM_chromvar_clust_dotplot.pdf'), width = 15, height = 7.5)

Idents(BM) <- 'BM_annot'
levels(BM) <- rev(clust_order)

DefaultAssay(BM) <- 'cisbpchromvar'

motifsBM_chromv <- data.frame()
motifsBM_chromv_list <- list()

for(i in levels(BM)) {
  
  da_motifs <- FindMarkers(
    object = BM,
    assay = 'cisbpchromvar',
    ident.1 = i,
    only.pos = TRUE,
    mean.fxn = rowMeans,
    fc.name = 'avg_diff'
  )
  
  if(i == 'Ba/M'){
    write.csv(da_motifs, paste0(csv_path_ch, '/vs_all/BM_', 'BaM', '_vsall.csv'))
  } else {
    write.csv(da_motifs, paste0(csv_path_ch, '/vs_all/BM_', i, '_vsall.csv'))
  }  
  
  if (min(da_motifs$p_val_adj) < 0.05 &&
      max(da_motifs$avg_diff) > 1.5) {
    top_motifs <- da_motifs[da_motifs$p_val_adj < 0.05, ]
    top_motifs <- da_motifs[da_motifs$avg_diff > 1.5, ]
    
    motifsBM_chromv_list[[i]] <- top_motifs
    
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
    
    motifsBM_chromv <- rbind(motifsBM_chromv, df)
  }
}

DotPlot(BM, features = rev(unique(motifsBM_chromv$motif)), cols = annot_colors$expr_cols_2) + RotatedAxis() + 
  theme(panel.background = element_rect(fill = 'white')) + 
  plot_annotation(title = 'DA TF motifs', theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
ggsave(paste0(fig_path, 'additional_plots/BM_chromvar_annot_dotplot.pdf'), width = 15, height = 7.5)

BM <- ScaleData(BM, assay = 'cisbpchromvar')

FeaturePlot(BM, features = c('EBF1', 'CEBPA', 'IRF8', 'GATA1', 'HOXB7', 'CTCF', 'SPI1', 'MECOM', 'MEF2C', 'ETV6', 'PAX5', 'TCF3'), 
            cols = annot_colors$expr_cols_1, slot = 'scale.data', min.cutoff = 'q50')

FeaturePlot(BM, features = 'EBF1', cols = annot_colors$expr_cols_1, slot = 'scale.data', min.cutoff = 'q50', pt.size = 0.5) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_chromvar_umap_EBF1_q50.pdf'))

FeaturePlot(BM, features = 'CEBPA', cols = annot_colors$expr_cols_1, slot = 'scale.data', min.cutoff = 'q50', pt.size = 0.5) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_chromvar_umap_CEBPA_q50.pdf'))

FeaturePlot(BM, features = 'GATA1', cols = annot_colors$expr_cols_1, slot = 'scale.data', min.cutoff = 'q50', pt.size = 0.5) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_chromvar_umap_GATA1_q50.pdf'))

FeaturePlot(BM, features = 'HOXB8', cols = annot_colors$expr_cols_1, slot = 'scale.data', min.cutoff = 'q50', pt.size = 0.5) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_chromvar_umap_HOXB8_q50.pdf'))

FeaturePlot(BM, features = 'TCF7', cols = annot_colors$expr_cols_1, slot = 'scale.data', min.cutoff = 'q50', pt.size = 0.5) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_chromvar_umap_TCF7_q50.pdf'))

