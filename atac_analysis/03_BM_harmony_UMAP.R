source('scripts/atac_analysis/00.R')

library(harmony)

BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

BM <- RunTFIDF(BM, assay = 'peaks')
BM <- FindTopFeatures(BM, min.cutoff = 'q0', assay = 'peaks')
BM <- RunSVD(BM, assay = 'peaks') 

DepthCor(BM)

# run Harmony again using the new peaks
BM <-
  RunHarmony(
    BM,
    group.by.vars = "sample_number",
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F,
    reduction.save = 'harmony'
  )

BM <- FindNeighbors(object = BM, reduction = 'harmony', dims = 2:20)
BM <- FindClusters(object = BM, verbose = FALSE, algorithm = 3, resolution = 0.7)
BM <- RunUMAP(BM, reduction = "harmony", dims = 2:20, return.model = TRUE, min.dist = 0.2,
              spread = 1, n.neighbors = 30)
DimPlot(BM, label = TRUE, pt.size = 0.85) + NoLegend() + coord_fixed()

BM[['umap']]@cell.embeddings[,'UMAP_2'] <- - (BM[['umap']]@cell.embeddings[,'UMAP_2'])

table(BM$seurat_clusters)

DimPlot(BM, label = TRUE, pt.size = 0.85) + coord_fixed() + NoLegend()

BM$seurat_clusters <- factor(BM$seurat_clusters, levels = 0:11)
BM$BM_annot <- factor(BM$BM_annot, levels = clust_order)
Idents(BM) <- 'BM_annot'

saveRDS(BM, paste0(rds_path, 'BM_only.rds'))

DimPlot(BM, label = TRUE, pt.size = 0.4, group.by = 'sample_number', label.size = 0, shuffle = T, cols = annot_colors$sample_numberBM) + coord_fixed()
ggsave(paste0(fig_path, 'sup1/BM_umap_harmony.pdf'), width = 6, height = 4)

DimPlot(BM, label = TRUE, pt.size = 0.6, group.by = 'BM_annot', cols = annot_colors$cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'sup1/BM_umap.pdf'), width = 6, height = 4)

DimPlot(BM, label = TRUE, pt.size = 0.6, group.by = 'seurat_clusters') + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/BM_umap.pdf'), width = 6, height = 4)

VlnPlot(BM, features = 'TSS.enrichment', group.by = 'sample_number', cols = annot_colors$sample_numberBM, pt.size = 0.02) + 
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'additional_plots/BM_TSS_vln.pdf'), width = 3)

VlnPlot(BM, features = 'nCount_peaks', group.by = 'BM_annot', cols = annot_colors$cluster_cols, pt.size = 0.02) + 
  NoLegend() +
  geom_hline(yintercept = median(BM@meta.data[WhichCells(BM, idents = 'HSC'), 'nCount_peaks']), color = "black")
ggsave(paste0(fig_path, 'sup1/BM_peak_vln.pdf'), width = 7.5, height = 4)

sup1 <- data.frame(matrix(table(BM$seurat_clusters, BM$sample_number), ncol = 2))
colnames(sup1) <- c('BM1', 'BM2')
rownames(sup1) <- 0:11
write.csv(sup1, paste0(sup_table_path, '/BM_clust_dist.csv'))

sup1 <- data.frame(matrix(table(BM$BM_annot, BM$sample_number), ncol = 2))
colnames(sup1) <- c('BM1', 'BM2')
rownames(sup1) <- clust_order
write.csv(sup1, paste0(sup_table_path, '/BM_annot_dist.csv'))

# statistics for peak counts per cluster
df <- BM@meta.data[, c('nCount_peaks', 'BM_annot')]
colnames(df)
table(df$BM_annot)
df$BM_annot <- factor(df$BM_annot, levels = clust_order)


kruskal.test(df$nCount_peaks ~ df$BM_annot, data = df)

res <- pairwise.wilcox.test(df$nCount_peaks, df$BM_annot, p.adjust.method = "bonferroni")
res

resdf <- data.frame(res$p.value)
write.csv(resdf, paste0(sup_table_path, '/BM_peak_counts_stats.csv'))
