source('scripts/atac_analysis/00.R')

library(harmony)

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

FL <- RunTFIDF(FL, assay = 'peaks')
FL <- FindTopFeatures(FL, min.cutoff = 'q0', assay = 'peaks')
FL <- RunSVD(FL, assay = 'peaks') 

DepthCor(FL)

# run Harmony again using the new peaks
FL <-
  RunHarmony(
    FL,
    group.by.vars = "sample_number",
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F,
    reduction.save = 'harmony'
  )

FL <- FindNeighbors(object = FL, reduction = 'harmony', dims = 2:20)
FL <- FindClusters(object = FL, verbose = FALSE, algorithm = 3, resolution = 0.8)
FL <- RunUMAP(FL, reduction = "harmony", dims = 2:20, return.model = TRUE, min.dist = 0.3,
              spread = 1, n.neighbors = 30)

DimPlot(FL, label = TRUE, pt.size = 0.85) + coord_fixed() + NoLegend()

# invert UMAP coordinates
FL[['umap']]@cell.embeddings[,'UMAP_1'] <- - (FL[['umap']]@cell.embeddings[,'UMAP_1'])

table(FL$seurat_clusters)

# renumber clusters
FL <- RenameIdents(FL, 
                   '4' = '0',
                   '1' = '1',
                   '10' = '2',
                   '11' = '3',
                   '7' = '4',
                   '5' = '7',
                   '8' = '8',
                   '6' = '6',
                   '2' = '5',
                   '12' = '9',
                   '0' = '10',
                   '3' = '12',
                   '13' = '11',
                   '9' = '13')

DimPlot(FL, label = TRUE, pt.size = 0.85) + coord_fixed() + NoLegend()

FL$seurat_clusters <- Idents(FL)
FL$seurat_clusters <- factor(FL$seurat_clusters, levels = 0:13)
Idents(FL) <- 'seurat_clusters'

saveRDS(FL, paste0(rds_path, 'FL_only.rds'))

DimPlot(FL, label = TRUE, pt.size = 0.4, group.by = 'sample_number', label.size = 0, shuffle = T, cols = annot_colors$sample_numberFL) + coord_fixed()
ggsave(paste0(fig_path, 'sup1/FL_umap_harmony.pdf'), width = 6, height = 4)

DimPlot(FL, label = TRUE, pt.size = 0.6, group.by = 'seurat_clusters', cols = annot_colors$seurat_cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'fig1/FL_umap.pdf'), width = 6, height = 4)

VlnPlot(FL, features = 'TSS.enrichment', group.by = 'sample_number', cols = annot_colors$sample_numberFL, pt.size = 0.02) + 
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
ggsave(paste0(fig_path, 'additional_plots/FL_TSS_vln.pdf'), width = 3)

VlnPlot(FL, features = 'nCount_peaks', group.by = 'seurat_clusters', cols = annot_colors$seurat_cluster_cols, idents = 0:13, pt.size = 0.02) + 
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0)) +
  geom_hline(yintercept = median(FL@meta.data[WhichCells(FL, idents = 0), 'nCount_peaks']), color = "black")
ggsave(paste0(fig_path, 'sup1/FL_peak_vln.pdf'), width = 7.5, height = 4)

sup1 <- data.frame(matrix(table(FL$seurat_clusters, FL$sample_number), ncol = 2))
colnames(sup1) <- c('FL1', 'FL2')
rownames(sup1) <- 0:13
write.csv(sup1, paste0(sup_table_path, '/tableS1_FL_clust_dist.csv'))

# statistics for peak counts per cluster
df <- FL@meta.data[, c('nCount_peaks', 'seurat_clusters')]
colnames(df)
table(df$seurat_clusters)
df$seurat_clusters <- factor(df$seurat_clusters, levels = 0:13)


kruskal.test(df$nCount_peaks ~ df$seurat_clusters, data = df)

res <- pairwise.wilcox.test(df$nCount_peaks, df$seurat_clusters, p.adjust.method = "bonferroni")
res

resdf <- data.frame(res$p.value)
write.csv(resdf, paste0(sup_table_path, '/FL_peak_counts_stats.csv'))
