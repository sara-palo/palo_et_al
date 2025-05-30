source('scripts/atac_analysis/00.R')

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
FL_CITE <- readRDS(paste0(rds_path, 'FL_CITE.rds'))
gene.activities.FL <- readRDS(paste0(rds_path, 'FL_geneact.rds'))
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

levels(FL_CITE) <- clust_order_CITE
DimPlot(FL_CITE, pt.size = 0.6, cols = annot_colors$cite_cols, reduction = 'umap') + coord_fixed()
ggsave(paste0(fig_path, 'sup4/FL_CITE_umap.pdf'))

FL[['RNA']] <- CreateAssayObject(counts = gene.activities.FL)
FL <- NormalizeData(
  object = FL,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(FL$nCount_RNA)
)

# project FL atac onto FL CITE-seq
DefaultAssay(FL) <- 'RNA'
DefaultAssay(FL_CITE) <- 'RNA'

FL_CITE <- FindVariableFeatures(FL_CITE, nfeatures = 2000)
genes.use <- VariableFeatures(FL_CITE)

FL <- ScaleData(FL, assay = 'RNA', features = genes.use)
FL_CITE <- ScaleData(FL_CITE, assay = 'RNA', features = genes.use)

transfer.anchors <- FindTransferAnchors(reference = FL_CITE, query = FL, features = genes.use, reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

FL <- MapQuery(
  anchorset = transfer.anchors,
  reference = FL_CITE,
  query = FL,
  refdata = FL_CITE$clust.names,
  reference.reduction = "cca",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

FL[['ref.umap']]@cell.embeddings[,'refUMAP_2'] <- - FL[['ref.umap']]@cell.embeddings[,'refUMAP_2']

FL$predicted.id <- factor(FL$predicted.id, levels = clust_order_CITE)

DimPlot(FL, reduction = 'ref.umap', group.by = 'predicted.id')

pred.id.CITE <- FL$predicted.id
pred.id.score.CITE <- FL$predicted.id.score
CITEumap <- FL[['ref.umap']]


# map CS22 cite on CS22 atac
FL <- FindVariableFeatures(FL, assay = 'RNA', nfeatures = 2000)
FL <- ScaleData(FL, assay = 'RNA', features = VariableFeatures(FL, assay = 'RNA'))
FL_CITE <- ScaleData(FL_CITE, assay = 'RNA', features = VariableFeatures(FL, assay = 'RNA'))

transfer.anchors <- FindTransferAnchors(reference = FL, query = FL_CITE, features = VariableFeatures(FL, assay = 'RNA'), 
                                        reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

FL_CITE <- MapQuery(
  anchorset = transfer.anchors,
  reference = FL,
  query = FL_CITE,
  refdata = FL$seurat_clusters,
  reference.reduction = "cca",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

FL_CITE[['ref.umap']]@cell.embeddings[,'refUMAP_1'] <- - FL_CITE[['ref.umap']]@cell.embeddings[,'refUMAP_1']

FL_CITE$predicted.id <- factor(FL_CITE$predicted.id, levels = 0:13)

pred.id.CS22 <- FL_CITE$predicted.id
pred.id.score.CS22 <- FL_CITE$predicted.id.score
CS22umap <- FL_CITE[['ref.umap']]

FL_CITE <- readRDS(paste0(rds_path, 'FL_CITE.rds'))

FL_CITE$pred.id.CS22 <- pred.id.CS22
FL_CITE$pred.id.score.CS22 <- pred.id.score.CS22
FL_CITE[['CS22umap']] <- CS22umap

saveRDS(FL_CITE, paste0(rds_path, 'FL_CITE.rds'))

DimPlot(FL_CITE, reduction = 'CS22umap', group.by = 'pred.id.CS22', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed() |
  DimPlot(FL_CITE, reduction = 'umap', group.by = 'pred.id.CS22', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/CITE_proj_onFL.pdf'), width = 8)


# project FL atac onto BM atac
DefaultAssay(FL) <- 'peaks'
DefaultAssay(BM) <- 'peaks'

# create a shared peak set
gr_merge <- reduce(c(granges(FL[['peaks']]), granges(BM[['peaks']])))
gr_merge

# Filter out possible too short/long peaks
peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]
gr_merge

FLcounts <- FeatureMatrix(Fragments(FL), 
                          features = gr_merge, 
                          cells = Cells(FL))

BMcounts <- FeatureMatrix(Fragments(BM), 
                          features = gr_merge, 
                          cells = Cells(BM))

FL[['mergedpeaks']] <- CreateChromatinAssay(counts = FLcounts, fragments = Fragments(FL))
BM[['mergedpeaks']] <- CreateChromatinAssay(counts = BMcounts, fragments = Fragments(BM))

DefaultAssay(FL) <- 'mergedpeaks'
DefaultAssay(BM) <- 'mergedpeaks'

FL <- RunTFIDF(FL, assay = 'mergedpeaks')
FL <- FindTopFeatures(FL, assay = 'mergedpeaks')
FL <- RunSVD(FL, assay = 'mergedpeaks')

BM <- RunTFIDF(BM, assay = 'mergedpeaks')
BM <- FindTopFeatures(BM, assay = 'mergedpeaks')
BM <- RunSVD(BM, assay = 'mergedpeaks')

#project FL onto BM
transfer.anchors <- FindTransferAnchors(
  reference = BM,
  reference.assay = 'mergedpeaks',
  query = FL,
  query.assay = 'mergedpeaks',
  reference.reduction = 'lsi',
  reduction = 'lsiproject',
  dims = 2:30
)

# map FL cells onto the BM UMAP
FL <- MapQuery(
  anchorset = transfer.anchors,
  reference = BM,
  query = FL,
  refdata = BM$BM_annot,
  reference.reduction = 'lsi',
  new.reduction.name = 'ref.lsi',
  reduction.model = 'umap'
)

FL[['ref.umap']]@cell.embeddings[,'refUMAP_2'] <- - FL[['ref.umap']]@cell.embeddings[,'refUMAP_2'] 

FL$predicted.id <- factor(FL$predicted.id, levels = clust_order)

preds <- FL$predicted.id
predscores <- FL$predicted.id.score
refumap <- FL[['ref.umap']]



#project BM onto FL
transfer.anchors <- FindTransferAnchors(
  reference = FL,
  reference.assay = 'mergedpeaks',
  query = BM,
  query.assay = 'mergedpeaks',
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

# map BM cells onto the FL UMAP
BM <- MapQuery(
  anchorset = transfer.anchors,
  reference = FL,
  query = BM,
  refdata = FL$seurat_clusters,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

BM[['ref.umap']]@cell.embeddings[,'refUMAP_1'] <- - BM[['ref.umap']]@cell.embeddings[,'refUMAP_1']

DimPlot(BM, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, pt.size = 0.8, cols = annot_colors$seurat_cluster_cols) + 
  coord_fixed() + NoLegend() + 
  labs(title = 'Bone marrow projected onto FL') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) 
ggsave(paste0(fig_path, 'additional_plots/BM_projected_on_FL.pdf'), width = 6.2, height = 4.65)

DimPlot(BM, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE, pt.size = 0.8, cols = annot_colors$seurat_cluster_cols) + 
  coord_fixed() + NoLegend() + 
  labs(title = 'Bone marrow - predicted IDs') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) 
ggsave(paste0(fig_path, 'additional_plots/BM_umap_predictedIDs.pdf'), width = 6.2, height = 4.65)


DimPlot(FL, reduction = 'ref.umap')

# reload the saved FL object and save only the predicted IDs, prediction scores and UMAPs
FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

FL$pred.id.BM <- preds
FL$pred.id.score.BM <- predscores
FL[['BMumap']] <- refumap

FL$pred.id.CITE <- pred.id.CITE
FL$pred.id.score.CITE <- pred.id.score.CITE
FL[['CITEumap']] <- CITEumap

DimPlot(FL, reduction = 'BMumap', group.by = 'pred.id.BM') | DimPlot(FL, reduction = 'CITEumap', group.by = 'pred.id.CITE')

saveRDS(FL, paste0(rds_path, 'FL_only.rds'))

# make plots
DimPlot(FL, reduction = 'CITEumap', group.by = 'pred.id.CITE', pt.size = 0.6, cols = annot_colors$cite_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/FL_proj_CITE.pdf'))

DimPlot(FL, reduction = 'CITEumap', group.by = 'seurat_clusters', pt.size = 0.5, cols = annot_colors$seurat_cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/FL_proj_CITE_seuratclust.pdf'))

DimPlot(FL, reduction = 'umap', group.by = 'pred.id.CITE', pt.size = 0.5, cols = annot_colors$cite_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/FL_predIDs_CITE.pdf'))

FeaturePlot(FL, reduction = 'umap', 'pred.id.score.CITE', cols = annot_colors$expr_cols_1, pt.size = 0.5) + coord_fixed()
ggsave(paste0(fig_path, 'sup3/FL_cite_predicted_ID_scores.pdf'))

VlnPlot(FL, 'pred.id.score.CITE', group.by = 'pred.id.CITE', cols = annot_colors$cite_cols, pt.size = 0.01) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/FL_pred_scores_vln_CITE.pdf'))

FL_CITE$clust.names <- factor(FL_CITE$clust.names, levels = clust_order_CITE)

predictions_clusters_CITE <- table(FL$seurat_clusters, FL$pred.id.CITE)
predictions_clusters_CITE <- predictions_clusters_CITE/rowSums(predictions_clusters_CITE)  # divide by number of cells for each cell type
predictions_clusters_CITE <- data.frame(predictions_clusters_CITE)

predictions_clusters_CITE$Var1 <- factor(predictions_clusters_CITE$Var1, levels = rev(c(0:13)))
predictions_clusters_CITE$Var2 <- factor(predictions_clusters_CITE$Var2, levels = clust_order_CITE)

ggplot(predictions_clusters_CITE, aes(x = Var2, y = Var1, fill = Freq)) + geom_tile() + 
  scale_fill_gradient(name = "Fraction of cells", low = "white", high = "navy") +  
  xlab("Predicted cell type - label transfer from CITE") +
  ylab("Seurat clusters") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(fig_path, 'additional_plots/FL_CITE_predictions_tile.pdf'))

predictions_clusters_CITE$Var1 <- factor(predictions_clusters_CITE$Var1, levels = 0:13)

ggplot(predictions_clusters_CITE, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = annot_colors$cite_cols, name = 'Predicted cell type') +
  xlab("Predicted cell type - label transfer from CITE") +
  theme_bw() +
  ylab("Frequency")
ggsave(paste0(fig_path, 'sup3/FL_CITE_predictions_bar.pdf'))


cite_df <- data.frame(table(FL_CITE$clust.names))
cite_df$cluster <- factor(cite_df$Var1, levels = clust_order_CITE)

ggplot(cite_df, aes(x = "", y = Freq, fill = cluster)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0, direction = -1) + scale_fill_manual(values = annot_colors$cite_cols) +
  theme_void() 
ggsave(paste0(fig_path, 'additional_plots/CITE_distr_pie.pdf'))


pred_df <- data.frame(table(FL$pred.id.CITE))
pred_df$cluster <- factor(pred_df$Var1, levels = clust_order_CITE)

ggplot(pred_df, aes(x = "", y = Freq, fill = cluster)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0, direction = -1) + scale_fill_manual(values = annot_colors$cite_cols) +
  theme_void() 
ggsave(paste0(fig_path, 'additional_plots/CITE_predictions_distr_pie.pdf'))



DimPlot(FL, reduction = 'BMumap', group.by = 'pred.id.BM', label = TRUE, repel = TRUE, pt.size = 0.7, 
        cols = annot_colors$cluster_cols) + coord_fixed() + 
  labs(title = 'Fetal liver projected onto BM') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) 
ggsave(paste0(fig_path, 'sup1/FL_projected_on_BM.pdf'), width = 6.2, height = 4.65)

FL$seurat_clusters <- factor(FL$seurat_clusters, levels = 0:13)

DimPlot(FL, reduction = 'umap', group.by = 'pred.id.BM', label = TRUE, repel = TRUE, pt.size = 0.5, cols = annot_colors$cluster_cols) + 
  coord_fixed() + 
  labs(title = 'Fetal liver - predicted IDs') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12))
ggsave(paste0(fig_path, 'fig1/FL_umap_predictedIDs_BM.pdf'), width = 6.2, height = 4.65)

FeaturePlot(FL, reduction = 'umap', 'pred.id.score.BM', cols = annot_colors$expr_cols_1, pt.size = 0.5) + 
  coord_fixed()
ggsave(paste0(fig_path, 'sup1/FL_umap_predscores_BM.pdf'), width = 6.2, height = 4.65)

FeaturePlot(FL, reduction = 'BMumap', 'pred.id.score.BM', cols = annot_colors$expr_cols_1, pt.size = 0.5) + 
  coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/FL_BMumap_predscores_BM.pdf'), width = 6.2, height = 4.65)

DimPlot(FL, reduction = 'BMumap', group.by = 'seurat_clusters', label = F, repel = TRUE, pt.size = 0.5, cols = annot_colors$seurat_cluster_cols) + 
  coord_fixed() + 
  labs(title = 'Fetal liver on BM UMAP') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12))
ggsave(paste0(fig_path, 'additional_plots/FL_BMrefumap.pdf'), width = 6.2, height = 4.65)


VlnPlot(FL, 'pred.id.score.BM', group.by = 'pred.id.BM', cols = annot_colors$cluster_cols, pt.size = 0.01) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/FL_pred_scores_vln_BM.pdf'))

predictions_clusters_BM <- table(FL$seurat_clusters, FL$pred.id.BM)
predictions_clusters_BM <- predictions_clusters_BM/rowSums(predictions_clusters_BM)  # divide by number of cells for each cell type
predictions_clusters_BM <- data.frame(predictions_clusters_BM)

predictions_clusters_BM$Var1 <- factor(predictions_clusters_BM$Var1, levels = rev(c(0:13)))
predictions_clusters_BM$Var2 <- factor(predictions_clusters_BM$Var2, levels = clust_order)

ggplot(predictions_clusters_BM, aes(x = Var2, y = Var1, fill = Freq)) + geom_tile() + 
  scale_fill_gradient(name = "Fraction of cells", low = "white", high = "navy") +  
  xlab("Predicted cell type - label transfer from BM") +
  ylab("Seurat clusters") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(fig_path, 'additional_plots/FL_predictions_tile_BM.pdf'))

predictions_clusters_BM$Var1 <- factor(predictions_clusters_BM$Var1, levels = 0:13)

ggplot(predictions_clusters_BM, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = annot_colors$cluster_cols, name = 'Predicted cell type') +
  xlab("Predicted cell type - label transfer from BM") +
  theme_bw() +
  ylab("Frequency")
ggsave(paste0(fig_path, 'fig1/FL_predictions_bar_BM.pdf'))


predictionsFL <- data.frame(table(FL$pred.id.BM, FL$sample_number))
colnames(predictionsFL) <- c('cell_type', 'sample_number', 'count')

#calculate which percentage of the sample each cluster/predicted cell type contains
predictionsFL$pct <- ifelse(predictionsFL$sample_number == 'FL1', (predictionsFL$count/(sum(FL$sample_number == 'FL1'))) * 100, 
                            (predictionsFL$count/(sum(FL$sample_number == 'FL2'))) * 100)

predictionsFL$sample_type <- 'Fetal liver'

clustersBM <- data.frame(table(BM$BM_annot, BM$sample_number))
colnames(clustersBM) <- c('cell_type', 'sample_number', 'count')

clustersBM$pct <- ifelse(clustersBM$sample_number == 'BM1', (clustersBM$count/(sum(BM$sample_number == 'BM1'))) * 100, 
                         (clustersBM$count/(sum(BM$sample_number == 'BM2'))) * 100)

clustersBM$sample_type <- 'Adult bone marrow'

plot_data1 <- rbind(predictionsFL, clustersBM)

plot_data1$sample_number <- factor(plot_data1$sample_number, levels = c('FL1', 'FL2', 'BM1', 'BM2'))
plot_data1$sample_type <- factor(plot_data1$sample_type, levels =  c('Fetal liver', 'Adult bone marrow'))

plot_data1$cell_type <- factor(plot_data1$cell_type, levels = c(clust_order))

ggplot(data = plot_data1, aes(x = cell_type, y = pct)) + theme_bw() +
  geom_bar(stat = "identity", aes(fill = sample_number), width = 0.7, position = position_dodge()) + 
  scale_fill_manual('sample number', values = c(annot_colors$sample_numberFL, annot_colors$sample_numberBM)) + 
  scale_y_continuous(limits = c(0, 25)) +
  labs(title = 'Cell type distribution - FL projected on BM') + xlab('Cell type') + ylab('% of total') +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.7))
ggsave(filename = paste0(fig_path, 'sup1/labeltransf_clust_distr.pdf'), width = 8, height = 3.2)




predictions_clusters_FL <- table(BM$BM_annot, BM$predicted.id)
predictions_clusters_FL <- predictions_clusters_FL/rowSums(predictions_clusters_FL)  # divide by number of cells for each cell type
predictions_clusters_FL <- data.frame(predictions_clusters_FL)

predictions_clusters_FL$Var1 <- factor(predictions_clusters_FL$Var1, levels = rev(clust_order))
predictions_clusters_FL$Var2 <- factor(predictions_clusters_FL$Var2, levels = 0:13)

ggplot(predictions_clusters_FL, aes(x = Var2, y = Var1, fill = Freq)) + geom_tile() + 
  scale_fill_gradient(name = "Fraction of cells", low = "white", high = "navy") +  
  xlab("Predicted cell type - label transfer from FL") +
  ylab("BM clusters") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0(fig_path, 'additional_plots/BM_predictions_tile_FL.pdf'))

predictions_clusters_FL$Var1 <- factor(predictions_clusters_FL$Var1, levels = clust_order)

ggplot(predictions_clusters_FL, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = annot_colors$seurat_cluster_cols, name = 'Predicted cell type') +
  xlab("Predicted cell type - label transfer from FL") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Frequency")
ggsave(paste0(fig_path, 'additional_plots/BM_predictions_bar_FL.pdf'))


predictionsBM <- data.frame(table(BM$predicted.id, BM$sample_number))
colnames(predictionsBM) <- c('cell_type', 'sample_number', 'count')

#calculate which percentage of the sample each cluster/predicted cell type contains
predictionsBM$pct <- ifelse(predictionsBM$sample_number == 'BM1', (predictionsBM$count/(sum(BM$sample_number == 'BM1'))) * 100, 
                            (predictionsBM$count/(sum(BM$sample_number == 'BM2'))) * 100)

predictionsBM$sample_type <- 'Adult bone marrow'

clustersFL <- data.frame(table(FL$seurat_clusters, FL$sample_number))
colnames(clustersFL) <- c('cell_type', 'sample_number', 'count')

clustersFL$pct <- ifelse(clustersFL$sample_number == 'FL1', (clustersFL$count/(sum(FL$sample_number == 'FL1'))) * 100, 
                         (clustersFL$count/(sum(FL$sample_number == 'FL2'))) * 100)

clustersFL$sample_type <- 'Fetal liver'

plot_data1 <- rbind(clustersFL, predictionsBM)

plot_data1$sample_number <- factor(plot_data1$sample_number, levels = c('FL1', 'FL2', 'BM1', 'BM2'))
plot_data1$sample_type <- factor(plot_data1$sample_type, levels =  c('Fetal liver', 'Adult bone marrow'))

plot_data1$cell_type <- factor(plot_data1$cell_type, levels = 0:13)

ggplot(data = plot_data1, aes(x = cell_type, y = pct)) + theme_bw() +
  geom_bar(stat = "identity", aes(fill = sample_number), width = 0.7, position = position_dodge()) + 
  scale_fill_manual('sample number', values = c(annot_colors$sample_numberFL, annot_colors$sample_numberBM)) + 
  scale_y_continuous(limits = c(0, 25)) +
  labs(title = 'Cell type distribution - BM projected on FL') + xlab('Cell type') + ylab('% of total') +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12))
ggsave(filename = paste0(fig_path, 'additional_plots/labeltransf_fromFL_clust_distr.pdf'), width = 8, height = 3.2)


