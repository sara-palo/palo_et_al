source('scripts/atac_analysis/00.R')

ranz <- readRDS(paste0(rds_path, 'ranzoni.rds'))
FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

DefaultAssay(FL) <- 'peaks'
DefaultAssay(BM) <- 'peaks'

gr_merge <- reduce(c(granges(FL[['peaks']]), granges(ranz[['peaks']])))

# Filter out possible too short/long peaks
peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]

FLcounts <- FeatureMatrix(Fragments(FL),
                          features = gr_merge,
                          cells = Cells(FL))

ccounts <- FeatureMatrix(Fragments(ranz), 
                         features = gr_merge, 
                         cells = Cells(ranz))

FL[['peaks']] <- CreateChromatinAssay(counts = FLcounts, fragments = Fragments(FL))
ranz[['peaks']] <- CreateChromatinAssay(counts = ccounts, fragments = Fragments(ranz))

DefaultAssay(FL) <- 'peaks'
DefaultAssay(ranz) <- 'peaks'

FL <- RunTFIDF(FL, assay = 'peaks')
FL <- FindTopFeatures(FL)
FL <- RunSVD(FL)

ranz <- RunTFIDF(ranz, assay = 'peaks')
ranz <- FindTopFeatures(ranz)
ranz <- RunSVD(ranz)

transfer.anchors <- FindTransferAnchors(
  reference = FL,
  reference.assay = 'peaks',
  query = ranz,
  query.assay = 'peaks',
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:10 
)

ranz <- MapQuery(
  anchorset = transfer.anchors,
  reference = FL,
  query = ranz,
  refdata = FL$seurat_clusters,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)


ranz[['ref.umap']]@cell.embeddings[,1] <- - ranz[['ref.umap']]@cell.embeddings[,1]

ranz$predicted.id <- factor(ranz$predicted.id, levels = 0:13)

DimPlot(ranz, reduction = 'ref.umap', group.by = 'predicted.id', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed() 

pred.id.CS22 <- ranz$predicted.id
pred.id.score.CS22 <- ranz$predicted.id.score
CS22umap <- ranz[['ref.umap']]

ranz <- readRDS(paste0(rds_path, 'ranzoni.rds'))
ranz$pred.id.CS22 <- pred.id.CS22
ranz$pred.id.score.CS22 <- pred.id.score.CS22
ranz[['CS22umap']] <- CS22umap
saveRDS(ranz, paste0(rds_path, 'ranzoni.rds'))


VlnPlot(ranz, 'pred.id.score.CS22', group.by = 'seurat_clusters', pt.size = 0.01, cols = annot_colors$seurat_cluster_cols) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/ranz_proj_on_CS22_scores.pdf'), width = 8) 

DimPlot(ranz, reduction = 'CS22umap', group.by = 'pred.id.CS22', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed() 
ggsave(paste0(fig_path, 'additional_plots/ranz_proj_onCS22umap.pdf'))


VlnPlot(ranz, 'nCount_peaks', group.by = 'pred.id.CS22', cols = annot_colors$seurat_cluster_cols) + NoLegend() 
ggsave(paste0(fig_path, 'additional_plots/ranz_peak_counts_predIDs_CS22.pdf'), width = 8)



# now project to BM
gr_merge <- reduce(c(granges(BM[['peaks']]), granges(ranz[['peaks']])))
gr_merge

peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]
gr_merge

BMcounts <- FeatureMatrix(Fragments(BM),
                          features = gr_merge,
                          cells = Cells(BM))

ccounts <- FeatureMatrix(Fragments(ranz), 
                         features = gr_merge, 
                         cells = Cells(ranz))

BM[['peaks']] <- CreateChromatinAssay(counts = BMcounts, fragments = Fragments(BM))
ranz[['peaks']] <- CreateChromatinAssay(counts = ccounts, fragments = Fragments(ranz))

DefaultAssay(BM) <- 'peaks'
DefaultAssay(ranz) <- 'peaks'

BM <- RunTFIDF(BM, assay = 'peaks')
BM <- FindTopFeatures(BM)
BM <- RunSVD(BM)

ranz <- RunTFIDF(ranz, assay = 'peaks')
ranz <- FindTopFeatures(ranz) 
ranz <- RunSVD(ranz)

transfer.anchors <- FindTransferAnchors(
  reference = BM,
  reference.assay = 'peaks',
  query = ranz,
  query.assay = 'peaks',
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:10 
)

ranz <- MapQuery(
  anchorset = transfer.anchors,
  reference = BM,
  query = ranz,
  refdata = BM$BM_annot,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

ranz[['ref.umap']]@cell.embeddings[,2] <- - ranz[['ref.umap']]@cell.embeddings[,2]
ranz$predicted.id <- factor(ranz$predicted.id, levels = clust_order)

DimPlot(ranz, reduction = 'ref.umap', group.by = 'predicted.id', pt.size = 0.6, cols = annot_colors$cluster_cols) + coord_fixed() 

pred.id.BM <- ranz$predicted.id
pred.id.score.BM <- ranz$predicted.id.score
BMumap <- ranz[['ref.umap']]

ranz <- readRDS(paste0(rds_path, 'ranzoni.rds'))
ranz$pred.id.BM <- pred.id.BM
ranz$pred.id.score.BM <- pred.id.score.BM
ranz[['BMumap']] <- BMumap

saveRDS(ranz, paste0(rds_path, 'ranzoni.rds'))

VlnPlot(ranz, 'pred.id.score.BM', group.by = 'paper_cluster', pt.size = 0.01) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/ranz_proj_on_BM_scores_paperclust.pdf'), width = 6) 

VlnPlot(ranz, 'pred.id.score.BM', group.by = 'pred.id.BM', pt.size = 0.01, cols = annot_colors$cluster_cols) + NoLegend() + ylim(0, 1.05)
ggsave(paste0(fig_path, 'sup3/ranz_proj_on_BM_scores_predIDs_BM.pdf'), width = 8) 

DimPlot(ranz, reduction = 'paper_UMAP', group.by = 'pred.id.BM', pt.size = 0.6, cols = annot_colors$cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/ranz_BM_predIDs.pdf'))

VlnPlot(ranz, 'pred.id.score.BM', group.by = 'seurat_clusters', pt.size = 0.01, cols = annot_colors$seurat_cluster_cols) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/ranz_proj_on_BM_scores.pdf'), width = 8) 

VlnPlot(ranz, 'nCount_peaks', group.by = 'pred.id.BM', cols = annot_colors$cluster_cols) + NoLegend() 
ggsave(paste0(fig_path, 'additional_plots/ranz_peak_counts_predIDs_BM.pdf'), width = 8)


