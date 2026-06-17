source('scripts/atac_analysis/00.R')

library(harmony)

cvejic <- readRDS(paste0(rds_path, 'cvejic_merge.rds'))
FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

gr_merge <- reduce(c(granges(FL[['peaks']]), granges(cvejic[['peaks']])))
gr_merge

# Filter out possible too short/long peaks
peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]
gr_merge

DefaultAssay(FL) <- 'peaks'

FLcounts <- FeatureMatrix(Fragments(FL),
                          features = gr_merge,
                          cells = Cells(FL))

ccounts <- FeatureMatrix(Fragments(cvejic), 
                         features = gr_merge, 
                         cells = Cells(cvejic))

FL[['peaks']] <- CreateChromatinAssay(counts = FLcounts, fragments = Fragments(FL))
cvejic[['peaks']] <- CreateChromatinAssay(counts = ccounts, fragments = Fragments(cvejic))

DefaultAssay(FL) <- 'peaks'
DefaultAssay(cvejic) <- 'peaks'

FL <- RunTFIDF(FL, assay = 'peaks')
FL <- FindTopFeatures(FL, min.cutoff = 'q75')
FL <- RunSVD(FL)

cvejic <- RunTFIDF(cvejic, assay = 'peaks')
cvejic <- FindTopFeatures(cvejic, min.cutoff = 'q75')
cvejic <- RunSVD(cvejic)

transfer.anchors <- FindTransferAnchors(
  reference = cvejic,
  reference.assay = 'peaks',
  query = FL,
  query.assay = 'peaks',
  reference.reduction = 'lsi',
  reduction = 'lsiproject',
  dims = 2:20 
)

FL <- MapQuery(
  anchorset = transfer.anchors,
  reference = cvejic,
  query = FL,
  refdata = cvejic$seurat_clusters,
  reference.reduction = 'lsi',
  new.reduction.name = 'ref.lsi',
  reduction.model = 'umap'
)

FL$predicted.id <- factor(FL$predicted.id, levels = 0:13)

table(FL$predicted.id)
table(FL$predicted.id, FL$seurat_clusters)

(DimPlot(cvejic, reduction = 'umap', pt.size = 0.7) ) | 
  (DimPlot(FL, reduction = 'ref.umap', group.by = 'predicted.id', pt.size = 0.7)) |
  (DimPlot(FL, reduction = 'ref.umap', group.by = 'seurat_clusters', pt.size = 0.7) )

# save ref umap and predicted IDs and later reload object - don't overwrite the peaks
cvejic.umap <- FL[['ref.umap']]
pred.id.cvejic <- FL$predicted.id
pred.id.score.cvejic <- FL$predicted.id.score

table(cvejic$seurat_clusters)

# subset cvejic object to contain clusters that have a significant proportion of CS22 cells projected onto them
preds <- table(FL$predicted.id)
clust.keep <- names(preds[preds > 100])

# reload cvejic object with the correct peaks
cvejic <- readRDS(paste0(rds_path, 'cvejic_merge.rds'))

cvejic <- subset(cvejic, idents = clust.keep)

DimPlot(cvejic)
table(Idents(cvejic))

table(rowSums(cvejic[['peaks']]@counts)<1)

# remove features with 0 counts
cvejic <- cvejic[names(which(rowSums(GetAssayData(object = cvejic, slot = "counts")) > 0))]

cvejic <- RunTFIDF(cvejic, assay = 'peaks') %>% FindTopFeatures(min.cutoff = 'q75') %>% RunSVD()

cvejic <-
  RunHarmony(
    cvejic,
    group.by.vars = 'sample',
    reduction.use = 'lsi',
    assay.use = 'peaks',
    plot_convergence = TRUE,
    project.dim = F,
    reduction.save = 'harmony'
  )

cvejic <- FindNeighbors(object = cvejic, dims = 2:20, reduction = 'harmony')
cvejic <- FindClusters(object = cvejic, verbose = FALSE, algorithm = 3, resolution = 0.7)
cvejic <- RunUMAP(cvejic, reduction = 'harmony', dims = 2:20, return.model = TRUE, min.dist = 0.3,
                  spread = 1, n.neighbors = 30)

DimPlot(cvejic, label = TRUE, pt.size = 0.6) + coord_fixed() + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset.pdf'), width = 10)

DimPlot(cvejic, label = F, pt.size = 0.4, group.by = 'sample', shuffle = T, cols = annot_colors$sample_number2nd) + coord_fixed() 
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_samp.pdf'), width = 10)

DimPlot(cvejic, label = F, pt.size = 0.6, split.by = 'sample') + coord_fixed() + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_persample.pdf'), width = 10)

saveRDS(cvejic, paste0(rds_path, 'cvejic_subset.rds'))


# now, project the 2nd trimester data onto the CS22 FL
gr_merge <- reduce(c(granges(FL[['peaks']]), granges(cvejic[['peaks']])))

# Filter out possible too short/long peaks
peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]


FLcounts <- FeatureMatrix(Fragments(FL),
                          features = gr_merge,
                          cells = Cells(FL))

ccounts <- FeatureMatrix(Fragments(cvejic), 
                         features = gr_merge, 
                         cells = Cells(cvejic))

FL[['peaks']] <- CreateChromatinAssay(counts = FLcounts, fragments = Fragments(FL))
cvejic[['peaks']] <- CreateChromatinAssay(counts = ccounts, fragments = Fragments(cvejic))

DefaultAssay(FL) <- 'peaks'
DefaultAssay(cvejic) <- 'peaks'

FL <- RunTFIDF(FL, assay = 'peaks')
FL <- FindTopFeatures(FL, min.cutoff = 'q75')
FL <- RunSVD(FL)

cvejic <- RunTFIDF(cvejic, assay = 'peaks')
cvejic <- FindTopFeatures(cvejic, min.cutoff = 'q75')
cvejic <- RunSVD(cvejic)

transfer.anchors <- FindTransferAnchors(
  reference = FL,
  reference.assay = 'peaks',
  query = cvejic,
  query.assay = 'peaks',
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:20
)

cvejic <- MapQuery(
  anchorset = transfer.anchors,
  reference = FL,
  query = cvejic,
  refdata = FL$seurat_clusters,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)


cvejic[['ref.umap']]@cell.embeddings[,'refUMAP_1'] <- - cvejic[['ref.umap']]@cell.embeddings[,'refUMAP_1']

cvejic$predicted.id <- factor(cvejic$predicted.id, levels = 0:13)

pred.id.CS22 <- cvejic$predicted.id
pred.id.score.CS22 <- cvejic$predicted.id.score
CS22umap <- cvejic[['ref.umap']]

cvejic <- readRDS(paste0(rds_path, 'cvejic_subset.rds'))
cvejic$pred.id.CS22 <- pred.id.CS22
cvejic$pred.id.score.CS22 <- pred.id.score.CS22
cvejic[['CS22umap']] <- CS22umap
saveRDS(cvejic, paste0(rds_path, 'cvejic_subset.rds'))

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))
FL[['cvejic_umap']] <- cvejic.umap
FL$pred.id.cvejic <- pred.id.cvejic
FL$pred.id.score.cvejic <- pred.id.score.cvejic
saveRDS(FL, paste0(rds_path, 'FL_only.rds'))

DimPlot(cvejic, reduction = 'umap', group.by = 'pred.id.CS22', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_proj_onFL.pdf'))

VlnPlot(FL, 'pred.id.score.cvejic', group.by = 'seurat_clusters', pt.size = 0.01, cols = annot_colors$seurat_cluster_cols) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/CS22_proj_on_cvejic_scores.pdf'), width = 8) 

VlnPlot(cvejic, 'pred.id.score.CS22', group.by = 'seurat_clusters', pt.size = 0.01, cols = annot_colors$seurat_cluster_cols) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/cvejic_proj_on_CS22_scores.pdf'), width = 8) 

DimPlot(FL, reduction = 'umap', group.by = 'pred.id.cvejic', label = TRUE, repel = TRUE, pt.size = 0.7) + 
  coord_fixed() + 
  labs(title = 'CS22 trimester - predicted IDs') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12))

(DimPlot(FL, reduction = 'cvejic_umap', group.by = 'pred.id.cvejic', pt.size = 0.6) + coord_fixed() + xlim(c(-12, 12)) + ylim(c(-12, 8))) |
  (DimPlot(FL, reduction = 'cvejic_umap', group.by = 'seurat_clusters', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed() + xlim(c(-12, 12)) + ylim(c(-12, 8)))
ggsave(paste0(fig_path, 'additional_plots/CS22_proj_on_cvejic_presubset.pdf'), width = 10) 

DimPlot(cvejic, reduction = 'CS22umap', group.by = 'pred.id.CS22', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed() 
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_proj_onCS22umap.pdf'))

VlnPlot(cvejic, 'nCount_peaks', group.by = 'seurat_clusters') + NoLegend() + ylim(0, max(cvejic$nCount_peaks))
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_peak_counts_clust.pdf'), width = 8)

VlnPlot(cvejic, 'nCount_peaks', group.by = 'sample', cols = annot_colors$sample_number2nd) + NoLegend() + ylim(0, max(cvejic$nCount_peaks))
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_peak_counts_samp.pdf'), width = 4.5)

VlnPlot(cvejic, 'nCount_peaks', group.by = 'pred.id.CS22', cols = annot_colors$seurat_cluster_cols) + NoLegend() + ylim(0, max(cvejic$nCount_peaks))
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_peak_counts_predIDs.pdf'), width = 8)

# now project to BM
BM <- readRDS(paste0(rds_path, 'BM_only.rds'))

DefaultAssay(BM) <- 'peaks'

cvejic

gr_merge <- reduce(c(granges(BM[['peaks']]), granges(cvejic[['peaks']])))
gr_merge

# Filter out possible too short/long peaks
peakwidths <- width(gr_merge)
hist(peakwidths)
gr_merge <- gr_merge[peakwidths < 10000 & peakwidths > 20]
gr_merge

BMcounts <- FeatureMatrix(Fragments(BM),
                          features = gr_merge,
                          cells = Cells(BM))

ccounts <- FeatureMatrix(Fragments(cvejic), 
                         features = gr_merge, 
                         cells = Cells(cvejic))

BM[['peaks']] <- CreateChromatinAssay(counts = BMcounts, fragments = Fragments(BM))
cvejic[['peaks']] <- CreateChromatinAssay(counts = ccounts, fragments = Fragments(cvejic))

DefaultAssay(BM) <- 'peaks'
DefaultAssay(cvejic) <- 'peaks'

BM <- RunTFIDF(BM, assay = 'peaks')
BM <- FindTopFeatures(BM, min.cutoff = 'q75')
BM <- RunSVD(BM)

cvejic <- RunTFIDF(cvejic, assay = 'peaks')
cvejic <- FindTopFeatures(cvejic, min.cutoff = 'q75')
cvejic <- RunSVD(cvejic)

transfer.anchors <- FindTransferAnchors(
  reference = BM,
  reference.assay = 'peaks',
  query = cvejic,
  query.assay = 'peaks',
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:20
)

cvejic <- MapQuery(
  anchorset = transfer.anchors,
  reference = BM,
  query = cvejic,
  refdata = BM$BM_annot,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

cvejic[['ref.umap']]@cell.embeddings[,'ref.umap_2'] <- - cvejic[['ref.umap']]@cell.embeddings[,'ref.umap_2']
cvejic$predicted.id <- factor(cvejic$predicted.id, levels = clust_order)

pred.id.BM <- cvejic$predicted.id
pred.id.score.BM <- cvejic$predicted.id.score
BMumap <- cvejic[['ref.umap']]

cvejic <- readRDS(paste0(rds_path, 'cvejic_subset.rds'))
cvejic$pred.id.BM <- pred.id.BM
cvejic$pred.id.score.BM <- pred.id.score.BM
cvejic[['BMumap']] <- BMumap

saveRDS(cvejic, paste0(rds_path, 'cvejic_subset.rds'))

DimPlot(cvejic, reduction = 'BMumap', group.by = 'pred.id.BM', pt.size = 0.6, cols = annot_colors$cluster_cols) + coord_fixed() |
  DimPlot(cvejic, reduction = 'umap', group.by = 'pred.id.BM', pt.size = 0.6, cols = annot_colors$cluster_cols) + coord_fixed()
ggsave(paste0(fig_path, 'additional_plots/cvejic_subset_proj_onBM.pdf'), width = 10)


VlnPlot(cvejic, 'pred.id.score.BM', group.by = 'pred.id.BM', cols = annot_colors$cluster_cols, pt.size = 0.01) + NoLegend() + ylim(0, 1.05)
ggsave(paste0(fig_path, 'sup4/mard_proj_on_BM_scores_predIDs_BM.pdf'), width = 8) 

# load Ranzoni et al data and plot BM UMAP using the same x and y axis ranges for all 4 datasets
ranz <- readRDS(paste0(rds_path, 'ranzoni.rds'))

xlims <- range(c(BM[['umap']]@cell.embeddings[,1], FL[['BMumap']]@cell.embeddings[,1], ranz[['BMumap']]@cell.embeddings[,1], cvejic[['BMumap']]@cell.embeddings[,1]))
ylims <- range(c(BM[['umap']]@cell.embeddings[,2], FL[['BMumap']]@cell.embeddings[,2], ranz[['BMumap']]@cell.embeddings[,2], cvejic[['BMumap']]@cell.embeddings[,2]))

DimPlot(ranz, reduction = 'BMumap', group.by = 'pred.id.BM', pt.size = 0.6, cols = annot_colors$cluster_cols) + coord_fixed() + xlim(xlims) + ylim(ylims)
ggsave(paste0(fig_path, 'sup3/ranz_proj_onBM.pdf'))

DimPlot(FL, reduction = 'BMumap', group.by = 'pred.id.BM', pt.size = 0.6, cols = annot_colors$cluster_cols) + coord_fixed() + xlim(xlims) + ylim(ylims)
ggsave(paste0(fig_path, 'fig2/FL_projected_on_BM.pdf'))

DimPlot(BM, reduction = 'umap', group.by = 'BM_annot', pt.size = 0.6, cols = annot_colors$cluster_cols, label = T) + coord_fixed() + xlim(xlims) + ylim(ylims)
ggsave(paste0(fig_path, 'fig2/BM_umap.pdf'))

DimPlot(BM, reduction = 'umap', group.by = 'sample_number', pt.size = 0.5, cols = annot_colors$sample_numberBM, shuffle = T) + coord_fixed() + xlim(xlims) + ylim(ylims)
ggsave(paste0(fig_path, 'sup2/BM_umap_harmony.pdf'))

DimPlot(cvejic, reduction = 'BMumap', group.by = 'pred.id.BM', pt.size = 0.6, cols = annot_colors$cluster_cols) + coord_fixed() + xlim(xlims) + ylim(ylims)
ggsave(paste0(fig_path, 'sup4/mard_proj_onBM.pdf'))


