FL_CITE <- readRDS('../CITE-seq_FL/FL.combined.rds')

DefaultAssay(FL_CITE) <- 'RNA'
FL_CITE <- FindVariableFeatures(FL_CITE, nfeatures = 1500)

#re-run umap since the model wasn't saved and can't be used for label transfer (parameters from Blood Adv paper)
FL_CITE <- RunUMAP(FL_CITE,
                   dims = 1:9,
                   n.neighbors = 30,
                   reduction = 'pca',
                   assay = 'RNA',
                   umap.method = "umap-learn",
                   verbose = TRUE,
                   return.model = TRUE
)

DimPlot(FL_CITE)

#invert UMAP y axis to make it the same as before
FL_CITE@reductions$umap@cell.embeddings[, 2] <- -(FL_CITE@reductions$umap@cell.embeddings[, 2])

DimPlot(FL_CITE) + coord_fixed()
ggsave(paste0(fig_path, 'FL_CITE_new_umap.pdf'))

saveRDS(FL_CITE, paste0(rds_path, 'FL_CITE.rds'))

