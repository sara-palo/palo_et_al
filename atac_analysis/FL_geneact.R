source('scripts/atac_analysis/00.R')

FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

DefaultAssay(FL) <- 'peaks'

gene.activities.FL <- GeneActivity(FL)

saveRDS(gene.activities.FL, paste0(rds_path, 'FL_geneact.rds'))


FL[['RNA']] <- CreateAssayObject(counts = gene.activities.FL)

FL <- NormalizeData(
  object = FL,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(FL$nCount_RNA)
)

Idents(FL) <- 'seurat_clusters'

VlnPlot(FL, features = 'nCount_RNA', group.by = 'seurat_clusters', split.by = 'sample_number', cols = annot_colors$sample_numberFL, pt.size = 0.01, idents = c(0,4,8,13))
ggsave(paste0(fig_path, 'additional_plots/FL_ncount_RNA.pdf'))

VlnPlot(FL, features = 'nFeature_RNA', group.by = 'seurat_clusters', split.by = 'sample_number', cols = annot_colors$sample_numberFL, pt.size = 0.01, idents = c(0,4,8,13))
ggsave(paste0(fig_path, 'additional_plots/FL_nfeat_RNA.pdf'))
