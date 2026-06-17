source('scripts/atac_analysis/00.R')

cvejic <- readRDS(paste0(rds_path, 'cvejic_merge.rds'))
FL <- readRDS(paste0(rds_path, 'FL_only.rds'))

# add cvejic cell annotations from correspondence with authors
md <- cvejic@meta.data
# sample information from ENA - link between the two naming conventions
ENA <- read.table(paste0(ext_data_path, 'marderstein/E-MTAB-13070.sdrf.txt'), header = T, sep = '\t') 
# metadata from authors
tab <- read.table(paste0(ext_data_path, 'marderstein/multi/Multiome.RNA_ATAC.meta.txt'), header = T, sep = '\t')

table(tab$dataset)
table(ENA$Characteristics.original.source.name., ENA$Source.Name)

# rename ENA 156836A/B to 15836A/B based on communication with Cvejic lab
ENA$Characteristics.original.source.name.[str_detect(ENA$Characteristics.original.source.name., pattern = '156836')]
ENA$Characteristics.original.source.name. <- 
  str_replace(ENA$Characteristics.original.source.name., pattern = '156836', replacement = '15836')

tab$dataset2 <- str_remove(tab$dataset, pattern = 'nuclei') %>% str_remove_all(pattern = ' ') %>% str_remove_all(pattern = '_')
table(tab$dataset2)

ENA$sample <- str_remove_all(ENA$Source.Name, pattern = '_')
table(ENA$sample)

# subset for atac only and make new metadata column with corresponding names
ENA <- subset(ENA, Comment.LIBRARY_STRATEGY. == 'ATAC-seq')

id_to_name <- table(ENA$Characteristics.original.source.name., ENA$Source.Name)
colnames(id_to_name) <- str_remove_all(colnames(id_to_name), pattern = '_') %>% 
  str_remove_all(pattern = 'ATAC') %>% str_remove_all(pattern = 'TS21')

rownames(id_to_name)[! rownames(id_to_name) %in% tab$dataset2]
# ENA sample 16171B is metadata sample 16171?
rownames(id_to_name) <- str_replace(rownames(id_to_name), pattern = '16171B', replacement = '16171')

id_to_name

match <- apply(id_to_name, 1, function(x) which(x>0))
match

tab$sample2 <- colnames(id_to_name)[match[tab[, 'dataset2']]]

head(Cells(cvejic))

tab$cell2 <- paste0(tolower(tab$sample2), '_', tab$cell)
head(tab$cell2)

table(Cells(cvejic) %in% tab$cell2)

subs <- tab[tab$sample2 %in% c('H5A', 'H5B', 'H6A'), ]
table(subs$cell2 %in% Cells(cvejic))

rownames(subs) <- subs$cell2
cvejic <- AddMetaData(cvejic, metadata = subs['subclust_v6'])

clust_order_cvejic <- c('HSCs', 'B cells', 'NK cells', 'pDCs', 'Granulocyte progenitors', 'Neutrophils', 'Pro-inflammatory macrophages', 'Kupffer cells',
                        'MEMPs', 'Early erythroid', 'Late erythroid', 'Mast cells', 'Megakaryocytes', 'Stroma')

cvejic$subclust_v6 <- factor(cvejic$subclust_v6, levels = clust_order_cvejic)
table(cvejic$subclust_v6)

cols <- colorRampPalette(brewer.pal('Paired', n = 10))(14)
plot_cols <- setNames(cols, clust_order_cvejic)

DimPlot(cvejic, group.by = 'subclust_v6', label = F, pt.size = 0.6, cols = plot_cols)
ggsave(paste0(fig_path, 'additional_plots/cvejic_paperlabels.pdf'))

DimPlot(cvejic, group.by = 'subclust_v6', label = F, pt.size = 0.6, split.by = 'sample', cols = plot_cols)
ggsave(paste0(fig_path, 'additional_plots/cvejic_paperlabels_persamp.pdf'), width = 8.5, height = 4)


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
  refdata = cvejic$subclust_v6,
  reference.reduction = 'lsi',
  new.reduction.name = 'ref.lsi',
  reduction.model = 'umap'
)

FL$predicted.id <- factor(FL$predicted.id, levels = clust_order_cvejic)

table(FL$predicted.id)
table(FL$predicted.id, FL$seurat_clusters)

FL$pred.id.cvejic <- FL$predicted.id
FL$pred.id.score.cvejic <- FL$predicted.id.score

(DimPlot(cvejic, reduction = 'umap', pt.size = 0.7) ) | 
  (DimPlot(FL, reduction = 'ref.umap', group.by = 'predicted.id', pt.size = 0.7)) |
  (DimPlot(FL, reduction = 'ref.umap', group.by = 'seurat_clusters', pt.size = 0.7) )


VlnPlot(FL, 'pred.id.score.cvejic', group.by = 'pred.id.cvejic', pt.size = 0.01, cols = plot_cols) + NoLegend()
ggsave(paste0(fig_path, 'additional_plots/CS22_proj_on_cvejic_scores.pdf'), width = 8, height = 5) 

(DimPlot(FL, reduction = 'ref.umap', group.by = 'pred.id.cvejic', pt.size = 0.6, cols = plot_cols) + coord_fixed() + xlim(c(-12, 12)) + ylim(c(-12, 8))) |
  (DimPlot(FL, reduction = 'ref.umap', group.by = 'seurat_clusters', pt.size = 0.6, cols = annot_colors$seurat_cluster_cols) + coord_fixed() + xlim(c(-12, 12)) + ylim(c(-12, 8)))
ggsave(paste0(fig_path, 'additional_plots/CS22_proj_on_cvejic_presubset.pdf'), width = 13, height = 5) 

predictions_clusters <- table(FL$seurat_clusters, FL$pred.id.cvejic)
predictions_clusters <- predictions_clusters/rowSums(predictions_clusters)  # divide by number of cells for each cell type
predictions_clusters <- data.frame(predictions_clusters)

predictions_clusters$Var1 <- factor(predictions_clusters$Var1, levels = 0:13)
predictions_clusters$Var2 <- factor(predictions_clusters$Var2, levels = clust_order_cvejic)

ggplot(predictions_clusters, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = plot_cols) +
  xlab("Predicted cell type - label transfer from Marderstein et al") +
  theme_bw() +
  ylab("Frequency")
ggsave(paste0(fig_path, 'sup2/FL_cvejic_predictions_bar_seuratclust.pdf'), width = 8.7, height = 5)

predictions_clusters <- table(FL$pred.id.BM, FL$pred.id.cvejic)
predictions_clusters <- predictions_clusters/rowSums(predictions_clusters)  # divide by number of cells for each cell type
predictions_clusters <- data.frame(predictions_clusters)

predictions_clusters$Var1 <- factor(predictions_clusters$Var1, levels = clust_order)
predictions_clusters$Var2 <- factor(predictions_clusters$Var2, levels = clust_order_cvejic)

ggplot(predictions_clusters, aes(fill = Var2, y = Freq, x = Var1)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = plot_cols) +
  xlab("Predicted cell type - label transfer from Marderstein et al") +
  theme_bw() +
  ylab("Frequency")
ggsave(paste0(fig_path, 'sup2/FL_cvejic_predictions_bar_BMannot.pdf'), width = 9, height = 5)



