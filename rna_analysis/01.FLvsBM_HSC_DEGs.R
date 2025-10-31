# Load libraries
library(tidyverse)  # dplyr, readr, etc.
library(Seurat)
library(stringr)
library(jsonlite)
library(plotly)   

set.seed(123)

PRJ <- "C:/Users/keiki/OneDrive/デスクトップ/Sara_project"
setwd(PRJ)
out_dir <- sprintf("output_%s", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir

# Define reference files
DATAFILES <- list(
  "human" = "cellradar_datasets/human.json",
  "human_celltypes" = "cellradar_datasets/human_celltypes.txt",
  "mouse" = "cellradar_datasets/mouse.json",
  "mouse_celltypes" = "cellradar_datasets/mouse_celltypes.txt"
)

# Normalize and fetch radar data for selected genes
prep_data <- function(dataset, genes) {
  normalize <- function(a) {
    a <- a - min(a)
    a <- a / max(a)
    return(c(a, a[1]))
  }
  
  if (length(genes) == 0) return(list(msg = "Please enter at least one gene name"))
  if (!(dataset %in% names(DATAFILES))) return(list(msg = "Selected dataset is invalid"))
  
  genes <- toupper(gsub("[^0-9a-zA-Z]+-", "", genes))
  ref_data <- fromJSON(DATAFILES[[dataset]])
  celltypes <- read.csv(str_c("cellradar_datasets/", dataset, "_celltypes.txt"), header = FALSE)[[1]]
  
  data_mean <- data.frame(matrix(ncol = length(celltypes)))[-1,]
  colnames(data_mean) <- celltypes
  data_std <- data.frame(matrix(ncol = length(celltypes)))[-1,]
  colnames(data_std) <- celltypes
  
  gene_order <- c()
  lookup <- list()
  nonexisting_genes <- c()
  
  for (i in genes) {
    if (!(i %in% names(ref_data))) {
      nonexisting_genes <- c(nonexisting_genes, i)
    } else {
      values <- ref_data[[i]]
      data_mean[i, ] <- values[1, ]
      data_std[i, ] <- values[2, ]
      lookup[[i]] <- NULL
      gene_order <- c(gene_order, i)
    }
  }
  
  if (length(gene_order) == 0) return(list(msg = "None of the entered gene names is valid"))
  
  y1 <- normalize(apply(data_mean, 2, median))
  y2 <- normalize(apply(data_mean - data_std, 2, median))
  y3 <- normalize(apply(data_mean + data_std, 2, median))
  
  return(list(
    mean_df = data_mean,
    std_df = data_std,
    values = list(y1, y2, y3),
    nonexisting_genes = nonexisting_genes
  ))
}

# Plot radar chart using plotly
plot_cellradar <- function(res_list = NULL, colors = NULL, names = NULL, opac = NULL, print_negs = TRUE) {
  fig <- plotly::plot_ly(
    type = 'scatterpolar',
    fill = 'toself',
    mode = 'lines',
    line = list(shape = 'spline')
  )
  
  if (is.null(opac)) opac <- rep(1, length(res_list))
  
  neg_list <- list()
  for (i in seq_along(res_list)) {
    res <- res_list[[i]]
    neg_list[[i]] <- res$nonexisting_genes
    fig <- plotly::add_trace(
      fig,
      r = res$values[[1]],
      theta = names(res$values[[1]]),
      name = names[i],
      line = list(color = colors[i]),
      fillcolor = colors[i],  
      opacity = opac[i]
    )
  }
  
  fig <- plotly::layout(fig, polar = list(
    angularaxis = list(direction = 'clockwise'),
    radialaxis = list(visible = FALSE, range = c(0, 1))
  ))
  
  if (print_negs) print(neg_list)
  return(fig)
}

# Load Seurat object and subset HSCs
seurat <- readRDS("rds/FL_BM.combined0426.rds")
hsc <- subset(seurat, combined_clusts %in% c("FL_HSC", "BM_HSC-I", "BM_HSC-II", "BM_HSC-III"))
Idents(hsc) <- "dataset"

# FL vs BM differential expression
FL_dif <- FindMarkers(hsc, ident.1 = "FL", ident.2 = "BM", only.pos = TRUE) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1)
write.csv(FL_dif, paste0(out_dir,"/scRNA_DEGs_RNA_FL_HSC_cellrader.csv"))
res_fl <- prep_data(dataset = 'human', genes = rownames(FL_dif))
fig_fl <- plot_cellradar(list(res_fl), colors = c('#1c6494'), names = c("FL"))

# BM vs FL differential expression
BM_dif <- FindMarkers(hsc, ident.1 = "BM", ident.2 = "FL", only.pos = TRUE) %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1)
write.csv(BM_dif, paste0(out_dir,"/scRNA_DEGs_RNA_BM_HSC_cellrader.csv"))
res_bm <- prep_data(dataset = 'human', genes = rownames(BM_dif))
fig_bm <- plot_cellradar(list(res_bm), colors = c('#1b4b1f'), names = c("BM"))

# Combine plots
fig_combined <- plot_cellradar(
  res_list = list(res_fl, res_bm),
  colors = c('#1c6494', '#1b4b1f'),
  names = c("FL", "BM")
)

# Save missing genes
missinggenes <- data.frame(
  FL = str_flatten(res_fl$nonexisting_genes, collapse = ','),
  BM = str_flatten(res_bm$nonexisting_genes, collapse = ',')
)

write.csv(missinggenes, paste0(out_dir,"/scRNA_missing_genes_FL_BM_HSC_cellrader.csv"), row.names = FALSE)
fig_combined # export as pdf 
