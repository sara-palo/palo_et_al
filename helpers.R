helper_packages <- c('RColorBrewer', 'cowplot', 'tictoc', 'parallel', 'assertthat', 
                     'tidyr', 'colorRamp2', 'ComplexHeatmap', 'patchwork', 'jsonlite', 'plotly')
suppressPackageStartupMessages(invisible(lapply(helper_packages, library, character.only = TRUE)))

tf_dotplot <- function(sobj = NULL,
                       end_cluster = NULL,
                       assay = 'chromvar',
                       #can also be 'metadata' meaning module scores
                       features = NULL,
                       slot_to_use = 'scale.data',
                       colors = NULL,
                       palette = 'Tableau 10',
                       xmax = NULL,
                       include.dots = TRUE,
                       scale.pseudotime = FALSE,
                       mark_range = NULL,
                       ymax = NULL,
                       ymin = NULL,
                       scale.modscore = FALSE) {
  end_cluster <- as.character(end_cluster)
  
  slingshot_df <- sobj@meta.data
  
  slingshot_df <- as.data.frame(
    list(
      "barcodes" = rownames(slingshot_df),
      "clust.names" = slingshot_df[, 'clust.names'],
      "pseudotime" = slingshot_df[, end_cluster]
    )
  ) %>%
    drop_na()
  
  if (assay == 'metadata') {
    plot_df <-
      subset(sobj@meta.data, select = features)
    colnames(plot_df) <- features
  } else {
    plot_df <- GetAssayData(object = sobj,
                            assay = assay,
                            slot = slot_to_use)
  }
  
  if (assay == 'chromvar') {
    motif.names <-
      ConvertMotifID(sobj@assays$peaks@motifs, id = rownames(plot_df))
  }
  
  if (assay == 'metadata') {
    plot_df <- data.frame(plot_df[slingshot_df$barcodes,])
    colnames(plot_df) <- features
    
    if (scale.modscore){
      plot_df <- scale(plot_df)
    }
    
  } else {
    plot_df <-
      plot_df[, slingshot_df$barcodes] %>% t() %>% data.frame()
  }
  
  if (assay == 'chromvar') {
    colnames(plot_df) <- motif.names
  }
  
  plot_df <- data.frame(cbind(plot_df, slingshot_df$pseudotime))
  
  colnames(plot_df)[ncol(plot_df)] <- 'pseudotime'
  
  if (scale.pseudotime) {
    plot_df$pseudotime <-
      scale(plot_df$pseudotime,
            center = F,
            scale = max(plot_df$pseudotime))[, 1]
  }
  
  if (is.null(colors)) {
    colors <- palette.colors(n = length(features),
                             palette = palette,
                             recycle = TRUE)
    names(colors) <- features
  }
  
  if (is.null(xmax)) {
    xmax <- max(plot_df$pseudotime)
  }
  
  ylab_titles <-
    c('chromvar' = 'Chromvar activity',
      'RNA' = 'Gene expression',
      'metadata' = 'Score')
  guide_titles <-
    c('chromvar' = 'TF',
      'RNA' = 'Gene',
      'metadata' = 'Module score')
  
  if(is.null(ymax)){ 
  ymax <- max(plot_df[, 1:(ncol(plot_df) - 1)])
  }
  
  if(is.null(ymin)){
  ymin <- min(plot_df[, 1:(ncol(plot_df) - 1)])
  }
  
  p <- ggplot() +
    theme_bw() +
    ylab(ylab_titles[assay]) +
    xlab(ifelse(scale.pseudotime, 'Scaled pseudotime', 'Pseudotime')) +
    guides(alpha = 'none', color = guide_legend(guide_titles[assay])) +
    ggtitle(end_cluster) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) +
    xlim(-0.01, xmax) + 
    ylim(ymin, ymax)
  
  geom_list <- list()
  
  if (!is.null(mark_range)) {
    geom_list <- c(geom_list,
                   geom_rect(
                     mapping = aes(
                       xmin = mark_range[1],
                       xmax = mark_range[2],
                       ymin = ymin,
                       ymax = ymax
                     ),
                     alpha = 0.3,
                     color = 'lightgrey'
                   ))
  }
  
  if (include.dots) {
    for (feature in features) {
      geom_list <- c(geom_list, geom_point(
        data = plot_df,
        aes(
          x = pseudotime,
          y = .data[[feature]],
          colour = !!as.character(feature),
          alpha = 0.01,
          stroke = 0
        )
      ))
    }
  }
  
  for (feature in features) {
    geom_list <- c(geom_list, geom_smooth(
      method = 'gam',
      formula = y ~ s(x, bs = "cs"), # ??
      data = plot_df,
      aes(
        x = pseudotime,
        y = .data[[feature]],
        colour = !!as.character(feature)
      )
    ))
  }
  
  p$layers <- geom_list
  return(p)
}

# Cellradar
DATAFILES <- list(
  "human" = paste0(ext_data_path, "cellradar_datasets/human.json"),
  "human_celltypes" = paste0(ext_data_path, "cellradar_datasets/human_celltypes.txt"),
  "mouse" = paste0(ext_data_path, "cellradar_datasets/mouse.json"),
  "mouse_celltypes" = paste0(ext_data_path, "cellradar_datasets/mouse_celltypes.txt")
)

prep_data <- function(dataset, genes) {
  normalize <- function(a) {
    a <- a - min(a)
    a <- a / max(a)
    return(c(a, a[1]))
  }
  
  blank_return <- list()
  
  if (length(genes) == 0) {
    return(c(blank_return, list(msg = "Please enter at least one gene name")))
  }
  if (!(dataset %in% names(DATAFILES))) {
    return(list(msg = "Selected dataset is invalid"))
  }
  genes <- toupper(gsub("[^0-9a-zA-Z]+-", "", genes))
  
  ref_data <- fromJSON(DATAFILES[[dataset]])
  celltypes <- read.csv(str_c(ext_data_path, "cellradar_datasets/", dataset, "_celltypes.txt"), header = F)[[1]]
  data_mean <- data.frame(matrix(ncol = length(celltypes)))[-1,]
  colnames(data_mean) = celltypes
  data_std <- data.frame(matrix(ncol = length(celltypes)))[-1,]
  colnames(data_std) = celltypes
  gene_order <- c()
  lookup <- list()
  nonexisting_genes <- c()
  
  for (i in genes) {
    if (i %in% names(lookup)) {
      next
    } else if (!(i %in% names(ref_data))) {
      nonexisting_genes <- c(nonexisting_genes, i)
    } else {
      values <- ref_data[[i]]
      data_mean[i, ] <- values[1,]
      data_std[i, ] <- values[2,]
      lookup[[i]] <- NULL
      gene_order <- c(gene_order, i)
    }
  }
  
  if (length(gene_order) == 0) {
    return(list(msg = "None of the entered gene names is valid"))
  }
  
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

plot_cellradar <- function(res_list = NULL,
                           colors = NULL,
                           names = NULL,
                           opac = NULL,
                           print_negs = TRUE,
                           dash = NULL){
  neg_list <- list()
  
  fig <- plot_ly(
    type = 'scatterpolar',
    fill = 'toself',
    mode = 'lines'
  )
  
  if (is.null(opac)){
    opac <- rep(1, length(res_list))
  }
  
  if (is.null(dash)){
    dash <- rep('solid', length(res_list))
  }
  
  for (i in 1:length(res_list)) {
    res <- res_list[[i]]
    
    neg_list[[i]] <- res$nonexisting_genes
    
    fig <- add_trace(
      fig,
      r = res$values[[1]],
      theta = names(res$values[[1]]),
      name = names[i],
      color = I(colors[i]),
      opacity = opac[i],
      line = list(shape = 'spline', 
                  dash = dash[i])
    )
  }
  
  fig <- layout(fig, polar = list(
    angularaxis = list(direction = 'clockwise'),
    radialaxis = list(visible = F, range = c(0, 1))
  ))
  if(print_negs){print(neg_list)}
  return(fig)
}

