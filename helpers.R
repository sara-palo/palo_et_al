helper_packages <- c('RColorBrewer', 'cowplot', 'tictoc', 'parallel', 'assertthat', 
                     'tidyr', 'colorRamp2', 'ComplexHeatmap', 'patchwork', 'jsonlite', 'plotly')
suppressPackageStartupMessages(invisible(lapply(helper_packages, library, character.only = TRUE)))

do_heatmap <- function(sobj = NULL,
                       features = NULL,
                       slingshot_df,
                       slot_to_use = "data",
                       show_row_names = TRUE,
                       cluster_features = FALSE,
                       supplied_data = NULL) {
  if (!is.null(supplied_data)) {
    mat_to_use <- supplied_data
  } else {
    mat_to_use <- slot(sobj@assays$RNA, slot_to_use)
  }
  subset_matrix <- mat_to_use[
    features,
    slingshot_df$barcodes
  ] %>%
    as.matrix()
  col_fun <- colorRamp2(
    c(0, quantile(subset_matrix, probs = c(0.99))[[1]]),
    c("darkblue", "green")
  )
  present_clusters <- unique(slingshot_df$clust.names)
  n_clusters <- (sobj[["clust.names"]] %>% unique() %>% dim())[1]
  
  #FL_colours <- brewer.pal(9, "Set1")
  #FL_colours <- colorRampPalette(FL_colours)(n_clusters)
  # FL_colours <- palette.colors(
  #   n = n_clusters, palette = "Classic Tableau"
  # )
  #names(FL_colours) <- as.character(unique(sobj[["clust.names"]])[,1])
  FL_colours <- palette.colors(
    n = n_clusters, palette = "Classic Tableau", recycle = TRUE
  )
  names(FL_colours) <- as.character(unique(sobj[["clust.names"]])[,1])
  cluster_annotation <- HeatmapAnnotation(
    df = slingshot_df["clust.names"],
    name = "clust.names",
    #col = list(clust.names = FL_colours[present_clusters]),
    annotation_legend_param = list(
      legend_direction = "horizontal",
      title_position = "lefttop-rot"
    )
  )
  return(
    Heatmap(
      mat = subset_matrix,
      show_column_names = F,
      show_row_names = show_row_names,
      cluster_rows = cluster_features,
      cluster_columns = F,
      col = col_fun,
      column_title = "Cells ordered by pseudotime -->",
      row_title = "Features",
      heatmap_legend_param = list(
        title = "Expression",
        title_position = "lefttop-rot"
      ),
      border_gp = gpar(col = "black", lty = "solid"),
      column_order = slingshot_df$barcodes,
      top_annotation = cluster_annotation,
      use_raster = FALSE
    )
  )
}


correlation_heatmap <- function(sobj = NULL,
                                end_cluster = NULL,
                                quant = .95,
                                slot_to_use = "data",
                                cluster_features = FALSE,
                                supplied_data = NULL) {
  end_cluster <- as.character(end_cluster)
  quantil <- quantile(
    abs(correlations[[end_cluster]]),
    probs = quant, na.rm = T
  )
  highest_features <-
    VariableFeatures(sobj)[which(abs(correlations[[end_cluster]]) > quantil)]
  # slingshot_df <- colData(sce[[traj_dict[[end_cluster]][["lvl1"]]]])
  # slingshot_df[["slingshot"]] <- NULL
  # slingshot_df <- slingshot_df %>% as.data.frame()
  slingshot_df <- sobj@meta.data
  slingshot_df <- as.data.frame(list(
    "barcodes" = rownames(slingshot_df),
    "clust.names" = slingshot_df$clust.names,
    "pseudotime" = sobj@meta.data[, end_cluster]
    # "pseudotime" = slingshot_df[[str_c("slingPseudotime_",
    # traj_dict[[end_cluster]][["lvl2"]])]]
    # "pseudotime" = pseudotimes_of_interest[[end_cluster]]
  )) %>%
    arrange(pseudotime) %>%
    drop_na()
  do_heatmap(sobj = sobj,
             features = highest_features,
             slingshot_df = slingshot_df,
             slot_to_use = slot_to_use,
             cluster_features = cluster_features,
             supplied_data = supplied_data,
  ) +
    labs(caption = end_cluster)
  
}


mod_heatmap <- function(sobj = NULL,
                        end_cluster = NULL,
                        cluster_features = FALSE,
                        modules = NULL,
                        slot_to_use = 'data') {
  end_cluster <- as.character(end_cluster)
  slingshot_df <- sobj@meta.data
  slingshot_df <- as.data.frame(list(
    "barcodes" = rownames(slingshot_df),
    "clust.names" = slingshot_df$clust.names,
    "pseudotime" = sobj@meta.data[, end_cluster]
  )) %>%
    arrange(pseudotime) %>%
    drop_na()
  do_mod_heatmap(sobj = sobj,
                 modules = modules,
                 slingshot_df = slingshot_df,
                 slot_to_use = slot_to_use,
                 cluster_features = cluster_features,
  ) +
    labs(caption = end_cluster)
  
}


do_mod_heatmap <- function(sobj = NULL,
                           slingshot_df,
                           slot_to_use = slot_to_use,
                           show_row_names = TRUE,
                           cluster_features = FALSE,
                           modules = NULL) {
  
  mat_to_use <- sobj@meta.data[, modules]
  mat_to_use <- t(mat_to_use)
  
  subset_matrix <- mat_to_use[
    ,
    slingshot_df$barcodes
  ] %>%
    as.matrix()
  col_fun <- colorRamp2(
    c(0, quantile(subset_matrix, probs = c(0.99))[[1]]),
    c("darkblue", "green")
  )
  present_clusters <- unique(slingshot_df$clust.names)
  n_clusters <- (sobj[["clust.names"]] %>% unique() %>% dim())[1]
  FL_colours <- palette.colors(
    n = n_clusters, palette = "Classic Tableau", recycle = TRUE
  )
  # names(FL_colours) <-
  #   as.character(unique(sobj[["clust.names"]])[, 1])
  names(FL_colours) <-
    as.character(levels(sobj$clust.names))
  cluster_annotation <- HeatmapAnnotation(
    df = slingshot_df["clust.names"],
    name = "clust.names",
    col = list(clust.names = FL_colours[as.character(present_clusters)]),
    annotation_legend_param = list(legend_direction = "horizontal",
                                   title_position = "lefttop-rot")
  )
  return(
    Heatmap(
      mat = subset_matrix,
      show_column_names = F,
      show_row_names = show_row_names,
      cluster_rows = cluster_features,
      cluster_columns = F,
      col = col_fun,
      column_title = "Cells ordered by pseudotime -->",
      row_title = "Module scores",
      heatmap_legend_param = list(
        title = "Expression",
        title_position = "lefttop-rot"
      ),
      border_gp = gpar(col = "black", lty = "solid"),
      column_order = slingshot_df$barcodes,
      top_annotation = cluster_annotation,
      use_raster = FALSE
    )
  )
}


#not only for motifs, also enhancers etc
correlation_heatmap_motifs <- function(sobj = NULL,
                                       end_cluster = NULL,
                                       quant = .95,
                                       slot_to_use = "data",
                                       cluster_features = FALSE,
                                       supplied_data = NULL,
                                       assay = NULL,
                                       show_row_names = FALSE,
                                       use_raster = FALSE) {
  end_cluster <- as.character(end_cluster)
  quantil <- quantile(abs(correlations[[end_cluster]]),
                      probs = quant, na.rm = T)
  highest_features <-
    rownames(sobj[[assay]]@data)[which(abs(correlations[[end_cluster]]) > quantil)]
  slingshot_df <- sobj@meta.data
  slingshot_df <- as.data.frame(
    list(
      "barcodes" = rownames(slingshot_df),
      "clusters" = slingshot_df$clusters,
      "pseudotime" = sobj@meta.data[, end_cluster]
    )
  ) %>%
    arrange(pseudotime) %>%
    drop_na()
  do_motif_heatmap(
    sobj = sobj,
    features = highest_features,
    slingshot_df = slingshot_df,
    slot_to_use = slot_to_use,
    cluster_features = cluster_features,
    supplied_data = supplied_data,
    assay = assay,
    show_row_names = show_row_names,
    use_raster = use_raster
  ) +
    labs(caption = end_cluster)
}



do_motif_heatmap <- function(sobj = NULL,
                             features = NULL,
                             slingshot_df,
                             slot_to_use = "data",
                             show_row_names = show_row_names,
                             cluster_features = FALSE,
                             supplied_data = NULL,
                             assay = NULL,
                             use_raster = use_raster) {
  if (!is.null(supplied_data)) {
    mat_to_use <- supplied_data
  } else {
    mat_to_use <- slot(sobj[[assay]], slot_to_use)
    #mat_to_use <- sobj[[assay]]@data
  }
  subset_matrix <- mat_to_use[features,
                              slingshot_df$barcodes] %>%
    as.matrix()
  col_fun <- colorRamp2(c(0, quantile(subset_matrix, probs = c(0.99))[[1]]),
                        c("darkblue", "green"))
  present_clusters <- unique(slingshot_df$clust.names)
  n_clusters <- (sobj[["clust.names"]] %>% unique() %>% dim())[1]
  FL_colours <- palette.colors(
    n = n_clusters, palette = "Classic Tableau", recycle = TRUE
  )
  # names(FL_colours) <-
  #   as.character(unique(sobj[["clust.names"]])[, 1])
  names(FL_colours) <-
    as.character(levels(sobj$clust.names))
  cluster_annotation <- HeatmapAnnotation(
    df = slingshot_df["clust.names"],
    name = "clust.names",
    col = list(clust.names = FL_colours[as.character(present_clusters)]),
    annotation_legend_param = list(legend_direction = "horizontal",
                                   title_position = "lefttop-rot")
  )
  return(
    Heatmap(
      mat = subset_matrix,
      show_column_names = F,
      show_row_names = show_row_names,
      cluster_rows = cluster_features,
      cluster_columns = F,
      col = col_fun,
      column_title = "Cells ordered by pseudotime -->",
      row_title = "Features",
      heatmap_legend_param = list(title = "Activity",
                                  title_position = "lefttop-rot"),
      border_gp = gpar(col = "black", lty = "solid"),
      column_order = slingshot_df$barcodes,
      top_annotation = cluster_annotation,
      use_raster = use_raster
    )
  )
}


#added 2023-07-19
da_enhancer_heatmap <- function(sobj = NULL,
                                end_cluster = NULL,
                                #quant = .95,
                                slot_to_use = "data",
                                cluster_features = FALSE,
                                supplied_data = NULL,
                                assay = NULL,
                                show_row_names = FALSE,
                                features = features,
                                use_raster = FALSE) {
  end_cluster <- as.character(end_cluster)
  slingshot_df <- sobj@meta.data
  slingshot_df <- as.data.frame(
    list(
      "barcodes" = rownames(slingshot_df),
      "clusters" = slingshot_df$clusters,
      "pseudotime" = sobj@meta.data[, end_cluster]
    )
  ) %>%
    arrange(pseudotime) %>%
    drop_na()
  do_motif_heatmap(
    sobj = sobj,
    features = features,
    slingshot_df = slingshot_df,
    slot_to_use = slot_to_use,
    cluster_features = cluster_features,
    supplied_data = supplied_data,
    assay = assay,
    show_row_names = show_row_names,
    use_raster = use_raster
  ) +
    labs(caption = end_cluster)
}


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
      "clust.names" = slingshot_df[, 'clust.names'], #slingshot_df[, 'clusters']
      "pseudotime" = slingshot_df[, end_cluster]
    )
  ) %>%
    drop_na()
  
  if (assay == 'metadata') {
    plot_df <-
      subset(sobj@meta.data, select = features) #change to include rownames in dataframe as column
    colnames(plot_df) <- features
    # rownames(plot_df) <- rownames(sobj@meta.data)
  } else {
    plot_df <- GetAssayData(object = sobj,
                            assay = assay,
                            slot = slot_to_use)
  }
  
  if (assay == 'chromvar') {
    motif.names <-
      ConvertMotifID(sobj@assays$peaks@motifs, id = rownames(plot_df))
    #rownames(plot_df) <- motif.names #unnecessary since they will get messed up after transposing anyway
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

make_pseudobulk <- function(sobj = NULL, assay = NULL, slot = 'data', avg.by = NULL){
  
  pseudobulk <- AggregateExpression(sobj, assays = assay, slot = slot) 
  ps <- data.frame(pseudobulk[[assay]], check.names = F)
  
  if(! is.null(avg.by)){
    cellcounts <- data.frame(table(sobj[[avg.by]]))
    rownames(cellcounts) <- cellcounts[, avg.by]
    cellcounts <- cellcounts[colnames(ps), ] # reorder according to pseudobulk df
    ps <- t(t(ps)/cellcounts$Freq)
  }
  
  return(ps)
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

