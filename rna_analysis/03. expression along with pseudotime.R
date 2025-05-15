library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)


# === Prepare expression data ===
FL <- readRDS(file="FL_CITE_trajectory.rds")
DefaultAssay(FL)<-"RNA"
ly_data_fl <- FetchData(FL, vars = c("l2_ps", "EBF1","HOXB7","SPIB", "TCF3"))
ly_data_fl <- na.omit(ly_data_fl)
ly_data_fl <- ly_data_fl %>%
  pivot_longer(cols = c("EBF1", "HOXB7", "SPIB",  "TCF3"), names_to = "Gene", values_to = "Expression") %>%
  mutate(Expression = scale(Expression))


tf_dotplot_l2ps <- function(data = NULL,
                            features = NULL,
                            colors  = c("EBF1" = "#E15759", 
                                        "HOXB7" = "#4E79A7", 
                                        "SPIB" = "#76B7B2",  
                                        "TCF3" = "#F28E2B"),
                            palette = 'Tableau 10',
                            xmax = NULL,
                            ymin = NULL,
                            ymax = NULL,
                            mark_range=NULL) {
  
  if (is.null(data) || is.null(features)) {
    stop("Data and features must be provided.")
  }
  features <- as.character(features)
  if (is.null(colors)) {
    colors <- palette.colors(n = length(features),
                             palette = palette,
                             recycle = TRUE)
    names(colors) <- features
  }
  if (is.null(xmax)) {
    xmax <- max(data$l2_ps, na.rm = TRUE)
  }
  if (is.null(ymax)) { 
    ymax <- max(data$Expression, na.rm = TRUE)
  }
  if (is.null(ymin)) {
    ymin <- min(data$Expression, na.rm = TRUE)
  }
  
  p <- ggplot(ly_data_fl, aes(x = l2_ps, y = Expression, color = Gene, group = Gene)) +
    #geom_point(alpha = 0.3, stroke = 0) +
    geom_smooth(method = 'gam', se = F, span = 0.1) +
    scale_y_continuous(limits = c(-1, 5)) + 
    scale_color_manual(values = c("EBF1" = "#E15759", 
                                  "HOXB7" = "#4E79A7", 
                                  "SPIB" = "#76B7B2", 
                                  "TCF3" = "#F28E2B",
                                  "PAX5"="darkred")) +
    theme_classic()
  if (!is.null(mark_range)) {
    p <- p + annotate("rect", xmin = mark_range[1], xmax = mark_range[2], 
                      ymin = ymin, ymax = ymax, fill = "lightgrey", alpha = 0.3)
  }
  return(p)
}

# === Generate and save plot ===
plot_ly <- tf_dotplot_l2ps(data = ly_data_fl, features = c("EBF1", "HOXB7", "SPIB",  "TCF3","PAX5"))
print(plot_ly)

pdf(file="250311_Ly2_trajectory_4TFs.pdf", height = 4, width = 8)
print(plot_ly)
dev.off()




# === Prepare expression data ===
gmp_data_fl <- FetchData(FL, vars = c("gmp_ps", "EBF1","HOXB7","SPIB", "CEBPA"))
gmp_data_fl <- na.omit(gmp_data_fl)
gmp_data_fl <- gmp_data_fl %>%
  pivot_longer(cols = c("EBF1", "HOXB7", "SPIB",  "CEBPA"), names_to = "Gene", values_to = "Expression") %>%
  mutate(Expression = scale(Expression))


tf_dotplot_gmp_ps <- function(data = NULL,
                              features = NULL,
                              xmax = NULL,
                              ymin = NULL,
                              ymax = NULL,
                              mark_range = NULL) {  
  
  if (is.null(data) || is.null(features)) {
    stop("Data and features must be provided.")
  }
  
  features <- as.character(features)
  if (is.null(xmax)) {
    xmax <- max(data$gmp_ps, na.rm = TRUE)
  }
  if (is.null(ymax)) { 
    ymax <- max(data$Expression, na.rm = TRUE)
  }
  if (is.null(ymin)) {
    ymin <- min(data$Expression, na.rm = TRUE)
  }
  colors <- c("EBF1" = "#E15759", 
              "HOXB7" = "#4E79A7",
              "SPIB" = "#76B7B2", 
              "CEBPA" = "#59A14F")
  p <- ggplot(data, aes(x = gmp_ps, y = Expression, color = Gene, group = Gene)) +
    #geom_point(alpha = 0.3, stroke = 0) +
    geom_smooth(method = 'gam', se = F, span = 0.1) + 
    scale_y_continuous(limits = c(-1, 5)) + # 適切な範囲に設定
    scale_color_manual(values = colors) + 
    theme_classic()
  if (!is.null(mark_range)) {
    p <- p + annotate("rect", xmin = mark_range[1], xmax = mark_range[2], 
                      ymin = ymin, ymax = ymax, fill = "lightgrey", alpha = 0.3)
  }
  return(p)
}

# === Generate and save plot ===
plot_gmp <- tf_dotplot_gmp_ps(data = gmp_data_fl, 
                              features = c("EBF1", "HOXB7", "SPIB", "CEBPA"))

pdf(file="250311_gmp_trajectory_4TFs.pdf", height = 4, width = 8)
print(plot_gmp)
dev.off()

