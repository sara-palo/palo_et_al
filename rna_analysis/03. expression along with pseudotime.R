library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)
library(tidyverse)
set.seed(123)

PRJ <- "C:/Users/keiki/OneDrive/デスクトップ/Sara_project"
setwd(PRJ)
out_dir <- sprintf("output_%s", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir


FL <- readRDS(file="rds/FL_CITE_trajectory_gmp2.rds")
DefaultAssay(FL)<-"RNA"
  # Lmphoid
  ly_data_fl <- FetchData(FL, vars = c("l2_ps", "EBF1","HOXB7","SPIB", "TCF3"))
  ly_data_fl <- na.omit(ly_data_fl)
  ly_data_fl <- ly_data_fl %>%
    pivot_longer(cols = c("EBF1", "HOXB7", "SPIB",  "TCF3"), names_to = "Gene", values_to = "Expression") %>%
    mutate(Expression = scale(Expression))
  
  colors  = c("EBF1" = "#E15759", "HOXB7" = "#4E79A7", "SPIB" = "#76B7B2",  "TCF3" = "#F28E2B")
  
  p_ly <- ggplot(ly_data_fl, aes(x = l2_ps, y = Expression, color = Gene, group = Gene)) +
    #geom_point(alpha = 0.3, stroke = 0) +
    geom_smooth(method = 'gam',formula = y ~ s(x), se = F, span = 0.1) +
    coord_cartesian(ylim = c(-1, 5)) +
    #scale_y_continuous(limits = c(-1, 5))+
    scale_color_manual(values = colors) +
    theme_classic()
  
  pdf(file=paste0(out_dir,"/scRNA_Ly2_trajectory_4TFs.pdf"), height = 4, width = 8)
  p_ly 
  dev.off()

      # GMP
      gmp_data_fl <- FetchData(FL, vars = c("gmp_ps", "EBF1","HOXB7","SPIB", "CEBPA"))
      gmp_data_fl <- na.omit(gmp_data_fl)
      gmp_data_fl <- gmp_data_fl %>%
        pivot_longer(cols = c("EBF1", "HOXB7", "SPIB",  "CEBPA"), names_to = "Gene", values_to = "Expression") %>%
        mutate(Expression = scale(Expression))
      colors <- c("EBF1" = "#E15759", "HOXB7" = "#4E79A7", "SPIB" = "#76B7B2", "CEBPA" = "#59A14F")
      
      p_gmp <- ggplot(gmp_data_fl, aes(x = gmp_ps, y = Expression, color = Gene, group = Gene)) +
        #geom_point(alpha = 0.3, stroke = 0) +
        geom_smooth(method = 'gam',formula = y ~ s(x, bs = "cs"), se = F, span = 0.1) +
        coord_cartesian(ylim = c(-1, 5)) +
        #scale_y_continuous(limits = c(-1, 5))+
        scale_color_manual(values = colors) + 
        theme_classic()
      
      pdf(file=paste0(out_dir,"/scRNA_gmp_trajectory_4TFs.pdf"), height = 4, width = 8)
      p_gmp
      dev.off()
    
    
    # GMP2
    gmp2_data_fl <- FetchData(FL, vars = c("gmp2_ps", "EBF1","HOXB7","SPIB", "CEBPA"))
    gmp2_data_fl <- na.omit(gmp2_data_fl)
    gmp2_data_fl <- gmp2_data_fl %>%
      pivot_longer(cols = c("EBF1", "HOXB7", "SPIB",  "CEBPA"), names_to = "Gene", values_to = "Expression") %>%
      mutate(Expression = scale(Expression))
    colors <- c("EBF1" = "#E15759", "HOXB7" = "#4E79A7", "SPIB" = "#76B7B2", "CEBPA" = "#59A14F")
    p_gmp2 <- ggplot(gmp2_data_fl, aes(x = gmp2_ps, y = Expression, color = Gene, group = Gene)) +
      #geom_point(alpha = 0.3, stroke = 0) +
      geom_smooth(method = 'gam',formula = y ~ s(x), se = F, span = 0.1) +
      coord_cartesian(ylim = c(-1, 5)) +
      #scale_y_continuous(limits = c(-1, 5))+
      scale_color_manual(values = colors) +
      theme_classic()
    p_gmp2
  
    FeaturePlot(FL,features=c("gmp2_ps","CEBPA"),reduction = "umap",pt.size = 2)
    

    pdf(file=paste0(out_dir,"/(gmp2)scRNA_gmp2_trajectory_4TFs.pdf"), height = 4, width = 8)
    p_gmp2
    dev.off()

    
