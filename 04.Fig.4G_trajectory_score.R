library(Seurat)
library(UCell)
library(ggplot2)
set.seed(123)


FL <- readRDS(file="FL_CITE_trajectory.rds")
Idents(FL)<-"clust.names"
###make a dataset 
lygenes_in_FLdata <- FindMarkers(FL,ident.1 = "Ly-II", ident.2 = "HSC",  only.pos=F)
lygenes_in_FLdata_filtered_pos <- lygenes_in_FLdata[lygenes_in_FLdata$avg_log2FC > 1, ]
lygenes_in_FLdata_filtered_neg <- lygenes_in_FLdata[lygenes_in_FLdata$avg_log2FC < -1, ]
rownames(lygenes_in_FLdata_filtered_neg) <- paste(rownames(lygenes_in_FLdata_filtered_neg), "-", sep="")
lygenes_in_FLdata_filtered <- list(ly_in_FL=c(rownames(lygenes_in_FLdata_filtered_pos),
                                              rownames(lygenes_in_FLdata_filtered_neg)))


gmpgenes_in_FLdata <- FindMarkers(FL,ident.1 = "GMP",  ident.2 = "HSC",only.pos=F)
gmpgenes_in_FLdata_filtered_pos <- gmpgenes_in_FLdata[gmpgenes_in_FLdata$avg_log2FC > 1, ]
gmpgenes_in_FLdata_filtered_neg <- gmpgenes_in_FLdata[gmpgenes_in_FLdata$avg_log2FC < -1, ]
rownames(gmpgenes_in_FLdata_filtered_neg) <- paste(rownames(gmpgenes_in_FLdata_filtered_neg), "-", sep="")
gmpgenes_in_FLdata_filtered <- list(gmp_in_FL=c(rownames(gmpgenes_in_FLdata_filtered_pos),
                                                rownames(gmpgenes_in_FLdata_filtered_neg)))


###Addmetadata
FL_data<-FL@assays$RNA@counts

ly2_score_in_FLdata.scores <- ScoreSignatures_UCell(FL_data,features=lygenes_in_FLdata_filtered)
FL <- AddMetaData(FL,metadata=ly2_score_in_FLdata.scores)

gmp_score_in_FLdata.scores <- ScoreSignatures_UCell(FL_data,features=gmpgenes_in_FLdata_filtered)
FL <- AddMetaData(FL,metadata=gmp_score_in_FLdata.scores)

FL <- SetIdent(FL, value = factor(Idents(FL), levels = c("HSC", "Ly-II", "GMP")))
pdf("violinplot_FLsubset_GMPscore_1.0.pdf")
      VlnPlot(
        FL, 
        features = "gmp_in_FL_UCell", 
        idents = c("HSC", "Ly-II", "GMP"), 
        cols = c('#E9CB53', "#793885", '#56ace5'),
        y.max =0.75,pt.size = 1)
      
dev.off()

pdf("violinplot_FLsubset_ly2score_1.0.pdf") 
      VlnPlot(
        FL, 
        features = "ly_in_FL_UCell", 
        idents = c("HSC", "Ly-II", "GMP"), 
        cols = c('#E9CB53', '#793885', '#56ace5'),
        y.max =0.75,pt.size = 1)
dev.off()      

            
            DefaultAssay(FL)<-"RNA"
            ly_data_fl <- FetchData(FL, vars = c("l2_ps", 
                                                 "ly_in_FL_UCell",
                                                 "EBF1","HOXB7","SPIB","TCF3"))
            gmp_data_fl <- FetchData(FL, vars = c("l2_ps", "gmp_in_FL_UCell",
                                                  "EBF1","HOXB7","SPIB","TCF3"))
            
            
            ly_data_fl <- na.omit(ly_data_fl)
            gmp_data_fl <- na.omit(gmp_data_fl)
            
            
            ly_data_fl <- ly_data_fl %>%
              mutate(group = "FL: Ly2",
                     UCell_Score = ly_in_FL_UCell) %>%
              select(-ly_in_FL_UCell)
            
            gmp_data_fl <- gmp_data_fl %>%
              mutate(group = "FL: GMP",
                     UCell_Score = gmp_in_FL_UCell) %>%
              select(-gmp_in_FL_UCell)
            combined_data <- bind_rows(ly_data_fl, gmp_data_fl)
            
            pdf("Ly2 GMP scoring in Ly2 Trajectories_wo_dot.pdf",height=4,width=8)
            
            ggplot(combined_data, aes(x = l2_ps, y = UCell_Score, color = group)) +
               geom_smooth(method = "gam", formula = y ~ s(x), se = T, size = 1) +
              geom_point(alpha = 0.3, stroke = 0) +
              labs(
                x = "l2_ps",
                y = "UCell Score",
                title = "Ly2/GMP scoring in Ly2 Trajectories (FLonly)"
              ) +
              scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
              theme_classic()

            dev.off()


ly_data_fl <- FetchData(FL, vars = c("gmp_ps", "ly_in_FL_UCell",
                                     "EBF1","HOXB7","SPIB","TCF3"))
gmp_data_fl <- FetchData(FL, vars = c("gmp_ps", "gmp_in_FL_UCell",
                                      "EBF1","HOXB7","SPIB","TCF3"))
ly_data_fl <- na.omit(ly_data_fl)
gmp_data_fl <- na.omit(gmp_data_fl)

ly_data_fl <- ly_data_fl %>%
  mutate(group = "FL: Ly2",
         UCell_Score = ly_in_FL_UCell) %>%
  select(-ly_in_FL_UCell)

gmp_data_fl <- gmp_data_fl %>%
  mutate(group = "FL: GMP",
         UCell_Score = gmp_in_FL_UCell) %>%
  select(-gmp_in_FL_UCell)
combined_data <- bind_rows(ly_data_fl, gmp_data_fl)
            
            pdf(file="Ly2 GMP scoring in gmp Trajectories_wo_dot.pdf",height = 4,width=8)
            ggplot(combined_data, aes(x = gmp_ps, y = UCell_Score, color = group)) +
              geom_smooth(method = "gam", formula = y ~ s(x), se = T, size = 1) +
              geom_point(alpha = 0.3, stroke = 0) +
              labs(
                x = "gmp_ps",
                y = "UCell Score",
                title = "Ly2/GMP scoring in GMP Trajectories (FLonly)"
              ) +
              scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
              theme_classic()
            
            dev.off()
            


            



