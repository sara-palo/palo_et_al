# Load required libraries
library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(123)

# Load data
FL <- readRDS(file="C:/Users/keiki/OneDrive/デスクトップ/Sara_project/FL_CITE_trajectory.rds")
Idents(FL) <- "clust.names"

# Identify DEGs between Ly-II vs HSC, and GMP vs HSC
  ly_markers <- FindMarkers(FL, ident.1 = "Ly-II", ident.2 = "HSC", only.pos = FALSE)
  ly_pos <- rownames(subset(ly_markers, avg_log2FC > 1, p_val_adj <0.05))
  ly_neg <- paste0(rownames(subset(ly_markers, avg_log2FC < -1, p_val_adj <0.05)), "-")
  ly_signature <- list(ly_in_FL = c(ly_pos, ly_neg))
  
  gmp_markers <- FindMarkers(FL, ident.1 = "GMP", ident.2 = "HSC", only.pos = FALSE)
  gmp_pos <- rownames(subset(gmp_markers, avg_log2FC > 1, p_val_adj <0.05))
  gmp_neg <- paste0(rownames(subset(gmp_markers, avg_log2FC < -1, p_val_adj <0.05)), "-")
  gmp_signature <- list(gmp_in_FL = c(gmp_pos, gmp_neg))

# Score with UCell
expr_data <- FL@assays$RNA@counts
  FL <- AddMetaData(FL, metadata = ScoreSignatures_UCell(expr_data, features = ly_signature),col.name = "ly_in_FL_UCell")
  FL <- AddMetaData(FL, metadata = ScoreSignatures_UCell(expr_data, features = gmp_signature),col.name = "gmp_in_FL_UCell")

# Set cluster identity order
FL <- SetIdent(FL, value = factor(Idents(FL), levels = c("HSC", "Ly-II", "GMP")))

# Plot violin plots
pdf("violinplot_FLsubset_GMPscore_1.0.pdf")
  VlnPlot(FL, features = "gmp_in_FL_UCell", idents = c("HSC", "Ly-II", "GMP"),
          cols = c('#E9CB53', "#793885", '#56ace5'), y.max = 1, pt.size = 1)
dev.off()

pdf("violinplot_FLsubset_ly2score_1.0.pdf")
  VlnPlot(FL, features = "ly_in_FL_UCell", idents = c("HSC", "Ly-II", "GMP"),
          cols = c('#E9CB53', "#793885", '#56ace5'), y.max = 1, pt.size = 1)
dev.off()

# Pseudotime scoring along Ly2 and GMP trajectories
DefaultAssay(FL) <- "RNA"
# Compare differences in last 10% bins
FL$l2_ps_cut <- cut(FL$l2_ps, breaks = 10)
print(unique(FL$l2_ps_cut)) #20.5-22.8 for l2_ps
FL$gmp_ps_cut <- cut(FL$gmp_ps, breaks = 10)
print(unique(FL$gmp_ps_cut)) #15.8-17.6 for gmp_ps


# Ly2 trajectory
      ly_df_lytrajectory <- FetchData(FL, vars = c("l2_ps", "ly_in_FL_UCell")) %>%
        mutate(group = "FL: Ly2", UCell_Score = ly_in_FL_UCell) %>%
        select(-ly_in_FL_UCell)
      
      gmp_df_lytrajectory <- FetchData(FL, vars = c("l2_ps", "gmp_in_FL_UCell")) %>%
        mutate(group = "FL: GMP", UCell_Score = gmp_in_FL_UCell) %>%
        select(-gmp_in_FL_UCell)
      
      combined_lytrajectory <- bind_rows(ly_df_lytrajectory, gmp_df_lytrajectory)
      
            pdf("Ly2_GMP_scoring_Ly2_Trajectory_highlight_band.pdf", height = 4, width = 8)
            ggplot(combined_lytrajectory, aes(x = l2_ps, y = UCell_Score, color = group)) +
              annotate("rect", xmin = 20.5, xmax = 22.8, ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.5) +
              geom_point(alpha = 0.3, stroke = 0) +
              geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, size = 1) +
              scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
              theme_classic() +
              labs(title = "Ly2/GMP scoring in Ly2 Trajectories (last 10% highlighted)",
                   x = "l2_ps", y = "UCell Score")
            dev.off()

# GMP trajectory
      ly_df_gmtrajectory <- FetchData(FL, vars = c("gmp_ps", "ly_in_FL_UCell")) %>%
        mutate(group = "FL: Ly2", UCell_Score = ly_in_FL_UCell) %>%
        select(-ly_in_FL_UCell)
      
      gmp_df_gmtrajectory <- FetchData(FL, vars = c("gmp_ps", "gmp_in_FL_UCell")) %>%
        mutate(group = "FL: GMP", UCell_Score = gmp_in_FL_UCell) %>%
        select(-gmp_in_FL_UCell)
      
      combined_gmtrajectory <- bind_rows(ly_df_gmtrajectory, gmp_df_gmtrajectory)

            pdf("Ly2_GMP_scoring_GMP_Trajectory_highlight_band.pdf", height = 4, width = 8)
            ggplot(combined_gmtrajectory, aes(x = gmp_ps, y = UCell_Score, color = group)) +
              annotate("rect", xmin = 15.8, xmax = 17.6, ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.5) +
              geom_point(alpha = 0.3, stroke = 0) +
              geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, size = 1) +
              scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
              theme_classic() +
              labs(title = "Ly2/GMP scoring in GMP Trajectories (last 10% highlighted)",
                   x = "gmp_ps", y = "UCell Score")
            dev.off()



l2_df <- FL@meta.data[FL$l2_ps_cut == levels(FL$l2_ps_cut)[10], c("ly_in_FL_UCell", "gmp_in_FL_UCell","l2_ps")]
l2_df$diff <- l2_df$ly_in_FL_UCell - l2_df$gmp_in_FL_UCell
l2_df <- l2_df[!is.na(l2_df$diff), ]

gmp_df <- FL@meta.data[FL$gmp_ps_cut == levels(FL$gmp_ps_cut)[10], c("ly_in_FL_UCell", "gmp_in_FL_UCell")]
gmp_df$diff <- gmp_df$gmp_in_FL_UCell - gmp_df$ly_in_FL_UCell
gmp_df <- gmp_df[!is.na(gmp_df$diff), ]


# Wilcoxon test
wc <- wilcox.test(l2_df$diff, gmp_df$diff)

# Summary stats
stats_df <- data.frame(
  ly_ps_mean_diff = mean(l2_df$diff),
  my_ps_mean_diff = mean(gmp_df$diff),
  ly_ps_median_diff = median(l2_df$diff),
  my_ps_median_diff = median(gmp_df$diff),
  wilcox_pval = wc$p.value
)
stats_df
write.csv(stats_df, file='pseudotime_modscores_FL_Ly_my_stat.csv')

# Boxplot
l2_df$group <- "Lymphoid-trajectory (Ly)"
gmp_df$group <- "Myeloid-trajectory (GMP)"
plot_df <- rbind(l2_df[, c("diff", "group")], gmp_df[, c("diff", "group")])

pdf("difference_boxplot_lineage_specificity.pdf", width = 4, height = 8)
ggplot(plot_df, aes(x = group, y = diff, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = c("Lymphoid-trajectory (Ly)" = "#793885", "Myeloid-trajectory (GMP)" = "#56ace5")) +
  labs(title = "Difference in UCell Scores", y = "Score difference (Ly - GMP)", x = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")
dev.off()
library(monocle3)
sessionInfo()

