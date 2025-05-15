# Load required libraries
library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(123)

# Load data
FL <- readRDS(file="D:/Sara_RNAseq/FL_CITE_trajectory.rds")
Idents(FL) <- "clust.names"

# Identify DEGs between Ly-II vs HSC, and GMP vs HSC
ly_markers <- FindMarkers(FL, ident.1 = "Ly-II", ident.2 = "HSC", only.pos = FALSE)
ly_pos <- rownames(subset(ly_markers, avg_log2FC > 1))
ly_neg <- paste0(rownames(subset(ly_markers, avg_log2FC < -1)), "-")
ly_signature <- list(ly_in_FL = c(ly_pos, ly_neg))

gmp_markers <- FindMarkers(FL, ident.1 = "GMP", ident.2 = "HSC", only.pos = FALSE)
gmp_pos <- rownames(subset(gmp_markers, avg_log2FC > 1))
gmp_neg <- paste0(rownames(subset(gmp_markers, avg_log2FC < -1)), "-")
gmp_signature <- list(gmp_in_FL = c(gmp_pos, gmp_neg))

# Score with UCell
expr_data <- FL@assays$RNA@counts
FL <- AddMetaData(FL, metadata = ScoreSignatures_UCell(expr_data, features = ly_signature), col.name = "ly_in_FL_UCell")
FL <- AddMetaData(FL, metadata = ScoreSignatures_UCell(expr_data, features = gmp_signature), col.name = "gmp_in_FL_UCell")
# Set cluster identity order
FL <- SetIdent(FL, value = factor(Idents(FL), levels = c("HSC", "Ly-II", "GMP")))

# Plot violin plots
pdf("violinplot_FLsubset_GMPscore_1.0.pdf")
VlnPlot(FL, features = "gmp_in_FL_UCell", idents = c("HSC", "Ly-II", "GMP"),
        cols = c('#E9CB53', "#793885", '#56ace5'), y.max = 0.75, pt.size = 1)
dev.off()

pdf("violinplot_FLsubset_ly2score_1.0.pdf")
VlnPlot(FL, features = "ly_in_FL_UCell", idents = c("HSC", "Ly-II", "GMP"),
        cols = c('#E9CB53', "#793885", '#56ace5'), y.max = 0.75, pt.size = 1)
dev.off()

# Pseudotime scoring along Ly2 and GMP trajectories
DefaultAssay(FL) <- "RNA"

# Ly2 trajectory
ly_df <- FetchData(FL, vars = c("l2_ps", "ly_in_FL_UCell")) %>%
  mutate(group = "FL: Ly2", UCell_Score = ly_in_FL_UCell) %>%
  select(-ly_in_FL_UCell)

gmp_df <- FetchData(FL, vars = c("l2_ps", "gmp_in_FL_UCell")) %>%
  mutate(group = "FL: GMP", UCell_Score = gmp_in_FL_UCell) %>%
  select(-gmp_in_FL_UCell)

combined <- bind_rows(ly_df, gmp_df)

# Highlight last 10% in pseudotime
l2_breaks <- quantile(combined$l2_ps, probs = seq(0, 1, 0.1), na.rm = TRUE)
x_start <- l2_breaks[10]
x_end <- l2_breaks[11]

pdf("Ly2_GMP_scoring_Ly2_Trajectory_highlight_band.pdf", height = 4, width = 8)
ggplot(combined, aes(x = l2_ps, y = UCell_Score, color = group)) +
  annotate("rect", xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.5) +
  geom_point(alpha = 0.3, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, size = 1) +
  scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
  theme_classic() +
  labs(title = "Ly2/GMP scoring in Ly2 Trajectories (last 10% highlighted)",
       x = "l2_ps", y = "UCell Score")
dev.off()

pdf("Ly2_GMP_scoring_Ly2_Trajectory_wo_band.pdf", height = 4, width = 8)
ggplot(combined, aes(x = l2_ps, y = UCell_Score, color = group)) +
  geom_point(alpha = 0.3, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, size = 1) +
  scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
  theme_classic() +
  labs(title = "Ly2/GMP scoring in Ly2 Trajectories (last 10% highlighted)",
       x = "l2_ps", y = "UCell Score")
dev.off()
# GMP trajectory
ly_df <- FetchData(FL, vars = c("gmp_ps", "ly_in_FL_UCell")) %>%
  mutate(group = "FL: Ly2", UCell_Score = ly_in_FL_UCell) %>%
  select(-ly_in_FL_UCell)

gmp_df <- FetchData(FL, vars = c("gmp_ps", "gmp_in_FL_UCell")) %>%
  mutate(group = "FL: GMP", UCell_Score = gmp_in_FL_UCell) %>%
  select(-gmp_in_FL_UCell)

combined <- bind_rows(ly_df, gmp_df)

gmp_breaks <- quantile(combined$gmp_ps, probs = seq(0, 1, 0.1), na.rm = TRUE)
x_start <- gmp_breaks[10]
x_end <- gmp_breaks[11]

pdf("Ly2_GMP_scoring_GMP_Trajectory_highlight_band.pdf", height = 4, width = 8)
ggplot(combined, aes(x = gmp_ps, y = UCell_Score, color = group)) +
  annotate("rect", xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.5) +
  geom_point(alpha = 0.3, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, size = 1) +
  scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
  theme_classic() +
  labs(title = "Ly2/GMP scoring in GMP Trajectories (last 10% highlighted)",
       x = "gmp_ps", y = "UCell Score")
dev.off()

pdf("Ly2_GMP_scoring_GMP_Trajectory_wo_band.pdf", height = 4, width = 8)
ggplot(combined, aes(x = gmp_ps, y = UCell_Score, color = group)) +
  geom_point(alpha = 0.3, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x), se = TRUE, size = 1) +
  scale_color_manual(values = c("FL: Ly2" = '#793885', "FL: GMP" = "#56ace5")) +
  theme_classic() +
  labs(title = "Ly2/GMP scoring in GMP Trajectories (last 10% highlighted)",
       x = "gmp_ps", y = "UCell Score")
dev.off()



# Compare differences in last 10% bins
FL$l2_ps_cut <- cut(FL$l2_ps, breaks = 10)
FL$gmp_ps_cut <- cut(FL$gmp_ps, breaks = 10)

l2_df <- FL@meta.data[FL$l2_ps_cut == levels(FL$l2_ps_cut)[10], c("ly_in_FL_UCell", "gmp_in_FL_UCell","l2_ps")]
l2_df <- na.omit(l2_df)
l2_df$diff <- l2_df$ly_in_FL_UCell - l2_df$gmp_in_FL_UCell

gmp_df <- FL@meta.data[FL$gmp_ps_cut == levels(FL$gmp_ps_cut)[10], c("ly_in_FL_UCell", "gmp_in_FL_UCell")]
gmp_df <- na.omit(gmp_df)
gmp_df$diff <- gmp_df$gmp_in_FL_UCell - gmp_df$ly_in_FL_UCell

# Wilcoxon test
wc <- wilcox.test(l2_df$diff, gmp_df$diff)

# Summary stats
stats_df <- data.frame(
  l2_ps_mean_diff = mean(l2_df$diff),
  my_ps_mean_diff = mean(gmp_df$diff),
  l2_ps_median_diff = median(l2_df$diff),
  my_ps_median_diff = median(gmp_df$diff),
  wilcox_pval = wc$p.value
)
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

