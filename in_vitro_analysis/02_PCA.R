library(tidyverse)
library(ggplot2)
library(patchwork)
set.seed(123)

setwd("~/work/FL_adult_ATAC/2024_FL_ATAC")
out_dir <- "figs/"

lineage_colors <- c(
  "B" = "#785EF0",
  "NK" = "lightgreen",
  "Mono" = "#648FFF",
  "Mk" = "lightpink2",
  "RBC" = "red"
)

tissue_colors <- c(
  "FL" = "#1c6494",
  "CB" = "#8A0F37",
  "BM" = "#1b4b1f"
)

df_all <- read.csv("csv_files/sup/in_vitro_processed_data.csv", row.names = 1)
rownames(df_all) <- paste0(df_all$source, "_", df_all$sample)

counts <- df_all[, c("RBC", "Mono", "Mk", "B.inclCD11b.", "NKp46..inclCD11b.")]
counts <- rename(counts, "NK" = "NKp46..inclCD11b.", "B" = "B.inclCD11b.")

# normalize counts of each lineage by dividing by number of live cells
norm_counts <- counts/df_all$Alive

#PCA
pca_result <- prcomp(norm_counts, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca_result$x)
pca_df <- cbind(pca_df, norm_counts[, c("RBC", "Mono", "Mk", "B", "NK")], df_all[, c("source", "Group")])

# reverse PC2 to give plots the same orientation as before
pca_df$PC2 <- -pca_df$PC2

xlims <- c(-3, 4)
ylims <- c(-4, 7)

B_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = B)) +
  geom_point(size = 5, shape = 21, color = "black", stroke = 0.1) +
  scale_fill_gradient(low = "white", high = lineage_colors["B"]) +
  theme_classic() + coord_fixed() + xlim(xlims) + ylim(ylims) +
  labs(title = "B cell ratio", fill = "B/live cells")

NK_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = NK)) +
  geom_point(size = 5, shape = 21, color =  "black", stroke = 0.1) +
  scale_fill_gradient(low = "white", high = lineage_colors["NK"]) +
  theme_classic() + coord_fixed() + xlim(xlims) + ylim(ylims) +
  labs(title = "NK cell ratio", fill = "NK/live cells")

Mono_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Mono)) +
  geom_point(size = 5, shape = 21, color =  "black", stroke = 0.1) +
  scale_fill_gradient(low = "white", high = lineage_colors["Mono"]) +
  theme_classic() + coord_fixed() + xlim(xlims) + ylim(ylims) +
  labs(title = "Mono cell ratio", fill = "Mono/live cells")

Ery_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = RBC)) +
  geom_point(size = 5, shape = 21, color = "black", stroke = 0.1) +
  scale_fill_gradient(low = "white", high = lineage_colors["RBC"]) +
  theme_classic() + coord_fixed() + xlim(xlims) + ylim(ylims) +
  labs(title = "Ery cell ratio", fill = "Ery/live cells")

Mk_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Mk)) +
  geom_point(size = 5, shape = 21, color =  "black", stroke = 0.1) +
  scale_fill_gradient(low = "white", high = lineage_colors["Mk"]) +
  theme_classic() + coord_fixed() + xlim(xlims) + ylim(ylims) +
  labs(title = "Mk cell ratio", fill = "Mk/live cells")

pdf(paste0(out_dir, "sup1/PCA_plot_all_lineage_FL_CB_BM.pdf"), width = 10, height = 6)
print(B_plot + Mono_plot + NK_plot + Ery_plot + Mk_plot)
dev.off()

pdf(paste0(out_dir, "additional_plots/PCA_plot_all_lineage_FL_CB_BM_ntot.pdf"), width = 10, height = 6)
print(B_plot + Mono_plot + NK_plot + Ery_plot + Mk_plot) + labs(caption = paste0("N clones = ", nrow(pca_df)))
dev.off()

FL_plot <- ggplot() +
  geom_point(data = pca_df, aes(PC1, PC2),color = "grey80", size = 5, alpha = 0.6, stroke = 0) +
  stat_density_2d(data = subset(pca_df, Group == "FL"), aes(PC1, PC2),
                  geom = "polygon", fill = tissue_colors["FL"], alpha = 0.20, color = NA,
                  bins = 20, contour_var = "ndensity",h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6)) +
  geom_density_2d(data = subset(pca_df, Group == "FL"), aes(PC1, PC2),
                  color = "grey20", linewidth = 0.25, bins = 10) +
  geom_point(data = subset(pca_df, Group == "FL"),aes(PC1, PC2),
             size = 2, shape = 21, color = "black", fill = tissue_colors["FL"], stroke = 0.1) +
  coord_fixed() + xlim(xlims) + ylim(ylims) +
  theme_classic()+labs(title = "FL", x = "PC1", y = "PC2")

CB_plot <- ggplot() +
  geom_point(data = pca_df, aes(PC1, PC2),color = "grey80", size = 5, alpha = 0.6, stroke = 0) +
  stat_density_2d(data = subset(pca_df, Group == "CB"), aes(PC1, PC2),
                  geom = "polygon", fill = tissue_colors["CB"], alpha = 0.20, color = NA,
                  bins = 20, contour_var = "ndensity",h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6)) +
  geom_density_2d(data = subset(pca_df, Group == "CB"), aes(PC1, PC2),
                  color = "grey20", linewidth = 0.25, bins = 10) +
  geom_point(data = subset(pca_df, Group == "CB"),aes(PC1, PC2),
             size = 2, shape = 21, color = "black", fill = tissue_colors["CB"], stroke = 0.1) +
  coord_fixed() + xlim(xlims) + ylim(ylims) +
  theme_classic()+labs(title = "CB", x = "PC1", y = "PC2")

BM_plot <- ggplot() +
  geom_point(data = pca_df, aes(PC1, PC2),color = "grey80", size = 5, alpha = 0.6, stroke = 0) +
  stat_density_2d(data = subset(pca_df, Group == "BM"), aes(PC1, PC2),
                  geom = "polygon", fill = tissue_colors["BM"], alpha = 0.20, color = NA,
                  bins = 20, contour_var = "ndensity",h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6)) +
  geom_density_2d(data = subset(pca_df, Group == "BM"), aes(PC1, PC2),
                  color = "grey20", linewidth = 0.25, bins = 10) +
  geom_point(data = subset(pca_df, Group == "BM"),aes(PC1, PC2),
             size = 2, shape = 21, color = "black", fill = tissue_colors["BM"], stroke = 0.1) +
  coord_fixed() + xlim(xlims) + ylim(ylims) +
  theme_classic()+labs(title = "BM", x = "PC1", y = "PC2")

pdf(paste0(out_dir, "sup1/PCA_plot_per_tissue.pdf"), width = 11, height = 4)
print(FL_plot + CB_plot + BM_plot)
dev.off()

pdf(paste0(out_dir, "additional_plots/PCA_plot_per_tissue.pdf"), width = 10, height = 6)
(FL_plot + labs(subtitle = paste0("N clones = ", length(which(pca_df$Group == "FL"))))) + 
  (CB_plot + labs(subtitle = paste0("N clones = ", length(which(pca_df$Group == "CB"))))) +
  (BM_plot + labs(subtitle = paste0("N clones = ", length(which(pca_df$Group == "BM")))))
dev.off()


#bi plot for B and mono
plot_df <- df_all[, c("B.inclCD11b.", "Mono", "Group")] %>%
  mutate(B.inclCD11b. = ifelse(B.inclCD11b. == 0, 0.5, B.inclCD11b.),
         Mono = ifelse(Mono == 0, 0.5, Mono))
plot_df$Group <- factor(plot_df$Group, levels = c("FL", "CB", "BM"))


ggplot(plot_df, aes(x = Mono, y = B.inclCD11b.)) +
  stat_density_2d(aes(fill = Group, alpha = after_stat(level)),
                  geom = "polygon", color = NA,contour = TRUE,bins =10) +
  geom_point(aes(color = Group),size = 1,stroke = 1) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  scale_x_log10(limits = c(0.1, 250000)) +
  scale_y_log10(limits = c(0.1, 250000)) +
  coord_fixed()+
  theme_classic(base_size = 14) +
  labs(x = "Mono", y = "B",
       fill = "Group", color = "Group") +
  guides(alpha = "none") + 
  scale_color_manual(values = tissue_colors) +
  scale_fill_manual(values = tissue_colors) 
ggsave(paste0(out_dir, "additional_plots/B_mono_output_combined.pdf"))


# B vs Mono by group
ggplot(plot_df, aes(x = Mono, y = B.inclCD11b.)) +
  facet_wrap(~ Group, nrow = 1, scales = "fixed", axes = "all") +
  geom_point(aes(color = Group), size = 1, stroke = 1) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  scale_x_log10(limits = c(0.1, 250000)) +
  scale_y_log10(limits = c(0.1, 250000)) +
  coord_fixed() +
  scale_color_manual(values = tissue_colors) +
  theme_classic(base_size = 14) +
  labs(x = "Mono", y = "B",
       fill = "Group", color = "Group")
ggsave(paste0(out_dir, "fig1/B_mono_output_sep.pdf"), width = 10, height = 6)

ggplot(plot_df, aes(x = Mono, y = B.inclCD11b.)) +
  facet_wrap(~ Group, nrow = 1, scales = "fixed", axes = "all") +
  geom_point(aes(color = Group), size = 1, stroke = 1) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  scale_x_log10(limits = c(0.1, 250000)) +
  scale_y_log10(limits = c(0.1, 250000)) +
  coord_fixed() +
  scale_color_manual(values = tissue_colors) +
  theme_classic(base_size = 14) +
  labs(x = "Mono", y = "B",
       fill = "Group", color = "Group", 
       title = paste0("FL clones = ", length(which(plot_df$Group == "FL")), 
                      ", CB clones = ", length(which(plot_df$Group == "CB")), 
                      ", BM clones = ", length(which(plot_df$Group == "BM"))))
ggsave(paste0(out_dir, "additional_plots/B_mono_output_sep_counts.pdf"), width = 10, height = 6)

