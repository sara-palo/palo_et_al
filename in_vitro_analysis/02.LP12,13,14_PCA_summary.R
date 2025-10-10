library(tidyverse)
library(readxl)
library(ggplot2)
library(viridis)
library(patchwork)
set.seed(123)


PRJ <- "C:/Users/keiki/OneDrive/デスクトップ/Sara_project"
setwd(PRJ)
out_dir <- sprintf("output_%s", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir

lineage_colors <- c("B" = "#785EF0","NK" = "lightgreen","Mono" = "#648FFF",
                    "Mk" = "lightpink2","RBC" = "red")
origin_colors <- c("BM" = '#1b4b1f',"CB" = "#DE4968FF","FL" = '#1c6494')


df_combined <- as.data.frame(df_all) #%>%
  #filter(Group %in% c("FL", "BM")) # filter CB sample if necessary

df_combined <- mutate(df_combined,rowname = paste(source, Sample, sep = "_"))
rownames(df_combined) <- as.vector(df_combined$rowname)

keep_cols <- setdiff(colnames(df_combined), 
                     c("B", "NK", "CD71", "Sample", "source", "PositivityCombo", 
                       "Group_source", "n_positive", "PositivityCategory",
                       "B_pos", "NK_pos", "Mono_pos", "RBC_pos", "Mk_pos","...1","Group"))
combined_for_pca <- df_combined[, keep_cols]
alive_values <- combined_for_pca$Alive
    normalized <- combined_for_pca
    rownames(normalized) <-combined_for_pca$rowname
    normalized <- normalized[, colnames(normalized) != "rowname"]
    normalized <- normalized %>%
      mutate(across(-Alive, ~ .x / Alive)) %>%  # Aliveで割る
      select(-Alive)                            # Alive列を削除


#PCA
pca_result <- prcomp(normalized, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$origin <- case_when(grepl("_(BM)",rownames(pca_df)) ~ "BM",
                           grepl("_(CB)",rownames(pca_df)) ~ "CB",
                           grepl("_(FL)|LP14A_|LP14B_", rownames(pca_df)) ~ "FL")
pca_df <- pca_df %>% filter(!is.na(origin))
    pca_df$B_ratio   <- normalized[rownames(pca_df), "B(inclCD11b)"]
    pca_df$NK_ratio   <- normalized[rownames(pca_df), "NKp46+(inclCD11b)"]
    pca_df$Mono_ratio <- normalized[rownames(pca_df), "Mono"]
    pca_df$RBC_ratio  <- normalized[rownames(pca_df), "RBC"]
    pca_df$Mk_ratio  <- normalized[rownames(pca_df), "Mk"]


xlims <- c(-3,4)# xlims <- range(pca_df$PC1, na.rm = TRUE) #-2.56 to 3.61
ylims <- c(-4,7)# ylims <- range(pca_df$PC2, na.rm = TRUE) #-3.32 to 6.63

    B_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = B_ratio)) +
              geom_point(size = 5, shape = 21, color = "black", stroke = 0.1) +
              scale_fill_gradient(low = "white", high = lineage_colors["B"]) +
              theme_classic()+coord_fixed()+xlim(xlims)+ylim(ylims)+
              labs(title = "B cell ratio", fill = "B / Alive")
              
    NK_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = NK_ratio)) +
               geom_point(size = 5, shape = 21, color =  "black", stroke = 0.1) +
               scale_fill_gradient(low = "white", high = lineage_colors["NK"]) +
               theme_classic()+coord_fixed()+xlim(xlims) +ylim(ylims) +
              labs(title = "NK cell ratio", fill = "NK / Alive")
              
    Mono_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Mono_ratio)) +
                geom_point(size = 5, shape = 21, color =  "black", stroke = 0.1) +
                scale_fill_gradient(low = "white", high = lineage_colors["Mono"]) +
                theme_classic() +
                coord_fixed() +   xlim(xlims) + ylim(ylims) +
                labs(title = "Mono cell ratio", fill = "Mono / Alive")
              
     Ery_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = RBC_ratio)) +
                geom_point(size = 5, shape = 21, color = "black", stroke = 0.1) +
                scale_fill_gradient(low = "white", high = lineage_colors["RBC"]) +
                theme_classic() +coord_fixed() +xlim(xlims) +ylim(ylims) +
                labs(title = "Ery cell ratio", fill = "Ery / Alive")
              
     Mk_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = Mk_ratio)) +
                geom_point(size = 5, shape = 21, color =  "black", stroke = 0.1) +
                scale_fill_gradient(low = "white", high = lineage_colors["Mk"]) +
                theme_classic() +coord_fixed() +xlim(xlims)+ ylim(ylims) +
                labs(title = "Mk cell ratio", fill = "Mk / Alive")


pdf(paste0(out_dir,"/PCA_plot_all_lineage_FL_CB_BM_bigger_black.pdf"), width = 10, height = 6)
print(B_plot + NK_plot + Mono_plot + Ery_plot + Mk_plot)
dev.off()

              # p8_FL <- ggplot(pca_df %>% filter(origin == "FL"), aes(PC1, PC2)) +
              #   stat_density_2d(geom = "polygon", fill ='#1c6494', alpha = 0.20, color = NA,
              #     bins = 20,contour_var = "ndensity",h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6)) +
              #   geom_density_2d(color = "grey20", linewidth = 0.25, bins = 10) +
              #   coord_fixed() + xlim(xlims) + ylim(ylims) + theme_classic() +
              #   labs(title = "FL", x = "PC1", y = "PC2")
              # 
              # 
              # p8_BM <- ggplot(pca_df %>% filter(origin == "BM"), aes(PC1, PC2)) +
              #   stat_density_2d(
              #     geom = "polygon", fill ='#1b4b1f', alpha = 0.25, color = NA,
              #     bins = 20,
              #     contour_var = "ndensity",
              #     h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6)  # ← df ではなく pca_df
              #   ) +
              #   geom_density_2d(color = "grey20", linewidth = 0.3, bins = 10) +
              # #  geom_point(size = 4, shape = 21, color = "black", fill = '#1b4b1f', stroke = 0.1) +
              #   coord_fixed() + xlim(xlims) + ylim(ylims) + theme_classic() +
              #   labs(title = "BM", x = "PC1", y = "PC2")
              # 
              # 
              # 
              # p8_CB <- ggplot(pca_df %>% filter(origin == "CB"), aes(PC1, PC2)) +
              #   stat_density_2d(
              #     geom = "polygon", fill ="#DE4968FF", alpha = 0.25, color = NA,
              #     bins = 20,
              #     contour_var = "ndensity",
              #     h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6)  # ← df ではなく pca_df
              #   ) +
              #   geom_density_2d(color = "grey20", linewidth = 0.3, bins = 10) +
              # #geom_point(size = 4, shape = 21, color = "brown", fill = '#ff6db6', stroke = 0.1) +
              #   coord_fixed() + xlim(xlims) + ylim(ylims) + theme_classic() +
              #   labs(title = "CB", x = "PC1", y = "PC2")
              # 
              
              p8_FL2 <- ggplot() +
                geom_point(data = pca_df, aes(PC1, PC2),color = "grey80", size = 5, alpha = 0.6, stroke = 0)+
                stat_density_2d(data = dplyr::filter(pca_df, origin == "FL"),aes(PC1, PC2),
                  geom = "polygon", fill = "#1c6494", alpha = 0.20, color = NA,
                  bins = 20, contour_var = "ndensity",h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6))+
                geom_point(data = dplyr::filter(pca_df, origin == "FL"),aes(PC1, PC2),
                           size = 2, shape = 21, color = "black", fill = "#1c6494", stroke = 0.1) +
                geom_density_2d(data = dplyr::filter(pca_df, origin == "FL"), aes(PC1, PC2),
                                color = "grey20", linewidth = 0.25, bins = 10) +
                coord_fixed() + xlim(xlims) + ylim(ylims) +
                theme_classic()+labs(title = "FL", x = "PC1", y = "PC2")
              
              p8_BM2 <- ggplot() +
                geom_point(data = pca_df, aes(PC1, PC2),color = "grey80", size = 5, alpha = 0.6, stroke = 0)+
                stat_density_2d(data = dplyr::filter(pca_df, origin == "BM"),aes(PC1, PC2),
                  geom = "polygon", fill = '#1b4b1f', alpha = 0.1, color = NA,
                  bins = 20, contour_var = "ndensity",h= c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6)    ) +
                geom_point(data = dplyr::filter(pca_df, origin == "BM"),aes(PC1, PC2),
                           size = 2, shape = 21, color = "black", fill = '#1b4b1f', stroke = 0.1 ) +
                geom_density_2d(data = dplyr::filter(pca_df, origin == "BM"),aes(PC1, PC2),
                                color = "grey20", linewidth = 0.25, bins = 10) +
                coord_fixed() + xlim(xlims) + ylim(ylims) +theme_classic() +labs(title = "BM", x = "PC1", y = "PC2")
            
              p8_CB2 <- ggplot() +
                geom_point(data = pca_df, aes(PC1, PC2),color = "grey80", size = 5, alpha = 0.6, stroke = 0) +
                stat_density_2d(data = dplyr::filter(pca_df, origin == "CB"),aes(PC1, PC2),
                  geom = "polygon", fill = "#DE4968FF", alpha = 0.1, color = NA,bins = 20, contour_var = "ndensity",
                  h = c(sd(pca_df$PC1)*0.6, sd(pca_df$PC2)*0.6) ) +
                geom_point(data = dplyr::filter(pca_df, origin == "CB"),aes(PC1, PC2),
                           size = 2, shape = 21, color = "black", fill = "#DE4968FF", stroke = 0.1) +
                geom_density_2d(data = dplyr::filter(pca_df, origin == "CB"),aes(PC1, PC2),
                                color = "grey20", linewidth = 0.25, bins = 10 ) +
                coord_fixed() + xlim(xlims) + ylim(ylims)+ theme_classic() +
                labs(title = "CB", x = "PC1", y = "PC2")
              
              

# pdf(paste0(out_dir,"/PCA_marker_contour_plot_by_organ.pdf"), width = 11, height = 4)
# print(p8_FL+p8_CB+p8_BM)
# dev.off()
pdf(paste0(out_dir,"/PCA_marker_contour_plot_by_organ_v2.pdf"), width = 11, height = 4)
print(p8_FL2+p8_CB2+p8_BM2)
dev.off()


#bi plot for B and mono
combined_clean <- df_all[, c("B(inclCD11b)", "Mono", "Group")] %>%
  mutate(`B(inclCD11b)` = ifelse(`B(inclCD11b)` == 0, 0.5, `B(inclCD11b)`),
    Mono = ifelse(Mono == 0, 0.5, Mono) )
xlim_vals <- c(0.01, 500000)
ylim_vals <- c(0.01, 500000)

p <-ggplot(combined_clean, aes(x = `B(inclCD11b)`, y = Mono)) +
  stat_density_2d(aes(fill = Group, alpha = after_stat(level)),
    geom = "polygon", color = NA,contour = TRUE,bins =10) +
  geom_point(aes(color = Group),size = 1,stroke = 1) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  scale_x_log10(limits = xlim_vals) +
  scale_y_log10(limits = ylim_vals) +
  coord_fixed()+
  theme_classic(base_size = 14) +
  labs(title = "Filled contour per group (adjusted bins, fixed edges)",
    x = "B (log10), 0 was converted to 0.5", y = "Mono (log10), 0 was converted to 0.5",
    fill = "Group", color = "Group") +
  guides(alpha = "none") + 
  scale_color_manual(values = origin_colors) +
  scale_fill_manual(values = origin_colors) 
p

pdf(paste0(out_dir,"/B_vs_Mono_contour_fixed_incl0count.pdf"), width = 5, height = 5)
p +
  stat_density_2d(aes(color = Group),geom = "contour", bins = 10, linewidth = 0.5, show.legend = FALSE )
dev.off()

pdf(paste0(out_dir,"/B_vs_Mono_contour_fixed_incl0count_no line.pdf"), width = 5, height = 5)
p
dev.off()


# B vs Mono by group
p_by_group <- ggplot(combined_clean, aes(x = Mono, y = `B(inclCD11b)`)) +
  facet_wrap(~ Group, nrow = 1, scales = "fixed") +
  stat_density_2d(
    aes(fill = Group), alpha = 0.15,
    geom = "polygon", color = NA, bins = 10, contour = TRUE ) +
  geom_density_2d(aes(color = Group), bins = 5, linewidth = 0.5, show.legend = FALSE) +
  geom_point(aes(color = Group), size = 1, stroke = 1) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
  scale_x_log10(limits = xlim_vals) +
  scale_y_log10(limits = ylim_vals) +
  coord_fixed() +
  scale_fill_manual(values = origin_colors) +
  scale_color_manual(values = origin_colors) +
  theme_classic(base_size = 14) +
  labs(title = "B vs Mono by group", x = "Mono (log10), value 0 was replaced by 0.5",
    y = "B (log10), value 0 was replaced by 0.5")
p_by_group


ggsave(file.path(out_dir, "/B_vs_Mono_contour_by_group.pdf"),
       p_by_group, width = 12, height = 4, dpi = 300)




# p <- ggplot(combined_clean, aes(x = `B(inclCD11b)`, y = Mono)) +
#   annotate("rect", xmin = 10, xmax = 50, ymin = -Inf, ymax = Inf,
#            fill = "gray80", alpha = 0.5) +
#   annotate("rect", xmin = -Inf, xmax = Inf, ymin = 10, ymax = 50,
#            fill = "gray80", alpha = 0.5) +
#   geom_vline(xintercept = 50, linetype = "dashed", color = "black") +
#   geom_hline(yintercept = 50, linetype = "dashed", color = "black") +
#   geom_point(
#     aes(x = `B(inclCD11b)`, y = Mono, color = Group),
#     size = 4,
#     stroke = 1,
#     alpha = 1
#   )+
# 
#   stat_density_2d(
#     aes(fill = Group, alpha = after_stat(level)),
#     geom = "polygon", color = NA,
#     contour = TRUE,
#     bins =10,
#     xlim = xlim_vals,
#     ylim = ylim_vals
#   ) +
#   scale_fill_manual(values = origin_colors) +
#   scale_color_manual(values = origin_colors) +
#   scale_alpha(range = c(0.1, 0.5), guide = "none") +
#   scale_x_log10(limits = xlim_vals) +
#   scale_y_log10(limits = ylim_vals) +
#   coord_fixed()+
#   theme_classic(base_size = 14) +
#   labs(
#     title = "Filled contour per group (adjusted bins, fixed edges)",
#     x = "B (log10), 0 was converted to 0.1", y = "Mono (log10), 0 was converted to 0.1",
#     fill = "Group", color = "Group"
#   ) +
#   theme(
#     legend.position = "right",
#     axis.text = element_text(color = "black"),
#     axis.line = element_line(color = "black")
#   )

