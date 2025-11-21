library(tidyverse)
library(readxl)
library(ggplot2)
library(ggalluvial)
library(patchwork)

set.seed(123)

setwd("~/work/FL_adult_ATAC/2024_FL_ATAC")
out_dir <- "figs/"

cut_off <- 50

### Define plot colors
tissue_colors <- c(
  "FL" = "#1c6494",
  "CB" = "#8A0F37",
  "BM" = "#1b4b1f"
)

n_lin_colors <- c(
  "non-defined" = "#E0E9ED",
  "1" = "#A8C7DC",
  "2" = "#6F92AA",
  "3 or more" = "#3C5868")

B_mono_colors <- c(
  "B" = "#755D67",
  "B + Mono" = "#946562",
  "Mono" = "#C69999",
  "Neither" = "#E2E3E4"
)

### Read data
df_all1 <- read_excel("in_vitro_data/LP12_RO.xls") %>% mutate(source = "LP12")
df_all2 <- read_excel("in_vitro_data/LP13_RO.xls") %>% mutate(source = "LP13")
df_all3 <- read_excel("in_vitro_data/LP14_RO.xls") %>%
  mutate(source = case_when(str_detect(...1, "^FL[123]") ~ "LP14A",
                            str_detect(...1, "^FL[456]") ~ "LP14B",TRUE ~ "LP14"))

df_all <- (rbind(df_all1, df_all2, df_all3))
colnames(df_all)[1] <- "sample"

### exclude control, water, stroma samples.
df_all <- df_all[! df_all$sample %in% c("BM2_H9_H09_056.fcs", "CB2_G4_G04_044.fcs", "BM reanaysis_D3_D03_018.fcs", "Specimen_001_full stain CB_001.fcs"), ]

df_all <- mutate(
  df_all,
  Group = case_when(
    str_detect(sample, "^BM") ~ "BM",
    str_detect(sample, "^CB") ~ "CB",
    str_detect(sample, "^FL") ~ "FL"
  ),
  B_pos = `B(inclCD11b)` >= cut_off,
  NK_pos = `NKp46+(inclCD11b)` >= cut_off,
  Mono_pos = Mono >= cut_off,
  RBC_pos = RBC >= cut_off,
  Mk_pos = Mk >= cut_off
) %>% data.frame()

df_all[is.na(df_all$Group), "sample"]
df_all <-  filter(df_all,!is.na(Group))
df_all$Group <- factor(df_all$Group, levels = c("FL", "CB", "BM"))

df_all$npos <- apply(df_all[, grep("_pos", colnames(df_all))], 1, sum)
table(df_all$npos)

df_all$pos_cat <- df_all$npos

# calc number of positive lineages
df_all <- df_all %>%
  rowwise() %>%
  mutate(
    pos_cat = case_when(
      npos == 0 ~ "non-defined",
      npos == 1 ~ "1",
      npos == 2 ~ "2",
      npos >= 3 ~ "3 or more"
    )
  )

props <- table(df_all$Group, df_all$pos_cat)
props <- data.frame(props/rowSums(props))

plot_df <- data.frame(table(df_all$Group, df_all$pos_cat))
colnames(plot_df) <- c("group", "pos_cat", "n")
plot_df$prop <- round(props$Freq, digits = 2)

# version with text for double checking data
ggplot(plot_df, aes(x = group, y = n, fill = pos_cat)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = n), 
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual("Nr of lineages\nper clone", 
                    values = n_lin_colors) +
  labs(y = "# of clones", x = "Group") +
  theme_classic()
ggsave(paste0(out_dir, "additional_plots/lin_per_clone_tot.pdf"), width = 5)

ggplot(plot_df, aes(x = group, y = prop, fill = pos_cat)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = prop), 
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual("Nr of lineages\nper clone", 
                    values = n_lin_colors) +
  labs(y = "# of clones", x = "Group") +
  theme_classic()
ggsave(paste0(out_dir, "additional_plots/lin_per_clone_perc.pdf"), width = 5)


ggplot(plot_df, aes(x = group, y = n, fill = pos_cat)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, color = 'black', linewidth = 0.25) +
  scale_fill_manual("Nr of lineages\nper clone", values = n_lin_colors) +
  labs(y = "# of clones", x = "Group") +
  theme_classic()
ggsave(paste0(out_dir, "sup1/lin_per_clone_tot.pdf"), width = 5)

ggplot(plot_df, aes(x = group, y = prop, fill = pos_cat)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = 'black', linewidth = 0.25) +
  scale_fill_manual("Nr of lineages\nper clone", 
                    values = n_lin_colors) +
  labs(y = "# of clones", x = "Group") +
  theme_classic()
ggsave(paste0(out_dir, "fig1/lin_per_clone_perc.pdf"), width = 5)

df_all$Group_source <- paste0(df_all$Group, "_", df_all$source)

df_all <- data.frame(df_all)

# sum the number of clones positive for each lineage per sample
df_list <- split(df_all, df_all$Group_source, drop = T)
counts <- lapply(df_list, function(x) apply(x[, grep("_pos", colnames(x))], 2, sum))
count_df <- data.frame(counts)
count_df <- data.frame(t(count_df))

# divide by number of total wells per sample
count_df_perc <- count_df/unlist(lapply(df_list, nrow), use.names = F)

count_df$tissue <- str_sub(rownames(count_df), 1, 2) 
count_df$tissue <- factor(count_df$tissue, levels = c("FL", "CB", "BM"))

count_df_perc$tissue <- count_df$tissue 

ggplot(count_df_perc, aes(x = tissue, y = B_pos, color = tissue)) + geom_boxplot(fatten = NULL, width = 0.6) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.6) +
  geom_point(color = "black", size = 1.4) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1)
ggsave(paste0(out_dir, "fig1/B_perc_of_clones.pdf"), width = 4)

ggplot(count_df_perc, aes(x = tissue, y = Mono_pos, color = tissue)) + geom_boxplot(fatten = NULL, width = 0.6) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.6) +
  geom_point(color = "black", size = 1.4) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1)
ggsave(paste0(out_dir, "fig1/Mono_perc_of_clones.pdf"), width = 4)

ggplot(count_df_perc, aes(x = tissue, y = Mk_pos, color = tissue)) + geom_boxplot(fatten = NULL, width = 0.6) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.6) +
  geom_point(color = "black", size = 1.4) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1)
ggsave(paste0(out_dir, "fig1/Mk_perc_of_clones.pdf"), width = 4)

ggplot(count_df_perc, aes(x = tissue, y = RBC_pos, color = tissue)) + geom_boxplot(fatten = NULL, width = 0.6) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.6) +
  geom_point(color = "black", size = 1.4) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1)
ggsave(paste0(out_dir, "fig1/RBC_perc_of_clones.pdf"), width = 4)

ggplot(count_df_perc, aes(x = tissue, y = NK_pos, color = tissue)) + geom_boxplot(fatten = NULL, width = 0.6) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.6) +
  geom_point(color = "black", size = 1.4) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1)
ggsave(paste0(out_dir, "fig1/NK_perc_of_clones.pdf"), width = 4)

# plot lineage output % with numbers for data checking
b <- ggplot(count_df_perc, aes(x = tissue, y = B_pos, color = tissue)) +
  geom_point(color = "black", size = 0.7) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1) +
  geom_text(aes(label = round(B_pos, 3)), position = position_jitter(width = 0.4, height = 0.07)) 

m <- ggplot(count_df_perc, aes(x = tissue, y = Mono_pos, color = tissue)) +
  geom_point(color = "black", size = 0.7) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1) +
  geom_text(aes(label = round(Mono_pos, 3)), position = position_jitter(width = 0.4, height = 0.07)) 

mk <- ggplot(count_df_perc, aes(x = tissue, y = Mk_pos, color = tissue)) +
  geom_point(color = "black", size = 0.7) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1) +
  geom_text(aes(label = round(Mk_pos, 3)), position = position_jitter(width = 0.4, height = 0.07)) 

rbc <- ggplot(count_df_perc, aes(x = tissue, y = RBC_pos, color = tissue)) +
  geom_point(color = "black", size = 0.7) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1) +
  geom_text(aes(label = round(RBC_pos, 3)), position = position_jitter(width = 0.4, height = 0.07)) 

nk <- ggplot(count_df_perc, aes(x = tissue, y = NK_pos, color = tissue)) +
  geom_point(color = "black", size = 0.7) + scale_color_manual(values = tissue_colors) + theme_classic() + ylim(0,1) +
  geom_text(aes(label = round(NK_pos, 3)), position = position_jitter(width = 0.4, height = 0.07)) 

pdf(paste0(out_dir, "additional_plots/perc_of_clones_numbers_FL_CB_BM.pdf"), width = 12, height = 8)
print(b + m + mk + rbc + nk)
dev.off()


# plot B/Mono output
df_all$B_mono <- paste0(df_all$B_pos, df_all$Mono_pos)

vars <- c("TRUETRUE" = "B + Mono", "TRUEFALSE" = "B", "FALSETRUE" = "Mono", "FALSEFALSE" = "Neither")

df_all$B_mono <- vars[match(df_all$B_mono, names(vars))]

plot_df <- data.frame(table(df_all$Group, df_all$B_mono))
colnames(plot_df) <- c("group", "B_mono", "n")
plot_df$group <- factor(plot_df$group, levels = c("FL", "CB", "BM"))

ggplot(plot_df, aes(x = group, y = n, fill = B_mono)) + geom_col(width = 0.5) + theme_classic() + 
  scale_fill_manual(values = B_mono_colors)
ggsave(paste0(out_dir, "sup1/B_mono_counts.pdf"), width = 4.5)

ggplot(data = plot_df,
       aes(x = group, y = n, alluvium = B_mono, 
           fill = B_mono, stratum = B_mono)) +
  geom_alluvium() +
  geom_stratum(width = 0.5) + theme_classic() + 
  scale_fill_manual(values = B_mono_colors)
ggsave(paste0(out_dir, "sup1/B_mono_counts_river.pdf"), width = 4.5)

ggplot(data = plot_df,
       aes(x = group, y = n, alluvium = B_mono, 
           fill = B_mono, stratum = B_mono)) +
  geom_alluvium() +
  geom_stratum(width = 0.5) + theme_classic() + 
  scale_fill_manual(values = B_mono_colors) + 
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5)
ggsave(paste0(out_dir, "additional_plots/B_mono_counts_river_numbers.pdf"), width = 4.5)


# divide by total number of colonies per sample group
plot_df$tot <- table(df_all$Group)[match(plot_df$group, names(table(df_all$Group)))]
plot_df$perc <- plot_df$n/plot_df$tot
plot_df$perc <- as.numeric(plot_df$perc)

ggplot(data = plot_df,
       aes(x = group, y = perc, alluvium = B_mono, 
           fill = B_mono, stratum = B_mono)) +
  geom_alluvium() +
  geom_stratum(width = 0.5) + theme_classic() + 
  scale_fill_manual(values = B_mono_colors)
ggsave(paste0(out_dir, "fig1/B_mono_perc_river.pdf"), width = 4.5)

ggplot(data = plot_df,
       aes(x = group, y = perc, alluvium = B_mono, 
           fill = B_mono, stratum = B_mono)) +
  geom_alluvium() +
  geom_stratum(width = 0.5) + theme_classic() + 
  scale_fill_manual(values = B_mono_colors) + 
  geom_text(aes(label = round(perc,3)), position = position_stack(vjust = 0.5), size = 3.5)
ggsave(paste0(out_dir, "additional_plots/B_mono_perc_river_numbers.pdf"), width = 4.5)


plot_df2 <- data.frame(table(df_all$Group_source, df_all$B_mono))
colnames(plot_df2) <- c("group", "B_mono", "n")
plot_df2$group <- factor(plot_df2$group, levels = c("FL_LP13", "FL_LP14A", "FL_LP14B", "CB_LP12", "CB_LP13", "CB_LP14","BM_LP12", "BM_LP13"))

ggplot(plot_df2, aes(x = group, y = n, fill = B_mono)) + geom_col(width = 0.5) + theme_classic() + 
  scale_fill_manual(values = B_mono_colors) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+ 
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5)
ggsave(paste0(out_dir, "additional_plots/B_mono_counts_persample.pdf"), width = 6, height = 4)


group_labels <- c(
  "FL_LP14A" = "FL_7pcw",
  "FL_LP14B" = "FL_8pcw",
  "FL_LP12" = "FL_9pcw",
  "FL_LP13" = "FL_8.5pcw",
  "CB_LP12" = "CB_#1",
  "CB_LP13" = "CB_#2",
  "CB_LP14" = "CB_#3",
  "BM_LP12" = "BM_#1",
  "BM_LP13" = "BM_#2"
)

df_all$sample_info <- group_labels[match(df_all$Group_source, names(group_labels))]

write.csv(df_all, "csv_files/sup/in_vitro_processed_data.csv")






