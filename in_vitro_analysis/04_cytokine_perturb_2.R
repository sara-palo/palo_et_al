# Load required libraries
library(tidyverse)
library(readxl)
library(patchwork)
library(scales)
library(forcats)
library(ggplot2)
library(renv)
set.seed(123)

## set work directory
PRJ <- "H:/Lineageplasticity/LP_revise_singleMS5"
setwd(PRJ)
out_dir <- sprintf("output_%s", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cut_off <- 50


### Read data using read_excel
df_all2 <- read_excel("LP17_RO_regate3.xls") %>% mutate(Sample= ...1,
                                                     source = "LP17",
                                                     Group = case_when(str_detect(Sample, "BM") ~ paste0("BM_",condition),
                                                                       str_detect(Sample, "CB") ~ paste0("CB_", condition),
                                                                       str_detect(Sample, "FL") ~ paste0("FL_",condition)),
                                                     batch = case_when(str_detect(Sample, "BM") ~ paste0("BM_",condition),
                                                                       str_detect(Sample, "CB") ~ paste0("CB_",condition),
                                                                       str_detect(Sample, "FL") ~ paste0("FL_", condition)))
df_all3 <- read_excel("LP18.1_RO_regate3.xls") %>% mutate(Sample= ...1,
                                                       source = "LP18",
                                                       Group = case_when(str_detect(Sample, "BM") ~ paste0("BM_",condition),
                                                                         str_detect(Sample, "CB") ~ paste0("CB_", condition)),
                                                       batch = case_when(str_detect(Sample, "BM") ~ paste0("BM_",condition),
                                                                         str_detect(Sample, "CB") ~ paste0("CB_", condition)))
df_all4 <- read_excel("LP18.2_RO_regate3.xls") %>% mutate(Sample= ...1,
                                                       source = "LP18",
                                                       Group = case_when(str_detect(Sample, "FL1") ~ paste0("FL_",condition),
                                                                         str_detect(Sample, "FL2") ~ paste0("FL_", condition)),
                                                       batch = case_when(str_detect(Sample, "FL1") ~ paste0("FL1_",condition),
                                                                         str_detect(Sample, "FL2") ~ paste0("FL2_", condition)))%>%na.omit
df_all5 <- read_excel("LP19.1_RO_regate3.xls") %>% mutate(Sample= ...1,
                                                       source = "LP19",
                                                       Group = case_when(str_detect(Sample, "BM") ~ paste0("BM_",condition),
                                                                         str_detect(Sample, "FL1") ~ paste0("FL_",condition),
                                                                         str_detect(Sample, "FL2") ~ paste0("FL_", condition)),
                                                       batch = case_when(str_detect(Sample, "BM") ~ paste0("BM_",condition),
                                                                         str_detect(Sample, "FL1") ~ paste0("FL1_",condition),
                                                                         str_detect(Sample, "FL2") ~ paste0("FL2_", condition)))%>%  na.omit()
df_all6 <- read_excel("LP19.2_RO_regate3.xls") %>% mutate(Sample= ...1,
                                                       source = "LP19",
                                                       Group = case_when(str_detect(Sample, "CB1") ~ paste0("CB_",condition),
                                                                         str_detect(Sample, "CB2") ~ paste0("CB_", condition)),
                                                       batch = case_when(str_detect(Sample, "CB1") ~ paste0("CB1_",condition),
                                                                         str_detect(Sample, "CB2") ~ paste0("CB2_", condition)))%>%na.omit()
### combine them and exclude the control
df_all <- bind_rows(df_all2,df_all3,df_all4,df_all5,df_all6)
df_all <- df_all[df_all$...1 != "CB MNC_B1_B01.fcs", ] 

### Processing function
df_all<-mutate(df_all,
               B_pos =`B` >= cut_off,
               NK_pos = `NK` >= cut_off,
               Mono_pos = Mono >= cut_off,
               RBC_pos = RBC >= cut_off,
               Mk_pos = Mgk >= cut_off)
df_all <-  filter(df_all,!is.na(Group))

# Positivity category
df_all <- df_all %>%
  rowwise() %>%
  mutate(n_positive = sum(c_across(ends_with("_pos"))),
         PositivityCategory = case_when(n_positive == 0 ~ "non-defined",
                                        n_positive == 1 ~ "1 positive",
                                        n_positive == 2 ~ "2 positives",
                                        n_positive >= 3 ~ "3+ positives"))
df_all <- mutate(df_all,Group = factor(Group, levels = c("FL_st+15","FL_st-7", "CB_st+15", "CB_st-7", "BM_st+15", "BM_st-7")))
write.csv(df_all, file.path(out_dir, "Processed_Lineage_Positivity_LP17-19.csv"), row.names = FALSE)

lineage_pos_df <- df_all %>%
  select(source, batch, Group, B_pos, NK_pos, Mono_pos, RBC_pos, Mk_pos) %>%
  pivot_longer(
    cols = ends_with("_pos"),
    names_to = "Lineage",
    values_to = "positive"
  ) %>%
  mutate(
    Lineage = recode(
      Lineage,
      B_pos = "B",
      NK_pos = "NK",
      Mono_pos = "Mono",
      RBC_pos = "RBC",
      Mk_pos = "Mgk"
    )
  ) %>%
  group_by(source, batch, Group, Lineage) %>%
  summarise(
    n_total = n(),
    n_positive = sum(positive, na.rm = TRUE),
    positive_rate = n_positive / n_total * 100,
    .groups = "drop"
  )

write.csv(lineage_pos_df, file.path(out_dir, "Lineage_Positive_Rate_LP17-19.csv"), row.names = FALSE)

p_pos <- ggplot(lineage_pos_df, aes(x = Group, y = positive_rate, color = Group)) +
  geom_boxplot( width = 0.6,outlier.shape = NA,  linewidth =1) +
  geom_point(position = position_jitter(width = 0.08), size = 1.5,color = "black") +
  facet_wrap(~ Lineage, ncol = 5, scales = "free_y") +
  scale_color_manual( values = c(      "FL_st+15" = "#1c6494","FL_st-7"  = "#1c6494",
                                       "CB_st+15" = "darkred","CB_st-7"  = "darkred",
                                       "BM_st+15" = "#1b4b1f","BM_st-7"  = "#1b4b1f")) +
  theme_classic(base_size = 16) +
  labs( title = "Positive rate by batch",x = NULL,  y = "Positive rate (%)",fill = NULL,shape = "Experiment") +
  theme(  plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          axis.text.y = element_text(face = "bold", color = "black"),
          axis.title.y = element_text(face = "bold"), legend.title = element_blank()  )+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 100, by = 20))

ggsave(filename = file.path(out_dir, "Lineage_Positive_Rate_LP17-19_v3.pdf"), plot = p_pos, width = 18, height = 6, dpi = 300)
p_pos

p_pos <- ggplot(lineage_pos_df, aes(x = Group, y = positive_rate, color = Group)) +
  geom_boxplot( width = 0.6,outlier.shape = NA,  linewidth =1) +
  geom_text(aes(label = round(positive_rate, 1)),    position = position_jitter(width = 0.08, height = 0),    size = 3,    color = "black"  ) +
  facet_wrap(~ Lineage, ncol = 5, scales = "free_y") +
  scale_color_manual( values = c(      "FL_st+15" = "#1c6494","FL_st-7"  = "#1c6494",
                                       "CB_st+15" = "darkred","CB_st-7"  = "darkred",
                                       "BM_st+15" = "#1b4b1f","BM_st-7"  = "#1b4b1f")) +
  theme_classic(base_size = 16) +
  labs( title = "Positive rate by batch",x = NULL,  y = "Positive rate (%)",fill = NULL,shape = "Experiment") +
  theme(  plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          axis.text.y = element_text(face = "bold", color = "black"),
          axis.title.y = element_text(face = "bold"), legend.title = element_blank()  )+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 100, by = 20))
p_pos
ggsave(filename = file.path(out_dir, "Lineage_Positive_Rate_LP17-19_number_v3.pdf"), plot = p_pos, width = 18, height = 6, dpi = 300)


