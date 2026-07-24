# Load required libraries
library(tidyverse)
library(readxl)
library(patchwork)
library(scales)
library(forcats)
library(ggplot2)
library(renv)
library(rstatix)
library(ggpubr)
set.seed(123)
group_cols <- c(  "FL" = "#1c6494", "CB" = "darkred","BM" = "#1b4b1f")
### set work directory
PRJ <- "H:/Lineageplasticity/LP revise_Myeloid"
setwd(PRJ)
out_dir <- sprintf("output_%s", format(Sys.Date(), "%Y%m%d"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

### Read data using read_excel
mye_all <- read_excel("Myeloid_7AAD_modified.xls") %>%
  mutate(source = case_when(str_detect(...1, "#18") ~ "18_mye",
                            str_detect(...1, "#17") ~ "17_mye",
                            str_detect(...1, "#19") ~ "19_mye")) %>%
  mutate(Group = case_when(str_detect(...1, "BM") ~ "BM",
                           str_detect(...1, "CB") ~ "CB",
                           str_detect(...1, "FL") ~ "FL")) %>%
  mutate(amount=210,input=50) %>%  na.omit()

mye_all$Group <- factor(mye_all$Group, levels = c("FL","CB","BM"))

###Percentage
mye_all$pct_mono <- ((mye_all$mono) / (mye_all$CD45)*100)

summary_pct_mono <- mye_all %>%
  group_by(Group) %>%
  summarise(mean_val = mean(pct_mono, na.rm = TRUE),sd_val   = sd(pct_mono, na.rm = TRUE),.groups = "drop")


###Krusukal.test
kruskal_mye <- kruskal.test(pct_mono ~ Group, data = mye_all) %>% tidy() #unsisgnificant
write.csv(kruskal_mye, file.path(out_dir, "Krusukal_Mye_results.csv"), row.names = FALSE)


p <- ggplot(mye_all, aes(x = Group, y = pct_mono,color=Group)) +
  geom_boxplot(linewidth = 1.5, width = 0.6 ) +
  geom_point(data = mye_all,aes(x = Group, y = pct_mono),color = "black",
             position = position_jitter(width = 0.05)) +
  labs(title = "Mono(%) in CD235a-CD45+ cells", y = "% in CD45+", x = NULL) +
  theme_classic(base_size = 18) +
  scale_color_manual(values=group_cols)+
  theme(plot.title = element_text(face = "bold", hjust = 0.5),axis.text = element_text(face = "bold"),axis.title.y = element_text(face = "bold"))
p
ggsave(filename = file.path(out_dir, "Mono_pct(CD45)_v2.pdf"), plot = p, width = 4, height = 5, dpi = 300)





