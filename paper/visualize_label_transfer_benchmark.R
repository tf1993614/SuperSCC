library(tidyverse)
library(RColorBrewer)
library(tidytext)
library(ggdist)

setwd("/home/fengtang/jupyter_notebooks/working_script/label_transfer/2nd_submission")

################################### when cell label is finest #####################################################
# get the SuperSCC scores from old score file
scores = list.files(path = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/", pattern = ".+finest_cell_label.+", recursive = FALSE, full.names = TRUE)
scores = map_dfr(scores[map_lgl(scores, ~ str_detect(.x, "singler|superscc|scn"))][c(1,3,5)], read.csv)

# get the rename files
rename = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/dataset_statistics_summary.csv", fileEncoding = "GBK")
rename[15, 2] = "Braga_2021_A"
rename[16, 2] = "Braga_2021_B"
colnames(rename)[2] = "new"
scores= scores

# set the factor level of Dataset
scores$new = factor(scores$new, 
                    levels = c("Madissoon_2019","Morse_2019", "Reyfman_2019",
                               "Bharat_2020", "Deprez_2020", "Habermann_2020",
                               "Travaglini_2020", "Braga_2021_A", "Braga_2021_B"
                    ))


# read Seurat scores
seurat_score = read.csv("/home/fengtang/jupyter_notebooks/working_script/label_transfer/seurat_score_on_finest_cell_label_2025-04-14 20:59:00.csv")

# read scmap scores
scmap_score = read.csv("/home/fengtang/jupyter_notebooks/working_script/label_transfer/scmap_score_on_finest_cell_label_2025-04-14 20:51:44.csv")

# read scANVI scores
scANVI_score = read.csv("/home/fengtang/jupyter_notebooks/working_script/label_transfer/scANVI_score_on_finest_cell_label_2025-04-14 20:53:14.csv")

all_scores = bind_rows(scores, seurat_score, scmap_score, scANVI_score) %>% 
  left_join(rename[c("old", "new")], by = join_by(Dataset == old))

# output supplementary table 1 - label transfer 
all_scores = all_scores %>% 
  dplyr::select(Dataset = new, 
                Label_transfer_method = method, 
                Scoring_method = score,
                Score = value
                )

write.csv(all_scores, "supple_table_1_label_transfer.csv", row.names = F)


colors = c("#DC0000FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF", "#7E6148FF","#F39B7FFF")
names(colors) = c("SuperSCC", "Seurat", "Scanpy", "scANVI", "scmap", "SingleCellNet")

p1 = all_scores %>%
  subset(score %in% c("accuracy_score", "f1_score", "v_measure_score", "matthews_corrcoef")) %>% 
  mutate(method = case_when(
    method == "seurat" ~ "Seurat",
    method == "Scanpy" ~ "Scanpy",
    .default = method
  )) %>% 
  ggplot(aes(reorder_within(method, value, score), value, fill = method)) +
  geom_violin(alpha = 0.5, scale = "width", linewidth = 0.6, trim = TRUE) +
  geom_boxplot(color = "white", outlier.fill = "black", outlier.color = "black", width = 0.4,  fill = NA) +
  geom_point(size = 1.5, position = position_jitter(width = 0.1, height = 0), colour = "black", aes(fill = method), stroke = 0.5, shape = 21) +
  scale_fill_manual(values = colors) +
  scale_x_reordered(labels = tidytext::reorder_func) + 
  facet_wrap(~ score, scales = "free") +
  theme_bw() +
  ylab("Evulate score") + 
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_text(color = "black", family = "Times", size = 12),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        # strip.text = element_text(color = "black", family = "Times", size = 12),
        # strip.background = element_rect(fill = "lightgray"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        legend.position = "none",
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))

ggsave(plot = p1, "fig2D.pdf", width = 2500, height = 2500, units = "px")

p2 = all_scores %>%
  subset(! score %in% c("accuracy_score", "f1_score", "v_measure_score", "matthews_corrcoef")) %>% 
  mutate(method = case_when(
    method == "seurat" ~ "Seurat",
    method == "Scanpy" ~ "Scanpy",
    .default = method
  )) %>% 
  ggplot(aes(reorder_within(method, value, score), value, fill = method)) +
  geom_violin(alpha = 0.5, scale = "width", linewidth = 0.6, trim = TRUE) +
  geom_boxplot(color = "white", outlier.fill = "black", outlier.color = "black", width = 0.4,  fill = NA) +
  geom_point(size = 1.5, position = position_jitter(width = 0.1, height = 0), colour = "black", aes(fill = method), stroke = 0.5, shape = 21) +
  scale_fill_manual(values = colors) +
  scale_x_reordered(labels = tidytext::reorder_func) + 
  facet_wrap(~ score, scales = "free") +
  theme_bw() +
  ylab("Evulate score") + 
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_text(color = "black", family = "Times", size = 12),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        legend.position = "none",
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))

ggsave(plot = p2, "supp fig2-part1.pdf", width = 2500, height = 2500, units = "px")


all_scores2 = all_scores %>%
  group_by(score, method) %>% 
  summarise(avg = mean(value), se = sd(value)/sqrt(length(value)))
