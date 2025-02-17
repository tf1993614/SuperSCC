library(tidyverse)
library(RColorBrewer)
library(tidytext)
library(ggdist)

setwd("~/jupyter_notebooks/working_script/label_transfer")


################################### when cell label is finest #####################################################
# get the score file
scores = list.files(path = "/home/fengtang/jupyter_notebooks/working_script/label_transfer", pattern = ".+finest_cell_label.+", recursive = FALSE, full.names = TRUE)
scores = map_dfr(scores[map_lgl(scores, ~ str_detect(.x, "singler|superscc|scn"))][c(1,3,5)], read.csv)

# get the rename files
rename = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/dataset_statistics_summary.csv", fileEncoding = "GBK")
rename[15, 2] = "Braga_2021_A"
rename[16, 2] = "Braga_2021_B"
colnames(rename)[2] = "new"
scores= scores %>% left_join(rename[c("old", "new")], by = join_by(Dataset == old))

# set the factor level of Dataset
scores$new = factor(scores$new,
                    levels = c("Madissoon_2019","Morse_2019", "Reyfman_2019",
                               "Bharat_2020", "Deprez_2020", "Habermann_2020",
                               "Travaglini_2020", "Braga_2021_A", "Braga_2021_B"
                    ))


# draw the heatmap to show the score per method per dataset per scoring metric
p = scores %>%
  ggplot(aes(new, method)) +
  geom_tile(aes(fill = value), color = "black", linewidth = 0.5) +
  # geom_bar(stat = "identity", position = "dodge") +
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                  "RdYlBu")))(200), name = "Evulation score") +
  geom_text(aes(label = round(value, 2)), color = "black", family = "Times") +
  facet_wrap(~ score) +
  theme_bw() +
  coord_fixed() +
  # ylab("Evulate score") +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        # axis.title.y.left = element_text(color = "black", family = "Times"),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text = element_text(color = "black", family = "Times"),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5),
        # legend.position = "none"
  )


ggsave(filename = "finest_label_transfer_score_per_dataset.pdf", p, width = 4000, height = 2000, units = "px")


# draw the bar plot to show the average score per method per scoring metric
colors = c("#DC0000FF", "#4DBBD5FF", "#3C5488FF")
names(colors) = c("SuperSCC", "SingleCellNet", "SingleR")

p = scores %>%
  ggplot(aes(reorder_within(method, value, score), value, fill = method)) +
  geom_violin(alpha = 0.5, scale = "width", linewidth = 0.6, trim = TRUE) +
  geom_boxplot(color = "white", outlier.fill = "black", outlier.color = "black", width = 0.4,  fill = NA) +
  geom_point(size = 1.5, position = position_jitter(width = 0.1, height = 0), colour = "black", aes(fill = method), stroke = 0.5, shape = 21) +
  scale_fill_manual(values = colors, name = "Method") +
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
        strip.text = element_text(color = "black", family = "Times", size = 12),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))

ggsave(filename = "average_finest_label_transfer_score.pdf", width = 3500, height = 2500, units = "px")



################################### when train on multi datasets and cell label is finest #####################################################

# get the score file
scores = list.files(path = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/", pattern = ".+multi_datasets+", recursive = FALSE, full.names = TRUE)
scores = map_dfr(scores[3:5], read.csv)

# get the rename files
rename = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/dataset_statistics_summary.csv", fileEncoding = "GBK")
rename[15, 2] = "Braga_2021_A"
rename[16, 2] = "Braga_2021_B"
colnames(rename)[2] = "new"
scores= scores %>% left_join(rename[c("old", "new")], by = join_by(Dataset == old))

# set the factor level of Dataset
scores$new = factor(scores$new,
                    levels = c("Morse_2019", "Deprez_2020", "Habermann_2020",
                               "Travaglini_2020", "Braga_2021_A", "Braga_2021_B"
                    ))


# draw the heatmap to show the score per method per dataset per scoring metric
p = scores %>%
  ggplot(aes(new, method)) +
  geom_tile(aes(fill = value), color = "black", linewidth = 0.5) +
  # geom_bar(stat = "identity", position = "dodge") +
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                                  "RdYlBu")))(200), name = "Evulation score") +
  geom_text(aes(label = round(value, 2)), color = "black", family = "Times") +
  facet_wrap(~ score) +
  theme_bw() +
  coord_fixed() +
  # ylab("Evulate score") +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times"),
        axis.text.y.left = element_text(color = "black", family = "Times"),
        # axis.title.y.left = element_text(color = "black", family = "Times"),
        legend.text = element_text(color = "black", family = "Times"),
        legend.title = element_text(color = "black", family = "Times"),
        axis.title.y.left = element_blank(),
        axis.title.x.bottom = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text = element_text(color = "black", family = "Times"),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5),
        # legend.position = "none"
  )


ggsave(filename = "train_on_multi_datasets_finest_label_transfer_score_per_dataset.png", p, width = 4000, height = 2000, units = "px")


# draw the bar plot to show the average score per method per scoring metric
colors = c("#DC0000FF", "#4DBBD5FF", "#3C5488FF")
names(colors) = c("SuperSCC", "SingleCellNet", "SingleR")

p = scores %>%
  ggplot(aes(reorder_within(method, value, score), value, fill = method)) +
  geom_violin(alpha = 0.5, scale = "width", linewidth = 0.6, trim = TRUE) +
  geom_boxplot(color = "white", outlier.fill = "black", outlier.color = "black", width = 0.4,  fill = NA) +
  geom_point(size = 1.5, position = position_jitter(width = 0.1, height = 0), colour = "black", aes(fill = method), stroke = 0.5, shape = 21) +
  scale_fill_manual(values = colors, name = "Method") +
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
        strip.text = element_text(color = "black", family = "Times", size = 12),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))

p + theme(text = element_blank())
ggsave(filename = "train_on_multi_datasets_average_finest_label_transfer_score.png", width = 3500, height = 2500, units = "px")


output = scores %>%
  pivot_wider(
    everything(),
    names_from = "score",
    values_from = "value"
  )

write.csv(scores, paste0("combined_on_multi_datasets_and_finest_cell_label_2", date(), ".csv"))
