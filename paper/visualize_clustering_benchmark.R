# load packages
library(tidyverse)
library(ggsci)
library(scales)
library(tidytext)
library(data.table)

# set working dir
setwd("/home/fengtang/jupyter_notebooks/working_script/evulate_clustering/3rd_submission/")

# file = list.files(path = "/home/fengtang/jupyter_notebooks/working_script/evulate_clustering/3rd_submission/lv2_clustering_performance_3nd_submission_2025-08-17 21:05:58.csv", full.names = T)
file = "/home/fengtang/jupyter_notebooks/working_script/evulate_clustering/3rd_submission/lv2_clustering_performance_3nd_submission_2025-08-17 21:05:58.csv"
scores = read_csv(file) 

# read files with chinese
rename_files = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/dataset_statistics_summary.csv", fileEncoding = "GBK")
rename_files[3, 1] = "GSE160269_食管鳞癌"
rename_files[15, 2] = "Braga_2021_A"
rename_files[16, 2] = "Braga_2021_B"

scores = scores %>% left_join(
  rename_files %>% dplyr::select(old, Dataset),
  by = join_by(dataset == old)
)

scores = scores %>% 
  mutate(new_method = paste0(Dataset, "_", method, "_", score_method)) %>%
  mutate(new_method = factor(new_method, levels = new_method))

# output supplementary table 1 - clustering
scores = scores %>% 
  dplyr::select(
    Dataset,
    Clustering_method = method,
    Scoring_method = score_method,
    Score = score
  )

write.csv(scores, "supple_table_1_clustering.csv", row.names = F)


# calculate the mean score values across datasets
score_mean = scores%>% 
  group_by(method, score_method) %>% 
  dplyr::summarise(average = mean(score), se = sd(score)/sqrt(n()))

score_dot = scores
score_dot$score = if_else(score_dot$score < 0, 0, score_dot$score)
score_dot$average = 1



colors = c("#DC0000FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF", "#7E6148FF", "#F39B7FFF","grey" ,"#808080","pink")
names(colors) = c("SuperSCC", "Seurat", "Scanpy", "SC3", "CIDR", "Sccaf", "Scshc", "Sclca", "Monocle")

# show the average score across datasets per method
p1 = score_mean %>%
  ggplot(aes(reorder_within(method, average, score_method), average, fill = method)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  geom_point(aes(y = score), data = score_dot, size = 1, alpha = 0.5,
             position = position_jitter(width = 0.1, height = 0), show.legend = FALSE) +
  geom_errorbar(aes(ymin = average, ymax = average + se),
                width = 0.4,
                linewidth = 0.2,
                position = position_dodge2(width = 0.5, padding = 0.5)) +
  scale_fill_manual(values = colors, name = "method") + 
  scale_x_reordered(labels = tidytext::reorder_func) + 
  facet_wrap(~ score_method, scales = "free") +
  theme_bw() +
  ylab("Evulate score") + 
  theme(
    # axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
    axis.text.x.bottom = element_blank(),
    axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
    axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
    legend.text = element_text(color = "black", family = "Times", size = 12),
    legend.title = element_text(color = "black", family = "Times", size = 15),
    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    # axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    strip.text = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    legend.position = "none",
    plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))

ggsave(plot = p1, "fig2C.pdf", width = 2500, height = 2500, units = "px")



# show score per method per dataset
score_per_group = scores %>% 
  dplyr::arrange(score, .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(new_method = paste0(Dataset, "_", method, "_", score_method)) %>%
  mutate(new_method = factor(new_method, levels = new_method))


new_colors = case_when(str_detect(levels(score_per_group$new_method), "SC3") ~ "#00A087FF",
                       str_detect(levels(score_per_group$new_method), "SuperSCC") ~ "#DC0000FF",
                       str_detect(levels(score_per_group$new_method), "Seurat") ~ "#4DBBD5FF",
                       str_detect(levels(score_per_group$new_method), "Scanpy") ~ "#3C5488FF",
                       str_detect(levels(score_per_group$new_method), "CIDR") ~ "#7E6148FF",
                       str_detect(levels(score_per_group$new_method), "Sccaf") ~ "#F39B7FFF",
                       str_detect(levels(score_per_group$new_method), "Scshc") ~ "grey",
                       str_detect(levels(score_per_group$new_method), "Sclca") ~ "black",
                       str_detect(levels(score_per_group$new_method), "Monocle") ~ "pink",
)

p2 = score_per_group %>%
  ggplot(aes(Dataset, score)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = new_method)) +
  scale_fill_manual(values = new_colors) +
  scale_x_reordered() + 
  facet_wrap(~ score_method, scales = "free") +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times"),
        legend.title = element_text(color = "black", family = "Times"),
        axis.title.y.left = element_blank(), 
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        # strip.text = element_text(color = "black", family = "Times"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5),
        legend.position = "none"
  ) 

ggsave(filename = paste0("Scores_compare_between_dataset", ".pdf"), plot = p2, width = 3500, height = 2000, dpi = 300, unit = "px")




# compare SuperSCC with Scanpy and Seurat under 
# different resolutions (weighted geometric mean to aggregate scores under different resolutions)
file = "Aggregate_mean_scanpy_seurat_under_different_resolution_2025-08-16 16:12:03.csv"
scores = read_csv(file) #%>% dplyr::select(2:5)

scores = scores %>% left_join(
  rename_files %>% dplyr::select(old, Dataset),
  by = join_by(dataset == old)
)


# output the scores per method per dataset per resolution
scores = scores %>% 
  dplyr::select(c(-1, -2)) %>%
  dplyr::select(Dataset, Clustering_method = method,
                Scoring_method = score_method,
                Resolution = resolution,
                Score = score) %>% 
  mutate(Resolution = if_else(Clustering_method == "SuperSCC", "Null", as.character(Resolution)))
  
write.csv(scores, "clustering_scores_scnapy_seurat_per_resolution_per_dataset.csv",row.names = F)


scores = scores %>% 
  mutate(new_method = paste0(Dataset, "_", method, "_", score_method)) %>%
  mutate(new_method = factor(new_method, levels = new_method))


score_mean = scores%>% 
  group_by(method, score_method) %>% 
  dplyr::summarise(average = mean(score), se = sd(score)/sqrt(n()))

score_dot = scores
score_dot$average = 1


p3 = score_mean%>%
  ggplot(aes(reorder_within(method, average, score_method), average, fill = method)) +
  geom_bar(stat = "identity", color = "black") +
  geom_point(aes(y = score), data = score_dot, size = 1, alpha = 0.5,
             position = position_jitter(width = 0.1, height = 0), show.legend = FALSE) +
  geom_errorbar(aes(ymin = average, ymax = average + se),
                width = 0.4,
                position = position_dodge2(width = 0.5, padding = 0.5)) +
  scale_fill_manual(values = colors, name = "method") + 
  scale_x_reordered(labels = tidytext::reorder_func) + 
  facet_wrap(~ score_method, scales = "free") +
  theme_bw() +
  ylab("Evulate score") + 
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_text(color = "black", family = "Times", size = 15),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        # strip.text = element_text(color = "black", family = "Times"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        legend.position = "none",
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))

ggsave(plot = p3, "Aggregate_mean_scanpy_seurat_under_different_resolution.pdf", width = 2500, height = 2500, units = "px")



