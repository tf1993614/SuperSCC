library(tidyverse)
library(tidytext)
library(lmtest)
library(sandwich)
library(lmPerm)
library(broom)
library(dunn.test)
library(ComplexHeatmap)
library(effectsize)

setwd("/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/2nd_submssion/evaluation_res/")

# read score files
qwen_files = list.files(pattern = "^qwen_temp.+", 
                   path = "/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/2nd_submssion",
                   recursive = F,
                   full.names = T)

gpt_4_1_mini_files = list.files(pattern = "^gpt_4_1_mini_temp.+", 
                        path = "/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/2nd_submssion",
                        recursive = F,
                        full.names = T)

deepseek_v3_files = list.files(pattern = "^deepseek_v3_temp.+", 
                          path = "/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/2nd_submssion",
                          recursive = F,
                          full.names = T)

top_10_markers_files = list.files(pattern = "top10_markers", 
                                       path = "/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/2nd_submssion",
                                       recursive = F,
                                       full.names = T)



all_files = c(qwen_files, gpt_4_1_mini_files, deepseek_v3_files)

file_name = map_chr(all_files, ~ basename(.x) %>% 
            str_remove("_evaulation_res.csv")) %>%
            str_remove("top10_markers_")
        


# load score file
scores = map(all_files, read.csv)
names(scores) = file_name

scores = imap(
  scores,
  function(x, y) {
    rows = which(x$Dataset == "GSE136831_Kaminski_2020")
    remain_data = x[-c(rows), ]
    amend_data = x[c(rows), ]
    amend_data = amend_data %>% mutate(across(where(is.numeric), ~ if_else(.x == 0, 0.5, .x)))
    final_data = bind_rows(remain_data, amend_data)
  }
)


mymean = function(x) {
  return(mean(x, na.rm = T))
}

mymedian = function(x) {
  return(median(x, na.rm = T))
}

scores_dot = map(
  scores,
  function(x) {
    data = x %>% 
      group_by(GeneSet, Dataset) %>% 
      summarise(across(c(RelevantGeneRatio, BiologicalRelevanceScore), list(mean = mymean, median = mymedian)))
  }
)

names(scores_dot) = file_name

scores_dot = imap(scores_dot, ~ .x %>% mutate(model_and_temperatue = .y))


# aggregate scoring values 
score_average = map(
  scores_dot,
  function(x) {
    data = x %>% 
      group_by(GeneSet) %>% 
      summarise(
        mean_gr = mean(RelevantGeneRatio_mean, na.rm = T),
        mean_gr = mean(RelevantGeneRatio_mean, na.rm = T),
        se_gr = sd(RelevantGeneRatio_mean, na.rm = T)/sqrt(length(RelevantGeneRatio_mean)),
        median_gr = mean(RelevantGeneRatio_median, na.rm = T),
        mean_br = mean(BiologicalRelevanceScore_mean, na.rm = T),
        se_br = sd(BiologicalRelevanceScore_mean, na.rm = T)/sqrt(length(BiologicalRelevanceScore_mean)),
        median_br = mean(BiologicalRelevanceScore_median, na.rm = T)
      )
  }
)


score_average = imap(score_average, ~ .x %>% mutate(model_and_temperatue = .y))


all_scores_dot = bind_rows(scores_dot) %>%
  mutate(method = str_remove(GeneSet, "_gene_set")) 

# output supplementary table 1
supply_table_1 = all_scores_dot %>% 
  subset(! method %in% c("Scanpy_log_reg", "Seurat_roc")) %>%
  separate_wider_delim(cols = model_and_temperatue, 
                       delim = "_temp_", names = c("model", "temperature")) %>%
  rename("model" = "Model", "temperature" = "Temperature", "method" = "Method")

rename_files = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/dataset_statistics_summary.csv", fileEncoding = "GBK")
rename_files = rename_files %>% dplyr::select(old, new_dataset = Dataset)

supply_table_1 = supply_table_1 %>% 
  mutate(Dataset = case_when(
    Dataset == "49万分的pbmc" ~ "血液",
    Dataset == "49万分的大肠" ~ "大肠",
    Dataset == "肺数据集全部取的两万个" ~ "所有细胞分2万",
    .default = Dataset
  )) %>%
  mutate(
    Temperature = case_when(
      Temperature == "one_tenth" ~ "0.1",
      Temperature == "five_tenth" ~ "0.5",
      Temperature == "nine_tenth" ~ "0.9",
      .default = Temperature
    )
  ) %>%
  left_join(rename_files, by = join_by(Dataset == old)) %>% 
  dplyr::select(-Dataset) %>% 
  dplyr::select(GeneSet, 
                Dataset = new_dataset, 
                RelevantGeneRatio_mean, 
                BiologicalRelevanceScore_mean, 
                Model, 
                Temperature,
                Method) 

write.csv(supply_table_1, "supp_table_1_feature_selection.csv", row.names = F)


all_scores_average = bind_rows(score_average) %>%
  mutate(method = str_remove(GeneSet, "_gene_set")) 

colors = c("#DC0000FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF", "#7E6148FF")
names(colors) = c("SuperSCC", "Seurat", "Scanpy", "Seurat_roc", "Scanpy_log_reg")


# create bar plot function 
plot_setting =  function(x) {
  x +
  scale_fill_manual(values = colors) +
  scale_x_reordered(labels = tidytext::reorder_func) +
  facet_wrap(~ model_and_temperatue, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_blank(),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(32,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        legend.position = "none",
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))
}


bar_plot = function(data, x, y,  se_y, data2, x2, y2, ylab = NULL, order_column = "model_and_temperatue"){
  plot_setting(
    data %>% 
    ggplot(aes(reorder_within(!! sym(x), !! sym(y), !! sym(order_column)), !! sym(y), fill = !! sym(x))) +
    geom_bar(stat = "identity", color = "black") +
    geom_point(aes(reorder_within(!! sym(x2), !! sym(y2), !! sym(order_column)), y = !! sym(y2)),
               data = data2 ,
               size = 1,
               alpha = 0.5,
               position = position_jitter(width = 0.1, height = 0),
               show.legend = FALSE) +
    geom_errorbar(aes(ymin = !! sym(y), ymax = !! sym(y) + !! sym(se_y)), 
                  width = 0.4,
                  position = position_dodge2(width = 0.5, padding = 0.5)) +
    geom_hline(yintercept = 0.75, linetype =2) +
    labs(y = ylab)
  )
}

# split the aggregate scoring list
all_scores_average_list = all_scores_average %>% 
  mutate(group = str_extract(model_and_temperatue, "[^_]+")) %>% 
  subset(! method %in% c("Scanpy_log_reg", "Seurat_roc")) %>% 
  mutate(method = case_when(
    method == "Seurat_wilcox" ~ "Seurat",
    method == "Scanpy_t_test" ~ "Scanpy",
    .default = method
  )) %>% 
  group_by(group) %>% 
  group_split()

names(all_scores_average_list) = map_chr(all_scores_average_list, ~ unique(.x$group))


all_scores_average_list = imap(
  all_scores_average_list,
  ~ .x %>%
    mutate(model_and_temperatue = factor(model_and_temperatue, levels = paste0(.y, c("_temp_one_tenth", "_temp_five_tenth",  "_temp_nine_tenth"))))
)

all_scores_dot_list = all_scores_dot %>% 
  mutate(group = str_extract(model_and_temperatue, "[^_]+")) %>% 
  subset(! method %in% c("Scanpy_log_reg", "Seurat_roc")) %>% 
  mutate(method = case_when(
    method == "Seurat_wilcox" ~ "Seurat",
    method == "Scanpy_t_test" ~ "Scanpy",
    .default = method
  )) %>% 
  group_by(group) %>% 
  group_split()

names(all_scores_dot_list) = map_chr(all_scores_dot_list, ~ unique(.x$group))


all_scores_dot_list = imap(
  all_scores_dot_list,
  ~ .x %>%
    mutate(model_and_temperatue = factor(model_and_temperatue, levels = paste0(.y, c("_temp_one_tenth", "_temp_five_tenth",  "_temp_nine_tenth"))))
)

all(names(all_scores_average_list) == names(all_scores_dot_list))


# show relevant gene ratio bar plot
gr_plots = map2(
  all_scores_average_list,
  all_scores_dot_list,
  ~ bar_plot(
    data = .x,
    x = "method",
    y = "mean_gr",
    se_y = "se_gr",
    data2 = .y,
    x2 = "method",
    y2 = "RelevantGeneRatio_mean",
    ylab = "Relevant gene ratio"
  )
)


map2(
  gr_plots,
  c("deepseek", "gpt", "qwen"),
  ~ ggsave(filename = paste0("top10_markers_", .y, "_gr_plot.pdf"), plot = .x, width = 2500, height = 2000, units = "px")
)


# show biological relevance score bar plot
br_plots = map2(
  all_scores_average_list,
  all_scores_dot_list,
  ~ bar_plot(
    data = .x,
    x = "method",
    y = "mean_br",
    se_y = "se_br",
    data2 = .y,
    x2 = "method",
    y2 = "BiologicalRelevanceScore_mean",
    ylab = "Biological relevance Score"
  )
)


map2(
  br_plots,
  c("deepseek", "gpt", "qwen"),
  ~ ggsave(filename = paste0("top10_markers_", .y, "_br_plot.pdf"), plot = .x, width = 2500, height = 2000, units = "px")
)




