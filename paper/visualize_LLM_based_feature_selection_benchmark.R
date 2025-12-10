library(tidyverse)
library(tidytext)
library(lmtest)
library(sandwich)
library(lmPerm)
library(broom)
library(dunn.test)
library(ComplexHeatmap)
library(effectsize)

setwd("/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/5th_submission/")

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

deepseek_files_temperature_0 = "/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/5th_submission/top20_markers_deepseek_temp_zero_evaulation_res.csv"

qwen_files_temperature_0 = "/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/5th_submission/top20_markers_qwen_temp_zero_evaulation_res.csv"

gpt_4_1_files_temperature_0 = "/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/5th_submission/top20_markers_gpt_4_1_mini_temp_zero_evaulation_res.csv"

all_files = c(qwen_files, gpt_4_1_mini_files, deepseek_v3_files,
              deepseek_files_temperature_0, qwen_files_temperature_0, gpt_4_1_files_temperature_0)

file_name = map_chr(all_files, ~ basename(.x) %>% 
            str_remove("_evaulation_res.csv")) %>%
  str_remove("top10_markers_") %>%
  str_remove("top20_markers_") %>% 
  str_replace("_", "")

file_name[10] = "deepseekv3_temp_zero"        

# load score file
scores = map(all_files, read.csv)
names(scores) = file_name


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

colors = c("#DC0000FF", "#4DBBD5FF", "#3C5488FF")
names(colors) = c("SuperSCC", "Seurat", "Scanpy")


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
  mutate(method = case_when(
    method == "Seurat_wilcox" ~ "Seurat",
    method == "Scanpy_t_test" ~ "Scanpy",
    .default = method
  )) %>% 
  group_by(group) %>% 
  group_split()

names(all_scores_average_list) = map_chr(all_scores_average_list, ~ unique(.x$group))


# all_scores_average_list = imap(
#   all_scores_average_list,
#   ~ .x %>%
#     mutate(model_and_temperatue = factor(model_and_temperatue, levels = paste0(.y, c("_temp_one_tenth", "_temp_five_tenth",  "_temp_nine_tenth"))))
# )

all_scores_dot_list = all_scores_dot %>% 
  mutate(group = str_extract(model_and_temperatue, "[^_]+")) %>% 
  mutate(method = case_when(
    method == "Seurat_wilcox" ~ "Seurat",
    method == "Scanpy_t_test" ~ "Scanpy",
    .default = method
  )) %>% 
  group_by(group) %>% 
  group_split()

names(all_scores_dot_list) = map_chr(all_scores_dot_list, ~ unique(.x$group))


# all_scores_dot_list = imap(
#   all_scores_dot_list,
#   ~ .x %>%
#     mutate(model_and_temperatue = factor(model_and_temperatue, levels = paste0(.y, c("_temp_one_tenth", "_temp_five_tenth",  "_temp_nine_tenth"))))
# )

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

# aggregate score values 


# prare data for summarizing reuslt by line plot
summary_level_data_for_line_plot = all_scores_average %>%
  separate_wider_regex(cols = model_and_temperatue, patterns = c(model = "[^_]+", other = ".+temp[^a-z]+", temperature = ".*")) %>%
  mutate(
    model = case_when(model == "deepseek" ~ "DeepSeek-v3",
                      model == "gpt" ~ "GPT-4.1-mini",
                      model == "qwen" ~ "Qwen-max",
                      model == "gpt5" ~ "GPT-5",
                      .default = model),
    temperature = case_when(
      temperature == "five_tenth" ~ "Temperature_0.5",
      temperature == "one_tenth" ~ "Temperature_0.1",
      temperature == "nine_tenth" ~ "Temperature_0.9",
      temperature == "zero" ~ "Temperature_0",
      .default = temperature
    ),
    method = case_when(
      method == "Seurat_wilcox" ~ "Seurat",
      method == "Scanpy_t_test" ~ "Scanpy",
      .default = method
    )) %>%
  subset(method %in% c("Seurat", "Scanpy", "SuperSCC"))

p1 = summary_level_data_for_line_plot %>%
  ggplot(aes(temperature, mean_br, group = method, color = method)) +
  geom_pointrange(aes(ymin = mean_br - se_br, ymax = mean_br + se_br)) +
  geom_line() +
  facet_wrap(~ model) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_blank(),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none",
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5)) +
  ylab("Biological relevance score")


p2 = summary_level_data_for_line_plot %>%
  ggplot(aes(temperature, mean_gr, group = method, color = method)) +
  geom_pointrange(aes(ymin = mean_gr - se_gr, ymax = mean_gr + se_gr)) +
  geom_line() +
  facet_wrap(~ model) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_blank(),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "inside",
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5)) +
  ylab("Relevant gene ratio")

p1_p2 = patchwork::wrap_plots(p1, p2)
ggsave(plot = p1_p2, filename = "gr_br_summarize_level_plot.pdf", width = 3500, height = 2000, units = "px")


# multi anova to decide the main factor affecting relevant gene ratio score
# and biological relevance score
test_table = all_scores_dot %>% 
  subset(method %in% c("SuperSCC", "Scanpy_t_test", "Seurat_wilcox")) %>% 
  separate_wider_regex(cols = model_and_temperatue, 
                       patterns = c(model = "[^_]+", other = ".+temp[^a-z]+", temperature = ".*")) %>%
  filter(str_detect(other, "_4o_temp_") == FALSE) %>%
  dplyr::select(everything(), mean_gr = RelevantGeneRatio_mean, mean_br = BiologicalRelevanceScore_mean)

write.csv(test_table, "test_table.csv")


test_table_model = test_table %>% 
  group_by(model) %>%
  group_split()

names(test_table_model) = c("DeepSeek-v3", "GPT 4.1-mini", "Qwen-max")


test_table_temperature = test_table %>% 
  group_by(temperature) %>%
  group_split()

names(test_table_temperature) = c("Temperature_0.5", "Temperature_0.9", "Temperature_0.1",  "Temperature_0")

test_table_method = test_table %>% 
  group_by(method) %>%
  group_split()

names(test_table_method) = c("Scanpy", "Seurat", "SuperSCC")

# permutation test when model is steady
steady_model = map(
  test_table_model,
  function(x) {
    res1 = kruskal.test(mean_br ~ temperature, data = x)
    res2 = kruskal.test(mean_br ~ method, data = x)
    res3 = dunn.test(x$mean_br, x$temperature, method="bh", list = T)
    res4 = dunn.test(x$mean_br, x$method, method="bh", list = T)
    
    res5 = kruskal.test(mean_gr ~ temperature, data = x)
    res6 = kruskal.test(mean_gr ~ method, data = x)
    res7 = dunn.test(x$mean_gr, x$temperature, method="bh", list = T)
    res8 = dunn.test(x$mean_gr, x$method, method="bh", list = T)
    
    list(br_kw_test_1 = res1, 
         br_kw_test_2 = res2, 
         br_dunn_test_1 = res3, 
         br_dunn_test_2 = res4, 
         gr_kw_test_1 = res5, 
         gr_kw_test_2 = res6, 
         gr_dunn_test_1 = res7, 
         gr_dunn_test_2 = res8)
  }
)

# permutation test when temperature is steady
steady_temperature = map(
  test_table_temperature,
  function(x) {
    res1 = kruskal.test(mean_br ~ model, data = x)
    res2 = kruskal.test(mean_br ~ method, data = x)
    res3 = dunn.test(x$mean_br, x$model, method="bh", list = T)
    res4 = dunn.test(x$mean_br, x$method, method="bh", list = T)
    
    res5 = kruskal.test(mean_gr ~ model, data = x)
    res6 = kruskal.test(mean_gr ~ method, data = x)
    res7 = dunn.test(x$mean_gr, x$model, method="bh", list = T)
    res8 = dunn.test(x$mean_gr, x$method, method="bh", list = T)
    
    list(br_kw_test_1 = res1, 
         br_kw_test_2 = res2, 
         br_dunn_test_1 = res3, 
         br_dunn_test_2 = res4, 
         gr_kw_test_1 = res5, 
         gr_kw_test_2 = res6, 
         gr_dunn_test_1 = res7, 
         gr_dunn_test_2 = res8)
  }
)


# permutation test when method is steady
steady_method = map(
  test_table_method,
  function(x) {
    res1 = kruskal.test(mean_br ~ model, data = x)
    res2 = kruskal.test(mean_br ~ temperature, data = x)
    res3 = dunn.test(x$mean_br, x$model, method="bh", list = T)
    res4 = dunn.test(x$mean_br, x$temperature, method="bh", list = T)
    
    res5 = kruskal.test(mean_gr ~ model, data = x)
    res6 = kruskal.test(mean_gr ~ temperature, data = x)
    res7 = dunn.test(x$mean_gr, x$model, method="bh", list = T)
    res8 = dunn.test(x$mean_gr, x$temperature, method="bh", list = T)
    
    list(br_kw_test_1 = res1, 
         br_kw_test_2 = res2, 
         br_dunn_test_1 = res3, 
         br_dunn_test_2 = res4, 
         gr_kw_test_1 = res5, 
         gr_kw_test_2 = res6, 
         gr_dunn_test_1 = res7, 
         gr_dunn_test_2 = res8)
  }
)

# effect size evaluation
effect_size = list(
  br_method = rank_eta_squared(mean_br ~ method, data = test_table),
  gr_method = rank_eta_squared(mean_gr ~ method, data = test_table),
  br_model = rank_eta_squared(mean_br ~ model, data = test_table),
  gr_model = rank_eta_squared(mean_gr ~ model, data = test_table),
  br_temperature = rank_eta_squared(mean_br ~ temperature, data = test_table),
  gr_temperature = rank_eta_squared(mean_gr ~ temperature, data = test_table)
)

effect_size = map_dfr(
  effect_size,
  function(x) {
    data.frame(effect_size = x$rank_eta_squared)
  }
)

effect_size$category = rep(c("method", "model", "temperature"), each = 2)
effect_size$score = rep(c("br", "gr"), 3)

p3 = effect_size %>% 
  ggplot(aes(category, effect_size, fill = score)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#FBDD49", "#5F559BFF")) +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_blank(),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none",
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(color = "black", family = "Times", margin = margin(40,8,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))


# only focus on low-frequency cell type 
freq_table = read.csv("/home/fengtang/jupyter_notebooks/working_script/evulate_feature_selection/2nd_submssion/evaluation_res/cell_type_frequency_table.csv", row.names = 1)

all_scores = imap_dfr(
  scores,
  function(x, y) {
    data = x %>% 
      dplyr::select(-X) %>% 
      mutate(model_and_temperatue = y) %>% 
      separate_wider_regex(cols = model_and_temperatue, patterns = c(model = "[^_]+", other = ".+temp[^a-z]+", temperature = ".*")) %>%
      mutate(
        model = case_when(model == "deepseek" ~ "DeepSeek-v3",
                          model == "gpt" ~ "GPT-4.1-mini",
                          model == "qwen" ~ "Qwen-max",
                          model == "gpt5" ~ "GPT-5",
                          .default = model),
        temperature = case_when(
          temperature == "five_tenth" ~ "Temperature_0.5",
          temperature == "one_tenth" ~ "Temperature_0.1",
          temperature == "nine_tenth" ~ "Temperature_0.9",
          temperature == "zero" ~ "Temperature_0",
          .default = temperature
        )
      )
    
    
    if(any(colnames(data) == "cell_type")) {
      data = data %>% rename("CellType" = "cell_type")
    }
    return(data)
  }
)


all_scores = all_scores %>%
  left_join(
    freq_table,
    join_by(Dataset == dataset, CellType == cell_type)
  ) %>% 
  subset(GeneSet %in% c("SuperSCC_gene_set", "Seurat_wilcox_gene_set", "Scanpy_t_test_gene_set")) %>%
  mutate(
    GeneSet = case_when(
      GeneSet == "SuperSCC_gene_set" ~ "SuperSCC",
      GeneSet == "Seurat_wilcox_gene_set" ~ "Seurat",
      GeneSet == "Scanpy_t_test_gene_set" ~ "Scanpy",
    )
  )

low_freq_scores = all_scores %>% 
  subset(proportion < 0.001) %>% 
  arrange(RelevantGeneRatio)


p4 = low_freq_scores %>% 
  ggplot(aes(proportion, BiologicalRelevanceScore, color = GeneSet)) +
  geom_line( size = 1) +
  geom_point(colour = "black", alpha = 0.5) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = colors ) +
  facet_grid(model ~ temperature, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x.bottom = element_text(color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_text(color = "black", family = "Times", size = 15),
        axis.title.x.bottom = element_text(color = "black", family = "Times", size = 12),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text.x = element_text(color = "black", family = "Times", margin = margin(16,8,8,8, unit = "pt")),
        strip.text.y = element_text(color = "black", family = "Times", margin = margin(8,16,8,8, unit = "pt")),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        legend.position = "none",
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5),
  ) +
  labs(x = "Frequency", y = "Biological relevance score")

ggsave(plot = p4, filename = "low_frequency_cell_type_feature_selection_compare.pdf", width = 3500, height = 2000, units = "px")
