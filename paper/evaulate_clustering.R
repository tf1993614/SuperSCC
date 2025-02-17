# load packages
library(tidyverse)
library(ggsci)
library(scales)
library(tidytext)
library(data.table)

# set working dir
setwd("/home/fengtang/jupyter_notebooks/working_script/evulate_clustering")

# read files in
files = list.files(pattern = ".+csv$", path = "./clustering_performance//", recursive = FALSE, full.names = TRUE)
names = files %>% basename() %>% str_remove("\\.csv")
scores = map(files, ~ read_csv(.x, col_select = 2:4))
names(scores) = names

# read files with chinese
rename_files = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/dataset_statistics_summary.csv", fileEncoding = "GBK")
rename_files[3, 1] = "GSE160269_食管鳞癌"
rename_files[15, 2] = "Braga_2021_A"
rename_files[16, 2] = "Braga_2021_B"

# only focus level 2 and compare different seurat resolutions
new_names = names[str_detect(names, "lv2")]
scores1 = scores[str_detect(names(scores), "lv2")]

scores1 = map(names(scores1), ~ scores1[[.x]] %>% 
               mutate(score = str_to_upper(str_extract(.x, "ami|nmi|fmi|ari")),
                      resolution = str_extract(.x, "\\d+.\\d|\\d+")))

scores1 = bind_rows(scores1)

# rename dataset name
scores1 = scores1 %>% left_join(
  rename_files %>% dplyr::select(old, Dataset),
  by = join_by(variable == old)
)

scores1_group = scores1 %>% 
  filter(method %in% c("SuperSCC", "Seurat")) %>% 
  group_by(resolution) %>% 
  group_split()

colors = c("#DC0000FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF", "#7E6148FF", "#F39B7FFF")
names(colors) = c("SuperSCC", "Seurat", "Scanpy", "SC3", "CIDR", "DUBStepR")

map(
  scores1_group,
  function(x){
    p = x %>%
  ggplot(aes(Dataset, value, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = colors, name = "Method") + 
  facet_wrap(~ score, scales = "free_y") +
  ggtitle(paste0("Seurat_resolution_", unique(x$resolution))) +
  theme_bw() +
  ylab("Evulate score") + 
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times"),
        axis.text.y.left = element_text(color = "black", family = "Times"),
        axis.title.y.left = element_text(color = "black", family = "Times"),
        legend.text = element_text(color = "black", family = "Times"),
        legend.title = element_text(color = "black", family = "Times"),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text = element_text(color = "black", family = "Times"),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5),
        legend.position = "right")
    
  
   ggsave(filename = paste0("Seurat_resolution_", unique(x$resolution), ".png"), p, width = 3500, height = 2000, dpi = 300, unit = "px")
  }
)

# only focus level 2 and compare different methods 
# seurat resolution is 0.8
scores2 = scores1 %>% subset(resolution == 0.8) %>%
  group_by(Dataset, score) %>% 
  dplyr::arrange(value, .by_group = TRUE) %>% 
  ungroup() %>% 
  mutate(new_method = paste0(Dataset, "_", method, "_", score)) %>%
  mutate(new_method = factor(new_method, levels = new_method))

new_colors = case_when(str_detect(levels(scores2$new_method), "SC3") ~ "#00A087FF",
          str_detect(levels(scores2$new_method), "SuperSCC") ~ "#DC0000FF",
          str_detect(levels(scores2$new_method), "Seurat") ~ "#4DBBD5FF",
          str_detect(levels(scores2$new_method), "Scanpy") ~ "#3C5488FF",
          str_detect(levels(scores2$new_method), "CIDR") ~ "#7E6148FF",
          str_detect(levels(scores2$new_method), "DUBStepR") ~ "#F39B7FFF",
          )

names(new_colors) = levels(scores2$method)

scores2$dataset = factor(scores2$Dataset, 
                         levels = c("Madissoon_2019","Morse_2019", "Reyfman_2019",
                                    "Adams_2020",  "Bharat_2020", "Deprez_2020",
                                    "Habermann_2020", "Travaglini_2020", "Wang_2020", 
                                    "Braga_2021_A", "Braga_2021_B", "Zhang_2021",
                                    "Jones_2022_Gut", "Jones_2022_PBMC", "Perez_2022",
                                    "Glasner_2023", "Kumar_2023", "Sikkema_2023"
                                    ))

p = scores2 %>%
      ggplot(aes(Dataset, value)) +
      # geom_tile(aes(fill = value), color = "black", linewidth = 0.5) +
      geom_bar(stat = "identity", position = "dodge", aes(fill = new_method)) +
      scale_fill_manual(values = new_colors) +
      scale_x_reordered() + 
      # # scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name =
      #                                                                 "RdYlBu")))(200), name = "Evulation score") + 
      # geom_text(aes(label = round(value, 2)), color = "black", family = "Times") +
      facet_wrap(~ score, scales = "free") +
      theme_bw() +
      # coord_fixed() +
      # ylab("Evulate score") + 
      theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
            axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
            # axis.title.y.left = element_text(color = "black", family = "Times"),
            legend.text = element_text(color = "black", family = "Times"),
            legend.title = element_text(color = "black", family = "Times"),
            axis.title.y.left = element_blank(), 
            axis.title.x.bottom = element_blank(),
            # panel.grid = element_blank(),
            panel.border = element_rect(colour = "black"),
            strip.text = element_text(color = "black", family = "Times"),
            strip.background = element_rect(fill = "lightgray"),
            plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5),
            legend.position = "none"
           ) 
    
ggsave(filename = paste0("Scores_compare_between_dataset", ".pdf"), p, width = 3500, height = 2000, dpi = 300, unit = "px")


# calculate the mean score values across datasets
score_mean = scores2 %>% 
  group_by(method, score) %>% 
  dplyr::summarise(average = mean(value), se = sd(value)/sqrt(n()))

# prepare dot dataset that can be 
# used to draw the dots above the bar plot
score_dot = scores2
score_dot$average = score_dot$value

p2 = score_mean%>%
  ggplot(aes(reorder_within(method, average, score), average, fill = method)) +
  geom_bar(stat = "identity", color = "black") +
  geom_point(aes(y = value), data = score_dot, size = 1, alpha = 0.5, 
             position = position_jitter(width = 0.1, height = 0), show.legend = FALSE) +
  geom_errorbar(aes(ymin = average, ymax = average + se),
                width = 0.4,
                position = position_dodge2(width = 0.5, padding = 0.5)) +
  scale_fill_manual(values = colors, name = "Method") + 
  scale_x_reordered(labels = tidytext::reorder_func) + 
  facet_wrap(~ score, scales = "free") +
  theme_bw() +
  ylab("Evulate score") + 
  theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", family = "Times", size = 12),
        axis.text.y.left = element_text(color = "black", family = "Times", size = 12),
        axis.title.y.left = element_text(color = "black", family = "Times", size = 12),
        legend.text = element_text(color = "black", family = "Times", size = 12),
        legend.title = element_text(color = "black", family = "Times", size = 15),
        axis.title.x.bottom = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text = element_text(color = "black", family = "Times"),
        strip.background = element_rect(fill = "lightgray"),
        plot.title = element_text(color = "black", family = "Times", size = 15, hjust = 0.5))

ggsave(filename = paste0("Score_summary_plot", ".pdf"), p2, width = 3000, height = 2000, dpi = 300, unit = "px")




# evulate feature selection 
score_df = map(c("SuperSCC", "SingleR", "SingleCellNet", "scMAP"),function(x){
  if(x == "SuperSCC"){
    tibble(Method = x, Accuracy = 83.81, homogeneity = 85.35, completeness = 90.26, V_measure = 86.74, Precision = 86.95, Recall = 83.81, F1 = 88.71)
  }
  else if(x == "SingleR"){
    tibble(Method = x, Accuracy = 82.95, homogeneity = 85.15, completeness = 90.71, V_measure = 87.84, Precision = 81.47, Recall = 82.96, F1 = 79.51)
  }
  else if(x == "SingleCellNet"){
    tibble(Method = x, Accuracy = 81.03, homogeneity = 81.45, completeness = 88.70, V_measure = 84.92, Precision = 84.33, Recall = 81.02, F1 = 84.72)
  }
  else if(x == "scMAP"){
    tibble(Method = x, Accuracy = 53.82, homogeneity = 67.24, completeness = 80.16, V_measure = 73.14, Precision = 74.70, Recall = 53.82, F1 = 58.60)
  }
})

score_df = bind_rows(score_df) %>% dplyr::select(Method, Accuracy,F1,Precision,Recall)

score_df = score_df %>% pivot_longer(-c(Method), values_to = "value", names_to = "score")

fig = score_df %>%
  ggplot(aes(score, value, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = rev(pal_npg("nrc")(10)[1:4])) + 
  theme_bw() +
  scale_y_continuous(labels = ~ .x/100, limits = c(NA, 100)) +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", family = "Times"),
        axis.text.y.left = element_text(color = "black", family = "Times"),
        axis.title.y.left = element_text(color = "black", family = "Times"),
        legend.text = element_text(color = "black", family = "Times"),
        legend.title = element_text(color = "black", family = "Times"),
        axis.title.x.bottom = element_blank()) +
  labs(y = "Evulation value") 

ggsave(filename = paste0("autoannotation", ".png")q, fig, width = 2000, height = 1000, dpi = 300, unit = "px")

# 