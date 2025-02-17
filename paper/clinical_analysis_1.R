library(tidyverse)
library(data.table)
library(Seurat)
library(ggpubr)
library(broom)

setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/F_cluster")

files = list.files(path = "/mnt/disk5/zhongmin/superscc/结果位置/应激_免疫内皮metadata", pattern = ".+csv$", full.names = T)
file_names = basename(files) %>% str_extract("[^_]+_\\d+")
names(files) = file_names
meta = map(files, ~ read.csv(.x))

# focus on laughney_2020 dataset
data = meta[[2]]

# calculate the frequency
patient_Fclustr_freq = data %>% group_by(patient, Fcluster) %>% count() %>% ungroup()
prop = patient_Fclustr_freq %>% group_by(Fcluster) %>% mutate(prop = n / sum(n))
prop_lists = map(
  unique(data$Fcluster)[unique(data$Fcluster) != 6],
  function(x) {
    x1 = prop %>% subset(Fcluster == 6) %>% dplyr::select(patient, Fcluster1 = Fcluster, prop1 = prop)
    x2 = prop %>% subset(Fcluster == x) %>% dplyr::select(patient, Fcluster2 = Fcluster, prop2 = prop)
    x3 = x1 %>% left_join(x2, by = "patient")
  }
)

prop_lists2 = prop_lists[map_lgl(prop_lists, ~ all(!is.na(.x$prop2)))]


# do correlation analysis between prop of stress B cell in F cluster 6 and
# each other cluster
p1 = map(
  prop_lists2,
  function(x) {
    fig = ggscatter(x,
              x = "prop1",
              y = "prop2",
              fill = "grey",
              alpha = 0.5,
              color = "black",
              shape = 21,
              size = 2, # Points color, shape and size
              add = "reg.line",  # Add regressin line
              add.params = list(color = "#E64B35FF", fill = "grey"), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
              cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor,
              cor.coef.size = 4)
      fig = fig + labs(x = "Fcluster 6", y = paste0("Fcluster ", unique(x$Fcluster2)))

  }
)

# get the pearson correlation test result
pearson_res = map_dfr(
  prop_lists2,
  ~ cor.test(.x$prop1, .x$prop2) %>%
    broom::tidy() %>%
    mutate(compare = unique(.x$Fcluster2))
)

# organize the result
pearson_res2 = pearson_res %>%
  dplyr::select(estimate, p.value, compare) %>%
  mutate(init = "F_cluster_6") %>%
  mutate(compare = paste0("F_cluster_", compare)) %>%
  mutate(sign = if_else(estimate > 0.5 & p.value < 0.05, "*", ""))

# visualize the correlation analysis in heatmap
heatmap = pearson_res2 %>%
ggplot(aes(compare, init, fill = estimate)) +
geom_tile() +
geom_text(aes(compare, init, label = sign)) +
coord_fixed() +
scale_fill_gradientn(colors = c("#3C5488FF",  "lightgrey", "#E64B35FF")) +
theme(
  axis.ticks = element_blank(),
  axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 1, colour = "black"),
  axis.text.y.left = element_text(colour = "black"),
  axis.title = element_blank(),
  panel.background = element_rect(fill = "white", colour = "white"),
  legend.position = "top"
)

ggsave(heatmap, filename = "correlation_heatmap.pdf", width = 2000, height = 1000, units = "px")
