library(tidyverse)
library(openxlsx)

setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/F_cluster")

data = openxlsx::read.xlsx("/mnt/disk5/zhongmin/superscc/结果位置/addmodulescore/师兄算的结果/addmodule/F_gene_module/addmodulescore_f_gene_module.xlsx")
disease_info = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/作者年_表达矩阵_metadata_pkl/作者年_表达矩阵_metadata_pkl.csv", fileEncoding = "GBK")
data = data %>% left_join(disease_info %>% dplyr::select(new, organization, diseases), by = "new")
data = data %>%
  pivot_longer(
    starts_with("F_gene_module"),
    names_to = "module",
    values_to = "value"
  )

immune_response_module = c(2, 11, 13, 15, 16, 21, 25, 27, 33, 38, 39, 40, 43, 46, 47, 51, 52, 53, 55, 56, 58, 61, 63, 65)
angiogenesis_module = c(1, 30, 36)

all_stat = data %>%
  subset(value != "NO") %>%
  group_by(module) %>%
  count(value) %>%
  group_by(module) %>%
  summarise(module = module, value = value, ratio = n / sum(n)) %>%
  ungroup() %>%
  mutate(module_num = as.numeric(str_extract(module, "\\d+")))

immune_all_stat = all_stat %>% subset(module_num %in% immune_response_module)
angiogenesis_all_stat = all_stat %>% subset(module_num %in% angiogenesis_module)
combined_all_stat = all_stat %>% subset(module_num %in% c(angiogenesis_module, immune_response_module))

tissue_stat = data %>%
  subset(value != "NO") %>%
  group_by(module, organization) %>%
  count(value) %>%
  group_by(module, organization) %>%
  summarise(module = module, tissue = organization, value = value, ratio = n / sum(n)) %>%
  ungroup() %>%
  mutate(module_num = as.numeric(str_extract(module, "\\d+")))

immune_tissue_stat = tissue_stat %>% subset(module_num %in% immune_response_module)
angiogenesis_tissue_stat = tissue_stat %>% subset(module_num %in% angiogenesis_module)
combined_tissue_stat = tissue_stat %>% subset(module_num %in% c(angiogenesis_module, immune_response_module))

F_module_rename = data.frame(
  original = c(33, 21, 27, 56, 38, 25, 65, 46,
               58, 51, 55, 52, 61, 16, 15, 53,
               13, 40, 2, 43, 11, 47, 39, 63,
               14, 6, 59, 28, 32, 17, 4, 30,
               1, 36, 8, 5, 64, 24, 12, 34, 50,
               26, 57, 3, 54, 49, 10, 42, 31,
               60, 35, 44, 9, 45, 23, 48, 18,
               19, 41, 7, 37, 62, 20, 29, 22),
  rename = seq(65)
)

all_stat = map(
  list(immune_all_stat,
       angiogenesis_all_stat,
       combined_all_stat),
  function(x) { x %>%
      left_join(F_module_rename, by = join_by(module_num == original)) %>%
      mutate(rename = paste0("F_GM", rename))
  }
)

tissue_stat = map(
  list(immune_tissue_stat,
       angiogenesis_tissue_stat,
       combined_tissue_stat),
  function(x) { x %>%
      left_join(F_module_rename, by = join_by(module_num == original)) %>%
      mutate(rename = paste0("F_GM", rename)) %>%
      mutate(tissue = str_to_title(tissue))
  }
)

all_stat[[3]]$rename = factor(all_stat[[3]]$rename,
                              levels = c("F_GM10", "F_GM15", "F_GM19", "F_GM13", "F_GM18", "F_GM16",
                                         "F_GM1", "F_GM17", "F_GM8", "F_GM12", "F_GM2", "F_GM14",
                                         "F_GM5", "F_GM3", "F_GM6", "F_GM4", "F_GM23", "F_GM7",
                                         "F_GM9", "F_GM24", "F_GM22", "F_GM11", "F_GM20", "F_GM21",
                                         "F_GM33", "F_GM34", "F_GM32"))

# plot all stat
p1 = imap(all_stat,
         function(x, y) {

  if(y < 3) {
    p = x %>%  subset(value == "TRUE") %>% ggplot(aes(reorder(rename, ratio), ratio))
  }
  else{
    p = x %>%  subset(value == "TRUE") %>% ggplot(aes(rename, ratio))
  }
  p +
  geom_bar(stat = "identity", fill = "#C4DFE6", color = "black", width = 0.5) +
  geom_hline(yintercept = 0.7, colour = "black", linetype = "dotted") +
  # geom_vline(xintercept = 24.5, colour = "black", size = 0.5) +
  scale_y_continuous(labels = ~ scales::percent(.x)) +
  ylab("Percentage") +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", family = "serif", size = 10),
    axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.y.left = element_text(family = "serif", size = 15),
    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_line(colour = "black", linewidth = 1),
    axis.ticks.y.left = element_line(colour = "black", linewidth = 1),
    panel.border = element_rect(color = "black", linewidth  = 1)
  )
         }
)


ggsave(plot = p1[[1]], filename = "The perenctage of auccurate each F immune module scoring.pdf", width = 3000, height = 800, units = "px")
ggsave(plot = p1[[2]], filename = "The perenctage of auccurate each F angiogenesis module scoring.pdf", width = 3000, height = 800, units = "px")
ggsave(plot = p1[[3]], filename = "The perenctage of auccurate each combined module scoring.pdf", width = 3000, height = 800, units = "px")


new_colors = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A", "#FF4E50")
new_colors = colorRampPalette(new_colors)(100)

all_tissues = map(tissue_stat, ~ unique(.x %>% subset(value == "TRUE") %>% pull(tissue)))

tissue_stat2 = map2(tissue_stat,
                   all_tissues,
                   function(x, y) {
                     x %>%
  subset(value == "TRUE") %>%
  group_by(module) %>%
  group_split() %>%
  map_dfr(
    function(x) {
      excluded_tissues = y[! y %in% unique(x$tissue)]
      if(length(excluded_tissues) > 0) {
        new_df = map_dfr(
          excluded_tissues,
          function(z) {
            data.frame(
              module = unique(x$module),
              organization = z,
              tissue = z,
              value = "TRUE",
              ratio = 0,
              module_num = NA,
              rename = unique(x$rename)
            )
          }
        )
        bind_rows(x, new_df)
      }
      else {
        return(x)
      }
    }
  )
                   }
)

tissue_stat2[[3]]$rename = factor(
  tissue_stat2[[3]]$rename,
  levels = c("F_GM10", "F_GM15", "F_GM19", "F_GM13", "F_GM18", "F_GM16",
             "F_GM1", "F_GM17", "F_GM8", "F_GM12", "F_GM2", "F_GM14",
             "F_GM5", "F_GM3", "F_GM6", "F_GM4", "F_GM23", "F_GM7",
             "F_GM9", "F_GM24", "F_GM22", "F_GM11", "F_GM20", "F_GM21",
             "F_GM33", "F_GM34", "F_GM32")
)

tissue_stat3 =  imap(
  tissue_stat2,
  function(x, y) {
    x$tissue =  factor(x$tissue, level = str_sort(unique(x$tissue), numeric = TRUE))
    if(length(unique(x$rename)) > 3 & length(unique(x$rename)) < 25) {
      x$rename = factor(
        x$rename,
        levels = c("F_GM10", "F_GM15", "F_GM19", "F_GM13", "F_GM18", "F_GM16",
                   "F_GM1", "F_GM17", "F_GM8", "F_GM12", "F_GM2", "F_GM14",
                   "F_GM5", "F_GM3", "F_GM6", "F_GM4", "F_GM23", "F_GM7",
                   "F_GM9", "F_GM24", "F_GM22", "F_GM11", "F_GM20", "F_GM21")
      )
    }
    else{
      x$rename = factor(x$rename, levels = c("F_GM32", "F_GM34", "F_GM33"))
    }
    if(y == 3){
      x$rename = factor(
        x$rename,
        levels = c("F_GM10", "F_GM15", "F_GM19", "F_GM13", "F_GM18", "F_GM16",
                   "F_GM1", "F_GM17", "F_GM8", "F_GM12", "F_GM2", "F_GM14",
                   "F_GM5", "F_GM3", "F_GM6", "F_GM4", "F_GM23", "F_GM7",
                   "F_GM9", "F_GM24", "F_GM22", "F_GM11", "F_GM20", "F_GM21",
                   "F_GM33", "F_GM34", "F_GM32")
      )
    }

    return(x)
  }
)


p2 = map(tissue_stat2,
         function(x) {
           x %>%
  subset(value == "TRUE") %>%
  ggplot(aes(rename, tissue, fill = ratio)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradientn(colors = new_colors, labels = ~ scales::percent(.x), name = "Percentage") +
  # coord_fixed() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(family = "serif", colour = "black"),
    axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 1),
    legend.ticks = element_line(linewidth = 0.5, colour = "white"),
    legend.text = element_text(family = "serif", colour = "black"),
    legend.title = element_text(family = "serif", colour = "black")
  )
         }
)

ggsave(plot = p2[[1]], filename = "The perenctage of auccurate each F immune module scoring per tissue.pdf", width = 3000, height = 2000, units = "px")
ggsave(plot = p2[[2]], filename = "The perenctage of auccurate each F angiogenesis module scoring per tissue.pdf", width = 3000, height = 2000, units = "px")
ggsave(plot = p2[[3]], filename = "The perenctage of auccurate each combined module scoring per tissue.pdf", width = 3000, height = 2000, units = "px")
