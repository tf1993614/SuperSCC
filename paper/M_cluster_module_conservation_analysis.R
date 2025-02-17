library(tidyverse)
library(openxlsx)

setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/M_cluster")

immune_data = openxlsx::read.xlsx("/mnt/disk5/zhongmin/superscc/结果位置/addmodulescore/师兄算的结果/addmodule/M_gene_module/addmodulescore_m_gene_module.xlsx")
angiogenesis_data = openxlsx::read.xlsx("/mnt/disk5/zhongmin/superscc/结果位置/addmodulescore/师兄算的结果/addmodule/M_gene_module/addmodulescore_m_gene_module_endothelial.xlsx")
combine_data = list(immune = immune_data, angiogenesis = angiogenesis_data)

disease_info = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/作者年_表达矩阵_metadata_pkl/作者年_表达矩阵_metadata_pkl.csv", fileEncoding = "GBK")

combine_data = map(
  combine_data,
  ~ .x %>% left_join(disease_info %>% dplyr::select(new, organization, diseases), by = "new")
)


combine_data = map(
  combine_data,
  ~ .x %>%
  pivot_longer(
      starts_with("M_gene_module"),
      names_to = "module",
      values_to = "value"
  )
)


all_stat = map(
  combine_data,
  ~ .x %>%
    subset(value != "NO") %>%
    group_by(module) %>%
    count(value) %>%
    group_by(module) %>%
    summarise(module = module, value = value, ratio = n / sum(n)) %>%
    ungroup() %>%
    mutate(module_num = as.numeric(str_extract(module, "\\d+")))
)


tissue_stat = map(
  combine_data,
  ~ .x %>%
    subset(value != "NO") %>%
    group_by(module, organization) %>%
    count(value) %>%
    group_by(module, organization) %>%
    summarise(module = module, tissue = organization, value = value, ratio = n / sum(n)) %>%
    ungroup() %>%
    mutate(module_num = as.numeric(str_extract(module, "\\d+")))
)

M_module_rename = data.frame(
  original = c(15, 13, 7, 3, 14, 5, 6, 11,
               12, 8, 9, 1, 10,  2, 4),
  rename = seq(15)
)

all_stat = map(
  all_stat,
  ~ .x %>%
    left_join(M_module_rename, by = join_by(module_num == original)) %>%
    mutate(rename = paste0("F_GM", rename))
)


tissue_stat = map(
  tissue_stat,
  ~ .x %>%
    left_join(M_module_rename, by = join_by(module_num == original)) %>%
    mutate(rename = paste0("F_GM", rename)) %>%
    mutate(tissue = str_to_title(tissue))
)

combine_all_stat = bind_rows(all_stat)
combine_all_stat$rename = factor(combine_all_stat$rename,
                                 levels = c("F_GM7", "F_GM1", "F_GM8", "F_GM3",
                                            "F_GM2", "F_GM4", "F_GM5", "F_GM6",
                                            "F_GM12", "F_GM11"))

combine_tissue_stat = bind_rows(tissue_stat)
combine_tissue_stat$rename = factor(combine_tissue_stat$rename,
                                 levels = c("F_GM7", "F_GM1", "F_GM8", "F_GM3",
                                            "F_GM2", "F_GM4", "F_GM5", "F_GM6",
                                            "F_GM12", "F_GM11"))

# plot all stat
p1 = combine_all_stat %>%
    subset(value == "TRUE") %>%
    ggplot(aes(rename,ratio)) +
    geom_bar(stat = "identity", fill = "#C4DFE6", color = "black") +
    geom_hline(yintercept = 0.7, colour = "black", linetype = "dotted") +
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
)


ggsave(plot = p1, filename = "The perenctage of auccurate each M module scoring.pdf", width = 3000, height = 2000, units = "px")


# plot tissue state
colors = c("#FBDD49", "#FF8103", "#FF1C6A", "#E200A3", "#9B04DB", "#6D1DC6", "white")
colors = rev(colorRampPalette(colors)(100))

new_colors = c("#E1F5C4", "#EDE574", "#F9D423", "#FC913A", "#FF4E50")
new_colors = colorRampPalette(new_colors)(100)

all_tissues = unique(combine_tissue_stat$tissue)

combine_tissue_stat = combine_tissue_stat %>%
  subset(value == "TRUE") %>%
  group_by(module) %>%
  group_split() %>%
  map_dfr(
    function(x) {
      excluded_tissues = all_tissues[! all_tissues %in% x$tissue]
      if(length(excluded_tissues) > 0) {
        new_df = map_dfr(
          excluded_tissues,
          function(y) {
            data.frame(
              module = unique(x$module),
              organization = y,
              tissue = y,
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

p2 = combine_tissue_stat %>%
  subset(value == "TRUE") %>%
  ggplot(aes(rename, tissue, fill = ratio)) +
  geom_tile(color = "black", size = 0.5) +
  scale_fill_gradientn(colors = new_colors, labels = ~ scales::percent(.x), name = "Percentage") +
  coord_fixed() +
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

ggsave(plot = p2, filename = "The perenctage of auccurate each M module scoring per tissue.pdf", width = 3000, height = 2000, units = "px")
