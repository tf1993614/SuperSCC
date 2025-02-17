library(tidyverse)

setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/F_cluster/")


F_module_scores = readRDS("F_cluster_all_cell_scores_update.rds")
M_module_scores = readRDS("../M_cluster/M_cluster_all_gene_modules_cell_scores.rds")

F_immune_response_module = c(2, 11, 13, 15, 16, 21, 25, 27, 33, 38, 39, 40, 43, 46, 47, 51, 52, 53, 55, 56, 58, 61, 63, 65)
F_angiogenesis_module = c(1, 30, 36)

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


immune_module_score = map(
  immune_response_module,
  function(x) {
    map(
      module_scores,
      function(y) {
       dataset = y
       score = as.data.frame(dataset[[x]][[1]])
       df = data.frame(V1 = score[["/1"]])
       colnames(df) = paste0("GM_", x)
       df
      }
    )
  }
)


angiogenesis_module_score = map(
  angiogenesis_module,
  function(x) {
    map(
      module_scores,
      function(y) {
        dataset = y
        score = as.data.frame(dataset[[x]][[1]])
        df = data.frame(V1 = score[["/1"]])
        colnames(df) = paste0("GM_", x)
        df
      }
    )
  }
)

dataset_name = names(module_scores)

pearson_res = map_dfr(
  seq(99),
  function(num) {
      imap_dfr(

      immune_response_module,
      function(x, y) {
        imap_dfr(
          angiogenesis_module,
          function(z, index) {

            compare1 = immune_module_score[[y]][[num]][[paste0("GM_", x)]]
            compare2 = angiogenesis_module_score[[index]][[num]][[paste0("GM_", z)]]

            immune_GM = paste0("GM_", x)
            angiogenesis_GM = paste0("GM_", z)

            res = broom::tidy(cor.test(compare1, compare2))

            res$immune_GM = immune_GM
            res$angiogenesis_GM = angiogenesis_GM
            res$dataset = dataset_name[num]
            res = res %>% dplyr::select(dataset, immune_GM, angiogenesis_GM, everything())
          }
        )
      }
    )
  }
)


pearson_res_per_angiogenesis_module = pearson_res %>%
  group_by(angiogenesis_GM) %>%
  group_split()


x1 = immune_module_score[[7]][[67]][[1]]
x2 = angiogenesis_module_score[[1]][[67]][[1]]

df = data.frame(x1 = x1, x2 = x2)


ggscatter(df, x = "x1", y = "x2",
          fill = "lightgrey",
          color = "black", shape = 21, size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "#E64B35FF", fill = "grey"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson",  label.sep = "\n",
                                label.x.npc = "center", label.y.npc = "top"
                                )
) +
  labs(x = "immune_GM_6", y = "angiogenesis_GM_33") +
  ggtitle("Dong_2020_2")
