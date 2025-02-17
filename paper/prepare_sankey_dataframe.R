library(tidyverse)
library(openxlsx)

m2f_info = map(seq(65), ~ openxlsx::read.xlsx("/mnt/disk5/zhongmin/superscc/结果位置/gene_module/F_M_gene_module/M_F_gene_module_modified.xlsx", sheet = .x))

# cluster F modules based on
# enrichment result
F_immune_response_module = c(2, 11, 13, 15, 16, 21, 25, 27, 33, 38, 39, 40, 43, 46, 47, 51, 52, 53, 55, 56, 58, 61, 63, 65)
F_cell_cycle_module = c(4, 17, 28, 32, 59)
F_angiogenesis_module = c(1, 30, 36)
F_cytoskeleton_module = c(5, 8, 12, 24, 26, 34, 50, 64)
F_energy_synthesis_module = c(18, 19, 23, 45, 48)
F_antioxidant_stress_module = c(7, 20, 37, 41, 62)
F_cell_adhesion_module = c(6, 14)
F_mesenchymal_module = c(3, 10, 49, 54, 57)
F_protein_synthesis_membrane_raft_module = c(9, 35, 44, 60)
F_proteolysis_module = c(31, 42)
F_unassign_module = c(22, 29)

F_modules = list(
  F_immune_response_module,
  F_cell_cycle_module,
  F_angiogenesis_module,
  F_cytoskeleton_module,
  F_energy_synthesis_module,
  F_antioxidant_stress_module,
  F_cell_adhesion_module,
  F_mesenchymal_module,
  F_protein_synthesis_membrane_raft_module,
  F_proteolysis_module,
  F_unassign_module
)

# cluster M modules based on
# enrichment result
M_immune_response_module = c(3, 5, 6, 7, 11, 13, 14, 15)
M_cell_development_module = c(8, 12)
M_angiogenesis_module = c(1, 9)
M_cytoskeleton_module = c(10)
M_defense_response_module = c(2)
M_protein_target_module = c(4)

M_modules = list(
  M_immune_response_module,
  M_cell_development_module,
  M_angiogenesis_module,
  M_cytoskeleton_module,
  M_defense_response_module,
  M_protein_target_module
)

sankey_df = map_dfr(
  m2f_info,
  function(x) {
    df = x %>% count(M_gene_module)
    df$F_gene_module = unique(x$F_gene_module)
    df = df %>% dplyr::select(M_gene_module, F_gene_module, n)
    df$M_num = str_extract(df$M_gene_module, "\\d+")
    df$F_num = str_extract(df$F_gene_module, "\\d+")
    df
  }
)

sankey_df = sankey_df %>%
  mutate(
    F_func = case_when(
      F_num %in% F_immune_response_module ~ "F_immune_response",
      F_num %in% F_cell_cycle_module ~ "F_cell_cycle",
      F_num %in% F_cytoskeleton_module ~ "F_cytoskeleton",
      F_num %in% F_energy_synthesis_module ~ "F_energy_synthesis",
      F_num %in% F_antioxidant_stress_module ~ "F_antioxidant/stress",
      F_num %in% F_cell_adhesion_module ~ "F_cell_adhesion",
      F_num %in% F_mesenchymal_module ~ "F_mesenchymal",
      F_num %in% F_protein_synthesis_membrane_raft_module ~ "F_protein_synthesis/membrane_raft",
      F_num %in% F_proteolysis_module ~ "F_proteolysis",
      F_num %in% F_angiogenesis_module ~ "F_angiogenesis",
      F_num %in% F_unassign_module ~ "F_unassigned",
      .default = F_num
  )
)

sankey_df = sankey_df %>%
  mutate(
    M_func = case_when(
      M_num %in% M_immune_response_module ~ "M_immune_response",
      M_num %in% M_cell_development_module ~ "M_cell_development",
      M_num %in% M_cytoskeleton_module ~ "M_cytoskeleton",
      M_num %in% M_angiogenesis_module ~ "M_angiogenesis",
      M_num %in% M_defense_response_module ~ "M_defense_response",
      M_num %in% M_protein_target_module ~ "M_protein_target",
      .default = M_num
    )
  )

sankey_df$M_func = if_else(is.na(sankey_df$M_func), "Unknown source", sankey_df$M_func)

# renmae the F modules based on
# the order of each module in the
# diagonal of the heatmap of Fig3C
F_module_rename = data.frame(
  original = c(33, 21, 27, 56, 38, 25, 65, 46,
               58, 51, 55, 52, 61, 16, 15, 53,
               13, 40, 2, 43, 11, 47, 39, 63,
               14, 6, 59, 28, 32, 17, 4, 30,
               1, 36, 8, 5, 64, 24, 12, 34, 50,
               26, 57, 3, 54, 49, 10, 42, 31,
               60, 35, 44, 9, 45, 23, 48, 18,
               19, 41, 7, 37, 62, 20, 29, 22),
  F_rename = seq(65)
)

F_module_rename = F_module_rename %>% mutate(across(everything(), as.character))

sankey_df = sankey_df %>% left_join(F_module_rename, by = join_by(F_num == original))

# renmae the M modules based on
# the order of each module in the
# diagonal of the heatmap of Fig3C
M_module_rename = data.frame(
  original = c(15, 13, 7, 3, 14, 5, 6, 11, 12, 8,
               9, 1, 10, 2, 4),
  M_rename = seq(15)
)


M_module_rename = M_module_rename %>% mutate(across(everything(), as.character))

sankey_df = sankey_df %>% left_join(M_module_rename, by = join_by(M_num == original))

write.csv(sankey_df, "sankey_df.csv")
