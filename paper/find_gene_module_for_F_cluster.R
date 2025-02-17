library(tidyverse)
library(snow)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(VennDiagram)
library(ggsci)
library(dendextend)
library(geneModule)

setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/F_cluster/")

############################## finding gene module on Fcluster ##############################
#                                                                                           #
#              markers of each Fcluster were retrieve by default setting                    #
#                                                                                           #
#############################################################################################

# read top 50 makers per F cluster per study
data = read.csv("Fcluster_gene_lists_with_default_retrieve_markers_settings.csv", fileEncoding = "GBK")
data = data[colnames(data)!= "X"]
data = data %>% mutate(across(everything(), as.character))

genes = data %>% mutate(across(everything(), ~ str_extract(.x, "[^/]+")))
values = data %>% mutate(across(everything(), ~ str_remove(str_extract(.x, "(?=/).+"), "/")))

# read numbers2ids file in
num2id = read.csv("gene_name_list.csv", header = FALSE)
genes_for_trans = genes %>% dplyr::select(where(~ all(str_detect(.x, "^\\d+"))))
values_for_trans = values[, colnames(values) %in% colnames(genes_for_trans)]
genes_for_trans2 = genes_for_trans %>% mutate(across(everything(), ~ num2id[as.integer(.x)+1, "V2"]))

new_genes = genes[, ! colnames(genes) %in% colnames(genes_for_trans2)]
new_values = values[, ! colnames(values) %in% colnames(genes_for_trans2)]

new_data1 = map_dfc(
  seq(dim(new_genes)[2]),
  function(x) {
    name = colnames(new_genes)[x]
    gene = new_genes[, x, drop = TRUE]
    value = new_values[, x, drop = TRUE]
    df = data.frame(X1 = paste0(gene, "/", value))
    colnames(df) = name
    df
  }
)

new_data2 = map_dfc(
  seq(dim(genes_for_trans2)[2]),
  function(x) {
    name = colnames(genes_for_trans2)[x]
    gene = genes_for_trans2[, x, drop = TRUE]
    value = values_for_trans[, x, drop = TRUE]
    df = data.frame(X1 = paste0(gene, "/", value))
    colnames(df) = name
    df
  }
)

final_data = bind_cols(new_data1, new_data2)[, colnames(data)]

# get the entrez, ensemble ids and the corresponding gene names
id2symbol = AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = c("ALIAS", "ENSEMBL"))

# convert gene symbol to ensembl id
ensemble_id_files = map(
  final_data,
  function(z) {
    x = str_remove(z, "/.+")
    y = str_remove(z, "[^/]+/")
    if (all(str_like(x, "^ENSG\\d+"))) {
      return(paste0(x, "/", y))
    }
    else if (any(str_like(x, "^\\d+"))) {
      last = function(x){x[[length(x)]]}
      return(paste0(mapIds(org.Hs.eg.db, keys = x, column = "ENSEMBL", keytype = "ENTREZID", multiVals = last), "/", y))
    }
    else {
      id2symbol = AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = c("ALIAS", "ENSEMBL")) %>%
        dplyr::select(-ENTREZID)
      intermediate = data.frame(x = x) %>% left_join(id2symbol, by = join_by(x == ALIAS), multiple = "first")
      return(paste0(intermediate$ENSEMBL, "/", y))

    }
  }
)

ensemble_id_files = bind_cols(ensemble_id_files)

# discard columns with NA
ensemble_id_files = ensemble_id_files[, map_lgl(ensemble_id_files, ~ ! any(str_like(.x, "NA.+")))]
saveRDS(ensemble_id_files, "Fcluster_gene_lists_with_default_retrieve_markers_settings_ensemble_id_files.rds")

# run get_gene_module
gene_module = get_gene_module(ensemble_id_files, parallel_num = 32)
saveRDS(gene_module, "F_cluster_with_default_retrieve_markers_settings_gene_module_first_multi.rds")

# read the gene module file
gene_module = readRDS("F_cluster_with_default_retrieve_markers_settings_gene_module_first_multi.rds")
ensemble_id_files = readRDS("Fcluster_gene_lists_with_default_retrieve_markers_settings_ensemble_id_files.rds")
ensemble_id_files = ensemble_id_files[, map_lgl(ensemble_id_files, ~ ! any(str_like(.x, "NA.+")))]

# calculate jaccard similarity matrix
data2 = ensemble_id_files %>% mutate(across(everything(), ~ str_remove(.x, "/.+")))

jaccard_matrix = do_jaccard_calculation(data2, parallel_num = 32)
jaccard_matrix = jaccard_matrix %>%
  pivot_wider(everything(), values_from = "jaccard_value", names_from = "compare_name") %>%
  column_to_rownames("base_name")

# output the gene module lists
gene_module_lists = gene_module[[1]][map_lgl(gene_module[[2]], ~ length(.x) > 2)]
gene_module_lists = map(gene_module_lists, ~ str_remove(.x, "/.+"))
names(gene_module_lists) = paste0("gene_module_", seq(length(gene_module_lists)))
gene_module_lists = as.data.frame(gene_module_lists)
write.csv(gene_module_lists, "F_cluster_with_default_retrieve_markers_settings_gene_module_list.csv", fileEncoding = "utf-8")

# do enrichment analysis for each gene module
enrich_data = map(
  gene_module_lists,
  ~ AnnotationDbi::mapIds(org.Hs.eg.db, keys = .x, column = c("ENTREZID"), keytype = "ENSEMBL")
)

enrichment = map(
  enrich_data,
  function(x) {
    enrichGO(x, org.Hs.eg.db, ont = "ALL")
  }
)

sig_go_enrichment = map(enrichment, ~ .x@result %>% subset(p.adjust < 0.05))
openxlsx::write.xlsx(sig_go_enrichment, "sig_go_enrichment.xlsx")


# cluster gene module based on GO enrichment results
immune_response_module = c(2, 11, 13, 15, 16, 21, 25, 27, 33, 38, 39, 40, 43, 46, 47, 51, 52, 53, 55, 56, 58, 61, 63, 65)
cell_cycle_module = c(4, 17, 28, 32, 59)
angiogenesis_module = c(1, 30, 36)
cytoskeleton_module = c(5, 8, 12, 24, 26, 34, 50, 64)
energy_synthesis_module = c(18, 19, 23, 45, 48)
antioxidant_stress_module = c(7, 20, 37, 41, 62)
cell_adhesion_module = c(6, 14)
mesenchymal_module = c(3, 10, 49, 54, 57)
protein_synthesis_membrane_raft_module = c(9, 35, 44, 60)
proteolysis_module = c(31, 42)
unassign_module = c(22, 29)

modules = list(immune_response_module, cell_cycle_module, angiogenesis_module, cytoskeleton_module,
            energy_synthesis_module, antioxidant_stress_module, cell_adhesion_module, mesenchymal_module,
            proteolysis_module, protein_synthesis_membrane_raft_module, unassign_module)
names(modules) = c("immune_response", "cell_cycle", "angiogenesis", "cytoskeleton",
                   "energy_synthesis", "antioxidant/stress", "cell_adhesion", "mesenchymal",
                   "protein_synthesis", "proteolysis", " unassigned")

modules = modules[c(1, 7, 2, 3,  4, 8, 9, 10, 5, 6, 11)]

# only retain the members contributing to gene modules with above 2 members
robust_candidates = gene_module[[2]][map_lgl(gene_module[[2]], ~ length(.x) > 2)] %>% unlist()
jaccard_matrix_sub = jaccard_matrix[, robust_candidates]
cond = map_int(robust_candidates, ~ which(rownames(jaccard_matrix_sub) %in% .x))
jaccard_matrix_sub = jaccard_matrix_sub[cond, robust_candidates]

# calculate the overlap between gene sets
intersection = find_overlap_singlet(gene_module_lists, full_search = T)

# recoder module names by its
# final order in the heatmap of Fig.3C
jaccard_similarity = do_jaccard_calculation(gene_module_lists)
jaccard_similarity$compare_num = as.numeric(str_extract(jaccard_similarity$compare_name, "\\d+"))
jaccard_similarity$base_num = as.numeric(str_extract(jaccard_similarity$base_name, "\\d+"))

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

recode_immune_resposne_module = c(1:24)
recode_cell_adhesion_module = c(25, 26)
recode_cell_cycle_module = c(27:31)
recode_angiogenesis_module = c(32:34)
recode_cytoskeleton_module = c(35:42)
recode_mesenchymal_module = c(43:47)
recode_proteolysis_module = c(48, 49)
recode_protein_synthesis_memebrane_raft_module = c(50:53)
recode_energy_synthesis_module = c(54:58)
recode_antioxidant_stress_module = c(59:63)
recode_unassigned_module = c(64, 65)

jaccard_similarity = jaccard_similarity %>% left_join(F_module_rename, by = join_by(compare_num == original))
jaccard_similarity = jaccard_similarity %>% left_join(F_module_rename, by = join_by(base_num == original))

jaccard_similarity = jaccard_similarity %>% dplyr::select(jaccard_value, rename.x, rename.y)
colnames(jaccard_similarity)[c(2, 3)] = c("compare_name", "base_name")

jaccard_similarity = jaccard_similarity %>%
  dplyr::select(compare_name, base_name, jaccard_value) %>%
  pivot_wider(everything(), names_from = "compare_name", values_from = "jaccard_value") %>%
  column_to_rownames("base_name")

label_order = hclust(dist(jaccard_similarity), method = "average")$order
label_color = map_chr(
  colnames(jaccard_similarity)[label_order],
  function(x) {
    if(x %in% recode_immune_resposne_module) {
      rgb(249/255, 54/255, 229/255, alpha = 0.5)
    }
    else if(x %in% recode_cell_cycle_module) {
      rgb(99/255, 202/255, 96/255, alpha = 0.5)
    }
    else if(x %in% recode_angiogenesis_module) {
      rgb(146/255, 29/255, 176/255, alpha = 0.5)
    }
    else if(x %in% recode_cytoskeleton_module) {
      rgb(54/255, 181/255, 196/255, alpha = 0.5)
    }
    else if(x %in% recode_energy_synthesis_module) {
      rgb(72/255, 224/255, 195/255, alpha = 0.5)
    }
    else if(x %in% recode_antioxidant_stress_module) {
      rgb(153/255, 91/255, 99/255, alpha = 0.5)
    }
    else if(x %in% recode_cell_adhesion_module) {
      rgb(68/255, 121/255, 3/255, alpha = 0.5)
    }
    else if(x %in% recode_mesenchymal_module) {
      rgb(156/255, 35/255, 50/255, alpha = 0.5)
    }
    else if(x %in% recode_protein_synthesis_memebrane_raft_module) {
      rgb(57/255, 159/255, 223/255, alpha = 0.5)
    }
    else if(x %in% recode_proteolysis_module) {
      rgb(165/255, 53/255, 224/255, alpha = 0.5)
    }
    else if(x %in% recode_unassigned_module) {
      rgb(110/255, 213/255, 210/255, alpha = 0.5)
    }
  }
)

dend = hclust(dist(jaccard_similarity), method = "average") %>% as.dendrogram()
dend = dend %>%
  set("labels", paste0(colnames(jaccard_similarity)[label_order])) %>%
  set("labels_colors", "black") %>%
  set("labels_cex", 0.5) %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 6) %>%  # node point size
  set("leaves_col", label_color) %>%
  as.ggdend()

p1 = ggplot(dend, labels = T) +
  scale_y_reverse(expand = c(0, 0.25)) +
  coord_polar(theta="x") +
  theme(
    text = element_text(family = "serif"),
    title = element_text(family = "serif")
  )

ggsave(plot = p1, filename = "F_clusters_heirarchy_dendrogram_20241218.pdf", width = 4000,
       height = 2000, units = "px")

# reorder gene module using heirarchy clustering
# based on the overlap between gene sets
intersection = intersection %>%
  pivot_wider(everything(), names_from = "compare_name", values_from = "number") %>%
  column_to_rownames("base_name")

between_module_reorder = map(
  modules,
  function(x) {
    sub_data = intersection[x, x]
    model = hclust(dist(sub_data), method = "average")
    as.integer(str_extract(colnames(sub_data)[rev(model$order)], "\\d+"))
  }
)

# reorder the members in each gene module using heriarchy
reorder_candidates = gene_module[[2]][map_lgl(gene_module[[2]], ~ length(.x) > 2)][unlist(between_module_reorder)]
within_module_reorder = map(
  reorder_candidates,
  function(x) {
    sub_data = jaccard_matrix_sub[, x]
    model = hclust(dist(t(sub_data)))
    colnames(sub_data)[rev(model$order)]
  }
)

cond2 = map_int(unlist(within_module_reorder), ~ which(rownames(jaccard_matrix_sub) %in% .x))

jaccard_matrix_sub_reorder = jaccard_matrix_sub[cond2, unlist(within_module_reorder)]

# generate the annotation side bar data frame for heatmap
robust_candidates = gene_module[[2]][map_lgl(gene_module[[2]], ~ length(.x) > 2)]
annotation_col = map2_dfr(modules, names(modules), function(x, y) {
  map_dfr(
    x,
    function(z) {
      intermediate = robust_candidates[z]
      res = unlist(intermediate)
      df = tibble(func = rep(y, length(res)))
      df$cell = res
      df$module = paste0("module_", z)
      df = df %>% column_to_rownames("cell")
      return(df)
    }
  )
})

new_colors = c(brewer.pal(11, "Set3"))
names(new_colors) = names(modules)
color_list = list(func = new_colors)


colors = rev(colorRampPalette(brewer.pal(11, "Spectral"))(100))
custom_magma = c(colorRampPalette(c("white", rev(magma(200, begin = 0.15))[1]))(10), rev(magma(200, begin = 0.18)))
png("Fcluster_gene_module_heatmap.png", width = 2000, height = 1500, res = 300)
pheatmap::pheatmap(jaccard_matrix_sub_reorder,
                   annotation_col = annotation_col,
                   annotation_colors = color_list,
                   show_colnames = F, show_rownames = F, color = colors, cluster_rows =
                     F, cluster_cols = F)
dev.off()

# better heatmap visualization by ComplexHeatMap
new_module = map(unique(annotation_col$module), function(z) {
  data = annotation_col %>% subset(module == z) %>% pull(module)
  return(c(paste0("/", str_extract(z, "\\d+"), "*"), rep(NA, (length(data) - 1))))
}
)

annotation_col$new_module = unlist(new_module)
annotation_col = annotation_col[map_int(rownames(jaccard_matrix_sub_reorder), ~ which(rownames(annotation_col) %in% .x)), ]

colors2 = ggsci::pal_npg("nrc")(10)
colors3 = annotation_col %>% mutate(new_color = case_when(
  annotation_col$func == "immune_response" ~ colors2[1],
  annotation_col$func == "cell_adhesion" ~ colors2[2],
  annotation_col$func == "cell_cycle" ~ colors2[3],
  annotation_col$func == "angiogenesis" ~ colors2[4],
  annotation_col$func == "cytoskeleton" ~ colors2[5],
  annotation_col$func == "protein_synthesis" ~ colors2[6],
  annotation_col$func == "energy_synthesis" ~ colors2[7],
  str_detect(annotation_col$func, "antioxidant.+") ~ colors2[9],
  annotation_col$func == "mesenchymal" ~ colors2[10],
  annotation_col$func == "proteolysis" ~ "grey60",
  .default =  "black"))

colors4 = colors3$new_color
names(colors4) = colors3$func

ha = HeatmapAnnotation(m = anno_simple(annotation_col$module,
                                       pch = annotation_col$new_module,
                                       pt_size = unit(0.2, "snpc")),
                       f = anno_simple(annotation_col$func, col = colors4))

row = rowAnnotation(m = anno_simple(annotation_col$module,
                                    pch = annotation_col$new_module,
                                    pt_size = unit(0.2, "snpc"),
), f = anno_simple(annotation_col$func))

Heatmap(jaccard_matrix_sub_reorder, name = "mat", col = colors,
        top_annotation = ha,
        left_annotation = row,
        cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F)


# calculate the ratio of health/disease samples in each gene module
disease_info = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/作者年_表达矩阵_metadata_pkl/作者年_表达矩阵_metadata_pkl.csv", fileEncoding = "GBK")

sample_ratio = map(robust_candidates, ~ str_remove(.x, "_cluster\\d+"))
sample_ratio = imap_dfr(sample_ratio, function(x, idx) data.frame(member = x, module = idx))
sample_ratio$member = ifelse(str_detect(sample_ratio$member, "K眉rten_2021"), "Kürten_2021", sample_ratio$member)
sample_ratio$member = str_replace(sample_ratio$member, "\\.", "-")

sample_ratio = sample_ratio %>%
  left_join(disease_info, by = join_by(member == new)) %>%
  dplyr::select(member, diseases, organization, module) %>%
  mutate(diseases = case_when(diseases != "Normal" ~ "disease",
                              .default = diseases))

sample_ratio = sample_ratio %>% mutate(
  module_type = case_when(
    module %in% immune_response_module ~ "immune_response",
  module %in% cell_cycle_module ~ "cell_cycle",
  module %in% angiogenesis_module ~ "angiogenesis",
  module %in% cytoskeleton_module ~ "cytoskeleton",
  module %in% energy_synthesis_module ~ "energy_synthesis",
  module %in% antioxidant_stress_module ~ "antioxidant_stress",
  module %in% cell_adhesion_module ~ "cell_adhesion",
  module %in% mesenchymal_module ~ "mesenchymal",
  module %in% protein_synthesis_membrane_raft_module ~ "protein_synthesis/membrane_raft",
  module %in% proteolysis_module ~ "proteolysis",
  .default = "unassigend"
  )
)

# calculate the ratio of each module type
module_ratio = sample_ratio %>% count(module_type) %>% mutate(percentage = n / sum(n))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size=14, face="bold")
  )

pie_colors = c(pal_npg("nrc")(10), "grey50")

new_colors = c("#E64B35FF", "#3C5488FF", "#F39B7FFF", "#00A087FF", "#91D1C2FF", "#8491B4FF",
               "#4DBBD5FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "black")
names(new_colors) = c("immune_response", "angiogenesis", "cytoskeleton", "cell_cycle",
                      "energy_synthesis", "proteolysis", "cell_adhesion",
                      "protein_synthesis/membrane_raft","antioxidant_stress",
                      "mesenchymal", "unassigend")


pie <- ggplot(module_ratio, aes(x="", y=percentage, fill= module_type)) +
  geom_bar(width = 1, stat = "identity", position = "stack") +
  # scale_fill_npg() +
  scale_fill_manual(values = new_colors) +
  coord_polar("y", start=0)+
  blank_theme #+
  # theme(axis.text.x=element_blank(), legend.position = "none")

ggsave(filename = "pie_all.png", pie, width = 1000, height = 1000, units = "px")

# calculate the ratio of each module type in
# lung , blood, breast, intestine, brain and airway
tissue_module_ratio =
  map(
    c("lung", "blood", "breast", "intestine", "airway", "brain"),
    function(x) {
    sample_ratio %>%
    subset(organization == x) %>%
    count(module_type) %>%
    mutate(percentage = n / sum(n))
  }
)


map2(
  tissue_module_ratio,
  c("lung", "blood", "breast", "intestine", "airway", "brain"),
  function(x, y) {
   fig =  ggplot(x, aes(x="", y=percentage, fill= module_type)) +
    geom_bar(width = 1, stat = "identity", position = "stack") +
    # scale_fill_npg() +
    scale_fill_manual(values = new_colors) +
    coord_polar("y", start=0)+
    blank_theme +
    theme(axis.text.x=element_blank(), legend.position = "none")

   ggsave(filename = paste0("pie_", y, ".png"), fig, width = 1000, height = 1000, units = "px")

   }
)
