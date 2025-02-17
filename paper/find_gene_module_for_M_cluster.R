library(readxl)
library(tidyverse)
library(snow)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(RColorBrewer)
library(cluster)
library(viridis)
library(geneModule)

# set working directory
setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/M_cluster/")

############################## finding gene module on Mcluster ##############################
#                                                                                           #
#              markers of each Fcluster were retrieve by default setting                    #
#                                                                                           #
#############################################################################################

# read data
data = read.csv("~/jupyter_notebooks/working_script/gene_module/M_cluster/Mcluster_gene_lists_with_default_retrieve_markers_settings.csv", fileEncoding = "GBK")
data = data[colnames(data)!= "X"]
data = data %>% mutate(across(everything(), as.character))

genes = data %>% mutate(across(everything(), ~ str_extract(.x, "[^/]+")))
values = data %>% mutate(across(everything(), ~ str_remove(str_extract(.x, "(?=/).+"), "/")))

#  read numbers2ids file in
num2id = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/gene_name/gene_name_list.csv", header = FALSE)
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
ensemble_id_files = ensemble_id_files[, map_lgl(ensemble_id_files, ~ ! any(str_like(.x, "NA.+")))]

# run get_gene_module
gene_module = get_gene_module(ensemble_id_files, parallel_num = 32)
saveRDS(gene_module, "Mcluster_gene_lists_with_default_retrieve_markers_settings_gene_module_first_multi_99datasets.rds")

# read the gene module file
gene_module = readRDS("Mcluster_gene_lists_with_default_retrieve_markers_settings_gene_module_first_multi_99datasets.rds")

# calculate jaccard similarity matrix
data2 = ensemble_id_files %>% mutate(across(everything(), ~ str_remove(.x, "/.+")))

jaccard_matrix = map_dfr(
  seq(dim(data2)[2]),
  function(x) {
    base = data2[, x, drop= TRUE]
    base_name = colnames(data2)[x]
    map_dfc(colnames(data2) ,
            function(y) {
              compare = data2[, y, drop = TRUE]
              compare_name = y
              jaccard = jaccard_score(base,compare) #length(dplyr::intersect(base, compare))
              df = data.frame(value = jaccard)
              colnames(df) = compare_name
              return(df)
            }
          )
  }
)

rownames(jaccard_matrix) = colnames(jaccard_matrix)

# output the gene module lists
gene_module_lists = gene_module[[1]][map_lgl(gene_module[[2]], ~ length(.x) > 2)]
gene_module_lists = map(gene_module_lists, ~ str_remove(.x, "/.+"))
names(gene_module_lists) = paste0("gene_module_", seq(length(gene_module_lists)))
gene_module_lists = as.data.frame(gene_module_lists)
write.csv(gene_module_lists, "M_cluster_with_default_retrieve_markers_settings_gene_module_list.csv", fileEncoding = "utf-8")

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
immune_response_module = c(3, 5, 6, 7, 11, 13, 14, 15)
cell_development_module = c(8, 12)
angiogenesis_module = c(1, 9)
cytoskeleton_module = c(10)
defense_response_module = c(2)
protein_target_module = c(4)
# cell_adhesion_module = c(15)

modules = list(immune_response_module, cell_development_module,
               angiogenesis_module, cytoskeleton_module,
               defense_response_module, protein_target_module)

names(modules) = c("immune_response", "cell_cycle", "angiogenesis", "cytoskeleton",
                   "defense_response", "protein_target")

# only focus on gene modules with above 2 members
robust_candidates = gene_module[[2]][map_lgl(gene_module[[2]], ~ length(.x) > 2)] %>% unlist()
jaccard_matrix_sub = jaccard_matrix[, robust_candidates]
cond = map_int(robust_candidates, ~ which(rownames(jaccard_matrix_sub) %in% .x))
jaccard_matrix_sub = jaccard_matrix_sub[cond, robust_candidates]

# reorder gene module using heirarchy clustering
intersection = find_overlap_singlet(gene_module_lists, full_search = T)
intersection = intersection %>%
  pivot_wider(everything(), names_from = "compare_name", values_from = "number") %>%
  column_to_rownames("base_name")

between_module_reorder = map(
  modules,
  function(x) {
    if (length(x) > 1) {
      sub_data = intersection[x, x]
      model = hclust(dist(sub_data), method = "average")
      as.integer(str_extract(colnames(sub_data)[rev(model$order)], "\\d+"))
    }
    else{
      x
    }
  }
)

between_module_reorder[["immune_response"]] = c(15, 13, 7, 3, 14, 5, 6, 11)

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

# heatmap showing the jaccard similarity between members in each gene module
colors = rev(colorRampPalette(brewer.pal(11, "Spectral"))(100))
# custom_magma <- c(colorRampPalette(c("white", rev(magma(200, begin = 0.15))[1]))(10), rev(magma(200, begin = 0.18)))

png("Mcluster_gene_module_heatmap.png", width = 2000, height = 1500, res = 300)

pheatmap::pheatmap(jaccard_matrix_sub_reorder,
                   annotation_col = annotation_col,
                   border_color = NA,
                   show_colnames = F, show_rownames = F, color = colors, cluster_rows = F, cluster_cols = F)
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
  annotation_col$func == "cell_cycle" ~ colors2[3],
  annotation_col$func == "angiogenesis" ~ colors2[4],
  annotation_col$func == "cytoskeleton" ~ colors2[5],
  annotation_col$func == "protein_target" ~ colors2[6],
  annotation_col$func == "defense_response" ~ colors2[7],
  .default =  "black"))
)

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
    module %in% cell_development_module ~ "cell_development",
    module %in% angiogenesis_module ~ "angiogenesis",
    module %in% cytoskeleton_module ~ "cytoskeleton",
    module %in% defense_response_module ~ "defense_response",
    module %in% protein_target_module ~ "protein_target",
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


new_colors = c("#E64B35FF", "#3C5488FF", "#F39B7FFF", "#00A087FF", "#91D1C2FF", "#8491B4FF")
names(new_colors) = c("immune_response", "angiogenesis", "cytoskeleton", "cell_development", "defense_response", "protein_target")

pie <- ggplot(module_ratio, aes(x="", y=percentage, fill= module_type)) +
  geom_bar(width = 1, stat = "identity", position = "stack") +
  # scale_fill_npg() +
  scale_fill_manual(values = new_colors) +
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank(), legend.position = "none")

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
