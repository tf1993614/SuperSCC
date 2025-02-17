library(tidyverse)
library(snow)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(VennDiagram)
library(ggsci)
library(Seurat)
library(slider)
library(data.table)
library(geneModule)

setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/F_cluster/")

############################## correlation of gene modules on Fcluster ######################
#                                                                                           #
#                                 only focus on shared modules                              #
#                                                                                           #
#############################################################################################

# read files
health_disease_state_ratio = read.csv("F_cluster_health_disease_state_ratio.csv")
gene_module = readRDS("F_cluster_with_default_retrieve_markers_settings_gene_module_first_multi.rds")
disease_info = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/作者年_表达矩阵_metadata_pkl/作者年_表达矩阵_metadata_pkl.csv", fileEncoding = "GBK")
num2id = read.csv("gene_name_list.csv", header = FALSE)

# get the robust gene modules
# and the corresponding contributors
robust_candidates = gene_module[[2]][map_lgl(gene_module[[2]], ~ length(.x) > 2)]
names(robust_candidates) = paste0("GM_", seq(65))
robust_modules = gene_module[[1]][map_lgl(gene_module[[2]], ~ length(.x) > 2)]
robust_modules = map(robust_modules, ~ str_remove(.x, "/.+"))
names(robust_modules) = paste0("GM_", seq(65))

max_length = max(map_int(robust_candidates, ~ length(.x)))
new_robust_candidates = map(
  robust_candidates,
  function(x) {
    if(length(x) < max_length) {
      remain = max_length - length(x)
      c(x, rep(NA, remain))
    }
    else {
      x
    }
  }
)

write.csv(as.data.frame(new_robust_candidates), "F_modules_candidates.csv")

modules_by_ensembl = robust_modules

# convert ensembl id to gene symbol
modules_by_symbol = map(
  robust_modules,
  function(x) {
    first = function(y){y[[1]]}
    mapIds(org.Hs.eg.db, keys = x, keytype = "ENSEMBL", column = "SYMBOL", multiVals = first)
  }
)


# load the counts matrix to be scored
# do normalization and use the gene module
# score the cells
loop_range = slide(seq(99), ~ .x, .step = 20, .after = 19)
loop_range = loop_range[map_lgl(loop_range, ~ ! is.null(.x))]

cell_scoring = function(input, index) {

    message(paste0(date(), ": load data"))
    data = fread(input[["counts"]][index], data.table = F)

    if("V1" %in% colnames(data)) {
      data = data %>% column_to_rownames("V1")
    }
    else if("cell_name" %in% colnames(data)) {
      data = data %>% column_to_rownames("cell_name")
    }
    else {
      data = data %>% column_to_rownames("index")
    }

    # check duplicated names
    # if yes, only keep the first observation
    duplicated_names = colnames(data) %>% table()
    duplicated_names = duplicated_names[duplicated_names > 1]
    if(length(duplicated_names) != 0) {
      duplicated_names = names(duplicated_names)
      remove_duplicated_index = map(duplicated_names, ~ which(colnames(data) %in% .x)[-1]) %>% unlist() %>% sort()
      data = data %>% as.data.frame() %>% .[, - remove_duplicated_index]
    }

    data = data[, ! colnames(data) %in% c("cluster", "cell_type", "Cluster")]
    data = data %>% dplyr::select(where(is.numeric))

    data = t(data)
    message(paste0(date(), ": finish loading data"))

    message(paste0(date(), ": check feature type"))
    is_ensembl = all(str_detect(rownames(data), "^ENSG.+"))
    is_entrez = all(str_detect(rownames(data), "^V\\d+"))
    message(paste0(date(), ": finish checking feature type"))
    if(is_ensembl) {
      message(paste0(date(), ": feature type is ENSEMBL ID"))
    }
    else if(is_entrez) {
      message(paste0(date(), ": feature type is ENTREZ ID"))
      df = data.frame(row.names = rownames(data))
      df = df %>% left_join(num2id, by = join_by(row.names == V1))
      rownames(data) = df$V1
    }
    else {
      message(paste0(date(), ": feature type is gene symbol"))
    }

    if(is_entrez) {
      df = data.frame(V1 = rownames(data))
      df = df %>% left_join(num2id, by = V1)
      rownames(data) = df$V2
    }

    message(paste0(date(), ": check normalization"))
    is_normalization = any(map_lgl(data[data > 0], ~ .x %% 1 != 0))
    message(paste0(date(), ": finish checking normalization - normalization is ", is_normalization))

    message(paste0(date(), ": start contruction of Seurat object"))
    if(is_normalization) {
      obj = CreateAssayObject(data = data)
      obj = CreateSeuratObject(obj)
    }
    else {
      obj = CreateAssayObject(counts = data)
      obj = CreateSeuratObject(obj)
      obj = NormalizeData(obj)
    }
    message(paste0(date(), ": finish contruction of Seurat object"))

    message(paste0(date(), ": start cell scoring"))
    res = map(
      seq(length(modules_by_ensembl)),
      # seq(length(shared_modules_by_ensembl)),
      function(idx) {
        if(is_entrez) {
          # obj = AddModuleScore(obj, features = list(c(shared_modules_by_symbol[[idx]])), name = paste0(names(shared_modules_by_symbol)[idx], "/"))
          # remained_features = rownames(obj)[rownames(obj) %in% shared_modules_by_symbol[[idx]]]
          obj = AddModuleScore(obj, features = list(c(modules_by_symbol[[idx]])), name = paste0(names(modules_by_symbol)[idx], "/"))
          remained_features = rownames(obj)[rownames(obj) %in% modules_by_symbol[[idx]]]
        }
        else if(is_ensembl) {
          # obj = AddModuleScore(obj, features = list(c(shared_modules_by_ensembl[[idx]])), name = paste0(names(shared_modules_by_symbol)[idx], "/"))
          # remained_features = rownames(obj)[rownames(obj) %in% shared_modules_by_ensembl[[idx]]]
          obj = AddModuleScore(obj, features = list(c(modules_by_ensembl[[idx]])), name = paste0(names(modules_by_ensembl)[idx], "/"))
          remained_features = rownames(obj)[rownames(obj) %in% modules_by_ensembl[[idx]]]
        }
        else {
          # obj = AddModuleScore(obj, features = list(c(shared_modules_by_symbol[[idx]])), name = paste0(names(shared_modules_by_symbol)[idx], "/"))
          # remained_features = rownames(obj)[rownames(obj) %in% shared_modules_by_symbol[[idx]]]
          obj = AddModuleScore(obj, features = list(c(modules_by_symbol[[idx]])), name = paste0(names(modules_by_symbol)[idx], "/"))
          remained_features = rownames(obj)[rownames(obj) %in% modules_by_symbol[[idx]]]
        }
        return(list(score = obj@meta.data, score_features = remained_features, normalization = is_normalization, dataset = input[["new"]][index]))
      },
      .progress = T
    )
    message(paste0(date(), ": finish cell scoring"))

    return(res)
}

cell_scoring_parallel = function(data, idx, parallel_num = 10) {
  cl = snow::makeCluster(get("parallel_num"))
  snow::clusterEvalQ(cl, {library(tidyverse)})
  snow::clusterEvalQ(cl, {library(Seurat)})
  snow::clusterEvalQ(cl, {library(data.table)})
  snow::clusterExport(cl, "num2id")
  snow::clusterExport(cl, "disease_info")
  # snow::clusterExport(cl, "shared_modules_by_ensembl")
  # snow::clusterExport(cl, "shared_modules_by_symbol")
  snow::clusterExport(cl, "modules_by_ensembl")
  snow::clusterExport(cl, "modules_by_symbol")
  snow::clusterExport(cl, "cell_scoring")
  snow::clusterExport(cl, "data")

  res = clusterApply(cl, idx, function(index) cell_scoring(data, index))
  stopCluster(cl)

  return(res)
}

cell_scores_1 = cell_scoring_parallel(disease_info, loop_range[[1]], 50)
cell_scores_2 = cell_scoring_parallel(disease_info, loop_range[[2]], 50)
cell_scores_3 = cell_scoring_parallel(disease_info, loop_range[[3]], 50)
cell_scores_4 = cell_scoring_parallel(disease_info, loop_range[[4]], 50)
cell_scores_5 = cell_scoring_parallel(disease_info, loop_range[[5]], 50)

all_cell_scores = c(cell_scores_1, cell_scores_2, cell_scores_3, cell_scores_4,
                    cell_scores_5)

names(all_cell_scores) = disease_info$new

saveRDS(all_cell_scores, "all_gene_modules_cell_scores.rds")
