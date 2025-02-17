library(tidyverse)
library(data.table)
library(Seurat)
library(snow)

setwd("/home/fengtang/jupyter_notebooks/working_script/gene_module/M_cluster/")

all_cell_scores = readRDS("/home/fengtang/jupyter_notebooks/working_script/gene_module/M_cluster/M_cluster_all_gene_modules_cell_scores.rds")
disease_info = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/作者年_表达矩阵_metadata_pkl/作者年_表达矩阵_metadata_pkl.csv", fileEncoding = "GBK")
disease_info = disease_info %>% subset(new != "Moreno_2022")

dataset_names = map_chr(
  all_cell_scores,
  function(x) {
    unique(map_chr(
      x,
      function(z) {
        z[[4]]
      }
    ))
  }
)

names(all_cell_scores) = dataset_names

all_cell_scores = map(
  all_cell_scores,
  function(x) {
    names(x) = paste0("GM_", seq(15))
    x
  }
)


# Min Max scaling function
min_max_scale = function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}


# basic cell scoring plot function
cell_scoring_plot = function(input, input2, index, index2) {

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


  meta = read.csv(input[["metadata"]][index])
  obj@meta.data = bind_cols(obj@meta.data, meta)
  cell_type = colnames(meta)[colnames(meta) %in%
                               c("Manuscript_Identity", "author_cell_type",
                                 "ann_finest_level", "cell_type",
                                 "celltype",
                                 "broad_cell_type", "annotated_type", "broad_celltype",
                                 "detail_celltype")]
  cell_type = cell_type[1]
  name = input[["new"]][index]

  scores = input2[input[["new"]][index]]
  obj[[index2]] = min_max_scale(scores[[1]][[index2]][[1]][[paste0(index2, "/1")]])

  obj = obj %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()  %>%
    RunUMAP(dims = 1:30)

  Idents(obj) = obj@meta.data[[cell_type]]
  fig1 = DimPlot(obj, pt.size = 1.5) + ggtitle(paste0(name, "_", index2))
  fig2 = FeaturePlot(obj, features = index2,
                     pt.size = 1.5,
                     cols = rev(colorRampPalette(brewer.pal(11, "Spectral"))(100))) + ggtitle(paste0(name, "_", index2))
  fig3 = DimPlot(obj, pt.size = 1.5, label = T) + ggtitle(paste0(name, "_", index2))

  ggsave(plot = fig1, paste0(name, "_", index2, "_A_.png"), width = 8000, height = 4000, units = "px")
  ggsave(plot = fig2, paste0(name, "_", index2, "_B_.png"), width = 8000, height = 4000, units = "px")
  ggsave(plot = fig3, paste0(name, "_", index2, "_C_.png"), width = 8000, height = 4000, units = "px")
}

# cell scoring plot supported by parallely running
cell_scoring_plot_parallel = function(input, input2, index, index2, parallel_num = 10) {

  cl = snow::makeCluster(get("parallel_num"))
  snow::clusterEvalQ(cl, {library(tidyverse)})
  snow::clusterEvalQ(cl, {library(Seurat)})
  snow::clusterEvalQ(cl, {library(data.table)})
  snow::clusterEvalQ(cl, {library(grDevices)})
  snow::clusterEvalQ(cl, {library(RColorBrewer)})
  snow::clusterExport(cl, "disease_info")
  snow::clusterExport(cl, "all_cell_scores")
  snow::clusterExport(cl, "min_max_scale")
  snow::clusterExport(cl, "cell_scoring_plot")
  snow::clusterEvalQ(cl, "input")
  snow::clusterEvalQ(cl, "input2")
  snow::clusterEvalQ(cl, "index2")

  res = snow::clusterApply(cl, index, function(i) cell_scoring_plot(input, input2, i, index2))
  stopCluster(cl)

  return(res)
}


# scoring the angiogenesis and immune response modules
home = getwd()
map(
  paste0("GM_", c(1, 9, 3, 5, 6, 7, 11, 13, 14, 15)),
  function(x) {
    dir.create(x)
    setwd(x)
    cell_scoring_plot_parallel(input = disease_info,
                               input2 = all_cell_scores,
                               index = 1:98,
                               index2 = x,
                               parallel_num = 30)
    setwd(home)
  }
)
