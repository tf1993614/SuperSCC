library(scmap)
library(SingleCellExperiment)
library(tidyverse)
library(snow)
library(data.table)

setwd("/home/fengtang/jupyter_notebooks/working_script/label_transfer/scmap")

# @ counts matrix should be features x cells 
# @ annotation should be a dataframe with cell type column; row index should the colnames of counts_matrix

run_scamp = function(count, ref_data, test_data, filename) {
  
  ref_cell_index = read.csv(ref_data)
  
  test_cell_index = read.csv(test_data)
  
  exp_mat = fread(count)
  
  ref = exp_mat[exp_mat$V1 %in% ref_cell_index$X, ]
  ref = as.data.frame(ref) %>% column_to_rownames("V1") %>% t() %>% as.data.frame() %>% dplyr::select(ref_cell_index$X)
  
  ref_sce = SingleCellExperiment(assays = list(normcounts = as.matrix(ref)), colData = ref_cell_index[, "cell_type", drop = FALSE])
  logcounts(ref_sce) = log2(normcounts(ref_sce) + 1)
  rowData(ref_sce)$feature_symbol = rownames(ref_sce) # use gene names as feature symbols
  ref_sce = ref_sce[!duplicated(rownames(ref_sce)), ] # remove features with duplicated names
  ref_sce = selectFeatures(ref_sce, suppress_plot = TRUE)
  ref_sce = indexCluster(ref_sce, cluster_col = "cell_type")
  
  test = exp_mat[exp_mat$V1 %in% test_cell_index$X, ]
  test = as.data.frame(test) %>% column_to_rownames("V1") %>% t()  %>% as.data.frame() %>% dplyr::select(test_cell_index$X)
  
  query_sce = SingleCellExperiment(assays = list(normcounts = as.matrix(test)))
  logcounts(query_sce) = log2(normcounts(query_sce) + 1)
  rowData(query_sce)$feature_symbol = rownames(query_sce) # use gene names as feature symbols
  query_sce = query_sce[!duplicated(rownames(query_sce)), ] # remove features with duplicated names
  query_sce = selectFeatures(query_sce, suppress_plot = TRUE)
  
  scmapCluster_results = scmapCluster(
    projection = query_sce, 
    index_list = list(
      prediction = metadata(ref_sce)$scmap_cluster_index
    )
  )
  
  saveRDS(scmapCluster_results, paste0(filename, ".rds"))
  print(paste0("finish prediction on ", filename))
}



file = read.csv("/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleCellNet/label_transfer_evulate_data_loc.csv", row.names = 1)
file$ref_data = str_replace(file$ref_data, "代码/", "代码/SuperSCC/finest_cell_label_res/")
file$test_data = str_replace(file$test_data, "代码/", "代码/SuperSCC/finest_cell_label_res/")

# run
pmap(file, run_scamp)

# tidy up result
files = list.files(path = ".", pattern = ".+rds$", recursive = F, full.names = T)
file_name = basename(files) %>% str_remove(".rds")
scmap_pred = map(files, readRDS) 
walk2(
  scmap_pred,
  file_name,
  function(x, y) {
    data = x$scmap_cluster_labs %>% as.data.frame() %>% dplyr::select(y_pred = prediction)
    write.csv(data, paste0(y, "_scmap_prediction.csv"))
  }
)
