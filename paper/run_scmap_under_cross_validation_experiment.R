library(scmap)
library(SingleCellExperiment)
library(tidyverse)
library(snow)
library(data.table)

setwd("/home/fengtang/jupyter_notebooks/working_script/label_transfer/scmap")

# @ counts matrix should be features x cells 
# @ annotation should be a dataframe with cell type column; row index should the colnames of counts_matrix

run_scamp = function(ref_data, test_data, filename) {
  

  exp_mat = fread(ref_data)
  
  ref = exp_mat %>% 
    column_to_rownames("V1") %>% 
    dplyr::select(-cell_type) %>% 
    t() %>% 
    as.data.frame() %>% 
    dplyr::select(where(is.numeric))
  
  ref_sce = SingleCellExperiment(assays = list(normcounts = as.matrix(ref)), colData = exp_mat[, "cell_type", drop = FALSE])
  logcounts(ref_sce) = log2(normcounts(ref_sce) + 1)
  rowData(ref_sce)$feature_symbol = rownames(ref_sce) # use gene names as feature symbols
  ref_sce = ref_sce[!duplicated(rownames(ref_sce)), ] # remove features with duplicated names
  ref_sce = selectFeatures(ref_sce, suppress_plot = TRUE)
  ref_sce = indexCluster(ref_sce, cluster_col = "cell_type")
  
  test = fread(test_data, data.table = F)
  test = test %>% 
    column_to_rownames("V1") %>% 
    dplyr::select(is.numeric) %>% 
    t() %>% 
    as.data.frame() 
  
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
  
  saveRDS(scmapCluster_results, paste0(filename, "scmap_prediction_on_multi_datasets.rds"))
  print(paste0("finish prediction on ", filename))
}


# get the file loc
file = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/结果位置_3.csv", fileEncoding = "GBK")
file = file %>% subset(数据集 %in% c("Banovich_Kropski_2020", "Barbry_Leroy_2020", "Krasnow_2020", "Lafyatis_Rojas_2019", "Nawijn_2021", "Teichmann_Meyer_2019"))
file = file[c(1, 11)]
file$ref_data = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SuperSCC/train_datasets_from_different_study.csv"
names(file)[c(1, 2, 3)] = c("filename", "test_data", "ref_data")
file = file %>% dplyr::select(ref_data, test_data, filename)

pmap(file, run_scamp)


# tidy up prediction results
file = list.files(path = ".", pattern = ".+on_multi_datasets.rds$", full.names = T)
file_name = str_remove(file, "./") %>% str_remove("rds")

predictions = map(file, readRDS)
names(predictions) = file_name

imap(predictions, ~ .x$scmap_cluster_labs %>% as.data.frame() %>% write.csv(file = paste0(.y, "csv")))
