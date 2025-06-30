library(tidyverse)
library(data.table)
library(Seurat)

setwd("/home/fengtang/jupyter_notebooks/working_script/label_transfer/Seurat")

run_seurat= function(ref_data, test_data, filename) {
  
  exp_mat = fread(ref_data, data.table = F)
  
  ref = exp_mat %>% 
  column_to_rownames("V1") %>% 
  dplyr::select(-cell_type) %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::select(where(is.numeric))
  
  ref_obj = CreateAssayObject(counts = ref)
  ref_obj = CreateSeuratObject(ref_obj)
  ref_obj = NormalizeData(ref_obj)
  ref_obj = FindVariableFeatures(ref_obj)
  ref_obj = ScaleData(ref_obj)
  ref_obj = RunPCA(ref_obj)
  ref_obj$cell_type = exp_mat$cell_type 
  
  
  test = fread(test_data, data.table = F)
  test = test %>% 
    column_to_rownames("V1") %>% 
    dplyr::select(where(is.numeric)) %>% 
    t() %>% 
    as.data.frame() 
  
  query_obj = CreateAssayObject(counts = test)
  query_obj = CreateSeuratObject(query_obj)
  query_obj = NormalizeData(query_obj)
  
  anchors = FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:30, reference.reduction = "pca")
  
  predictions = TransferData(anchorset = anchors, refdata = ref_obj$cell_type, dims = 1:30)
  
  saveRDS(predictions, paste0(filename, "_seurat_prediction_on_multi_datasets.rds"))
  print(paste0("finish prediction on ", filename))
}

# get the file loc
file = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/结果位置_3.csv", fileEncoding = "GBK")
file = file %>% subset(数据集 %in% c("Banovich_Kropski_2020", "Barbry_Leroy_2020", "Krasnow_2020", "Lafyatis_Rojas_2019", "Nawijn_2021", "Teichmann_Meyer_2019"))
file = file[c(1, 11)]
file$ref_data = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SuperSCC/train_datasets_from_different_study.csv"
names(file)[c(1, 2, 3)] = c("filename", "test_data", "ref_data")
file = file %>% dplyr::select(ref_data, test_data, filename)

# recursive run
pmap(file, run_seurat)


# tidy up prediction results
file = list.files(path = ".", pattern = ".+on_multi_datasets.rds$", full.names = T)
file_name = str_remove(file, "./") %>% str_remove("rds")

predictions = map(file, readRDS)
names(predictions) = file_name

imap(predictions, ~ .x %>% dplyr::select(predicted.id) %>% write.csv(file = paste0(.y, "csv")))

