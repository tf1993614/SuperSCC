library(tidyverse)
library(data.table)
library(Seurat)

setwd("/home/fengtang/jupyter_notebooks/working_script/label_transfer/Seurat")

run_seurat= function(count, ref_data, test_data, filename) {
  
  ref_cell_index = read.csv(ref_data)
  
  test_cell_index = read.csv(test_data)
  
  exp_mat = fread(count)
  
  ref = exp_mat[exp_mat$V1 %in% ref_cell_index$X, ]
  ref = as.data.frame(ref) %>% column_to_rownames("V1") %>% t() %>% as.data.frame() %>% dplyr::select(ref_cell_index$X)
  
  ref_obj = CreateAssayObject(counts = ref)
  ref_obj = CreateSeuratObject(ref_obj)
  ref_obj = NormalizeData(ref_obj)
  ref_obj = FindVariableFeatures(ref_obj)
  ref_obj = ScaleData(ref_obj)
  ref_obj = RunPCA(ref_obj)
  ref_obj$cell_type = ref_cell_index$cell_type 
  
  
  test = exp_mat[exp_mat$V1 %in% test_cell_index$X, ]
  test = as.data.frame(test) %>% column_to_rownames("V1") %>% t()  %>% as.data.frame() %>% dplyr::select(test_cell_index$X)
  
  query_obj = CreateAssayObject(counts = test)
  query_obj = CreateSeuratObject(query_obj)
  query_obj = NormalizeData(query_obj)
  
  anchors = FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:30, reference.reduction = "pca")
  
  predictions = TransferData(anchorset = anchors, refdata = ref_obj$cell_type, dims = 1:30)
  
  saveRDS(predictions, paste0(filename, "_seurat_prediction.rds"))
  print(paste0("finish prediction on ", filename))
}

file = read.csv("/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleCellNet/label_transfer_evulate_data_loc.csv", row.names = 1)
file$ref_data = str_replace(file$ref_data, "代码/", "代码/SuperSCC/finest_cell_label_res/")
file$test_data = str_replace(file$test_data, "代码/", "代码/SuperSCC/finest_cell_label_res/")

pmap(file, run_seurat)

# tidy up results
files = list.files(path = ".", pattern = ".+rds$", recursive = F, full.names = T)
file_name = basename(files) %>% str_remove("_seurat_prediction.rds")
seurat_pred = map(files, readRDS) 
walk2(
  seurat_pred,
  file_name,
  function(x, y) {
    data = x %>% dplyr::select(y_pred = predicted.id)
    write.csv(data, paste0(y, "_seurat_prediction.csv"))
  }
)

