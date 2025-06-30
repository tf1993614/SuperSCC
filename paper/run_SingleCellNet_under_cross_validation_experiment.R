library(singleCellNet)
library(tidyverse)
library(data.table)

parent_dir = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleCellNet"
setwd(parent_dir)

run_SingleCellNet = function(ref_data, test_data, filename){
  
  if(dir.exists(filename) == FALSE) {
    dir.create(filename)
  }
  
  setwd(filename)
  
  
  ref = fread(ref_data, data.table = FALSE)
  ref_label = ref %>% as.data.frame() %>% dplyr::select(V1, cell_type) 
  names(ref_label) = c("x", "y")
  
  ref = as.data.frame(ref) %>% column_to_rownames("V1") %>% dplyr::select(where(is.numeric)) %>% t() %>% as.matrix() 
  test = fread(test_data, data.table = FALSE)
  test = as.data.frame(test) %>% column_to_rownames("V1") %>% dplyr::select(where(is.numeric)) %>% t() %>% as.matrix()
  
  if(len(list.files(pattern = ".+_SingleCellNet_Prediction.rds", path = parent_dir, recursive = TRUE)) == 0) {
    
    class_info = scn_train(
    stTrain = as.data.frame(ref_label),
    expTrain = ref,
    nTopGenes = 10,
    nRand = 70,
    nTrees = 1000,
    nTopGenePairs = 25,
    dLevel = "y",
    colName_samp = "x"
  )
   saveRDS(class_info, paste0(filename, "_SingleCellNet_Prediction.rds"))
  }
  else {
    class_info = readRDS(list.files(pattern = ".+_SingleCellNet_Prediction.rds", path = parent_dir, recursive = TRUE, full.names = TRUE)[1])
  }
  
  
  pred = scn_predict(cnProc=class_info[['cnProc']], expDat=test, nrand = 0)
  
  write.csv(t(as.data.frame(pred)), paste0(filename, "_SingleCellNet_Prediction_raw.csv"))
  
  # output label
  pred_cell_type = pred %>%
    as.data.frame() %>%
    map_dfr(function(column) {
      max_value = max(unlist(column))
      max_index = which(unlist(column) == max_value)
      cell_type = rownames(pred)[max_index[1]]
      tibble(cell_type = cell_type)
    })
  
  
  write_csv(data.frame(y_true = test_cell_index$cell_type, y_pred = pred_cell_type$cell_type), paste0(filename, "_SingleCellNet_Prediction_cell_type.csv"))
  
  setwd(parent_dir)
  
}


file = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/结果位置_3.csv")
file = file %>% subset(数据集 %in% c("Banovich_Kropski_2020", "Barbry_Leroy_2020", "Krasnow_2020", "Lafyatis_Rojas_2019", "Nawijn_2021", "Teichmann_Meyer_2019"))
file = file[c(1, 7)]
file$ref_data = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SuperSCC/train_datasets_from_different_study.csv"
names(file)[c(1,2)] = c("filename", "test_data")
file = file %>% dplyr::select(ref_data, test_data, filename)

# recursive run
pmap(file, run_SingleCellNet)

