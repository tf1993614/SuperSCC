library(SingleR)
library(tidyverse)
library(data.table)

parent_dir = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleR"
setwd(parent_dir)

run_singleR = function(ref_data, test_data, filename){
  
  if(dir.exists(filename) == FALSE) {
    dir.create(filename)
  }
  
  setwd(filename)
  
  ref = fread(ref_data, data.table = FALSE)
  ref_label = ref %>% as.data.frame() %>% dplyr::select(V1, cell_type) 
  names(ref_label) = c("x", "y")
  
  ref = as.data.frame(ref) %>% column_to_rownames("V1") %>% dplyr::select(where(is.numeric)) %>% t() %>% as.matrix() 
  ref = SingleCellExperiment::SingleCellExperiment(list(counts = ref))
  
  test = fread(test_data, data.table = FALSE)
  test_label = test$cell_type
  test = as.data.frame(test) %>% column_to_rownames("V1") %>% dplyr::select(where(is.numeric)) %>% t() %>% as.matrix()
  test = SingleCellExperiment::SingleCellExperiment(list(counts = test))
  
  # do normalization
  ref = scuttle::logNormCounts(ref)
  test = scuttle::logNormCounts(test)
  
  # do the annotation
  pred = SingleR(test=test, ref= ref, labels= ref_label$y, de.method="wilcox")
  saveRDS(pred, paste0(filename, "_SingleR_prediction.rds"))
  
  # write the prediction
  pred_table = tibble(y_true = test_label$cell_type, y_pred = pred$labels)
  write_csv(pred_table, paste0(filename, "_SingleR_prediction.csv"))
  
  setwd(parent_dir)
  
}

# get the file loc
file = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/结果位置_3.csv", fileEncoding = "GBK")
file = file %>% subset(数据集 %in% c("Banovich_Kropski_2020", "Barbry_Leroy_2020", "Krasnow_2020", "Lafyatis_Rojas_2019", "Nawijn_2021", "Teichmann_Meyer_2019"))
file = file[c(1, 11)]
file$ref_data = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SuperSCC/train_datasets_from_different_study.csv"
names(file)[c(1, 2, 3)] = c("filename", "test_data", "ref_data")
file = file %>% dplyr::select(ref_data, test_data, filename)

# recursive run 
pmap(file, run_singleR)

