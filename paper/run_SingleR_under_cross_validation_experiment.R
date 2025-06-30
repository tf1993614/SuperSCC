library(SingleR)
library(tidyverse)
library(data.table)

parent_dir = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleR"
setwd(parent_dir)

run_singleR = function(count, ref_data, test_data, filename){
  
  if(dir.exists(filename) == FALSE) {
    dir.create(filename)
  }
  
  setwd(filename)
  
  # read cell index and expression matrix in
  ref_cell_index = read.csv(ref_data)

  test_cell_index = read.csv(test_data)

  exp_mat = fread(count)

  ref = exp_mat[exp_mat$V1 %in% ref_cell_index$X, ]
  ref = as.data.frame(ref) %>% column_to_rownames("V1") %>% t() %>% as.data.frame() %>% dplyr::select(ref_cell_index$X)
  ref = SingleCellExperiment::SingleCellExperiment(list(counts = ref))

  test = exp_mat[exp_mat$V1 %in% test_cell_index$X, ]
  test = as.data.frame(test) %>% column_to_rownames("V1") %>% t() %>% as.data.frame() %>% dplyr::select(test_cell_index$X)
  test = SingleCellExperiment::SingleCellExperiment(list(counts = test))

  # do normalization
  ref = scuttle::logNormCounts(ref)
  test = scuttle::logNormCounts(test)

  # do the annotation
  pred = SingleR(test=test, ref= ref, labels= ref_cell_index$cell_type, de.method="wilcox")
  saveRDS(pred, paste0(filename, "_SingleR_prediction.rds"))

  # write the prediction
  pred_table = tibble(y_true = test_cell_index$cell_type, y_pred = pred$labels)
  write_csv(pred_table, paste0(filename, "_SingleR_prediction.csv"))

  setwd(parent_dir)
  
}

# get the file loc
file = read.csv("/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleR/file_loc.csv")
file2 = list.files(pattern = "(.+)_test_dataset.+", path = "/mnt/disk5/zhongmin/superscc/label_transfer/代码/SuperSCC/finest_cell_label_res", recursive = T, full.names = T)
file3 = list.files(pattern = "(.+)_train_dataset.+", path = "/mnt/disk5/zhongmin/superscc/label_transfer/代码/SuperSCC/finest_cell_label_res", recursive = T, full.names = T)
file$ref_data = file3
file$test_data = file2
file = file[c(2:5)]
file = file[c(1, 3, 4, 2)]

# recursive run 
pmap(file, run_singleR)

