library(singleCellNet)
library(tidyverse)
library(data.table)

parent_dir = "/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleCellNet"
setwd(parent_dir)

run_SingleCellNet = function(count, ref_data, test_data, filename){

  if(dir.exists(filename) == FALSE) {
    dir.create(filename)
  }

  setwd(filename)

  # read cell index and expression matrix in
  ref_cell_index = read.csv(ref_data)

  test_cell_index = read.csv(test_data)

  exp_mat = fread(count)

  ref = exp_mat[exp_mat$V1 %in% ref_cell_index$X, ]
  ref = as.data.frame(ref) %>% column_to_rownames("V1") %>% t() %>% as.data.frame() %>% dplyr::select(ref_cell_index$X) %>% as.matrix()

  test = exp_mat[exp_mat$V1 %in% test_cell_index$X, ]
  test = as.data.frame(test) %>% column_to_rownames("V1") %>% t() %>% as.data.frame() %>% dplyr::select(test_cell_index$X) %>% as.matrix()

  names(ref_cell_index) = c("x", "y")

  class_info = scn_train(
                        stTrain = as.data.frame(ref_cell_index),
                        expTrain = ref,
                        nTopGenes = 10,
                        nRand = 70,
                        nTrees = 1000,
                        nTopGenePairs = 25,
                        dLevel = "y",
                        colName_samp = "x"
                        )

  saveRDS(class_info, paste0(filename, "_SingleCellNet_Prediction.rds"))


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


  write_csv(data.frame(y_true = test_cell_index$cell_type, y_pred = pred_cell_type$cell_type), paste0(filename, "_SingleCellNet_Prediction_cell_type_", date(), ".csv"))

  setwd(parent_dir)

}


file = read.csv("/home/fengtang/jupyter_notebooks/working_script/label_transfer/SingleCellNet/label_transfer_evulate_data_loc.csv", row.names = 1)

# recursive run
pmap(file, run_SingleCellNet)
