
library(SC3)
library(foreach)
library(rngtools)
library(Seurat)

library(SC3)
library(tidyverse)
library(data.table)
library(SingleCellExperiment)
# SC3 manual
# https://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/SC3.html

# import expression matrix
file_path <- "/mnt/disk5/zhongmin/superscc/师兄整理的肺数据/未去批次效应couns数据/没有去除批次效应_Krasnow_2020数据.csv"
data <- fread(file_path)
data = as.data.frame(data)
data = data %>% column_to_rownames("V1")
data = t(data)

# import meta data
file_path <- "/mnt/disk5/zhongmin/superscc/师兄整理的肺数据/未去批次效应数据/没有去除批次效应_Krasnow_2020数据metadata.csv"
metadata <- fread(file_path)

# organize meta data
cell_type = data.frame(cell_type = metadata$ann_finest_level)
cell_type$cell = metadata$V1
cell_type2 = data.frame(cell = colnames(data))
cell_type3 = cell_type2 %>% left_join(cell_type, by = "cell")

obj <- CreateSeuratObject(counts = data,  meta.data = metadata)
obj <- NormalizeData(obj)
data2 = as.matrix(obj@assays$RNA$data)

# create a SingleCellExperiment object
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(data2),
    logcounts = as.matrix(data2)
  ),
  colData = cell_type3$cell_type
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)

# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# rub sc3
sce <- sc3(sce, ks = c(6), biology = TRUE) # the Ks is determined by the number of SuperSCC M clustering.
sce <- sc3_run_svm(sce, ks = c(6))

# get the SC3 clustering
col_data <- colData(sce)
col_data_df <- as.data.frame(col_data[, grep("sc3_", colnames(col_data))])

write.csv(col_data_df, "col_data.csv", row.names = TRUE)

saveRDS(sce, "Krasnow_2020_SC3.rds")
