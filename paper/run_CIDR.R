library(cidr)
library(tidyverse)
library(data.table)
library(Matrix)

# read expression matrix
data = readRDS("/mnt/disk5/zhongmin/superscc/师兄整理的肺数据/seurat结果/Krasnow_2020/2024-09-02/task1/Krasnow_2020数据_seurat.rds")
counts = as.matrix(data@assays$RNA$counts)

# run CIDR
scCIDR = scDataConstructor(as.matrix(counts))
scCIDR = determineDropoutCandidates(scCIDR)
scCIDR = wThreshold(scCIDR)
scCIDR = scDissim(scCIDR)
scCIDR = scPCA(scCIDR)
scCIDR = nPC(scCIDR)
nCluster(scCIDR)
scCIDR <- scCluster(scCIDR)

# write CIDR results
saveRDS(scCIDR, paste0("scCIDR-", today(), ".rds"))
