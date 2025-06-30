
library(CellChat)
library(patchwork)
library(Seurat)
library(data.table)
library(SeuratDisk)
library(rhdf5)

# load expression matrix and meta data
file_path <- "/mnt/disk1/zhongmin/superscc/D035_hallmarks_tumors_copy/Data_Neuroendocrine/Data_Dong2020_Neuroendocrine_1/未去批次效应counts数据/没有去除批次效应_Data_Dong2020_Neuroendocrine_1数据.csv"
data <- fread(file_path, data.table = FALSE)
rownames(data) <- data[,1]
data <- data[,-1]
head(data[1:10, 1:10])
data = as.data.frame(t(data))
data <- data[rownames(data) != "cell_type", ]
anyNA(data)
head(data[1:10, 1:10])

# load annotation file  
file_path <- "/mnt/disk1/zhongmin/superscc/D035_hallmarks_tumors_copy/Data_Neuroendocrine/Data_Dong2020_Neuroendocrine_1/未去批次效应counts数据/内皮细胞cellchat_obs.csv"
metadata <- fread(file_path, data.table = FALSE)
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]
head(metadata)

# do normalization and remove cells with unknown annotations
seurat <- CreateSeuratObject(counts = data,  meta.data = metadata)
zzm1 = seurat@meta.data
anyNA(seurat@assays$RNA$counts)
selected_cell_types <- c("Double")
seurat <- subset(seurat, subset = celltype_zzm %in% selected_cell_types,invert = T)
table(seurat$celltype_zzm)
seurat <- NormalizeData(seurat)
cellchat.matrix = seurat@assays$RNA$data
meta <- seurat@meta.data


# create cell chat object
cellchat=createCellChat(object = cellchat.matrix,meta = meta,group.by = "celltype_zzm")
cellchat <- setIdent(cellchat, ident.use = "celltype_zzm") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 1)
names(cellchat)

# calculate the cell-cell communication probabailites 
df.net <- subsetCommunication(cellchat)
head(df.net)
class(df.net)
cellchat <- computeCommunProbPathway(cellchat)
cellchat@netP$pathways

# aggregate communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat,"/mnt/disk5/zhongmin/superscc/结果位置/cellcahat/免疫内皮细胞.rds")
cellchat = readRDS("/mnt/disk5/zhongmin/superscc/结果位置/cellcahat/免疫内皮细胞.rds")

# visualization of cell-cell communications only focus on MHC-II signalling 
netVisual_bubble(
  object = cellchat, 
  sources.use = c("Immune_endothelial cells", "Non immune_endothelial cells"), 
  targets.use = c("B cells", "Macrophage"),
  signaling = c("MHC-II"),  
  remove.isolate = FALSE,
  thresh = 0.05
)

pathways.show <- c("MHC-II") 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP", thresh = 0.05)
netAnalysis_signalingRole_network(cellchat, signaling = "MHC-II", width = 8, height = 2.5, font.size = 10)

dev.off()


