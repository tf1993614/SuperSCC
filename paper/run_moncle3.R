# load packages
library(Seurat)
library(monocle3)

# create output path
input_dir <- "/data/beifen/zhongmin/superscc/18个数据集的rds文件/seurat_objects/"
output_dir_rds <- "/data/beifen/zhongmin/superscc/monocle3_results/rds"
output_dir_csv <- "/data/beifen/zhongmin/superscc/monocle3_results/csv"

dir.create(output_dir_rds, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_csv, recursive = TRUE, showWarnings = FALSE)

# get all rds files containing 18 benchmark datasets
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
file_names <- list.files(input_dir, pattern = "\\.rds$", full.names = FALSE)

# do monocle3 clustering
for (i in 1:length(rds_files)) {


  file_base <- sub("\\.rds$", "", file_names[i])

  cat("\n处理文件:", file_names[i], "\n")


  pbmc <- readRDS(rds_files[i])

  data = as.matrix(pbmc@assays$RNA$counts)
  cell_metadata <- pbmc@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
  cds <- preprocess_cds(cds, num_dim = 100)

  cds <- reduce_dimension(cds)

  plot_cells(cds, color_cells_by="scLCA_celltype")
  cds <- cluster_cells(cds, resolution=1e-5)

  save(cds, file = file.path(output_dir_rds, paste0(file_base, "_monocle3.RData")))

  result_df <- data.frame(
    cell_name = colnames(cds),
    scLCA_celltype = colData(cds)$scLCA_celltype,
    monocle3_cluster = clusters(cds)[colnames(cds)]
  )

  write.csv(result_df,
            file = file.path(output_dir_csv, paste0(file_base, "_clusters.csv")),
            row.names = FALSE)

  cat("完成文件:", file_names[i], "\n")
}
