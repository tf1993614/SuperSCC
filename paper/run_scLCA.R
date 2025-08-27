library(Seurat)
library(scLCA)
library(Matrix)

# Define input and output directories
input_dir <- "/mnt/disk5/zhongmin/superscc/results_location/seurat_objects"
output_dir <- "/mnt/disk5/zhongmin/superscc/results_location/scLCA_results"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get all RDS files
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "RDS files\n")

if (length(rds_files) == 0) {
  stop("No RDS files found in the specified directory")
}

# Create summary results table
summary_results <- data.frame(
  dataset_name = character(),
  file_name = character(),
  n_genes = integer(),
  n_cells = integer(),
  zero_percentage = numeric(),
  is_counts_matrix = logical(),
  n_true_clusters = integer(),
  n_estimated_clusters = integer(),
  ari_score = numeric(),
  status = character(),
  stringsAsFactors = FALSE
)

# Loop through each RDS file
for (i in 1:length(rds_files)) {
  rds_file <- rds_files[i]
  file_name <- basename(rds_file)
  dataset_name <- gsub("\\.rds$", "", file_name)  # Remove .rds suffix as dataset name

  cat("==========================================\n")
  cat("Processing file [", i, "/", length(rds_files), "]:", file_name, "\n")
  cat("Dataset name:", dataset_name, "\n")

  tryCatch({
    # Read Seurat object
    cat("  - Reading Seurat object...\n")
    seurat_obj <- readRDS(rds_file)

    # Check basic information of Seurat object
    cat("  - Seurat object information:\n")
    cat("    * Number of genes:", nrow(seurat_obj), "\n")
    cat("    * Number of cells:", ncol(seurat_obj), "\n")

    # Extract counts matrix
    datamatrix <- as.matrix(seurat_obj@assays$RNA$counts)

    # Extract true labels (directly use scLCA_celltype)
    if (!"scLCA_celltype" %in% colnames(seurat_obj@meta.data)) {
      stop("scLCA_celltype column not found")
    }

    truelabel <- seurat_obj$scLCA_celltype
    n_true_clusters <- length(unique(truelabel))

    cat("  - Number of true cell types:", n_true_clusters, "\n")
    cat("  - Cell types:", paste(unique(truelabel), collapse = ", "), "\n")

    # Calculate zero value proportion
    zero_percentage <- sum(datamatrix == 0) / length(datamatrix) * 100
    cat("  - Zero value proportion:", round(zero_percentage, 2), "%\n")

    # Create format required by scLCA
    myscExampleData <- list(
      datamatrix = datamatrix,
      truelabel = as.numeric(as.factor(truelabel))
    )

    # Print data information (similar to example)
    cat("  - scLCA data format:\n")
    cat("    * Data dimensions:", dim(myscExampleData$datamatrix), "\n")
    cat("    * Zero value proportion:", round(zero_percentage, 1), "%\n")
    cat("    * Label distribution:\n")
    print(table(myscExampleData$truelabel))

    # Run scLCA analysis
    cat("  - Starting scLCA clustering analysis...\n")
    myclust.res <- myscLCA(myscExampleData$datamatrix)

    # Extract clustering results
    cluster_labels <- myclust.res[[1]]
    n_estimated_clusters <- length(unique(cluster_labels))
    cat("  - Estimated number of clusters:", n_estimated_clusters, "\n")

    # Calculate confusion matrix
    confusion_matrix <- table(cluster_labels, myscExampleData$truelabel)
    cat("  - Confusion matrix:\n")
    print(confusion_matrix)

    # Calculate ARI
    ari_score <- NA
    if (requireNamespace("mclust", quietly = TRUE)) {
      ari_score <- mclust::adjustedRandIndex(cluster_labels, myscExampleData$truelabel)
      cat("  - Adjusted Rand Index (ARI):", round(ari_score, 4), "\n")
    }

    # Prepare result data frame
    result_df <- data.frame(
      cell_id = colnames(datamatrix),
      cluster = cluster_labels,
      true_label = truelabel,
      true_label_numeric = myscExampleData$truelabel,
      stringsAsFactors = FALSE
    )

    if (!is.na(ari_score)) {
      result_df$ARI <- ari_score
    }

    # Save detailed results
    output_file <- file.path(output_dir, paste0(dataset_name, "_scLCA_results.csv"))
    write.csv(result_df, output_file, row.names = FALSE)

    # Save confusion matrix
    confusion_file <- file.path(output_dir, paste0(dataset_name, "_confusion_matrix.csv"))
    write.csv(confusion_matrix, confusion_file)

    # Save all clustering solutions
    all_clusters_file <- file.path(output_dir, paste0(dataset_name, "_scLCA_all_solutions.rds"))
    saveRDS(myclust.res, all_clusters_file)

    # Save comparison results (similar to table output in example)
    comparison_table <- table(cluster_labels, myscExampleData$truelabel)
    comparison_file <- file.path(output_dir, paste0(dataset_name, "_cluster_comparison.csv"))
    write.csv(comparison_table, comparison_file)

    # Record successful results
    summary_results <- rbind(summary_results, data.frame(
      dataset_name = dataset_name, file_name = file_name,
      n_genes = nrow(datamatrix), n_cells = ncol(datamatrix),
      zero_percentage = round(zero_percentage, 2),
      is_counts_matrix = TRUE,
      n_true_clusters = n_true_clusters,
      n_estimated_clusters = n_estimated_clusters,
      ari_score = ifelse(is.na(ari_score), NA, round(ari_score, 4)),
      status = "Success"
    ))

    cat("  - Successfully processed:", dataset_name, "\n")
    cat("  - Results saved to:", output_file, "\n\n")

    # Clean up memory
    rm(seurat_obj, datamatrix, myclust.res, myscExampleData)
    gc()

  }, error = function(e) {
    cat("  - Error:", e$message, "\n")

    summary_results <<- rbind(summary_results, data.frame(
      dataset_name = dataset_name, file_name = file_name,
      n_genes = NA, n_cells = NA, zero_percentage = NA,
      is_counts_matrix = FALSE, n_true_clusters = NA,
      n_estimated_clusters = NA, ari_score = NA,
      status = paste("Error:", e$message)
    ))

    # Clean up memory
    gc()
  })
}

# Save summary results
summary_file <- file.path(output_dir, "scLCA_summary_results.csv")
write.csv(summary_results, summary_file, row.names = FALSE)

# Print final summary
cat("\n==========================================\n")
cat("Batch processing completed\n")
cat("==========================================\n")
cat("Total number of files:", nrow(summary_results), "\n")
cat("Successfully processed:", sum(summary_results$status == "Success"), "\n")
cat("Failed datasets:", sum(summary_results$status != "Success"), "\n")

# Display successful results
successful <- summary_results[summary_results$status == "Success", ]
if (nrow(successful) > 0) {
  cat("\nSuccessfully processed dataset results:\n")
  print(successful[, c("dataset_name", "n_true_clusters", "n_estimated_clusters", "ari_score")])
}

# Display failed datasets
failed <- summary_results[summary_results$status != "Success", ]
if (nrow(failed) > 0) {
  cat("\nFailed datasets:\n")
  print(failed[, c("dataset_name", "status")])
}

cat("\nAll results saved to:", output_dir, "\n")
cat("Summary results saved to:", summary_file, "\n")
