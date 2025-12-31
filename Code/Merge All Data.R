# ==================== 4. Merge All Data ====================
# 4.1 Identify common genes across all cancer types
common_genes <- Reduce(intersect, lapply(all_data, function(x) rownames(x$expression)))
cat(sprintf("Number of common genes: %d\n", length(common_genes)))

# 4.2 Extract expression matrices for common genes and merge
exp_matrices <- list()
sample_info_all <- data.frame()

for (cancer in names(all_data)) {
  # Extract expression data for common genes
  exp_data <- all_data[[cancer]]$expression[common_genes, , drop = FALSE]
  
  # Transpose: samples as rows, genes as columns
  exp_t <- as.data.frame(t(exp_data))  
  exp_t$Sample_ID <- rownames(exp_t)
  
  # Add sample information
  sample_info <- all_data[[cancer]]$sample_info
  exp_t <- exp_t %>%
    left_join(sample_info, by = "Sample_ID") %>%
    dplyr::select(Sample_ID, Cancer_Group, Sample_Type, Batch, everything())
  
  # Merge into overall data frame
  if (nrow(sample_info_all) == 0) {
    sample_info_all <- exp_t
  } else {
    sample_info_all <- bind_rows(sample_info_all, exp_t)
  }
}

sample_info_all$Sample_ID <- substr(sample_info_all$Sample_ID, 1, 16)  # Take the first part
# Save merged raw data
saveRDS(sample_info_all, file.path(main_dir, "01_processed", 
                                   "combined_raw_tpm_data.rds"))

# ==================== 5. Data Transformation and Batch Correction ====================
# 5.1 Separate expression matrix and sample information
expression_columns <- setdiff(colnames(sample_info_all), 
                              c("Sample_ID", "Cancer_Group", "Sample_Type", "Batch"))
expr_matrix <- as.matrix(sample_info_all[, expression_columns])
rownames(expr_matrix) <- sample_info_all$Sample_ID

# 5.2 log2 transformation (value + 1)
cat("Performing log2 transformation...\n")
expr_log2 <- round(log2(expr_matrix + 1), 2)
cat("Transformed data range:", range(expr_log2), "\n")
cat("First 5 rows and columns:\n")
print(expr_log2[1:5, 1:5])

# 5.3 Transpose to gene×sample matrix (required by ComBat)
expr_log2_t <- t(expr_log2)

# 5.4 Prepare batch information
batch_factor <- as.factor(sample_info_all$Cancer_Group)
cat("Batch information (by cancer type):\n")
print(table(batch_factor))