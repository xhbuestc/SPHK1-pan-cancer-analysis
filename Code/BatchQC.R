# 5.5 Perform batch correction using ComBat
cat("Performing ComBat batch correction...\n")
tryCatch({
  # Check number of batches
  if (length(unique(batch_factor)) > 1) {
    # Use mod parameter to preserve biological differences between cancer types
    mod <- model.matrix(~1, data = data.frame(batch = batch_factor))
    expr_combat <- ComBat(dat = expr_log2_t, 
                          batch = batch_factor, 
                          mod = mod,
                          par.prior = TRUE, 
                          prior.plots = FALSE)
    
    # Transpose back to sample×gene matrix
    expr_combat_t <- t(expr_combat)
    
    # Save batch-corrected data
    saveRDS(expr_combat_t, file.path(main_dir, "01_processed", 
                                     "combined_batch_corrected.rds"))
    
    # Save sample information
    saveRDS(sample_info_all, file.path(main_dir, "01_processed", 
                                       "combined_sample_info.rds"))
    
    cat("Batch correction completed and saved!\n")
  } else {
    cat("Only one batch, skipping batch correction\n")
  }
}, error = function(e) {
  cat("Batch correction failed:", e$message, "\n")
})

# ==================== 6. Batch Effect Evaluation ====================
cat("Performing batch effect evaluation...\n")
tryCatch({
  # Create evaluation data
  expr_for_batchqc <- expr_log2_t
  batch <- as.character(batch_factor)
  
  # Run BatchQC
  batchqc_results <- batchqc(expr_for_batchqc, batch = batch)
  
  # Save evaluation results
  saveRDS(batchqc_results, file.path(main_dir, "02_results", 
                                     "batchqc_results.rds"))
  
  # Generate report
  batchqc_report(batchqc_results, 
                 outfile = file.path(main_dir, "02_results", 
                                     "batchqc_report.html"))
  
  cat("Batch effect evaluation completed!\n")
}, error = function(e) {
  cat("Batch effect evaluation failed:", e$message, "\n")
})

# ==================== 7. Final Data Organization ====================
cat("Organizing final data...\n")

# Create final data frame
if (exists("expr_combat_t")) {
  final_expr <- expr_combat_t
  data_type <- "batch_corrected"
} else {
  final_expr <- expr_log2
  data_type <- "log2_transformed"
}

final_data <- data.frame(
  Sample_ID = rownames(final_expr),
  sample = substr(rownames(final_expr), 1, 15),
  final_expr,
  stringsAsFactors = FALSE
) %>%
  left_join(sample_info_all[, c("Sample_ID", "Cancer_Group", "Sample_Type", "Batch")], 
            by = "Sample_ID")  %>% rename(Cancer = Cancer_Group, Group = Sample_Type) %>% 
  dplyr::select(Sample_ID, sample, Cancer, Group, Batch, everything())
final_data$Cancer <- gsub("^TCGA-", "", final_data$Cancer)
final_data$Group <- factor(final_data$Group, levels = c("Normal", "Tumor"))

# Save final data
fwrite(final_data, 
       file.path(main_dir, "01_processed", 
                 paste0("TCGA_PanCancer_TPM_", data_type, ".csv")))

# Save in R data format
saveRDS(final_data, 
        file.path(main_dir, "01_processed", 
                  paste0("TCGA_PanCancer_TPM_", data_type, ".rds")))

cat("All processing completed!\n")
cat("Data saved in:", main_dir, "\n")