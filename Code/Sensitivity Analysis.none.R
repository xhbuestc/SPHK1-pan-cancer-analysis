# ==================== Sensitivity Analysis for SPHK1 Pan-Cancer Analysis ====================
# Based on TCGA raw count data to evaluate the impact of different preprocessing parameters
# Only compare different expression thresholds, with fixed batch correction method (ComBat)

# Load required packages
library(limma)
library(sva)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(SummarizedExperiment)
library(foreach)
library(doParallel)

# Set up parallel computing (optional, to speed up processing)
# registerDoParallel(cores = 4)

# 0. Load previously saved data
load("TCGA_PanCancer_Analysis/00_rawdata/TCGA_all_rawdata_tpm.Rdata")

# Check data structure
cat("Data structure check:\n")
cat(sprintf("Data contains %d cancer types\n", length(all_data_tpm)))
cat(sprintf("First cancer type (%s) information:\n", names(all_data_tpm)[1]))
cat(sprintf("  Number of genes: %d\n", nrow(all_data_tpm[[1]]$expression)))
cat(sprintf("  Number of samples: %d\n", ncol(all_data_tpm[[1]]$expression)))

# 1. Define sensitivity analysis function (threshold comparison only with fixed ComBat)
perform_sensitivity_analysis_for_cancer <- function(cancer_data, 
                                                    gene_symbol = "SPHK1",
                                                    filter_thresholds = c(5, 10, 20, 30)) {
  
  results <- list()
  
  # Extract count matrix and sample information
  count_matrix <- cancer_data$expression
  sample_info <- cancer_data$sample_info
  
  # Check if the gene exists
  if (!gene_symbol %in% rownames(count_matrix)) {
    cat(sprintf("  Warning: Gene %s not found in the data\n", gene_symbol))
    return(NULL)
  }
  
  # Keep only tumor and normal samples
  sample_info <- sample_info %>%
    filter(Sample_Type %in% c("Tumor", "Normal"))
  
  # Update count matrix to match filtered samples
  common_samples <- intersect(colnames(count_matrix), sample_info$Sample_ID)
  if (length(common_samples) == 0) {
    cat("  Warning: No matching samples\n")
    return(NULL)
  }
  
  count_matrix <- count_matrix[, common_samples, drop = FALSE]
  sample_info <- sample_info %>% filter(Sample_ID %in% common_samples)
  
  # Check if there are sufficient tumor and normal samples
  tumor_count <- sum(sample_info$Sample_Type == "Tumor")
  normal_count <- sum(sample_info$Sample_Type == "Normal")
  
  if (tumor_count < 3 || normal_count < 3) {
    cat(sprintf("  Warning: Insufficient samples (Tumor: %d, Normal: %d)\n", tumor_count, normal_count))
    return(NULL)
  }
  
  # Analyze each filter threshold
  for (threshold in filter_thresholds) {
    
    # Create unique analysis ID
    analysis_id <- paste("Threshold", threshold, sep = "_")
    
    # 1. Gene filtering
    # Keep genes with counts above threshold in at least 10% of samples
    min_samples <- ceiling(0.1 * ncol(count_matrix))
    keep_genes <- rowSums(count_matrix >= threshold) >= min_samples
    
    if (!gene_symbol %in% rownames(count_matrix)[keep_genes]) {
      cat(sprintf("    Gene %s filtered out at threshold %d\n", gene_symbol, threshold))
      next
    }
    
    filtered_counts <- count_matrix[keep_genes, , drop = FALSE]
    
    # 2. log2 transformation (counts + 1)
    expr_log2 <- log2(filtered_counts + 1)
    
    # 3. Fixed batch correction: ComBat (if batch information exists)
    expr_corrected <- expr_log2
    
    if ("Batch" %in% colnames(sample_info)) {
      # Use ComBat batch correction
      batch_factor <- as.factor(sample_info$Batch)
      
      if (length(unique(batch_factor)) > 1) {
        tryCatch({
          expr_corrected <- ComBat(dat = expr_log2, 
                                   batch = batch_factor,
                                   par.prior = TRUE,
                                   prior.plots = FALSE)
          cat(sprintf("    Applied ComBat batch correction for threshold %d\n", threshold))
        }, error = function(e) {
          cat(sprintf("    ComBat batch correction failed for threshold %d: %s\n", threshold, e$message))
          expr_corrected <- expr_log2
        })
      } else {
        cat(sprintf("    No batch correction needed for threshold %d (only 1 batch)\n", threshold))
      }
    } else {
      cat(sprintf("    No batch information available for threshold %d\n", threshold))
    }
    
    # 4. Differential expression analysis (Wilcoxon rank-sum test)
    tumor_samples <- sample_info$Sample_ID[sample_info$Sample_Type == "Tumor"]
    normal_samples <- sample_info$Sample_ID[sample_info$Sample_Type == "Normal"]
    
    tumor_expr <- expr_corrected[gene_symbol, tumor_samples]
    normal_expr <- expr_corrected[gene_symbol, normal_samples]
    
    # Wilcoxon test
    wilcox_test <- wilcox.test(tumor_expr, normal_expr, exact = FALSE)
    
    # Calculate log2 fold change
    log2fc <- median(tumor_expr, na.rm = TRUE) - median(normal_expr, na.rm = TRUE)
    
    # Save results
    results[[analysis_id]] <- data.frame(
      Filter_Threshold = threshold,
      Log2FC = log2fc,
      P_value = wilcox_test$p.value,
      N_Tumor = length(tumor_samples),
      N_Normal = length(normal_samples),
      Median_Tumor = median(tumor_expr, na.rm = TRUE),
      Median_Normal = median(normal_expr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  # Consolidate results
  if (length(results) == 0) {
    return(NULL)
  }
  
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  results_df$Analysis_ID <- paste("Threshold", results_df$Filter_Threshold, sep = "_")
  
  return(results_df)
}

# 2. Execute pan-cancer sensitivity analysis
cat("\nStarting pan-cancer sensitivity analysis (threshold comparison only with fixed ComBat)...\n")

# Store results for all cancer types
all_cancer_results <- list()

for (cancer_name in names(all_data_tpm)) {
  cat(sprintf("Analyzing cancer type: %s\n", cancer_name))
  
  cancer_data <- all_data_tpm[[cancer_name]]
  
  # Perform sensitivity analysis
  sensitivity_results <- perform_sensitivity_analysis_for_cancer(
    cancer_data = cancer_data,
    gene_symbol = "SPHK1",
    filter_thresholds = c(5, 10, 20, 30)
  )
  
  if (!is.null(sensitivity_results)) {
    sensitivity_results$Cancer_Type <- cancer_name
    all_cancer_results[[cancer_name]] <- sensitivity_results
    
    cat(sprintf("  Completed: %d threshold combinations\n", nrow(sensitivity_results)))
  } else {
    cat("  Skipped: No valid results\n")
  }
}

# 3. Consolidate all results
if (length(all_cancer_results) == 0) {
  stop("No analysis results obtained. Please check the data.")
}

combined_results <- do.call(rbind, all_cancer_results)

rownames(combined_results) <- NULL

cat(sprintf("\nAnalysis completed! Analyzed %d cancer types, obtained %d result records\n",
            length(all_cancer_results), nrow(combined_results)))

# 4. Calculate consistency metrics
calculate_consistency_metrics <- function(results_df) {
  
  # Group by cancer type
  cancer_groups <- split(results_df, results_df$Cancer_Type)
  
  consistency_list <- lapply(names(cancer_groups), function(cancer) {
    df <- cancer_groups[[cancer]]
    
    if (nrow(df) < 2) {
      return(data.frame(
        Cancer_Type = cancer,
        N_Analyses = nrow(df),
        Sign_Agreement = NA,
        Significance_Agreement = NA,
        Log2FC_CV = NA,
        Avg_Abs_Log2FC = NA,
        Min_P_Value = NA,
        Max_P_Value = NA,
        Consistent_Across_Thresholds = "Insufficient",
        stringsAsFactors = FALSE
      ))
    }
    
    # Calculate sign consistency (whether log2FC signs are consistent across all analyses)
    signs <- sign(df$Log2FC)
    sign_agreement <- ifelse(length(unique(signs)) == 1, 1, 0)
    
    # Calculate significance consistency (whether all analyses are at the same significance level, using 0.05 as threshold)
    sig_status <- ifelse(df$P_value < 0.05, "Significant", "Not_Significant")
    sig_agreement <- ifelse(length(unique(sig_status)) == 1, 1, 0)
    
    # Calculate coefficient of variation for log2FC
    log2fc_cv <- ifelse(mean(abs(df$Log2FC)) > 0,
                        sd(abs(df$Log2FC)) / mean(abs(df$Log2FC)),
                        NA)
    
    # Determine if consistent across all thresholds
    consistent <- ifelse(sign_agreement == 1 && sig_agreement == 1, 
                         "Yes", "No")
    
    data.frame(
      Cancer_Type = cancer,
      N_Analyses = nrow(df),
      Sign_Agreement = sign_agreement,
      Significance_Agreement = sig_agreement,
      Log2FC_CV = log2fc_cv,
      Avg_Abs_Log2FC = mean(abs(df$Log2FC), na.rm = TRUE),
      Min_P_Value = min(df$P_value, na.rm = TRUE),
      Max_P_Value = max(df$P_value, na.rm = TRUE),
      Consistent_Across_Thresholds = consistent,
      stringsAsFactors = FALSE
    )
  })
  
  consistency_df <- do.call(rbind, consistency_list)
  return(consistency_df)
}

# Calculate consistency metrics
consistency_metrics <- calculate_consistency_metrics(combined_results)

# 5. Create visualizations
# 5.1 Comparison of log2FC across different thresholds
p1 <- ggplot(combined_results, 
             aes(x = Cancer_Type, y = Log2FC, 
                 color = factor(Filter_Threshold))) +
  geom_point(size = 4, alpha = 0.8, 
             position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "gray", linewidth = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 10),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank()) +
  labs(title = "SPHK1 Differential Expression Sensitivity Analysis",
       subtitle = "Comparison of log2 Fold Change under different expression thresholds (ComBat batch correction)",
       x = "Cancer Type", y = "log2 Fold Change (Tumor vs Normal)",
       color = "Expression Threshold") +
  scale_color_manual(values = c("5" = "#FF6B6B", "10" = "#4ECDC4", 
                                "20" = "#45B7D1", "30" = "#96CEB4"))

# 5.2 Comparison of P-values across different thresholds
p2 <- ggplot(combined_results, 
             aes(x = Cancer_Type, y = -log10(P_value + 1e-10), 
                 color = factor(Filter_Threshold))) +
  geom_point(size = 4, alpha = 0.8, 
             position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue", linewidth = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank()) +
  labs(title = "Significance Level Comparison Across Thresholds",
       subtitle = "ComBat batch correction applied",
       x = "Cancer Type", y = "-log10(P-value)",
       color = "Expression Threshold") +
  scale_color_manual(values = c("5" = "#FF6B6B", "10" = "#4ECDC4", 
                                "20" = "#45B7D1", "30" = "#96CEB4"))

# 5.3 Threshold comparison boxplot for each cancer type
p3 <- ggplot(combined_results, 
             aes(x = factor(Filter_Threshold), y = Log2FC, 
                 fill = factor(Filter_Threshold))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_wrap(~ Cancer_Type, scales = "free_y", ncol = 5) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank()) +
  labs(title = "Threshold Sensitivity Analysis by Cancer Type",
       subtitle = "Effect of expression threshold on SPHK1 log2FC (ComBat batch correction)",
       x = "Expression Threshold", y = "log2 Fold Change (Tumor vs Normal)",
       fill = "Threshold") +
  scale_fill_manual(values = c("5" = "#FF6B6B", "10" = "#4ECDC4", 
                               "20" = "#45B7D1", "30" = "#96CEB4"))

# 5.4 Summary of consistency metrics
p4 <- consistency_metrics %>%
  filter(!is.na(Sign_Agreement)) %>%
  select(Cancer_Type, Sign_Agreement, Significance_Agreement) %>%
  pivot_longer(cols = c(Sign_Agreement, Significance_Agreement),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = ifelse(Metric == "Sign_Agreement", 
                         "Sign Consistency", "Significance Consistency")) %>%
  ggplot(aes(x = Cancer_Type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
           width = 0.6) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 0.5) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank()) +
  labs(title = "Result Consistency Across Different Expression Thresholds",
       subtitle = "ComBat batch correction applied",
       x = "Cancer Type", y = "Consistency Proportion", fill = "Consistency Type") +
  scale_fill_manual(values = c("Sign Consistency" = "#4DBBD5", 
                               "Significance Consistency" = "#E64B35"))

# 5.5 Identify cancer types with consistent results across all thresholds
consistent_cancers <- consistency_metrics %>%
  filter(Consistent_Across_Thresholds == "Yes")

if (nrow(consistent_cancers) > 0) {
  cat("\nCancer types with consistent results across all thresholds:\n")
  print(consistent_cancers$Cancer_Type)
  
  # Create expression pattern plot for consistent cancer types
  consistent_data <- combined_results %>%
    filter(Cancer_Type %in% consistent_cancers$Cancer_Type)
  
  p5 <- ggplot(consistent_data, 
               aes(x = Cancer_Type, y = Log2FC, 
                   fill = factor(Filter_Threshold))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "right",
          panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor = element_blank()) +
    labs(title = "Cancer Types with Consistent Results Across All Thresholds",
         subtitle = "SPHK1 Expression Patterns Remain Stable Under Different Expression Thresholds (ComBat batch correction)",
         x = "Cancer Type", y = "log2 Fold Change (Tumor vs Normal)",
         fill = "Expression Threshold") +
    scale_fill_manual(values = c("5" = "#FF6B6B", "10" = "#4ECDC4", 
                                 "20" = "#45B7D1", "30" = "#96CEB4"))
  
  # 5.6 Log2FC trend across thresholds for consistent cancers
  p6 <- consistent_data %>%
    ggplot(aes(x = factor(Filter_Threshold), y = Log2FC, group = Cancer_Type, color = Cancer_Type)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 3) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor = element_blank()) +
    labs(title = "Log2FC Trends Across Different Thresholds",
         subtitle = "For cancer types with consistent results (ComBat batch correction)",
         x = "Expression Threshold", y = "log2 Fold Change (Tumor vs Normal)",
         color = "Cancer Type")
}

# 6. Save results
output_dir <- "TCGA_PanCancer_Analysis/02_results/sensitivity_analysis_threshold"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save detailed results
write.csv(combined_results,
          file.path(output_dir, "SPHK1_sensitivity_analysis_threshold_comparison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(consistency_metrics,
          file.path(output_dir, "SPHK1_consistency_metrics_threshold.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Save visualizations
pdf(file.path(output_dir, "SPHK1_threshold_sensitivity_analysis_report.pdf"),
    width = 14, height = 10, onefile = TRUE)

print(p1)
print(p2)
print(p3)
print(p4)

if (nrow(consistent_cancers) > 0) {
  print(p5)
  print(p6)
}

dev.off()

# Also save individual plots
ggsave(file.path(output_dir, "SPHK1_log2FC_threshold_comparison.pdf"), 
       plot = p1, width = 12, height = 6)
ggsave(file.path(output_dir, "SPHK1_pvalue_threshold_comparison.pdf"), 
       plot = p2, width = 12, height = 6)
ggsave(file.path(output_dir, "SPHK1_threshold_sensitivity_by_cancer.pdf"), 
       plot = p3, width = 15, height = 10)
ggsave(file.path(output_dir, "SPHK1_consistency_metrics_threshold.pdf"), 
       plot = p4, width = 10, height = 6)

if (nrow(consistent_cancers) > 0) {
  ggsave(file.path(output_dir, "SPHK1_consistent_cancers_threshold.pdf"), 
         plot = p5, width = 10, height = 6)
  ggsave(file.path(output_dir, "SPHK1_log2fc_trends_threshold.pdf"), 
         plot = p6, width = 8, height = 6)
}

