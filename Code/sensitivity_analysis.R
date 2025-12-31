# ==================== Sensitivity Analysis for SPHK1 Pan-Cancer Analysis ====================
# Based on TCGA raw count data to evaluate the impact of different preprocessing parameters

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

# 1. Define sensitivity analysis function
perform_sensitivity_analysis_for_cancer <- function(cancer_data, 
                                                    gene_symbol = "SPHK1",
                                                    filter_thresholds = c(5,10, 20, 30),
                                                    batch_methods = c("ComBat", "limma_removeBatchEffect", "none")) {
  
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
  
  # Analyze each combination of filter threshold and batch correction method
  for (threshold in filter_thresholds) {
    for (batch_method in batch_methods) {
      
      # Create unique analysis ID
      analysis_id <- paste(threshold, batch_method, sep = "_")
      
      # 1. Gene filtering
      # Keep genes with counts above threshold in at least 10% of samples
      min_samples <- ceiling(0.1 * ncol(count_matrix))
      keep_genes <- rowSums(count_matrix >= threshold) >= min_samples
      
      if (!gene_symbol %in% rownames(count_matrix)[keep_genes]) {
        next
      }
      
      filtered_counts <- count_matrix[keep_genes, , drop = FALSE]
      
      # 2. log2 transformation (counts + 1)
      expr_log2 <- log2(filtered_counts + 1)
      
      # 3. Batch correction
      expr_corrected <- expr_log2
      
      if (batch_method == "ComBat" && "Batch" %in% colnames(sample_info)) {
        # Use ComBat batch correction
        batch_factor <- as.factor(sample_info$Batch)
        
        if (length(unique(batch_factor)) > 1) {
          tryCatch({
            expr_corrected <- ComBat(dat = expr_log2, 
                                     batch = batch_factor,
                                     par.prior = TRUE,
                                     prior.plots = FALSE)
          }, error = function(e) {
            cat(sprintf("    ComBat batch correction failed: %s\n", e$message))
            expr_corrected <- expr_log2
          })
        }
      } else if (batch_method == "limma_removeBatchEffect" && "Batch" %in% colnames(sample_info)) {
        # Use limma's removeBatchEffect
        batch_factor <- as.factor(sample_info$Batch)
        sample_type <- as.factor(ifelse(sample_info$Sample_Type == "Tumor", 1, 0))
        
        if (length(unique(batch_factor)) > 1) {
          design <- model.matrix(~ sample_type)
          expr_corrected <- removeBatchEffect(expr_log2, 
                                              batch = batch_factor,
                                              design = design)
        }
      }
      # If batch_method is "none", skip batch correction
      
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
        Batch_Method = batch_method,
        Log2FC = log2fc,
        P_value = wilcox_test$p.value,
        N_Tumor = length(tumor_samples),
        N_Normal = length(normal_samples),
        Median_Tumor = median(tumor_expr, na.rm = TRUE),
        Median_Normal = median(normal_expr, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Consolidate results
  if (length(results) == 0) {
    return(NULL)
  }
  
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  results_df$Analysis_ID <- paste(results_df$Filter_Threshold, 
                                  results_df$Batch_Method, sep = "_")
  
  return(results_df)
}

# 2. Execute pan-cancer sensitivity analysis
cat("\nStarting pan-cancer sensitivity analysis...\n")

# Store results for all cancer types
all_cancer_results <- list()

for (cancer_name in names(all_data_tpm)) {
  cat(sprintf("Analyzing cancer type: %s\n", cancer_name))
  
  cancer_data <- all_data_tpm[[cancer_name]]
  
  # Perform sensitivity analysis
  sensitivity_results <- perform_sensitivity_analysis_for_cancer(
    cancer_data = cancer_data,
    gene_symbol = "SPHK1",
    filter_thresholds = c(5,10, 20, 30),
    batch_methods = c("ComBat", "limma_removeBatchEffect", "none")
  )
  
  if (!is.null(sensitivity_results)) {
    sensitivity_results$Cancer_Type <- cancer_name
    all_cancer_results[[cancer_name]] <- sensitivity_results
    
    cat(sprintf("  Completed: %d parameter combinations\n", nrow(sensitivity_results)))
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
        Consistent_Across_Methods = "Insufficient",
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
    
    # Determine if consistent across all methods
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
      Consistent_Across_Methods = consistent,
      stringsAsFactors = FALSE
    )
  })
  
  consistency_df <- do.call(rbind, consistency_list)
  return(consistency_df)
}

# Calculate consistency metrics
consistency_metrics <- calculate_consistency_metrics(combined_results)

# 5. Create visualizations
# 5.1 Comparison of log2FC across different methods
p1 <- ggplot(combined_results, 
             aes(x = Cancer_Type, y = Log2FC, 
                 color =  interaction(Filter_Threshold, Batch_Method)) )+
  geom_point(size = 3, alpha = 0.7, 
             position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "gray", linewidth = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right",
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank()) +
  labs(title = "SPHK1 Differential Expression Sensitivity Analysis",
       subtitle = "Comparison of log2 Fold Change under different preprocessing parameters",
       x = "Cancer Type", y = "log2 Fold Change (Tumor vs Normal)",
       color = "Method (Threshold_BatchCorrection)") +
  scale_color_brewer(palette = "Set3", 
                     labels = function(x) gsub("\\.", "_", x))
ggsave(filename = "TCGA_PanCancer_Analysis/02_results/sensitivity_analysis/SPHK1 Differential Expression Sensitivity Analysis.pdf",plot = p1,width = 10,height = 5)
# 5.2 Comparison of P-values across different methods
p2 <- ggplot(combined_results, 
             aes(x = Cancer_Type, y = -log10(P_value + 1e-10), 
                 color =  interaction(Filter_Threshold, Batch_Method))) +
  geom_point(size = 3, alpha = 0.7, 
             position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue", linewidth = 0.3) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none",
        panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = element_blank()) +
  labs(title = "Significance Level Comparison",
       x = "Cancer Type", y = "-log10(P-value)",
       color = "Method (Threshold_BatchCorrection)") +
  scale_color_brewer(palette = "Set3")
ggsave(filename = "TCGA_PanCancer_Analysis/02_results/sensitivity_analysis/SPHK1 Significance Level Comparison.pdf",plot = p2,width = 10,height = 5)
# 5.3 Consistency heatmap
# Prepare data: calculate log2FC correlations between different methods for each cancer type
create_correlation_heatmap <- function(results_df) {
  # Convert results to wide format
  wide_data <- results_df %>%
    mutate(Method_ID = paste(Filter_Threshold, Batch_Method, sep = "_")) %>%
    select(Cancer_Type, Method_ID, Log2FC) %>%
    pivot_wider(names_from = Method_ID, values_from = Log2FC) %>%
    column_to_rownames("Cancer_Type")
  
  # Calculate correlation matrix (using complete pairwise observations)
  cor_matrix <- cor(wide_data, use = "pairwise.complete.obs", method = "spearman")
  
  # Convert to long format
  cor_long <- as.data.frame(cor_matrix) %>%
    rownames_to_column("Method1") %>%
    pivot_longer(cols = -Method1, names_to = "Method2", values_to = "Correlation") %>%
    filter(Method1 != Method2)  # Remove diagonal
  
  # Create heatmap
  p <- ggplot(cor_long, aes(x = Method1, y = Method2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", Correlation)), 
              color = "black", size = 3) +
    scale_fill_gradient2(low = "#2E8B57", mid = "white", high = "#E64B35", 
                         midpoint = 0, limits = c(-1, 1),
                         name = "Spearman\nCorrelation") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          panel.grid = element_blank(),
          legend.position = "right") +
    labs(title = "Correlation Between Different Analysis Methods",
         x = "Analysis Method", y = "Analysis Method")
  
  return(p)
}

p3 <- create_correlation_heatmap(combined_results)
ggsave(filename = "TCGA_PanCancer_Analysis/02_results/sensitivity_analysis/Correlation Between Different Analysis Methods.pdf",plot = p3,width = 8,height = 7)
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
  labs(title = "Summary of Result Consistency Metrics",
       x = "Cancer Type", y = "Consistency Proportion", fill = "Consistency Type") +
  scale_fill_manual(values = c("Sign Consistency" = "#4DBBD5", 
                               "Significance Consistency" = "#E64B35"))


ggsave(filename = "TCGA_PanCancer_Analysis/02_results/sensitivity_analysis/Summary of Result Consistency Metrics.pdf",plot = p4,width = 8,height = 4.5)
# 5.5 Identify cancer types with consistent results across all methods
consistent_cancers <- consistency_metrics %>%
  filter(Consistent_Across_Methods == "Yes")

if (nrow(consistent_cancers) > 0) {
  cat("\nCancer types with consistent results across all analysis methods:\n")
  print(consistent_cancers$Cancer_Type)
  
  # Create expression pattern plot for consistent cancer types
  consistent_data <- combined_results %>%
    filter(Cancer_Type %in% consistent_cancers$Cancer_Type)
  
  p5 <- ggplot(consistent_data, 
               aes(x = Cancer_Type, y = Log2FC, 
                   fill = interaction(Filter_Threshold, Batch_Method))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "bottom",
          legend.text = element_text(size = 8),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
          panel.grid.minor = element_blank()) +
    labs(title = "Cancer Types with Consistent Results Across All Methods",
         subtitle = "SPHK1 Expression Patterns Remain Stable Under Different Preprocessing Parameters",
         x = "Cancer Type", y = "log2 Fold Change (Tumor vs Normal)",
         fill = "Analysis Method") +
    scale_fill_brewer(palette = "Set3")
}

ggsave(filename = "TCGA_PanCancer_Analysis/02_results/sensitivity_analysis/SPHK1 Expression Patterns Remain Stable Under Different Preprocessing Parameters.pdf",plot = p5,width = 8.5,height = 4.5)
fwrite(combined_results,"combined_results.csv")

# 6. Save results
output_dir <- "TCGA_PanCancer_Analysis/02_results/sensitivity_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save detailed results
write.csv(combined_results,
          file.path(output_dir, "SPHK1_sensitivity_analysis_detailed.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

write.csv(consistency_metrics,
          file.path(output_dir, "SPHK1_consistency_metrics.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Save visualizations
pdf(file.path(output_dir, "SPHK1_sensitivity_analysis_report.pdf"),
    width = 16, height = 12, onefile = TRUE)

