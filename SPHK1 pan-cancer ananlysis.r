# ==============================================================================
# TCGA Pan-Cancer Transcriptome Analysis Pipeline
# Author: zgd
# Date: 2025-12-22
# Description: This script downloads, processes, and analyzes TCGA pan-cancer 
#              TPM data, performs batch correction, differential expression,
#              survival analysis, and visualization.
# ==============================================================================

# ========================== 1. LOAD REQUIRED PACKAGES =========================
library(TCGAbiolinks)    # For TCGA data download
library(SummarizedExperiment)
library(sva)             # For ComBat batch correction
library(tidyverse)       # Data manipulation
library(BatchQC)         # Batch effect evaluation
library(data.table)      # Fast reading/writing
library(ggplot2)         # Plotting
library(ggpubr)          # Publication-ready plots
library(pROC)            # ROC analysis
library(survival)        # Survival analysis
library(survminer)       # Survival plots
library(forestplot)      # Forest plots

# ======================== 2. SETUP WORKING DIRECTORIES ========================
main_dir <- "TCGA_PanCancer_Analysis"
sub_dirs <- c("00_rawdata", "01_processed", "02_results", "03_logs")

for (dir in file.path(main_dir, sub_dirs)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# ====================== 3. DEFINE CANCER TYPES TO ANALYZE =====================
tcga_cancer_types <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
  "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
  "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
  "UCEC", "UCS", "UVM"
)

# ================= 4. FUNCTION TO PROCESS EACH CANCER TYPE ====================
process_cancer_data <- function(cancer_type) {
  # Build query for TCGA data
  query <- GDCquery(
    project = paste0("TCGA-", cancer_type),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # Download data
  GDCdownload(query, method = "api", files.per.chunk = 5, directory = "GDCdata")
  data_se <- GDCprepare(query)
  
  # Extract TPM matrix
  tpm_matrix <- assay(data_se, "tpm_unstranded")
  
  # Extract gene names
  gene_info <- rowRanges(data_se)
  if ("gene_name" %in% names(mcols(gene_info))) {
    gene_names <- mcols(gene_info)$gene_name
    gene_names[is.na(gene_names)] <- rownames(tpm_matrix)[is.na(gene_names)]
  } else {
    gene_names <- rownames(tpm_matrix)
  }
  
  # Handle duplicate gene names by averaging
  unique_genes <- unique(gene_names)
  processed_matrix <- matrix(NA, nrow = length(unique_genes), 
                             ncol = ncol(tpm_matrix))
  rownames(processed_matrix) <- unique_genes
  colnames(processed_matrix) <- colnames(tpm_matrix)
  
  for (gene in unique_genes) {
    indices <- which(gene_names == gene)
    if (length(indices) > 1) {
      processed_matrix[gene, ] <- colMeans(tpm_matrix[indices, , drop = FALSE], 
                                           na.rm = TRUE)
    } else {
      processed_matrix[gene, ] <- tpm_matrix[indices, ]
    }
  }
  
  # Filter genes: keep genes with TPM > 1 in at least 10% of samples
  filter_threshold <- 0.1 * ncol(processed_matrix)
  keep_genes <- rowSums(processed_matrix > 1) >= filter_threshold
  processed_matrix <- processed_matrix[keep_genes, ]
  
  # Extract sample information
  sample_info <- colData(data_se) %>%
    as.data.frame() %>%
    mutate(
      Sample_ID = barcode,
      Cancer_Group = cancer_type,
      Sample_Type = case_when(
        substr(barcode, 14, 15) == "01" ~ "Tumor",
        substr(barcode, 14, 15) == "11" ~ "Normal",
        TRUE ~ "Other"
      ),
      # Use data submission center as batch indicator
      Batch = substr(barcode, 6, 7)
    ) %>%
    dplyr::select(Sample_ID, Cancer_Group, Sample_Type, Batch)
  
  return(list(
    expression = processed_matrix,
    sample_info = sample_info
  ))
}

# =========== 5. MAIN LOOP: PROCESS ALL CANCER TYPES & COMBINE DATA ============
all_data <- list()
batch_info <- data.frame()

for (i in seq_along(tcga_cancer_types)) {
  cancer <- tcga_cancer_types[i]
  cat(sprintf("\nProcessing cancer type %d/%d: %s\n", i, 
              length(tcga_cancer_types), cancer))
  
  tryCatch({
    result <- process_cancer_data(cancer)
    all_data[[cancer]] <- result
    
    # Record batch information
    batch_df <- result$sample_info %>%
      mutate(Cancer_Type = cancer) %>%
      dplyr::select(Sample_ID, Batch, Cancer_Type)
    batch_info <- rbind(batch_info, batch_df)
    
    # Save intermediate results
    saveRDS(result, file.path(main_dir, "00_rawdata", 
                              paste0(cancer, "_processed.rds")))
    
    cat(sprintf("  Completed: %d samples, %d genes\n", 
                ncol(result$expression), nrow(result$expression)))
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
  })
  
  # Pause to avoid API limits
  Sys.sleep(3)
}

# Identify common genes across all cancer types
common_genes <- Reduce(intersect, lapply(all_data, function(x) rownames(x$expression)))
cat(sprintf("Number of common genes: %d\n", length(common_genes)))

# Combine all data
exp_matrices <- list()
sample_info_all <- data.frame()

for (cancer in names(all_data)) {
  exp_data <- all_data[[cancer]]$expression[common_genes, , drop = FALSE]
  exp_t <- as.data.frame(t(exp_data))  # Transpose: samples as rows
  exp_t$Sample_ID <- rownames(exp_t)
  
  sample_info <- all_data[[cancer]]$sample_info
  exp_t <- exp_t %>%
    left_join(sample_info, by = "Sample_ID") %>%
    dplyr::select(Sample_ID, Cancer_Group, Sample_Type, Batch, everything())
  
  if (nrow(sample_info_all) == 0) {
    sample_info_all <- exp_t
  } else {
    sample_info_all <- bind_rows(sample_info_all, exp_t)
  }
}

# Save combined raw data
saveRDS(sample_info_all, file.path(main_dir, "01_processed", 
                                   "combined_raw_data.rds"))

# ====================== 6. DATA TRANSFORMATION & BATCH CORRECTION =============
# Separate expression matrix and sample info
expression_columns <- setdiff(colnames(sample_info_all),
                              c("Sample_ID", "Cancer_Group", "Sample_Type", "Batch"))
expr_matrix <- as.matrix(sample_info_all[, expression_columns])
rownames(expr_matrix) <- sample_info_all$Sample_ID

# Log2 transformation (value + 1)
cat("Performing log2 transformation...\n")
expr_log2 <- round(log2(expr_matrix + 1), 2)
expr_log2 <- t(expr_log2)  # Genes as rows, samples as columns

# Prepare batch information
batch_factor <- as.factor(sample_info_all$Cancer_Group)

# ComBat batch correction
tryCatch({
  expr_corrected <- ComBat(
    dat = expr_log2,
    batch = batch_factor,
    mod = NULL,
    par.prior = TRUE,
    prior.plots = FALSE
  )
}, error = function(e) {
  cat("Batch correction failed:", e$message, "\n")
})

# Evaluate batch correction effect with PCA
pca_before <- prcomp(t(expr_log2), scale. = TRUE, center = TRUE)
pca_after <- prcomp(t(expr_corrected), scale. = TRUE, center = TRUE)

# Save PCA results
pca_results <- list(
  before = pca_before,
  after = pca_after
)
saveRDS(pca_results, file.path(main_dir, "02_results", 
                               "pca_batch_analysis.rds"))

# ======================= 7. CREATE FINAL DATA FRAME (TPM_ALL) =================
expr_corrected_t <- t(expr_corrected) %>% as.data.frame()
if (is.null(expr_corrected_t$Sample_ID)) {
  expr_corrected_t$Sample_ID <- rownames(expr_corrected_t)
  cat("Converted row names to Sample_ID column.\n")
}

sample_info_only <- sample_info_all[, c("Sample_ID", "Cancer_Group", "Sample_Type"), 
                                    drop = FALSE]
tpm_all <- sample_info_only %>%
  left_join(expr_corrected_t, by = "Sample_ID")

# Save final TPM matrix
fwrite(tpm_all, "TCGA_PanCancer_Analysis/02_results/tpm_all.csv")

# ====================== 8. EXTRACT PAIRED SAMPLES DATA ========================
tpm_all$Patient_ID <- substr(tpm_all$Sample_ID, 1, 12)

# Identify patients with both tumor and normal samples
paired_patients <- tpm_all %>%
  group_by(Patient_ID, Cancer_Group) %>%
  summarise(
    has_tumor = any(Sample_Type == "Tumor"),
    has_normal = any(Sample_Type == "Normal"),
    .groups = 'drop'
  ) %>%
  filter(has_tumor & has_normal) %>%
  pull(Patient_ID)

# Extract paired sample data
paired_data <- tpm_all %>%
  filter(Patient_ID %in% paired_patients) %>%
  mutate(
    Group = ifelse(Sample_Type == "Tumor", "Tumor", "Normal"),
    ID = Patient_ID,
    Cancer = Cancer_Group
  ) %>%
  dplyr::select(Sample_ID, Cancer, Group, ID, everything()) %>%
  dplyr::select(-Cancer_Group, -Sample_Type, -Patient_ID)

# Save paired data
fwrite(paired_data, "TCGA_PanCancer_Analysis/02_results/paired_tpm.csv")
save(tpm_all, paired_data, 
     file = "TCGA_PanCancer_Analysis/02_results/paired_all_tpm.Rdata")

# ====================== 9. CLINICAL DATA DOWNLOAD & PROCESSING ================
meta <- fread("TCGA_PanCancer_Analysis/00_rawdata/Survival_SupplementalTable_S1_20171025_xena_sp.gz") %>%
  select(sample, `cancer type abbreviation`, age_at_initial_pathologic_diagnosis,
         gender, ajcc_pathologic_tumor_stage, OS:PFI.time) %>%
  rename(Cancer = `cancer type abbreviation`,
         event = OS,
         time_raw = OS.time,
         age_raw = `age_at_initial_pathologic_diagnosis`,
         stage_raw = `ajcc_pathologic_tumor_stage`) %>%
  mutate(patient_id = sub("(TCGA-[^-]+-[^-]+).*", "\\1", sample)) %>%
  mutate(
    time = round(time_raw / 365.25, 2),
    age = as.integer(age_raw),
    gender = substr(gender, 1, 1),
    stage = gsub("Stage ", "", stage_raw),
    DSS.time = round(DSS.time / 365.25, 2),
    DFI.time = round(DFI.time / 365.25, 2),
    PFI.time = round(PFI.time / 365.25, 2)
  ) %>%
  mutate(
    stage_trimmed = trimws(stage),
    clean_stage = case_when(
      stage_trimmed %in% c("I", "II", "III", "IV") ~ stage_trimmed,
      grepl("^I+[ABC]*$", stage_trimmed) & nchar(stage_trimmed) <= 4 ~ 
        gsub("[ABC]", "", stage_trimmed),
      stage_trimmed == "0" ~ "0",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(-stage_trimmed) %>%
  mutate(
    stage_numeric = case_when(
      clean_stage == "0" ~ 0,
      clean_stage == "I" ~ 1,
      clean_stage == "II" ~ 2,
      clean_stage == "III" ~ 3,
      clean_stage == "IV" ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  select(patient_id, Cancer, event, time, age, gender, stage_numeric, 
         DSS:PFI.time) %>%
  as.data.frame()

# Remove duplicate patient records
meta_dedup <- meta[!duplicated(meta$patient_id), ]
rownames(meta_dedup) <- meta_dedup$patient_id

cat(sprintf("Retained %d unique patients after deduplication (removed %d duplicates).\n",
            nrow(meta_dedup), nrow(meta)-nrow(meta_dedup)))

fwrite(meta_dedup, "TCGA_PanCancer_Analysis/02_results/meta_dedup.csv")

# ========================== 10. DOWNSTREAM ANALYSIS ===========================
# Define target gene for analysis
gene <- "SPHK1"

# Define TCGA cancer types for pan-cancer analysis
cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD",
  "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC",
  "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC",
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
  "UCEC", "UCS", "UVM"
)

# Define color palettes for visualization
color2 <- c("#4DBBD5", "#E64B35")  # For boxplots (Normal, Tumor)
color1 <- c("#E64B35", "#4DBBD5")  # For survival plots (High, Low)

# =================== 10.1 PAN-CANCER DIFFERENTIAL EXPRESSION =================
pan_boxplot <- function(gene, palette = color2, legend = "right", method = "wilcox.test") {
  # This function creates boxplots comparing gene expression between Tumor and Normal samples across all cancer types
  # Input: gene name, color palette, legend position, statistical test method
  # Output: ggplot object with boxplots and p-values
  
  p <- ggboxplot(dplyr::select(tpm, Cancer, Group, all_of(gene)),
                 x = "Cancer", y = gene, fill = "Group", xlab = "", color = "black",
                 palette = palette) +
    rotate_x_text(angle = 90) +
    grids(linetype = "dashed") +
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12)) +
    border("black") +
    theme(legend.position = legend)
  
  # Add statistical comparison between Tumor and Normal groups
  p + stat_compare_means(aes(group = Group), label = "p.signif",
                         method = method, label.y.npc = "top", hide.ns = TRUE)
}

# Generate pan-cancer differential expression boxplot for target gene
pan_boxplot(gene)

# ================ 10.2 PAIRED SAMPLE DIFFERENTIAL ANALYSIS ===================
pan_paired_boxplot <- function(gene, palette = color2, legend = "right", method = "wilcox.test") {
  # This function creates paired boxplots for tumor-normal pairs across cancer types with paired samples
  # Input: gene name, color palette, legend position, statistical test method
  # Output: ggplot object with paired boxplots and p-values
  
  # Cancer types with sufficient paired samples
  paired_cancers <- c("BLCA", "BRCA", "CESC","CHOL",, "COAD", "ESCA", "HNSC", "KICH",
                   "KIRC", "KIRP", "LIHC", "LUAD", "LUSC","PAAD","PCPG",
"PRAD", "PEAD","SARC","SKCM","STAD","THCA", "UCEC")
  
  # Select data for the target gene from paired samples
  df <- dplyr::select(paired_tpm, Cancer, Group, ID, all_of(gene))
  
  # Filter to include only cancer types with paired samples
  df <- df %>% filter(Cancer %in% paired_cancers)
  
  # Create paired boxplot
  p <- ggpaired(df, x = "Group", y = gene, id = "ID", color = "Group",
                palette = palette, add = "jitter", xlab = "", ylab = gene,
                line.color = "gray", line.size = 0.4, facet.by = "Cancer",
                nrow = 1) +
    theme_classic() +
    rotate_x_text(angle = 90) +
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12, colour = "black"),
          strip.text.x = element_text(size = 12, color = "black")) +
    theme(legend.position = legend)
  
  # Add statistical comparison for paired samples
  p + stat_compare_means(label = "p.signif", method = method,
                         paired = TRUE, label.x.npc = 0.4, label.y.npc = "top")
}

# ======================== 10.3 DIAGNOSTIC ROC CURVES ==========================
tcga_roc <- function(cancer, gene) {
  # This function creates ROC curves to evaluate diagnostic potential of a gene for specific cancer type
  # Input: cancer type abbreviation, gene name
  # Output: ggplot object with ROC curve, AUC, and 95% CI
  
  # Filter data for specific cancer type and select relevant columns
  df <- dplyr::filter(tpm, Cancer == cancer) %>%
    dplyr::select(Group, all_of(gene))
  
  # Calculate ROC curve
  res <- pROC::roc_(df, "Group", gene, aur = TRUE, ci = TRUE,
                    smooth = TRUE, levels = c("Normal", "Tumor"))
  
  # Create ROC plot
  p <- pROC::ggroc(res, color = "#4DBBD5", legacy.axes = TRUE) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                 color = "darkgrey", linetype = 4) +
    theme_bw() +
    ggtitle(paste0(cancer, " ROC Curve")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.title = element_text(size = 14, colour = "black"),
          legend.text = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 14, colour = "black"),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"))
  
  # Add AUC and 95% CI annotation
  p + ggplot2::annotate("text", x = 0.75, y = 0.25,
                        label = paste("AUC = ", round(res$auc, 3))) +
    ggplot2::annotate("text", x = 0.75, y = 0.2,
                      label = paste("95%CI: ", round(res$ci[1], 3),
                                    "-", round(res$ci[3], 3)))
}

# Generate ROC curves for all cancer types
for (cancer_type in cancers) {
  roc_plot <- tcga_roc(cancer_type, gene)
  
  # Save ROC plot as PDF
  pdf(file = paste0("00.fig/16.", cancer_type, gene, "tcga_roc.pdf"),
      width = 5, height = 5)
  print(roc_plot)
  dev.off()
}

# ====================== 10.4 BATCH SURVIVAL ANALYSIS =========================
tcga_km <- function(cancer, gene, palette = "jco") {
  # This function performs Kaplan-Meier survival analysis for a gene in specific cancer type
  # Input: cancer type abbreviation, gene name, color palette
  # Output: ggsurvplot object with survival curves and risk table
  
  # Prepare expression data for tumor samples
  exprSet <- subset(tpm, Group == "Tumor" & Cancer == cancer) %>%
    dplyr::select(all_of(gene)) %>%
    tibble::add_column(ID = stringr::str_sub(rownames(.), 1, 12)) %>%
    dplyr::filter(!duplicated(ID)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("ID") %>%
    dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == cancer))) %>%
    t() %>%
    as.matrix()
  
  # Get clinical data for samples with expression data
  cl <- meta[colnames(exprSet), ]
  
  # Dichotomize expression into high/low groups based on median
  cl$expression <- ifelse(exprSet[gene, ] > median(exprSet[gene, ]),
                          "high", "low")
  
  # Fit survival model
  sfit <- survival::survfit(survival::Surv(time, event) ~ expression,
                            data = cl)
  
  # Create survival plot
  survminer::ggsurvplot(sfit, pval = TRUE, palette = palette,
                        data = cl, legend = c(0.8, 0.8),
                        title = paste0("KMplot of ", gene, " in ", cancer),
                        risk.table = TRUE)
}

# Generate Kaplan-Meier plots for all cancer types
for (cancer_type in cancers) {
  km_plot <- tcga_km(cancer_type, gene, color1)
  
  # Save KM plot as PDF
  pdf(file = paste0("00.fig/20.", cancer_type, gene, "tcga_kmplot.pdf"),
      width = 5, height = 6)
  print(km_plot)
  dev.off()
}

# ====================== 10.5 BATCH COX REGRESSION ANALYSIS ====================
pan_forest <- function(gene, adjust = FALSE) {
  # This function performs Cox regression analysis across all cancer types
  # Input: gene name, whether to adjust for age
  # Output: forest plot of hazard ratios
  
  if (adjust == TRUE) {
    # Age-adjusted Cox regression
    cox_results <- list()
    
    for (cancer in cancers) {
      # Prepare expression data
      exprSet <- subset(tpm, Group == "Tumor" & Cancer == cancer) %>%
        tibble::add_column(ID = stringr::str_sub(rownames(.), 1, 12),
                           .before = "Cancer") %>%
        dplyr::filter(!duplicated(ID)) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("ID") %>%
        dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == cancer)))
      
      exprSet <- exprSet[, -(1:2)]
      exprSet <- as.matrix(t(exprSet))
      
      # Prepare clinical data
      cl <- meta[colnames(exprSet), ]
      cl$symbol <- exprSet[gene, ]
      
      # Fit Cox proportional hazards model adjusted for age
      m <- survival::coxph(survival::Surv(time, event) ~ symbol + age, data = cl)
      
      # Extract model coefficients
      beta <- coef(m)
      se <- sqrt(diag(vcov(m)))
      HR <- exp(beta)
      HRse <- HR * se
      
      # Compile results
      tmp <- round(cbind(coef = beta, se = se, z = beta/se,
                         p = 1 - pchisq((beta/se)^2, 1),
                         HR = HR, HRse = HRse,
                         HRz = (HR - 1)/HRse,
                         HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                         HRCILL = exp(beta - qnorm(0.975, 0, 1) * se),
                         HRCIUL = exp(beta + qnorm(0.975, 0, 1) * se)), 3)
      
      cox_results[[cancer]] <- tmp["symbol", ]
    }
    
    # Combine results from all cancer types
    cox_results <- do.call(rbind, cox_results)
    cox_results <- as.data.frame(cox_results[, c(5, 9:10, 4)])  # HR, CI, p-value
    
    # Format HR and 95% CI for display
    np <- paste0(cox_results$HR, " (", cox_results$HRCILL,
                 "-", cox_results$HRCIUL, ")")
    
    # Prepare table text for forest plot
    tabletext <- cbind(c("Cancer", rownames(cox_results)),
                       c("HR (95%CI)", np),
                       c("P Value", cox_results$p))
    
    # Create forest plot
    forestplot::forestplot(labeltext = tabletext, graph.pos = 3,
                           mean = c(NA, cox_results$HR),
                           lower = c(NA, cox_results$HRCILL),
                           upper = c(NA, cox_results$HRCIUL),
                           title = paste0("Hazard Ratio Plot of ", gene, " adjusted by age"),
                           hrzl_lines = list(`1` = grid::gpar(lwd = 2, col = "black"),
                                             `2` = grid::gpar(lwd = 2, col = "black"),
                                             `35` = grid::gpar(lwd = 2, col = "black")),
                           is.summary = c(TRUE, rep(FALSE, 33)),
                           col = forestplot::fpColors(box = "#4DBBD5",
                                                      lines = "#4DBBD5",
                                                      zero = "gray50"),
                           zero = 1, cex = 0.9, lineheight = "auto",
                           colgap = unit(8, "mm"),
                           txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 1)),
                           boxsize = 0.5, ci.vertices = TRUE,
                           ci.vertices.height = 0.3)
  } else {
    # Unadjusted Cox regression
    cox_results <- list()
    
    for (cancer in cancers) {
      # Prepare expression data
      exprSet <- subset(tpm, Group == "Tumor" & Cancer == cancer) %>%
        tibble::add_column(ID = stringr::str_sub(rownames(.), 1, 12),
                           .before = "Cancer") %>%
        dplyr::filter(!duplicated(ID)) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames("ID") %>%
        dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == cancer)))
      
      exprSet <- exprSet[, -(1:2)]
      exprSet <- as.matrix(t(exprSet))
      
      # Prepare clinical data
      cl <- meta[colnames(exprSet), ]
      cl$symbol <- exprSet[gene, ]
      
      # Fit unadjusted Cox proportional hazards model
      m <- survival::coxph(survival::Surv(time, event) ~ symbol, data = cl)
      
      # Extract model coefficients
      beta <- coef(m)
      se <- sqrt(diag(vcov(m)))
      HR <- exp(beta)
      HRse <- HR * se
      
      # Compile results
      tmp <- round(cbind(coef = beta, se = se, z = beta/se,
                         p = 1 - pchisq((beta/se)^2, 1),
                         HR = HR, HRse = HRse,
                         HRz = (HR - 1)/HRse,
                         HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                         HRCILL = exp(beta - qnorm(0.975, 0, 1) * se),
                         HRCIUL = exp(beta + qnorm(0.975, 0, 1) * se)), 3)
      
      cox_results[[cancer]] <- tmp["symbol", ]
    }
    
    # Combine results from all cancer types
    cox_results <- do.call(rbind, cox_results)
    cox_results <- as.data.frame(cox_results[, c(5, 9:10, 4)])  # HR, CI, p-value
    
    # Format HR and 95% CI for display
    np <- paste0(cox_results$HR, " (", cox_results$HRCILL,
                 "-", cox_results$HRCIUL, ")")
    
    # Prepare table text for forest plot
    tabletext <- cbind(c("Cancer", rownames(cox_results)),
                       c("HR (95%CI)", np),
                       c("P Value", cox_results$p))
    
    # Create forest plot
    forestplot::forestplot(labeltext = tabletext, graph.pos = 3,
                           mean = c(NA, cox_results$HR),
                           lower = c(NA, cox_results$HRCILL),
                           upper = c(NA, cox_results$HRCIUL),
                           title = paste0("Hazard Ratio Plot of ", gene),
                           hrzl_lines = list(`1` = grid::gpar(lwd = 2, col = "black"),
                                             `2` = grid::gpar(lwd = 2, col = "black"),
                                             `35` = grid::gpar(lwd = 2, col = "black")),
                           is.summary = c(TRUE, rep(FALSE, 33)),
                           col = forestplot::fpColors(box = "#4DBBD5",
                                                      lines = "#4DBBD5",
                                                      zero = "gray50"),
                           zero = 1, cex = 0.9, lineheight = "auto",
                           colgap = unit(8, "mm"),
                           txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 1)),
                           boxsize = 0.5, ci.vertices = TRUE,
                           ci.vertices.height = 0.3)
  }
}

# Generate forest plot for target gene (unadjusted)
pan_forest(gene)

# ========================== 11. TCGA PLOT ANALYSIS ============================
# Note: The following functions are assumed to be defined elsewhere in the analysis pipeline
# They perform various TCGA-specific analyses and visualizations
Library(TCGAPlot)
# Microsatellite instability (MSI) radar plot
gene_MSI_radar(gene)

# Tumor mutational burden (TMB) radar plot
gene_TMB_radar(gene)

# Immune cell infiltration heatmap
gene_immucell_heatmap(gene)  # immune cell analysis

# Immune checkpoint gene expression heatmap
gene_checkpoint_heatmap(gene)  # checkpoint analysis

# Chemokine gene expression heatmap
gene_chemokine_heatmap(gene)  # chemokine analysis

# Receptor gene expression correlation heatmap
gene_receptor_heatmap(gene, method = "spearman")  # receptor analysis

# Gene Set Enrichment Analysis (GSEA) for KEGG pathways
# Note: Requires specific cancer type; typically run in a loop for multiple cancers
 gene_gsea_kegg(cancer_type, gene)  # GSEA analysis
```
