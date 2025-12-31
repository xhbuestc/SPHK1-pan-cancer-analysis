# Load packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(sva)        # For ComBat batch correction
library(tidyverse)  # For data processing
library(BatchQC)    # For batch effect evaluation

# Set working directory
main_dir <- "TCGA_PanCancer_Analysis"
sub_dirs <- c("00_rawdata", "01_processed", "02_results", "03_logs")
for (dir in file.path(main_dir, sub_dirs)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# ==================== 1. Define Cancer Types ====================
tcga_cancer_types <- paste0("TCGA-", c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", 
  "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", 
  "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", 
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", 
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", 
  "UCEC", "UCS", "UVM"
))

# ==================== 2. Define Data Processing Function ====================
process_cancer_data <- function(cancer_type) {
  cat("Processing:", cancer_type, "\n")
  
  # Build query
  query <- GDCquery(
    project = cancer_type,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # Download data
  GDCdownload(query, method = "api", files.per.chunk = 5, directory = "GDCdata")
  
  # Prepare data
  data_se <- GDCprepare(query)
  
  # Extract TPM data
  tpm_matrix <- assay(data_se, "tpm_unstrand")
  
  # Extract gene names
  gene_info <- rowRanges(data_se)
  if ("gene_name" %in% names(mcols(gene_info))) {
    gene_names <- mcols(gene_info)$gene_name
    gene_names[is.na(gene_names)] <- rownames(tpm_matrix)[is.na(gene_names)]
  } else {
    gene_names <- rownames(tpm_matrix)
  }
  
  # Handle duplicate gene names: take the average
  unique_genes <- unique(gene_names)
  processed_matrix <- matrix(NA, nrow = length(unique_genes), 
                             ncol = ncol(tpm_matrix))
  rownames(processed_matrix) <- unique_genes
  colnames(processed_matrix) <- colnames(tpm_matrix)
  
  for (i in seq_along(unique_genes)) {
    gene <- unique_genes[i]
    indices <- which(gene_names == gene)
    if (length(indices) > 1) {
      processed_matrix[gene, ] <- colMeans(tpm_matrix[indices, , drop = FALSE], 
                                           na.rm = TRUE)
    } else {
      processed_matrix[gene, ] <- tpm_matrix[indices, ]
    }
  }
  
  # Strict filtering: Keep only genes with TPM > 1 in at least 10% of samples
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
      # Extract batch information (using data submission center as batch)
      Batch = substr(barcode, 6, 7)  # The XX part in TCGA-XX-XXXX-XXX
    ) %>%
    dplyr::select(Sample_ID, Cancer_Group, Sample_Type, Batch)
  
  return(list(
    expression = processed_matrix,
    sample_info = sample_info
  ))
}

# ==================== 3. Main Loop: Process All Cancer Types ====================
all_data <- list()
batch_info <- data.frame()

for (i in seq_along(tcga_cancer_types)) {
  cancer <- tcga_cancer_types[i]
  cat(sprintf("\nProcessing cancer type %d/%d: %s\n", i, length(tcga_cancer_types), cancer))
  
  tryCatch({
    # Process data
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

all_data_tpm <-  all_data
# Save data
save(all_data_conut, file = "TCGA_PanCancer_Analysis/00_rawdata/TCGA_all_rawdata_conut.Rdata")
save(all_data_tpm, file = "TCGA_PanCancer_Analysis/00_rawdata/TCGA_all_rawdata_tpm.Rdata")