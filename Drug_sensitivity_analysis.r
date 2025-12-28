# ==============================================================================
# DRUG SENSITIVITY ANALYSIS USING CELLMINER DATABASE
# This script analyzes the correlation between gene expression (SPHK1) and 
# drug sensitivity in NCI-60 cancer cell lines using CellMiner database data
# ==============================================================================

# -------------------------- 1. DATA DOWNLOAD INSTRUCTIONS ---------------------
# 1.1 Go to the CellMiner database homepage: https://discover.nci.nih.gov/cellminer/home.do
# 1.2 Click "Download Data Sets" to enter the data download interface
# 1.3 Under "Processed Data Set", download:
#     - RNA expression data (RNA: RNA-seq - composite expression): 
#       https://discover.nci.nih.gov/cellminer/download/processeddataset/nci60_RNA__RNA_seq_composite_expression.zip
#     - Drug sensitivity data (Compound activity: DTP NCI-60 - Average z score):
#       https://discover.nci.nih.gov/cellminer/download/processeddataset/DTP_NCI60_ZSCORE.zip

# ====================== 2. PREPARATION OF DRUG SENSITIVITY DATA ===============

# Clear workspace to start fresh analysis
rm(list = ls())

# Define target gene for analysis
target_gene <- "SPHK1"

# ---------------------- 2.1 LOAD REQUIRED PACKAGES ----------------------------
# Load readxl package for reading Excel files
library(readxl)

# ---------------------- 2.2 READ DRUG SENSITIVITY DATA ------------------------
# Read drug sensitivity data (skipping first 7 rows which contain metadata)
drug_data <- read_excel(path = "00.rawdata/00.CellMiner/DTP_NCI60_ZSCORE.xlsx", 
                        skip = 7)

# Set column names from the first row and remove the first row (now used as header)
colnames(drug_data) <- drug_data[1, ]
drug_data <- drug_data[-1, -c(67, 68)]  # Remove metadata columns

# Check FDA approval status distribution
# Note: 75 drugs in clinical trials, 188 drugs FDA approved
table(drug_data$`FDA status`)

# Filter drugs: keep only FDA approved or in clinical trial for reliability
drug_data <- drug_data[drug_data$`FDA status` %in% c("FDA approved", "Clinical trial"), ]

# Remove unnecessary columns (keeping only drug names and sensitivity data)
drug_data <- drug_data[, -c(1, 3:6)]

# Save filtered drug data
write.csv(drug_data, 
          file = "00.rawdata/00.CellMiner/00.CellMiner.drug.csv",
          row.names = FALSE)

# ====================== 3. PREPARATION OF GENE EXPRESSION DATA ===============

# Read RNA-seq composite expression data (skipping first 9 rows of metadata)
expression_data <- read_excel(
  path = "00.rawdata/00.CellMiner/RNA__RNA_seq_composite_expression.xls",
  skip = 9
)

# Set column names from the first row and remove the first row
colnames(expression_data) <- expression_data[1, ]
expression_data <- expression_data[-1, -c(2:6)]  # Remove metadata columns

# Save processed expression data
write.csv(expression_data,
          file = "00.rawdata/00.CellMiner/00.CellMiner.geneExp.csv",
          row.names = FALSE)

# ====================== 4. DRUG SENSITIVITY ANALYSIS =========================

# Clear workspace to prepare for correlation analysis
rm(list = ls())

# ---------------------- 4.1 LOAD REQUIRED PACKAGES ---------------------------
library(impute)  # For k-nearest neighbor imputation of missing data
library(limma)   # For averaging replicates
library(ggplot2) # For visualization
library(ggpubr)  # For arranging multiple plots

# ---------------------- 4.2 PROCESS DRUG SENSITIVITY DATA --------------------

# Read processed drug sensitivity data
drug_df <- read.csv("00.rawdata/00.CellMiner/00.CellMiner.drug.csv")

# Clean drug names: replace hyphens with underscores for consistency
drug_df$Drug.name <- gsub("-", "_", drug_df$Drug.name)

# Convert to matrix format for processing
drug_matrix <- as.matrix(drug_df)
rownames(drug_matrix) <- drug_matrix[, 1]  # Set drug names as row names
drug_sensitivity <- drug_matrix[, 2:ncol(drug_matrix)]  # Extract sensitivity data

# Create dimension names and convert to numeric matrix
dim_names <- list(rownames(drug_sensitivity), colnames(drug_sensitivity))
numeric_drug_data <- matrix(as.numeric(as.matrix(drug_sensitivity)),
                            nrow = nrow(drug_sensitivity),
                            dimnames = dim_names)

# Remove columns with >80% missing values
high_missing_cols <- colSums(is.na(numeric_drug_data)) / nrow(numeric_drug_data) > 0.8
filtered_drug_data <- numeric_drug_data[, !high_missing_cols]

# Impute missing values using k-nearest neighbors
imputed_data <- impute.knn(filtered_drug_data)
drug_imputed <- imputed_data$data

# Average replicates if present
drug_final <- avereps(drug_imputed)

# ---------------------- 4.3 PROCESS GENE EXPRESSION DATA ---------------------

# Read processed gene expression data
expression_df <- read.csv("00.rawdata/00.CellMiner/00.CellMiner.geneExp.csv",
                          row.names = 1)

# Check data dimensions
# Note: Contains 60 cancer cell lines and 23,805 genes
dim(expression_df)

# Preview first few rows and columns
expression_df[1:4, 1:4]

# ---------------------- 4.4 EXTRACT TARGET GENE EXPRESSION -------------------

# Clean target gene name (remove spaces)
target_gene <- gsub(" ", "", target_gene)

# Check if target gene exists in expression data
target_gene <- intersect(target_gene, rownames(expression_df))

# Extract expression of target gene for cell lines with drug sensitivity data
target_expression <- expression_df[target_gene, colnames(drug_final)]

# ---------------------- 4.5 CORRELATION ANALYSIS -----------------------------

# Initialize empty data frame for results
correlation_results <- data.frame()

# Loop through each target gene (typically one gene in this analysis)
for (gene in rownames(target_expression)) {
  # Get expression values for current gene
  gene_expr <- as.numeric(target_expression[gene, ])
  
  # Loop through each drug
  for (drug in rownames(drug_final)) {
    # Get drug sensitivity values
    drug_sens <- as.numeric(drug_final[drug, ])
    
    # Calculate Pearson correlation
    cor_test <- cor.test(gene_expr, drug_sens, method = "pearson")
    
    # Extract correlation coefficient and p-value
    cor_coef <- cor_test$estimate
    p_value <- cor_test$p.value
    
    # Store significant correlations (p < 0.01)
    if (p_value < 0.01) {
      result_row <- cbind(Gene = gene, Drug = drug, 
                          Correlation = cor_coef, P_value = p_value)
      correlation_results <- rbind(correlation_results, result_row)
    }
  }
}

# Sort results by p-value (most significant first)
correlation_results <- correlation_results[
  order(as.numeric(as.vector(correlation_results$P_value))),
]

# Save correlation results
write.table(correlation_results,
            file = "00.rawdata/00.CellMiner/00.CellMiner.drugCor.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------------- 4.6 VISUALIZATION OF RESULTS -------------------------

# Set number of top correlations to plot (maximum 16)
num_plots <- min(16, nrow(correlation_results))

# Initialize list to store individual plots
plot_list <- list()

# Generate scatter plots for top correlations
for (i in 1:num_plots) {
  gene_name <- correlation_results[i, 1]
  drug_name <- correlation_results[i, 2]
  
  # Extract expression and sensitivity data
  x_data <- as.numeric(expression_df[gene_name, ])
  y_data <- as.numeric(drug_final[drug_name, ])
  
  # Format correlation coefficient and p-value for display
  cor_formatted <- sprintf("%.03f", as.numeric(correlation_results[i, 3]))
  
  if (as.numeric(correlation_results[i, 4]) < 0.001) {
    p_display <- "p < 0.001"
  } else {
    p_display <- paste0("p = ", 
                       sprintf("%.03f", as.numeric(correlation_results[i, 4])))
  }
  
  # Create data frame for plotting
  plot_df <- data.frame(Expression = x_data, Sensitivity = y_data)
  
  # Create scatter plot with regression line
  plot <- ggplot(data = plot_df, aes(x = Expression, y = Sensitivity)) +
    geom_point(size = 1) +
    stat_smooth(method = "lm", se = FALSE, formula = y ~ x) +
    labs(x = "", y = "", 
         title = paste0(gene_name, ", ", drug_name),
         subtitle = paste0("Cor = ", cor_formatted, ", ", p_display)) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank()) +
    theme_bw()
  
  plot_list[[i]] <- plot
}

# Calculate grid dimensions for arranging plots
n_rows <- ceiling(sqrt(num_plots))
n_cols <- ceiling(num_plots / n_rows)

# Save combined plot as PDF
pdf(file = "00.rawdata/00.CellMiner/00.CellMiner.drugCor_SPHK1.pdf",
    width = 13, height = 9)
ggarrange(plotlist = plot_list, nrow = n_rows, ncol = n_cols)
dev.off()

