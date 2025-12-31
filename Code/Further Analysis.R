# ============================================
# Required R Packages
# ============================================
# Please ensure the following packages are installed:
# install.packages(c(
#   "ggplot2", "ggpubr", "survival", "survminer", 
#   "pROC", "forestplot", "dplyr", "stringr"
# ))
# ===========================================

# 1. Define TCGA cancer type list
tcga_cancer_types <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", 
  "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", 
  "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", 
  "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", 
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", 
  "UCEC", "UCS", "UVM"
)

# Define color schemes for visualization
plot_colors_2 <- c("#4DBBD5", "#E64B35")  # Normal-Tumor order
plot_colors_1 <- c("#E64B35", "#4DBBD5")  # High-Low expression order

# Extract Single Gene Expression  
SPHK1_tpm <- final_data[, c("sample", "Group", "SPHK1", "Cancer")] %>% na.omit(.)  
gene <- "SPHK1"  

# Concise Version – Pan-cancer Boxplot  
library(ggpubr)  

# Use your data frame SPHK1_tpm to draw a pan-cancer boxplot  
p <- ggboxplot(SPHK1_tpm,  
               x = "Cancer",   
               y = "SPHK1",   
               fill = "Group",   
               xlab = "",   
               color = "black",  
               palette = c("#4DBBD5", "#E64B35")) +   
  rotate_x_text(angle = 90) +   
  grids(linetype = "dashed") +   
  theme(  
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 12),   
    axis.title.x = element_text(size = 14),   
    axis.text.x = element_text(size = 12),   
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 12)  
  ) +   
  border("black") +   
  theme(legend.position = "right")  

# Add statistical test  
p <- p + stat_compare_means(  
  aes(group = Group),   
  label = "p.signif",   
  method = "wilcox.test",   
  label.y.npc = "top",   
  hide.ns = TRUE  
)  

# Display the plot  
print(p)  

# Save the plot  
ggsave("SPHK1_pan_cancer_boxplot_ggpubr.png",   
       p, width = 16, height = 8, dpi = 300)  

# Paired Sample Boxplot ---------------------------------------------------------  
SPHK1_pair_tpm <- paired_data %>% select(Sample_ID, Cancer, Group, ID, SPHK1) %>% na.omit(.)  
SPHK1_pair_tpm[1:5, ]  

# Concise Paired Boxplot  
# Direct plotting (no need to define a function)  
library(ggpubr)  

# Define cancer types  
pcancers <- c(  
  "BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH",   
  "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD",   
  "STAD", "THCA", "UCEC")  

# Filter data  
df <- SPHK1_pair_tpm %>%  
  filter(Cancer %in% pcancers) %>%  
  mutate(  
    Cancer = factor(Cancer, levels = pcancers),  
    Group = factor(Group, levels = c("Normal", "Tumor"))  
  )  

# Draw the plot  
p <- ggpaired(df,   
              x = "Group",   
              y = "SPHK1",   
              id = "ID",   
              color = "Group",   
              palette = c("#4DBBD5", "#E64B35"),   
              add = "jitter",   
              xlab = "",   
              ylab = "SPHK1",   
              line.color = "gray",   
              line.size = 0.4,   
              facet.by = "Cancer",   
              nrow = 1) +  
  theme_classic() +  
  rotate_x_text(angle = 90) +  
  theme(  
    legend.title = element_text(size = 14),   
    legend.text = element_text(size = 12),   
    axis.title.x = element_text(size = 14),   
    axis.text.x = element_text(size = 12, colour = "black"),   
    axis.title.y = element_text(size = 14),   
    axis.text.y = element_text(size = 12, colour = "black"),   
    strip.text.x = element_text(size = 12, color = "black")  
  )  

# Add statistical test  
p <- p + stat_compare_means(  
  label = "p.signif",   
  method = "wilcox.test",   
  paired = TRUE,   
  label.x.npc = 0.4,   
  label.y.npc = "top"  
)  

print(p)  
ggsave("SPHK1_paired_plot.png", p, width = 20, height = 6, dpi = 300)  

# Draw a forest plot --------------------------------  
library(survival)  
library(forestplot)  
library(dplyr)  
library(tibble)  

# Prepare data  
meta <- meta_dedup %>% select(patient_id, Cancer, event, time, age, gender) %>% rename(ID = patient_id)  
cancers <- unique(SPHK1_tpm$Cancer)  
cox_results <- list()  

for (cancer in cancers) {  
  exprSet <- subset(SPHK1_tpm, Group == "Tumor" & Cancer == cancer) %>%  
    tibble::add_column(ID = substr(.$sample, 1, 12), .before = "Cancer") %>%  
    dplyr::filter(!duplicated(ID)) %>%  
    tibble::column_to_rownames("ID") %>%  
    dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == cancer)))  
  
  cl <- meta[rownames(exprSet), ]  
  cl$symbol <- exprSet[["SPHK1"]]  
  
  tryCatch({  
    m <- coxph(Surv(time, event) ~ symbol, data = cl)  
    beta <- coef(m)  
    se <- sqrt(diag(vcov(m)))  
    
    tmp <- round(cbind(  
      HR = exp(beta),  
      HRCILL = exp(beta - qnorm(0.975, 0, 1) * se),  
      HRCIUL = exp(beta + qnorm(0.975, 0, 1) * se),  
      p = 1 - pchisq((beta/se)^2, 1)  
    ), 3)  
    
    cox_results[[cancer]] <- tmp  
  }, error = function(e) NULL)  
}  

name <- names(cox_results)  
cox_results <- do.call(rbind, cox_results)  
rownames(cox_results) <- name   
# Prepare forest plot data  
np <- paste0(  
  sprintf("%.2f", cox_results[,1]), " (",  
  sprintf("%.2f", cox_results[,2]), "-",  
  sprintf("%.2f", cox_results[,3]), ")"  
)  

tabletext <- cbind(  
  c("Cancer", rownames(cox_results)),  
  c("HR (95%CI)", np),  
  c("P Value", sprintf("%.3f", cox_results[,4]))  
)  

# Draw the forest plot  
forestplot(  
  labeltext = tabletext,  
  graph.pos = 3,  
  mean = c(NA, cox_results[,1]),  
  lower = c(NA, cox_results[,2]),  
  upper = c(NA, cox_results[,3]),  
  title = "Hazard Ratio Plot of SPHK1",  
  hrzl_lines = list(  
    `1` = gpar(lwd = 2, col = "black"),  
    `2` = gpar(lwd = 1, col = "black")  
  ),  
  is.summary = c(TRUE, rep(FALSE, nrow(cox_results))),  
  col = fpColors(box = "#1c61b6", lines = "#1c61b6", zero = "gray50"),  
  zero = 1,  
  boxsize = 0.5  
)  

# Batch Draw ROC Curves for Multiple Cancers ---------------------------------------------  
cancers_to_plot <- c("BLCA", "BRCA", "LUAD", "LUSC", "COAD")  

plots_list <- list()  

for (cancer in cancers_to_plot) {  
  # Filter data  
  df <- SPHK1_tpm %>%  
    filter(Cancer == cancer) %>%  
    select(Group, SPHK1) %>%  
    mutate(Group = factor(Group, levels = c("Normal", "Tumor")))  
  
  # Check if both normal and tumor samples exist  
  if (n_distinct(df$Group) == 2) {  
    # Calculate ROC  
    res <- roc(Group ~ SPHK1, data = df,   
               levels = c("Normal", "Tumor"),  
               direction = "<",  
               ci = TRUE)  
    
    # Create the plot  
    p <- ggroc(res, color = "#2E9FDF", legacy.axes = TRUE) +  
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),   
                   color = "darkgrey", linetype = 4) +  
      theme_bw() +  
      ggtitle(paste(cancer, "ROC Curve")) +  
      theme(  
        plot.title = element_text(hjust = 0.5, size = 14),  
        axis.title = element_text(size = 12),  
        axis.text = element_text(size = 10)  
      )  
    
    # Add AUC information  
    auc_value <- auc(res)  
    p <- p + annotate("text", x = 0.7, y = 0.3,   
                      label = paste("AUC =", round(auc_value, 3)),  
                      size = 4, fontface = "bold")  
    
    plots_list[[cancer]] <- p  
  }  
}  

# Combine all plots together  
library(patchwork)  
combined_plot <- wrap_plots(plots_list, ncol = 2)  
print(combined_plot)  

# Save the combined plot  
ggsave("SPHK1_multiple_cancers_ROC.png",   
       combined_plot, width = 12, height = 10, dpi = 300)  

# Batch Draw KM Curves for Multiple Cancer Types  
library(patchwork)  

genes_to_plot <- c("SPHK1")  
cancers_to_plot <- c("BLCA", "BRCA", "LUAD", "COAD")  

plots_list <- list()  

for (cancer in cancers_to_plot) {  
  for (gene in genes_to_plot) {  
    tryCatch({  
      # Prepare data  
      expr_data <- SPHK1_tpm %>%  
        filter(Cancer == cancer & Group == "Tumor") %>%  
        mutate(ID = substr(sample, 1, 12)) %>%  
        distinct(ID, .keep_all = TRUE) %>%  
        select(ID, all_of(gene)) %>%  
        filter(ID %in% meta$ID)  
      
      cl <- meta %>%  
        filter(ID %in% expr_data$ID) %>%  
        mutate(  
          expression = ifelse(  
            expr_data[[gene]][match(ID, expr_data$ID)] >   
              median(expr_data[[gene]], na.rm = TRUE),   
            "high", "low"  
          ),  
          expression = factor(expression, levels = c("low", "high"))  
        )  
      
      # Fit survival curve  
      surv_obj <- Surv(time = cl$time, event = cl$event)  
      fit <- survfit(surv_obj ~ expression, data = cl)  
      
      # Create simplified plot (space-saving)  
      p <- ggsurvplot(  
        fit,  
        data = cl,  
        pval = TRUE,  
        pval.coord = c(0.1, 0.1),  
        palette = c("#4DBBD5", "#E64B35"),  
        title = paste(cancer, "-", gene),  
        risk.table = FALSE,  
        xlab = "Time (years)",  
        ylab = "Survival",  
        legend = "top",  
        ggtheme = theme_minimal(),  
        size = 1  
      )$plot  
      
      plots_list[[paste(cancer, gene, sep = "_")]] <- p  
      
    }, error = function(e) {  
      message(paste("Error in", cancer, gene, ":", e$message))  
    })  
  }  
}  

# Combine all plots  
combined_plot <- wrap_plots(plots_list, ncol = 2)  
print(combined_plot)  

# Save the combined plot  
ggsave("SPHK1_multiple_cancers_KM.png",   
       combined_plot, width = 14, height = 10, dpi = 300)