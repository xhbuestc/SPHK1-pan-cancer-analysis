
library(TCGAPlot)


gene_MSI_radar(gene)
gene_TMB_radar(gene)
gene_immucell_heatmap(gene) # immucell
gene_checkpoint_heatmap(gene) # checkpoint
gene_chemokine_heatmap(gene)  # chemokine
gene_receptor_heatmap(gene, method = "spearman") # receptor
gene_gsea_kegg(cancertype, gene)  # GSEA  