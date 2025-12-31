# Extract paired sample data from tpm_all ===
# 1. Check data frame structure
# 2. Extract Patient ID (first 12 characters from Sample_ID)
tpm_all <- final_data
tpm_all$Patient_ID <- substr(tpm_all$Sample_ID, 1, 12)

# 3. Identify paired samples (patients with both tumor and normal samples)
# Note: Based on the example, sample type is judged by positions 14-15 of Sample_ID
paired_patients <- tpm_all %>%
  group_by(Patient_ID, Cancer) %>%
  summarise(
    has_tumor = any(Group == "Tumor"),
    has_normal = any(Group == "Normal"),
    .groups = 'drop'
  ) %>%
  filter(has_tumor & has_normal) %>%
  pull(Patient_ID)

# 4. Extract paired sample data
paired_data <- tpm_all %>%
  filter(Patient_ID %in% paired_patients) %>%
  mutate(
    # Create Group column: Tumor/Normal
    Group = ifelse(Group == "Tumor", "Tumor", "Normal"),
    # Create ID column: Patient ID
    ID = Patient_ID,
    # Create Cancer column: Cancer type
    Cancer = Cancer
  ) %>%
  # Select and rename needed columns
  dplyr::select(Sample_ID, Cancer, Group, ID, everything())

fwrite(paired_data, "TCGA_PanCancer_Analysis/02_results/tcga_paired_tpm.csv")

save(tpm_all, paired_data, file = "TCGA_PanCancer_Analysis/02_results/tcga_paired_tpm.Rdata")

save(final_data, file = "TCGA_PanCancer_Analysis/02_results/final_data.Rdata")