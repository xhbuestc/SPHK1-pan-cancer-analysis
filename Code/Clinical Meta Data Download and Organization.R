meta <-  fread("TCGA_PanCancer_Analysis/00_rawdata/Survival_SupplementalTable_S1_20171025_xena_sp.gz") %>% 
  select(sample, `cancer type abbreviation`, age_at_initial_pathologic_diagnosis, gender, ajcc_pathologic_tumor_stage, OS:PFI.time) %>%
  # 1. Rename core columns
  rename(Cancer = `cancer type abbreviation`,
         event = OS,
         time_raw = OS.time,
         age_raw = `age_at_initial_pathologic_diagnosis`,
         stage_raw = `ajcc_pathologic_tumor_stage`) %>%
  # 2. Extract patient ID from 'sample' column as new column and row name
  mutate(patient_id = sub("(TCGA-[^-]+-[^-]+).*", "\\1", sample)) %>%
  # 3. Clean and format other columns
  mutate(
    # Survival time: Convert from days to years, keep 2 decimal places
    time = round(time_raw / 365.25, 2),
    # Age: Usually in years, can be used directly, or convert as needed
    age = as.integer(age_raw),
    # Gender: Simplify to M/F
    gender = substr(gender, 1, 1),
    # Stage: Remove "Stage " prefix
    stage = gsub("Stage ", "", stage_raw),
    DSS.time = round(DSS.time / 365.25, 2),
    DFI.time = round(DFI.time / 365.25, 2),
    PFI.time = round(PFI.time / 365.25, 2)
  )  %>%
  mutate(
    # Step 1: Standardize Roman numeral format, remove extra spaces (e.g., "II " -> "II")
    stage_trimmed = trimws(stage),
    # Step 2: Core recoding logic
    clean_stage = case_when(
      # Keep clear standard Roman numeral stages I, II, III, IV
      stage_trimmed %in% c("I", "II", "III", "IV") ~ stage_trimmed,
      # Merge substages (e.g., IIA, IIIB) into main stage
      grepl("^I+[ABC]*$", stage_trimmed) & nchar(stage_trimmed) <= 4 ~ gsub("[ABC]", "", stage_trimmed),
      # Treat explicit numeric 0 as "0" (or separate category)
      stage_trimmed == "0" ~ "0",
      # Mark all other cases (including [Discrepancy], [Unknown], X, IS, NOS, etc.) as Unknown
      TRUE ~ "Unknown"
    )
  ) %>%
  # Remove temporary column
  select(-stage_trimmed) %>%
  mutate(
    # First use strategy A or B to obtain clear clean_stage
    # ... (Run code from strategy A or B) ...
    # Then convert clear stage to numeric
    stage_numeric = case_when(
      clean_stage == "0"   ~ 0,
      clean_stage == "I"   ~ 1,
      clean_stage == "II"  ~ 2,
      clean_stage == "III" ~ 3,
      clean_stage == "IV"  ~ 4,
      TRUE ~ NA_real_ # Set "Unknown", etc. to NA
    )
  ) %>%
  # 4. Filter and arrange final needed columns
  select(patient_id, sample, Cancer, event, time, age, gender, stage_numeric, DSS:PFI.time) %>%
  # 5. Set row names
  as.data.frame()

# Deduplicate by patient_id, keep the first record for each patient
meta_dedup <- meta[!duplicated(meta$sample), ]
# Now it's safe to set row names
rownames(meta_dedup) <- meta$sample
cat(sprintf("After deduplication: %d unique patients retained (removed %d duplicate records).\n", 
            nrow(meta_dedup), nrow(meta)-nrow(meta_dedup)))

fwrite(meta_dedup, "TCGA_PanCancer_Analysis/02_results/meta_dedup.csv")
meta <- meta_dedup
save(meta_dedup, file = "TCGA_PanCancer_Analysis/02_results/meta_dedup.Rdata")