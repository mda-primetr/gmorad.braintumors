# Use gtsummary package or similar package to summarize the metadata
source("src/libraries.R")
source("src/common_functions.R")

load("data/processed_data/physeq_samples.RData")
load("data/processed_data/physeq_scrubbed_plus_nc_rem.RData")



df_metadata <- sample_data(physeq_samples) %>%
    data.frame() %>%
    mutate(SampleID = rownames(.))

df_metadata_filtered <- sample_data(physeq_scrubbed_plus_nc_rem) %>%
    data.frame() %>%
    mutate(SampleID = rownames(.))

df_metadata_filtered_glioma <- sample_data(physeq_scrubbed_plus_nc_rem) %>%
    data.frame() %>%
    filter(tumor_category == "Glioma") %>%
    mutate(SampleID = rownames(.))


df_metadata_filtered_METS <- sample_data(physeq_scrubbed_plus_nc_rem) %>%
    data.frame() %>%
    filter(tumor_category == "Met") %>%
    mutate(SampleID = rownames(.))

# ────────────────────────────────────────────────────────────────────────────────────────────────────
#    Glioma ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────




# Stats for Glioma by progression for all available samples and MDACC1 cohort -----
# We have not taken normal brain tissue into account for this step
df_metadata %>%
    filter(microbiome_sample_type == "Tumor") %>%
    filter(tumor_category == "Glioma") %>%
    mutate(glioma_progression_overall = case_when(
        glioma_progression_overall == "Unknown" ~ NA,
        TRUE ~ glioma_progression_overall
    )) %>%
    mutate(have_bacterial_reads = ifelse(sample_id %in% df_metadata_filtered_glioma$Sample, "Yes", "No")) %>%
    mutate(glioma_fu_time = as.numeric(glioma_fu_time)) %>%
    dplyr::select(
        have_bacterial_reads,
        glioma_progression_overall,
        sex,
        age_at_surgery,
        race,
        steroids_at_surgery,
        grade_category,
        x16s_seq_batch,
        glioma_idh_mutation,
        glioma_fu_time,
        have_bacterial_reads,
        glioma_recurrent_status_at_surg
    ) %>%
    gtsummary::tbl_summary(by = have_bacterial_reads) %>%
    add_p() %>%
    gtsummary::modify_caption("**Glioma : Have detectable bacteria**") %>%
    as_gt() %>%
    gt::gtsave(., "output/tables/Glioma_have_bacterial_reads_ALL_16S.docx")



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# METS
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_metadata %>%
    filter(microbiome_sample_type == "Tumor") %>%
    filter(tumor_category == "Met") %>%
    mutate(br_m_progression_overall = case_when(
        br_m_progression_overall == "Unknown" ~ NA,
        TRUE ~ br_m_progression_overall
    )) %>%
    mutate(ecm_at_last_followup = stringr::str_to_title(ecm_at_last_followup)) %>%
    mutate(ecm_at_last_followup = case_when(
        ecm_at_last_followup == "Unknown" ~ NA,
        TRUE ~ ecm_at_last_followup
    )) %>%
    mutate(have_bacterial_reads = ifelse(sample_id %in% df_metadata_filtered_METS$Sample, "Yes", "No")) %>%
    mutate(br_m_total_follow_up_days = as.numeric(br_m_total_follow_up_days)) %>%
    dplyr::select(
        have_bacterial_reads,
        br_m_progression_overall,
        br_m_primary_category,
        sex,
        age_at_surgery,
        race,
        steroids_at_surgery,
        x16s_seq_batch,
        br_m_prior_treatment,
        prev_xrt,
        prev_it,
        prev_targeted,
        lmd_at_last_followup,
        ecm_at_last_followup,
        br_m_total_follow_up_days
    ) %>%
    gtsummary::tbl_summary(by = have_bacterial_reads) %>%
    add_p() %>%
    gtsummary::modify_caption("**Mets : Have detectable bacteria**") %>%
    as_gt() %>%
    gt::gtsave(., "output/tables/Mets_have_bacterial_reads_ALL_16S.docx")