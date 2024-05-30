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




















# Stats for Glioma by progression for all available samples and MDACC1 cohort -----
# We have not taken normal brain tissue into account for this step
read_csv("data/processed_data/metadata_cleaned.csv") %>%
    filter(microbiome_sample_type == "Tumor") %>%
    filter(tumor_category == "Glioma") %>%
    filter(dataset == "MDACC1") %>%
    filter(glioma_progression_overall != "Unknown") %>%
    mutate(have_bacterial_reads = ifelse(sample_id %in% df_metadata_cleaned$Sample, "Yes", "No")) %>%
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
    gtsummary::modify_caption("**Glioma (MDACC1) : Have detectable bacteria**") %>%
    as_gt() %>%
    gt::gtsave(., "output/figures/Glioma_have_bacterial_reads_MDACC1.docx")



# Stats for Glioma by progression for all available samples and all cohorts -----
# We have not taken normal brain tissue into account for this step
read_csv("data/processed_data/metadata_cleaned.csv") %>%
    filter(microbiome_sample_type == "Tumor") %>%
    filter(tumor_category == "Glioma") %>%
    filter(glioma_progression_overall != "Unknown") %>%
    mutate(have_bacterial_reads = ifelse(sample_id %in% df_metadata_cleaned$Sample, "Yes", "No")) %>%
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
    gtsummary::modify_caption("**Glioma All : Have detectable bacteria**") %>%
    as_gt() %>%
    gt::gtsave(., "output/figures/Glioma_have_bacterial_reads_all_cohorts.docx")


# Calculate the above but we are focused on the steroids at surgery so we will remove unknowns -----
read_csv("data/processed_data/metadata_cleaned.csv") %>%
    filter(microbiome_sample_type == "Tumor") %>%
    filter(tumor_category == "Glioma") %>%
    filter(steroids_at_surgery != "Unknown") %>%
    mutate(have_bacterial_reads = ifelse(sample_id %in% df_metadata_cleaned$Sample, "Yes", "No")) %>%
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
    gtsummary::modify_caption("**Glioma All (Steroids) : Have detectable bacteria**") %>%
    as_gt() %>%
    gt::gtsave(., "output/figures/Glioma_have_bacterial_reads_all_cohorts_STEROIDS.html")


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Glioma ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Test for progression in days and bacterial loads per gram of tissue ----

df_progression_reads <- physeq_all_filt %>%
    speedyseq::tax_glom(taxrank = "Genus", .) %>%
    psmelt() %>%
    filter(tumor_category == "Glioma") %>%
    filter(dataset == "MDACC1") %>%
    filter(!glioma_time_to_progression %in% c("Unknown", "NA")) %>%
    mutate(glioma_time_to_progression = as.numeric(glioma_time_to_progression)) %>%
    filter(!is.na(tumor_weight)) %>%
    group_by(Sample) %>%
    mutate(total_reads = sum(Abundance)) %>%
    mutate(tumor_weight = as.numeric(tumor_weight)) %>%
    ungroup() %>%
    distinct(Sample, .keep_all = TRUE) %>%
    rowwise() %>%
    mutate(total_reads_per_g = total_reads / tumor_weight) %>%
    dplyr::select(total_reads_per_g, Sample, glioma_time_to_progression)

df_cor_stat_glioma_progression <- cor.test(df_progression_reads$total_reads_per_g, df_progression_reads$glioma_time_to_progression,
    method = "spearman"
) %>%
    broom::tidy() %>%
    mutate(p.value = format.pval(p.value, eps = 0.01, digits = 2))

ggplot(df_progression_reads, aes(x = glioma_time_to_progression, y = total_reads_per_g)) +
    geom_point(size = 4) +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(title = "Glioma Progression vs Bacterial load per gram of tissue") +
    geom_text(
        data = df_cor_stat_glioma_progression, aes(
            x = Inf,
            y = Inf,
            label = paste0("Spearman's rho = ", round(estimate, 2))
        ),
        hjust = 2, vjust = 2, size = 5, color = "black"
    ) +
    geom_text(
        data = df_cor_stat_glioma_progression, aes(
            x = Inf,
            y = Inf,
            label = paste0(" p-value = ", p.value)
        ),
        hjust = 2.9, vjust = 3.5, size = 5, color = "black"
    )
ggsave("output/figures/Glioma_progression_vs_bacterial_load_per_g_tissue.pdf", width = 10, height = 10)




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# METS
# ────────────────────────────────────────────────────────────────────────────────────────────────────


# Stats for Mets by progression for all available samples -----
# We have not taken normal brain tissue into account for this step
read_csv("data/processed_data/metadata_cleaned.csv") %>%
    filter(microbiome_sample_type == "Tumor") %>%
    filter(tumor_category == "Met") %>%
    mutate(br_m_progression_overall = case_when(
        br_m_progression_overall == "Unknown" ~ NA,
        TRUE ~ br_m_progression_overall
    )) %>%
    mutate(have_bacterial_reads = ifelse(sample_id %in% df_metadata_cleaned$Sample, "Yes", "No")) %>%
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
    gt::gtsave(., "output/figures/Mets_have_bacterial_reads.docx")
