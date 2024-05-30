source("src/libraries.R")
source("src/common_functions.R")

# Load Unrarefied data
load("data/processed_data/physeq_SSO_clean.RData")
physeq_SSO_clean



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Glioma -----
# ────────────────────────────────────────────────────────────────────────────────────────────────────


# Stats for Glioma by IDH mutation -----
sample_data(physeq_SSO_clean) %>%
    data.frame() %>%
    filter(microbiome_sample_type == "Stool") %>%
    filter(tumor_category == "Glioma") %>%
    filter(glioma_idh_mutation != "Unknown") %>%
    dplyr::select(
        tumor_category,
        sex,
        age_at_surgery,
        race,
        steroids_at_surgery,
        grade_category,
        shotgun_seq_batch,
        glioma_idh_mutation
    ) %>%
    gtsummary::tbl_summary(by = glioma_idh_mutation) %>%
    gtsummary::modify_caption("**Glioma IDH Mutation**") %>%
    gtsummary::add_p() %>%
    as_gt() %>%
    gt::gtsave(., "output/tables/Glioma_by_IDH_mutation_SSO_WGS.docx")


# Stats for Glioma for Progression category ----
sample_data(physeq_SSO_clean) %>%
    data.frame() %>%
    filter(microbiome_sample_type == "Stool") %>%
    filter(tumor_category == "Glioma") %>%
    filter(glioma_idh_mutation != "Unknown") %>%
    filter(grade_category != "Unknown") %>%
    mutate(glioma_fu_time = as.numeric(glioma_fu_time)) %>%
    filter(glioma_progression_overall != "Unknown") %>%
    dplyr::select(
        glioma_progression_overall,
        sex,
        age_at_surgery,
        race,
        steroids_at_surgery,
        grade_category,
        shotgun_seq_batch,
        glioma_idh_mutation,
        glioma_fu_time,
        glioma_recurrent_status_at_surg
    ) %>%
    gtsummary::tbl_summary(by = glioma_progression_overall) %>%
    gtsummary::modify_caption("**Glioma Progression overall**") %>%
    gtsummary::add_p() %>%
    as_gt() %>%
    gt::gtsave(., "output/tables/Glioma_glioma_progression_overall_SSO_WGS.docx")



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Brm/Met -----
# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Stats for BM by Progression overall (Brain Mets) -----
sample_data(physeq_SSO_clean) %>%
    data.frame() %>%
    filter(microbiome_sample_type == "Stool") %>%
    filter(tumor_category == "Met") %>%
    filter(br_m_progression_overall != "Unknown") %>%
    mutate(BrM_total_follow_up.days. = as.numeric(br_m_total_follow_up_days)) %>%
    dplyr::select(
        BrM_progression_overall = br_m_progression_overall,
        BrM_primary_category = br_m_primary_category,
        sex,
        age_at_surgery,
        race,
        steroids_at_surgery,
        shotgun_seq_batch,
        BrM_prior_treatment = br_m_prior_treatment,
        prev_xrt,
        prev_it,
        prev_targeted,
        LMD_at_last_followup = lmd_at_last_followup,
        ECM_at_last_followup = ecm_at_last_followup,
        BrM_total_follow_up.days.
    ) %>%
    gtsummary::tbl_summary(by = BrM_progression_overall) %>%
    gtsummary::modify_caption("**Brain Mets Overall Progression **") %>%
    gtsummary::add_p() %>%
    as_gt() %>%
    gt::gtsave(., "output/tables/Brain_Mets_by_overall_progression_SSO_WGS.docx")




# Stats for BM  by BrM_prior_treatment (Brain Mets) -----
sample_data(physeq_SSO_clean) %>%
    data.frame() %>%
    filter(microbiome_sample_type == "Stool") %>%
    filter(tumor_category == "Met") %>%
    mutate(BrM_total_follow_up.days. = as.numeric(br_m_total_follow_up_days)) %>%
    dplyr::select(
        BrM_primary_category = br_m_primary_category,
        sex,
        age_at_surgery,
        race,
        steroids_at_surgery,
        shotgun_seq_batch,
        BrM_prior_treatment = br_m_prior_treatment,
    ) %>%
    gtsummary::tbl_summary(by = BrM_prior_treatment) %>%
    gtsummary::modify_caption("**Brain Mets Prior Treatment**") %>%
    gtsummary::add_p() %>%
    as_gt() %>%
    gt::gtsave(., "output/tables/Brain_Mets_by_prior_treatment_SSO_WGS.docx")
