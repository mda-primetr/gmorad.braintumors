source("src/libraries.R")
source("src/common_functions.R")

# Set seed for reproducibility ----
RNGkind(sample.kind = "Rounding")
set.seed(711)

# Load metadata file ----
# Read metadata ----
m1 <- read_csv("data/metadata/megametadata_combined_16S_WGS_2025_05_01.csv") %>%
    janitor::clean_names()


# Load 16s data after all filters ----
load("data/processed_data/16S/physeq_scrubbed_plus_nc_rem.RData")

# Sample data for 16S tumor samples after all filters ----
df_metadata_post_filtration_16S <- physeq_scrubbed_plus_nc_rem %>%
    tax_glom(taxrank = "Genus") %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    group_by(coded_id, Sample, tumor_category) %>%
    summarize(total_abundance = sum(Abundance)) %>%
    mutate(sequencing_technology = "16S") %>%
    filter(tumor_category %in% c("Met", "Glioma")) %>%
    ungroup() %>%
    dplyr::select(coded_id, Sample, tumor_category, total_abundance)


# Metadata for 16S tumor samples ----
m1_16S_tumor <- m1 %>%
    filter(microbiome_sample_type %in% c("Tumor")) %>%
    filter(sequencing_technology == "16S")



# Load WGS tumor data ---
# WGS post filter ----
load(file = "data/processed_data/WGS/physeq_bacmap_tumor_post_filter.Rdata")
physeq_bacmap_tumor_post_filter

# Get taxa that remained after all WGS tumor filters ----
list_of_wgs_remaining_taxa <- read_csv(
    "output/tables/WGS/Tumor_samples_wgs_after_all_filter_including_average_letters_filter.csv"
)

df_metadata_post_filtration_WGS <- physeq_bacmap_tumor_post_filter %>%
    tax_glom(taxrank = "Species") %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    filter(Species %in% list_of_wgs_remaining_taxa$Species) %>%
    group_by(coded_id, Sample, tumor_category) %>%
    summarize(total_abundance = sum(Abundance)) %>%
    filter(tumor_category %in% c("Met", "Glioma")) %>%
    ungroup() %>%
    dplyr::select(coded_id, Sample, tumor_category, total_abundance)



# Metadata for WGS tumor samples ----
m1_WGS_tumor <- m1 %>%
    filter(microbiome_sample_type %in% c("Tumor")) %>%
    filter(sequencing_technology == "WGS")


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Combine 16S and WGS data ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_combined_16S_WGS <- rbind(
    m1_WGS_tumor %>%
        left_join(df_metadata_post_filtration_WGS, by = c("coded_id", "sequencing_id" = "Sample", "tumor_category")),
    m1_16S_tumor %>%
        left_join(df_metadata_post_filtration_16S, by = c("coded_id", "sequencing_id" = "Sample", "tumor_category"))
)

df_list_coded_ids_with_bacterial_signal <- df_combined_16S_WGS %>%
    group_by(coded_id) %>%
    summarize(sum_abundance = sum(total_abundance, na.rm = T)) %>%
    filter(sum_abundance > 0)


df_combined_16S_WGS_unique <- df_combined_16S_WGS %>%
    mutate(bacterial_signal = case_when(
        coded_id %in% df_list_coded_ids_with_bacterial_signal$coded_id ~ "Yes",
        TRUE ~ "No"
    )) %>%
    distinct(coded_id, .keep_all = TRUE)


# ────────────────────────────────────────────────────────────────────────────────────────────────────
#  Compare bacterial signature with Multi lesion
# ────────────────────────────────────────────────────────────────────────────────────────────────────
df_combined_16S_WGS_unique %>%
    mutate(age_at_surgery = as.numeric(age_at_surgery)) %>%
    mutate(
        steroids_at_surgery = factor(
            steroids_at_surgery,
            levels = c("N", "Y")
        )
    ) %>%
    mutate(infection_days_before_surgery = as.numeric(infection_days_before_surgery)) %>%
    dplyr::select(
        age_at_surgery,
        sex,
        race,
        steroids_at_surgery,
        history_bacterial_infection,
        infection_days_before_surgery,
        antibiotic,
        multi_focal_lesion,
        resected_tumor_location,
        bacterial_signal
    ) %>%
    mutate(multi_focal_lesion = factor(
        multi_focal_lesion,
        levels = c("N", "Y")
    )) %>%
    gtsummary::tbl_summary(by = bacterial_signal) %>%
    gtsummary::modify_spanning_header(
        all_stat_cols() ~ paste0("**By Bacterial signature** : ")
    ) %>%
    gtsummary::add_p(
        test = list(
            all_categorical() ~ "chisq.test",
            all_continuous() ~ "wilcox.test"
        ),
        test.args = all_tests("chisq.test") ~ list(simulate.p.value = TRUE)
    ) %>%
    gtsummary::separate_p_footnotes() %>%
    as_gt() %>%
    gt::gtsave(
        .,
        paste0(
            "output/tables/BrM_Glioma_16S_WGS_combined_compare_bacteria_signal",
            ".docx"
        )
    )


# ────────────────────────────────────────────────────────────────────────────────────────────────────
#  Compare bacterial signature with No Multi lesion
# ────────────────────────────────────────────────────────────────────────────────────────────────────
df_combined_16S_WGS_unique %>%
    mutate(age_at_surgery = as.numeric(age_at_surgery)) %>%
    mutate(
        steroids_at_surgery = factor(
            steroids_at_surgery,
            levels = c("N", "Y")
        )
    ) %>%
    mutate(infection_days_before_surgery = as.numeric(infection_days_before_surgery)) %>%
    dplyr::select(
        age_at_surgery,
        sex,
        race,
        steroids_at_surgery,
        tumor_category,
        history_bacterial_infection,
        infection_days_before_surgery,
        antibiotic,
        resected_tumor_location,
        multi_focal_lesion,
        bacterial_signal
    ) %>%
    filter(multi_focal_lesion == "N") %>%
    gtsummary::tbl_summary(by = bacterial_signal) %>%
    gtsummary::modify_spanning_header(
        all_stat_cols() ~ paste0("**By Bacterial signature** : ")
    ) %>%
    gtsummary::add_p(
        test = list(
            all_categorical() ~ "chisq.test",
            all_continuous() ~ "wilcox.test"
        ),
        test.args = all_tests("chisq.test") ~ list(simulate.p.value = TRUE)
    ) %>%
    gtsummary::separate_p_footnotes() %>%
    as_gt() %>%
    gt::gtsave(
        .,
        paste0(
            "output/tables/BrM_Glioma_16S_WGS_combined_compare_bacteria_signal_with_no_multifocal",
            ".docx"
        )
    )
