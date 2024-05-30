source("src/libraries.R")
source("src/common_functions.R")


load("data/processed_data/physeq_SSO_clean.RData")
load("data/processed_data/physeq_WGS_HMP.RData")

ps_sso_clean_melt <- psmelt(physeq_SSO_clean) %>%
    filter(Abundance > 0) %>%
    dplyr::select(Sample, microbiome_sample_type, Abundance, OTU, Kingdom:Species) %>%
    mutate(cohort = "MDACC", .after = "microbiome_sample_type")




# Read the HMP metadata for Oral swab ----
m1_hmp_oral_swab <- read_tsv("data/metadata/HMP_SraRunTable.txt") %>%
    filter(isolation_source %in% c(
        "G_DNA_Buccal mucosa"
    )) %>%
    mutate(sample_id = Run) %>%
    distinct(submitted_subject_id, .keep_all = TRUE) %>% # So we can only select one sample from each subject
    column_to_rownames(var = "Run") %>%
    write_csv("data/processed_data/metadata_HMP_oral_swab.csv")


ps_melt_HMP_oral <- physeq_WGS_HMP %>%
    speedyseq::filter_sample_data(isolation_source %in% c(
        "G_DNA_Buccal mucosa"
    )) %>%
    speedyseq::filter_sample_data(sample_id %in% rownames(m1_hmp_oral_swab)) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    mutate(microbiome_sample_type = "Oral swab") %>%
    dplyr::select(Sample, microbiome_sample_type, Abundance, OTU, Kingdom:Species) %>%
    mutate(cohort = "HMP", .after = "microbiome_sample_type")

# Read the HMP metadata for Saliva ----
m1_hmp_saliva <- read_tsv("data/metadata/HMP_SraRunTable.txt") %>%
    filter(isolation_source %in% c(
        "G_DNA_Saliva"
    )) %>%
    mutate(sample_id = Run) %>%
    distinct(submitted_subject_id, .keep_all = TRUE) %>% # So we can only select one sample from each subject
    column_to_rownames(var = "Run") %>%
    write_csv("data/processed_data/metadata_HMP_saliva.csv")



ps_melt_HMP_saliva <- physeq_WGS_HMP %>%
    speedyseq::filter_sample_data(isolation_source %in% c(
        "G_DNA_Saliva"
    )) %>%
    speedyseq::filter_sample_data(sample_id %in% rownames(m1_hmp_saliva)) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    mutate(microbiome_sample_type = "Saliva") %>%
    dplyr::select(Sample, microbiome_sample_type, Abundance, OTU, Kingdom:Species) %>%
    mutate(cohort = "HMP", .after = "microbiome_sample_type")




# Read the HMP metadata for stool ----
m1_hmp_stool <- read_tsv("data/metadata/HMP_SraRunTable.txt") %>%
    filter(isolation_source %in% c(
        "G_DNA_Stool"
    )) %>%
    mutate(sample_id = Run) %>%
    distinct(submitted_subject_id, .keep_all = TRUE) %>% # So we can only select one sample from each subject
    column_to_rownames(var = "Run") %>%
    write_csv("data/processed_data/metadata_HMP_stool.csv")


ps_melt_HMP_stool <- physeq_WGS_HMP %>%
    speedyseq::filter_sample_data(isolation_source %in% c(
        "G_DNA_Stool"
    )) %>%
    speedyseq::filter_sample_data(sample_id %in% rownames(m1_hmp_stool)) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    mutate(microbiome_sample_type = "Stool") %>%
    dplyr::select(Sample, microbiome_sample_type, Abundance, OTU, Kingdom:Species) %>%
    mutate(cohort = "HMP", .after = "microbiome_sample_type")



# combine all samples from HMP and SSO ----
physeq_HMP_SSO_combo_species <- rbind(
    ps_sso_clean_melt,
    ps_melt_HMP_oral,
    ps_melt_HMP_saliva,
    ps_melt_HMP_stool
) %>%
    dplyr::rename("sample_id" = "Sample") %>%
    melt_to_physeq() %>%
    speedyseq::tax_glom(taxrank = "Species")



# Differential abundance between Glioma and Met combo vs HMP -----

# No need to control
for (i in c("Stool", "Saliva", "Oral swab")) {
    final_phyloseq_before_rarefaction_species_i <- physeq_HMP_SSO_combo_species %>%
        speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
        speedyseq::mutate_sample_data(cohort = factor(cohort, levels = c("HMP", "MDACC"))) %>%
        speedyseq::filter_sample_data(sample_sums(.) >= 100)


    if (i == "Stool") {
        final_phyloseq_before_rarefaction_species_i <- physeq_HMP_SSO_combo_species %>%
            speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
            speedyseq::mutate_sample_data(cohort = factor(cohort, levels = c("HMP", "MDACC"))) %>%
            speedyseq::filter_sample_data(sample_sums(.) >= 1000)
    }

    df_maaslin_all <- data.frame(otu_table(final_phyloseq_before_rarefaction_species_i) %>% t())

    df_maaslin_all_meta <- data.frame(sample_data(final_phyloseq_before_rarefaction_species_i))


    # Run Maaslin2
    maaslin_by_sample <- Maaslin2(
        input_data = df_maaslin_all,
        input_metadata = df_maaslin_all_meta,
        output = paste0("data/processed_data/maaslin_by_cohort_HMP_Vs_MDACC_all_", i),
        analysis_method = "NEGBIN",
        fixed_effects = c("cohort"),
        plot_scatter = FALSE,
        plot_heatmap = FALSE,
        min_abundance = 100,
        min_prevalence = 0.25,
        cores = 24,
        transform = "NONE",
        normalization = "TMM",
        reference = c("cohort,HMP"),
        standardize = TRUE
    )


    # Run AncomBC version 1
    ancombc_out <- ancombc(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0.25,
        formula = "cohort",
        struc_zero = FALSE,
        neg_lb = FALSE,
        n_cl = 24
    )

    # Get AncomBC1 results as a dataframe
    get_ancombc_df(ancombc_out$res, taxa_level = "Species", group_var = "cohort") %>%
        tibble() %>%
        write_tsv(paste0("data/processed_data/ancombc_by_cohort_HMP_vs_MDACC_all_", i, "_significant_results.tsv"))




    # Run AncomBC version 2
    ancombc_out2 <- ancombc2(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0.25,
        fix_formula = "cohort",
        struc_zero = FALSE,
        pseudo_sens = FALSE,
        neg_lb = FALSE,
        n_cl = 24
    )


    # Get AncomBC2 results as a dataframe
    get_ancombc2_df(ancombc_out2$res, taxa_level = "Species", group_var = "cohort") %>%
        write_tsv(paste0("data/processed_data/ancombc2_by_cohort_HMP_vs_MDACC_all_", i, "_significant_results.tsv"))
}




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Combo figures using Maaslin and AncomBC ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────


# Taxa table for Maaslin2 output look up
df_taxa_table <- data.frame(tax_table(physeq_HMP_SSO_combo_species)) %>%
    rownames_to_column(var = "tax_id") %>%
    dplyr::select(-c(Kingdom:Genus)) %>%
    mutate(Species = gsub("s__", "", Species)) %>%
    mutate(Species = gsub("_", " ", Species))


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Glioma Progression category
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_combo_res_hmp_sso <- data.frame()
for (i in c("Stool", "Saliva", "Oral swab")) {
    df_masslin_out <- read_tsv(paste0("data/processed_data/maaslin_by_cohort_HMP_Vs_MDACC_all_", i, "/significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        mutate(feature_copy = feature) %>%
        mutate(feature = gsub("X", "", feature)) %>%
        filter(metadata == "cohort") %>%
        mutate(log2fc = log2(exp(coef))) %>%
        mutate(groupings = case_when(
            coef > 0 ~ "MDACC",
            coef < 0 ~ "HMP"
        )) %>%
        filter(`N.not.0` >= N * 0.25) %>%
        filter(qval < 0.25) %>%
        filter(abs(log2fc) > 2) %>%
        inner_join(df_taxa_table, by = c("feature" = "tax_id")) %>%
        dplyr::select(Species = Species, log2fc, qval, groupings, sample_type) %>%
        mutate(method = "Maaslin2")


    df_ancombc_out <- read_tsv(paste0("data/processed_data/ancombc_by_cohort_HMP_vs_MDACC_all_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "cohort") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "MDACC",
            log2fc < 0 ~ "HMP"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC")


    df_ancombc2_out <- read_tsv(paste0("data/processed_data/ancombc2_by_cohort_HMP_vs_MDACC_all_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "cohort") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "MDACC",
            log2fc < 0 ~ "HMP"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC2")

    df_combo_res_hmp_sso <- bind_rows(df_combo_res_hmp_sso, df_masslin_out, df_ancombc_out, df_ancombc2_out)
}

df_combo_res_hmp_sso_plt <- df_combo_res_hmp_sso %>%
    filter(abs(log2fc) > 2) %>%
    filter(qval < 0.25) %>%
    write_csv("data/processed_data/df_HMP_MDACC_by_cohort_methods_combined_output_plt.csv")



# Filter for species that are significant with atleast 2 methods within each sample type
df_species_filt_combo_res_prog_cat <- df_combo_res_hmp_sso_plt %>%
    group_by(sample_type, Species) %>%
    summarize(count = n()) %>%
    filter(count > 1)

df_combo_plot_hmp_mdacc <- df_combo_res_hmp_sso_plt %>%
    inner_join(df_species_filt_combo_res_prog_cat, by = c("sample_type", "Species")) %>%
    mutate(sample_type = factor(sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
    mutate(log2fc = round(log2fc, 2))




ggplot(df_combo_plot_hmp_mdacc, aes(y = Species, x = method, color = groupings)) +
    geom_tile(aes(fill = log2fc), color = "white") +
    scale_fill_gradient2(midpoint = 0, mid = "white", high = "#e8a122", low = "#225ea8") +
    geom_text(aes(x = method, y = Species, label = log2fc),
        color = ifelse(abs(df_combo_plot_hmp_mdacc$log2fc) > 7, "white", "black"),
        fontface = "bold", size = 7
    ) +
    facet_grid(sample_type ~ ., scales = "free_y", space = "free_y") +
    theme_light(base_size = 20) +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    labs(title = "(Glioma + MET) vs HMP", x = "Method", y = "Species") +
    labs(caption = "Positive log2fc values indicate higher abundance in MDACC (Glioma + MET) samples")
ggsave("output/figures/MDACC_GL_PLUS_BM_vs_HMP_SSO.pdf", width = 13, height = 22, dpi = 300)
