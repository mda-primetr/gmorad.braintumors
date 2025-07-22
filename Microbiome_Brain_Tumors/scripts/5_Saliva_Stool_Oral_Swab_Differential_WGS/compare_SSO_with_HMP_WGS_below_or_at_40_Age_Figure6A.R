source("src/libraries.R")
source("src/common_functions.R")

# Load Unrarefied data
physeq_SSO_clean <- readRDS(
    "data/raw_data/WGS/combined_MDA_sub_cosmos_final.rds"
)

# Summarize the samples that went into this analysis ----
sample_data(physeq_SSO_clean) %>%
    data.frame() %>%
    mutate(age_at_surgery = as.numeric(age_at_surgery)) %>%
    filter(age_at_surgery <= 40) %>%
    dplyr::select(
        coded_id,
        age_at_surgery,
        tumor_category,
        microbiome_sample_type
    ) %>%
    group_by(tumor_category, microbiome_sample_type) %>%
    summarize(
        count = n()
    ) %>%
    write_csv(
        "output/tables/WGS/SSO_40_age_tumor_category_microbiome_sample_type.csv"
    )

# Load HMP data
physeq_WGS_HMP <- readRDS("data/raw_data/WGS/selected_SSO_HMP.rds") %>%
    speedyseq::mutate_sample_data(
        age_category = "less than 40"
    ) %>%
    speedyseq::mutate_sample_data(tumor_category = "NA") %>%
    speedyseq::mutate_sample_data(
        sex = case_when(
            host_sex == "male" ~ "M",
            host_sex == "female" ~ "F",
            TRUE ~ "Unknown"
        )
    )


ps_sso_clean_melt <- psmelt(physeq_SSO_clean) %>%
    filter(Abundance > 0) %>%
    mutate(age_at_surgery = as.numeric(age_at_surgery)) %>%
    filter(age_at_surgery <= 40) %>%
    dplyr::select(
        Sample,
        sex,
        tumor_category,
        microbiome_sample_type,
        Abundance,
        OTU,
        Kingdom:Species
    ) %>%
    mutate(cohort = "MDACC", .after = "microbiome_sample_type")


ps_melt_HMP_oral <- physeq_WGS_HMP %>%
    speedyseq::filter_sample_data(
        isolation_source %in%
            c(
                "G_DNA_Buccal mucosa"
            )
    ) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    mutate(microbiome_sample_type = "Oral swab") %>%
    dplyr::select(
        Sample,
        sex,
        tumor_category,
        microbiome_sample_type,
        Abundance,
        OTU,
        Kingdom:Species
    ) %>%
    mutate(cohort = "HMP", .after = "microbiome_sample_type")


ps_melt_HMP_saliva <- physeq_WGS_HMP %>%
    speedyseq::filter_sample_data(
        isolation_source %in%
            c(
                "G_DNA_Saliva"
            )
    ) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    mutate(microbiome_sample_type = "Saliva") %>%
    dplyr::select(
        Sample,
        sex,
        tumor_category,
        microbiome_sample_type,
        Abundance,
        OTU,
        Kingdom:Species
    ) %>%
    mutate(cohort = "HMP", .after = "microbiome_sample_type")



ps_melt_HMP_stool <- physeq_WGS_HMP %>%
    speedyseq::filter_sample_data(
        isolation_source %in%
            c(
                "G_DNA_Stool"
            )
    ) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    mutate(microbiome_sample_type = "Stool") %>%
    dplyr::select(
        Sample,
        sex,
        tumor_category,
        microbiome_sample_type,
        Abundance,
        OTU,
        Kingdom:Species
    ) %>%
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

saveRDS(
    physeq_HMP_SSO_combo_species,
    "data/processed_data/WGS/physeq_HMP_SSO_combo_species_age_40_and_below.rds"
)


for (i in c("Stool", "Saliva", "Oral swab")) {
    sample_data(physeq_HMP_SSO_combo_species) %>%
        data.frame() %>%
        filter(microbiome_sample_type == i) %>%
        dplyr::select(
            microbiome_sample_type,
            cohort,
            sex
        ) %>%
        gtsummary::tbl_summary(by = cohort) %>%
        gtsummary::modify_spanning_header(
            all_stat_cols() ~ paste0("**HMP vs SSO ** : ", i)
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
                "output/tables/WGS/HMP_SSO_WGS_40",
                i,
                ".docx"
            )
        )
}


# Differential abundance between Glioma and Met combo vs HMP -----

# No need to control
for (i in c("Stool", "Saliva", "Oral swab")) {
    if (i %in% c("Saliva", "Oral swab")) {
        final_phyloseq_before_rarefaction_species_i <- physeq_HMP_SSO_combo_species %>%
            speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
            speedyseq::mutate_sample_data(
                cohort = factor(cohort, levels = c("HMP", "MDACC"))
            ) %>%
            speedyseq::filter_taxa(
                .,
                function(x) sum(x >= 100) > (0.25 * length(x)),
                TRUE
            ) %>%
            speedyseq::filter_sample_data(sample_sums(.) >= 500)
    }

    if (i == "Stool") {
        final_phyloseq_before_rarefaction_species_i <- physeq_HMP_SSO_combo_species %>%
            speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
            speedyseq::mutate_sample_data(
                cohort = factor(cohort, levels = c("HMP", "MDACC"))
            ) %>%
            speedyseq::filter_taxa(
                .,
                function(x) sum(x >= 100) > (0.25 * length(x)),
                TRUE
            ) %>%
            speedyseq::filter_sample_data(sample_sums(.) >= 1000)
    }

    df_maaslin_all <- data.frame(
        otu_table(final_phyloseq_before_rarefaction_species_i) %>% t()
    )

    df_maaslin_all_meta <- data.frame(sample_data(
        final_phyloseq_before_rarefaction_species_i
    ))

    # Run Maaslin2

    maaslin_by_sample <- Maaslin2(
        input_data = df_maaslin_all,
        input_metadata = df_maaslin_all_meta,
        output = paste0(
            "data/processed_data/WGS/maaslin_by_cohort_HMP_Vs_MDACC_40_",
            i
        ),
        analysis_method = "NEGBIN",
        fixed_effects = c("cohort"),
        plot_scatter = FALSE,
        plot_heatmap = FALSE,
        min_abundance = 0,
        min_prevalence = 0,
        cores = 24,
        transform = "NONE",
        normalization = "TMM",
        reference = c("cohort,HMP"),
        save_models = TRUE,
        standardize = TRUE
    )

    # Run AncomBC version 1
    ancombc_out <- ancombc(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0,
        formula = "cohort",
        struc_zero = FALSE,
        neg_lb = FALSE,
        n_cl = 24
    )

    # Get AncomBC1 results as a dataframe
    get_ancombc_2_8_df(
        ancombc_out$res,
        taxa_level = "Species",
        group_var = "cohort",
        group_var_order = c("HMP", "MDACC"),
        phyloseq_obj = final_phyloseq_before_rarefaction_species_i
    ) %>%
        tibble() %>%
        write_tsv(paste0(
            "data/processed_data/WGS/ancombc_by_cohort_HMP_vs_MDACC_40_",
            i,
            "_significant_results.tsv"
        ))

    # Run AncomBC version 2
    ancombc_out2 <- ancombc2(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0,
        fix_formula = "cohort",
        struc_zero = FALSE,
        pseudo_sens = FALSE,
        neg_lb = FALSE,
        n_cl = 24
    )

    # Get AncomBC2 results as a dataframe
    get_ancombc2_2_8_df(
        ancombc_out2$res,
        taxa_level = "Species",
        group_var = "cohort",
        group_var_order = c("HMP", "MDACC"),
        phyloseq_obj = final_phyloseq_before_rarefaction_species_i
    ) %>%
        write_tsv(paste0(
            "data/processed_data/WGS/ancombc2_by_cohort_HMP_vs_MDACC_40_",
            i,
            "_significant_results.tsv"
        ))
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
    df_masslin_out <- read_tsv(paste0(
        "data/processed_data/WGS/maaslin_by_cohort_HMP_Vs_MDACC_40_",
        i,
        "/significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        mutate(feature_copy = feature) %>%
        mutate(feature = gsub("X", "", feature)) %>%
        filter(metadata == "cohort") %>%
        mutate(log2fc = log2(exp(coef))) %>%
        mutate(
            groupings = case_when(
                coef > 0 ~ "MDACC",
                coef < 0 ~ "HMP"
            )
        ) %>%
        filter(`N.not.0` >= N * 0.25) %>%
        filter(qval < 0.25) %>%
        filter(abs(log2fc) > 2) %>%
        inner_join(df_taxa_table, by = c("feature" = "tax_id")) %>%
        dplyr::select(
            Species = Species,
            log2fc,
            qval,
            groupings,
            sample_type
        ) %>%
        mutate(method = "Maaslin2")

    df_ancombc_out <- read_tsv(paste0(
        "data/processed_data/WGS/ancombc_by_cohort_HMP_vs_MDACC_40_",
        i,
        "_significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "cohort") %>%
        mutate(
            groupings = case_when(
                log2fc > 0 ~ "MDACC",
                log2fc < 0 ~ "HMP"
            )
        ) %>%
        dplyr::select(
            Species = feature,
            log2fc,
            qval,
            groupings,
            sample_type
        ) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC")

    df_ancombc2_out <- read_tsv(paste0(
        "data/processed_data/WGS/ancombc2_by_cohort_HMP_vs_MDACC_40_",
        i,
        "_significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "cohort") %>%
        mutate(
            groupings = case_when(
                log2fc > 0 ~ "MDACC",
                log2fc < 0 ~ "HMP"
            )
        ) %>%
        dplyr::select(
            Species = feature,
            log2fc,
            qval,
            groupings,
            sample_type
        ) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC2")

    df_combo_res_hmp_sso <- bind_rows(
        df_combo_res_hmp_sso,
        df_masslin_out,
        df_ancombc_out,
        df_ancombc2_out
    )
}

df_combo_res_hmp_sso_plt <- df_combo_res_hmp_sso %>%
    filter(abs(log2fc) > 2) %>%
    filter(qval < 0.25) %>%
    write_csv(
        "data/processed_data/WGS/df_HMP_MDACC_by_cohort_methods_combined_output_40_plt.csv"
    )


# Filter for species that are significant with atleast 2 methods within each sample type
df_species_filt_combo_res_prog_cat <- df_combo_res_hmp_sso_plt %>%
    group_by(sample_type, Species) %>%
    summarize(count = n()) %>%
    filter(count > 1)


# First, check consistency of groupings across methods for each species and sample type
species_grouping_consistency <- df_combo_res_hmp_sso_plt %>%
    group_by(Species, sample_type) %>%
    summarize(
        consistent_grouping = length(unique(groupings)) == 1,
        groupings_found = paste(unique(groupings), collapse = ","),
        .groups = "drop"
    )

# View species with inconsistent groupings
inconsistent_species <- species_grouping_consistency %>%
    filter(!consistent_grouping) %>%
    arrange(sample_type, Species)

df_combo_plot_hmp_mdacc <- df_combo_res_hmp_sso_plt %>%
    inner_join(
        df_species_filt_combo_res_prog_cat,
        by = c("sample_type", "Species")
    ) %>%
    # Join with consistency information
    left_join(
        species_grouping_consistency,
        by = c("Species", "sample_type")
    ) %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("Saliva", "Oral swab", "Stool")
        )
    ) %>%
    mutate(log2fc = round(log2fc, 2)) %>%
    # Filter out species with inconsistent directions across methods
    filter(consistent_grouping == TRUE)


ggplot(
    df_combo_plot_hmp_mdacc,
    aes(y = Species, x = method, color = groupings)
) +
    geom_tile(aes(fill = log2fc), color = "white") +
    scale_fill_gradient2(
        midpoint = 0,
        mid = "white",
        high = "#e8a122",
        low = "#225ea8"
    ) +
    geom_text(
        aes(x = method, y = Species, label = log2fc),
        color = ifelse(
            abs(df_combo_plot_hmp_mdacc$log2fc) > 7,
            "white",
            "black"
        ),
        fontface = "bold",
        size = 7
    ) +
    facet_grid(sample_type ~ ., scales = "free_y", space = "free_y") +
    theme_light(base_size = 20) +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(title = "(Glioma + MET) vs HMP (<=40)", x = "Method", y = "Species") +
    labs(
        caption = "Positive log2fc values indicate higher abundance in MDACC (Glioma + MET) samples"
    )
ggsave(
    "output/figures/WGS/MDACC_GL_PLUS_BM_vs_HMP_SSO_40.pdf",
    width = 13,
    height = 22,
    dpi = 300
)
