source("src/libraries.R")
source("src/common_functions.R")


# Load SSO data ----
physeq_SSO_clean <- readRDS(
    "data/raw_data/WGS/combined_MDA_sub_cosmos_final.rds"
)

# Metadata
df_metadata <- data.frame(sample_data(physeq_SSO_clean)) %>%
    rownames_to_column(var = "sequencer_output_id")


physeq_SSO_clean_species <- physeq_SSO_clean %>%
    speedyseq::tax_glom(taxrank = "Species")

sample_data(physeq_SSO_clean_species) %>% colnames()

# Load tumor size data ----
df_tumor_size <- read_csv(
    "data/metadata/Tumor_Size_metadata_2024_10_27_revised.csv"
)


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Brm/Met -----
# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Stats for BM by Progression overall (Brain Mets) -----
for (i in c("Stool", "Saliva", "Oral swab")) {
    sample_data(physeq_SSO_clean) %>%
        data.frame() %>%
        filter(microbiome_sample_type == i) %>%
        filter(tumor_category == "Met") %>%
        filter(br_m_progression_overall != "Unknown") %>%
        filter(br_m_progression_overall %in% c("Y", "N")) %>%
        mutate(age_at_surgery = as.numeric(age_at_surgery)) %>%
        mutate(
            BrM_total_follow_up.days. = as.numeric(br_m_total_follow_up_days)
        ) %>%
        left_join(df_tumor_size, by = c("coded_id" = "coded_id")) %>%
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
            BrM_total_follow_up.days.,
            "Tumor size" = tumor_size
        ) %>%
        gtsummary::tbl_summary(by = BrM_progression_overall) %>%
        gtsummary::modify_spanning_header(
            all_stat_cols() ~ paste0("**Brain Mets Overall Progression** : ", i)
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
                "output/tables/WGS/Brain_Mets_by_overall_progression_SSO_WGS_",
                i,
                ".docx"
            )
        )
}


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# MET by Progression caterogical differential abundance ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Maaslin ---

for (i in c("Stool", "Saliva", "Oral swab")) {
    final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
        speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
        speedyseq::filter_sample_data(tumor_category == "Met") %>%
        speedyseq::filter_sample_data(
            br_m_progression_overall %in% c("Y", "N")
        ) %>%
        speedyseq::filter_taxa(., function(x) sum(x >= 100) > (0.25 * length(x)), TRUE)

    df_maaslin_all <- data.frame(
        otu_table(final_phyloseq_before_rarefaction_species_i) %>% t()
    )

    df_maaslin_all_meta <- data.frame(sample_data(
        final_phyloseq_before_rarefaction_species_i
    )) %>%
        filter(br_m_progression_overall %in% c("Y", "N"))

    # Please refer output/tables/WGS
    # Since shotgun_seq_batch is a significant factor in the metadata for saliva and stool,
    # we will include it in the fixed effects
    fixed_effects_vector <- c("br_m_progression_overall")
    if (i %in% c("Oral swab", "Saliva")) {
        # For oral swab and saliva, we will include shotgun_seq_batch as a fixed effect
        fixed_effects_vector <- c(
            fixed_effects_vector,
            "shotgun_seq_batch"
        )
        reference_vector <- c(
            "br_m_progression_overall,N;shotgun_seq_batch,Batch 1"
        )
    } else {
        reference_vector <- c("br_m_progression_overall,N")
    }


    # Run Maaslin2
    maaslin_by_sample <- Maaslin2(
        input_data = df_maaslin_all,
        input_metadata = df_maaslin_all_meta,
        output = paste0(
            "data/processed_data/WGS/maaslin_phyloseq_by_Met_categorical_",
            "br_m_progression_overall",
            "_",
            i
        ),
        analysis_method = "NEGBIN",
        fixed_effects = fixed_effects_vector,
        random_effects = NULL,
        plot_scatter = FALSE,
        plot_heatmap = FALSE,
        min_abundance = 0,
        min_prevalence = 0,
        cores = 24,
        transform = "NONE",
        save_models = TRUE,
        reference = reference_vector,
        normalization = "TMM",
        standardize = TRUE
    )



    # Set the formula based on the covariates
    if (i %in% c("Oral swab", "Saliva")) {
        fixed_effects_formula <- "br_m_progression_overall + shotgun_seq_batch"
    } else {
        fixed_effects_formula <- "br_m_progression_overall"
    }

    # Run AncomBC
    ancombc_out <- ancombc(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0,
        formula = fixed_effects_formula,
        group = "br_m_progression_overall",
        struc_zero = TRUE,
        neg_lb = TRUE,
        n_cl = 24
    )



    get_ancombc_2_8_df(
        ancombc_out$res,
        taxa_level = "Species",
        group_var = "br_m_progression_overall",
        group_var_order = c("N", "Y"),
        phyloseq_obj = final_phyloseq_before_rarefaction_species_i
    ) %>%
        write_tsv(paste0(
            "data/processed_data/WGS/ancombc_phyloseq_by_Met_categorical_",
            "br_m_progression_overall",
            "_",
            i,
            "_significant_results.tsv"
        ))

    # Run AncomBC2
    ancombc_out2 <- ancombc2(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0,
        fix_formula = fixed_effects_formula,
        group = "br_m_progression_overall",
        struc_zero = TRUE,
        pseudo_sens = TRUE,
        neg_lb = TRUE,
        n_cl = 24
    )


    get_ancombc2_2_8_df(
        ancombc_out2$res,
        taxa_level = "Species",
        group_var = "br_m_progression_overall",
        group_var_order = c("N", "Y"),
        phyloseq_obj = final_phyloseq_before_rarefaction_species_i
    ) %>%
        write_tsv(paste0(
            "data/processed_data/WGS/ancombc2_phyloseq_by_Met_categorical_",
            "br_m_progression_overall",
            "_",
            i,
            "_significant_results.tsv"
        ))
}
# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Combo figures using Maaslin and AncomBC ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Taxa table for Maaslin2 output look up
df_taxa_table <- data.frame(tax_table(physeq_SSO_clean_species)) %>%
    rownames_to_column(var = "tax_id") %>%
    dplyr::select(-c(Kingdom:Genus)) %>%
    mutate(Species = gsub("s__", "", Species)) %>%
    mutate(Species = gsub("_", " ", Species))


df_combo_res_Met_prog_cat <- data.frame()
for (i in c("Stool", "Saliva", "Oral swab")) {
    df_masslin_out <- read_tsv(paste0(
        "data/processed_data/WGS/maaslin_phyloseq_by_Met_categorical_br_m_progression_overall_",
        i,
        "/significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        mutate(feature_copy = feature) %>%
        mutate(feature = gsub("X", "", feature)) %>%
        filter(metadata == "br_m_progression_overall") %>%
        mutate(log2fc = log2(exp(coef))) %>%
        mutate(
            groupings = case_when(
                coef > 0 ~ "Y",
                coef < 0 ~ "N"
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
        "data/processed_data/WGS/ancombc_phyloseq_by_Met_categorical_br_m_progression_overall_",
        i,
        "_significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "br_m_progression_overall") %>%
        mutate(
            groupings = case_when(
                log2fc > 0 ~ "Y",
                log2fc < 0 ~ "N"
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
        "data/processed_data/WGS/ancombc2_phyloseq_by_Met_categorical_br_m_progression_overall_",
        i,
        "_significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "br_m_progression_overall") %>%
        mutate(
            groupings = case_when(
                log2fc > 0 ~ "Y",
                log2fc < 0 ~ "N"
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

    df_combo_res_Met_prog_cat <- bind_rows(
        df_combo_res_Met_prog_cat,
        df_masslin_out,
        df_ancombc_out,
        df_ancombc2_out
    )
}

df_combo_res_MET_prog_cat_plt <- df_combo_res_Met_prog_cat %>%
    filter(abs(log2fc) > 2) %>%
    filter(qval < 0.25) %>%
    write_csv(
        "data/processed_data/WGS/df_combo_MET_by_progression_categorical_multiple_methods_combined_output_plt.csv"
    )


# Filter for species that are significant with atleast 2 methods within each sample type
df_species_filt_MET_combo_res_prog_cat <- df_combo_res_MET_prog_cat_plt %>%
    group_by(sample_type, Species) %>%
    summarize(count = n()) %>%
    filter(count > 1)

# First, check consistency of groupings across methods for each species and sample type
species_grouping_consistency <- df_combo_res_MET_prog_cat_plt %>%
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

# Now create the plot data, filtering out inconsistent species or keeping them with a note
df_combo_plot_prog_cat_met <- df_combo_res_MET_prog_cat_plt %>%
    inner_join(
        df_species_filt_MET_combo_res_prog_cat,
        by = c("Species", "sample_type")
    ) %>%
    # Join with consistency information
    left_join(species_grouping_consistency, by = c("Species", "sample_type")) %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("Saliva", "Oral swab", "Stool")
        )
    ) %>%
    mutate(log2fc = round(log2fc, 2)) %>%
    # Filter out species with inconsistent directions across methods
    filter(consistent_grouping == TRUE)


plot_met_progression_heatmap <- df_combo_plot_prog_cat_met %>%
    ggplot(., aes(y = Species, x = method, color = groupings)) +
    geom_tile(aes(fill = log2fc), color = "white") +
    scale_fill_gradient2(
        midpoint = 0,
        limits = c(-16, 16),
        mid = "white",
        low = "#e8a122",
        high = "#238443"
    ) +
    geom_text(
        aes(x = method, y = Species, label = log2fc),
        color = ifelse(
            abs(df_combo_plot_prog_cat_met$log2fc) > 7,
            "white",
            "black"
        ),
        fontface = "bold",
        size = 7
    ) +
    facet_grid(sample_type ~ ., scales = "free_y", space = "free_y") +
    theme_light(base_size = 20) +
    theme(
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(title = "MET Progression category", x = "Method", y = "Species") +
    labs(
        caption = "Positive log2fc values indicate higher abundance in Y progression group"
    )
ggsave(
    "output/figures/WGS/MET_Progression_category_SSO.pdf",
    width = 12,
    height = 14,
    dpi = 300
)


plot_met_progression_heatmap_data <- plot_met_progression_heatmap$data %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("Saliva", "Oral swab", "Stool")
        )
    ) %>%
    group_by(sample_type) %>%
    distinct(Species) %>%
    arrange(desc(Species), .by_group = T)


plot_met_progression_abundance_box <- physeq_SSO_clean_species %>%
    transform_sample_counts(., function(x) (x / sum(x)) * 100) %>%
    fn_clean_tax_table(.) %>%
    ps_filter(tumor_category == "Met") %>%
    ps_filter(br_m_progression_overall %in% c("Y", "N")) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    dplyr::select(
        microbiome_sample_type,
        Species,
        br_m_progression_overall,
        Abundance
    ) %>%
    mutate(
        microbiome_sample_type = factor(
            microbiome_sample_type,
            levels = c("Saliva", "Oral swab", "Stool")
        )
    ) %>%
    inner_join(
        plot_met_progression_heatmap_data,
        by = c("microbiome_sample_type" = "sample_type", "Species")
    ) %>%
    ggplot(
        .,
        aes(x = Abundance, y = Species, fill = br_m_progression_overall)
    ) +
    geom_boxplot(linewidth = 0.1) +
    geom_point(
        aes(group = br_m_progression_overall),
        position = position_dodge(width = .75)
    ) +
    scale_fill_manual(values = c("N" = "#e8a122", "Y" = "#238443")) +
    scale_x_log10(labels = scales::comma) +
    facet_grid(microbiome_sample_type ~ ., scales = "free_y") +
    theme_light(base_size = 18) +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    ylab(NULL) +
    ggh4x::force_panelsizes(rows = c(0.31, 0.060, 1)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )


plot_met_progression_heatmap +
    plot_met_progression_abundance_box +
    plot_layout(guides = "collect", width = c(0.2, 0.3), nrow = 1)
ggsave(
    "output/figures/WGS/BRM_MET_progression_SSO_combo_box.pdf",
    width = 17,
    height = 14,
    dpi = 300
)
