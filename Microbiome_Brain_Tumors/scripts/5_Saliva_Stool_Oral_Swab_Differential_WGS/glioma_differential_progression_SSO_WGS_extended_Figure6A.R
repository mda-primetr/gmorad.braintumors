source("src/libraries.R")
source("src/common_functions.R")


# Load Unrarefied data
physeq_SSO_clean <- readRDS(
    "data/raw_data/WGS/combined_MDA_sub_cosmos_final.rds"
)

# Metadata
df_metadata <- data.frame(sample_data(physeq_SSO_clean)) %>%
    rownames_to_column(var = "sequencer_output_id")


physeq_SSO_clean_species <- physeq_SSO_clean %>%
    speedyseq::tax_glom(taxrank = "Species") %>%
    microViz::ps_mutate(age_at_surgery = as.numeric(age_at_surgery))


# Load tumor size data ----
df_tumor_size <- read_csv(
    "data/metadata/Tumor_Size_metadata_2024_10_27_revised.csv"
)


# Stats for Glioma for Progression category ----
for (i in c("Stool", "Saliva", "Oral swab")) {
    sample_data(physeq_SSO_clean_species) %>%
        data.frame() %>%
        filter(microbiome_sample_type == i) %>%
        filter(tumor_category == "Glioma") %>%
        filter(glioma_idh_mutation != "Unknown") %>%
        filter(grade_category != "Unknown") %>%
        mutate(glioma_fu_time = as.numeric(glioma_fu_time)) %>%
        filter(glioma_progression_overall != "Unknown") %>%
        left_join(df_tumor_size, by = c("coded_id" = "coded_id")) %>%
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
            glioma_recurrent_status_at_surg,
            "Tumor size" = tumor_size
        ) %>%
        gtsummary::tbl_summary(by = glioma_progression_overall) %>%
        gtsummary::modify_spanning_header(
            all_stat_cols() ~ paste0("**Glioma Progression overall** : ", i)
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
                "output/tables/WGS/Glioma_glioma_progression_overall_SSO_WGS_",
                i,
                ".docx"
            )
        )
}


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# GLioma by Progression caterogical differential abundance ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Maaslin ---

for (i in c("Stool", "Saliva", "Oral swab")) {
    final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
        speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
        speedyseq::filter_sample_data(tumor_category == "Glioma") %>%
        speedyseq::filter_sample_data(
            glioma_progression_overall %in% c("Y", "N")
        ) %>%
        speedyseq::filter_sample_data(glioma_fu_time != "Unknown") %>%
        speedyseq::mutate_sample_data(
            glioma_fu_time = as.numeric(glioma_fu_time)
        ) %>%
        speedyseq::mutate_sample_data(
            steroids_at_surgery = factor(
                steroids_at_surgery,
                levels = c("N", "Y")
            )
        ) %>%
        speedyseq::mutate_sample_data(
            sex = factor(sex, levels = c("F", "M"))
        ) %>%
        speedyseq::mutate_sample_data(
            grade_category = factor(grade_category, levels = c("Low", "High"))
        ) %>%
        speedyseq::mutate_sample_data(
            glioma_idh_mutation = factor(
                glioma_idh_mutation,
                levels = c("Wild", "Mutant")
            )
        ) %>%
        speedyseq::filter_taxa(., function(x) sum(x >= 100) > (0.25 * length(x)), TRUE)
    # Filter for taxa that are present in at least 25% of the samples with abundance >= 100

    df_maaslin_all <- data.frame(
        otu_table(final_phyloseq_before_rarefaction_species_i) %>% t()
    )




    df_maaslin_all_meta <- data.frame(sample_data(
        final_phyloseq_before_rarefaction_species_i
    )) %>%
        filter(glioma_progression_overall %in% c("Y", "N"))

    # Please refer output/tables folder
    # Since shotgun_seq_batch is a significant factor in the metadata for saliva and stool,
    # we will include it in the fixed effects
    fixed_effects_vector <- c(
        "glioma_progression_overall",
        "age_at_surgery",
        "grade_category",
        "glioma_idh_mutation"
    )
    if (i %in% c("Stool")) {
        fixed_effects_vector <- c(fixed_effects_vector, "steroids_at_surgery")
        reference_vector <- c(
            "glioma_progression_overall,N;steroids_at_surgery,N;grade_category,Low;glioma_idh_mutation,Wild"
        )
    } else {
        reference_vector <- c(
            "glioma_progression_overall,N;grade_category,Low;glioma_idh_mutation,Wild"
        )
    }

    # Please refer output/tables folder
    maaslin_by_sample <- Maaslin2(
        input_data = df_maaslin_all,
        input_metadata = df_maaslin_all_meta,
        output = paste0(
            "data/processed_data/WGS/maaslin_phyloseq_by_Glioma_categorical_glioma_progression_overall_",
            i
        ),
        analysis_method = "NEGBIN",
        fixed_effects = fixed_effects_vector,
        plot_scatter = FALSE,
        plot_heatmap = FALSE,
        min_abundance = 0,
        min_prevalence = 0,
        cores = 24,
        transform = "NONE",
        reference = reference_vector,
        normalization = "TMM",
        save_models = TRUE,
        standardize = TRUE
    )

    # Conditional covariates adjustments based on the metadata stats -----
    fixed_effects_formula <- c(
        "glioma_progression_overall+ age_at_surgery+ grade_category+ glioma_idh_mutation"
    )
    if (i %in% c("Stool")) {
        fixed_effects_formula <- c(
            "glioma_progression_overall+ age_at_surgery+ steroids_at_surgery + grade_category+ glioma_idh_mutation"
        )
    }

    # AncomBC version 1----
    ancombc_out <- ancombc(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0,
        formula = fixed_effects_formula,
        group = "glioma_progression_overall",
        struc_zero = TRUE,
        neg_lb = TRUE,
        n_cl = 24
    )

    # store AncomBC1 results
    get_ancombc_2_8_df(
        ancombc_out$res,
        taxa_level = "Species",
        group_var = "glioma_progression_overall",
        group_var_order = c("N", "Y"),
        phyloseq_obj = final_phyloseq_before_rarefaction_species_i
    ) %>%
        tibble() %>%
        write_tsv(paste0(
            "data/processed_data/WGS/ancombc_phyloseq_by_Glioma_categorical_glioma_progression_overall_",
            i,
            "_significant_results.tsv"
        ))

    # AncomBC version 2----
    ancombc_out2 <- ancombc2(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0,
        fix_formula = fixed_effects_formula,
        group = "glioma_progression_overall",
        struc_zero = TRUE,
        pseudo_sens = TRUE,
        neg_lb = TRUE,
        n_cl = 24
    )

    # store AncomBC2 results
    get_ancombc2_2_8_df(
        ancombc_out2$res,
        taxa_level = "Species",
        group_var = "glioma_progression_overall",
        group_var_order = c("N", "Y"),
        phyloseq_obj = final_phyloseq_before_rarefaction_species_i
    ) %>%
        write_tsv(paste0(
            "data/processed_data/WGS/ancombc2_phyloseq_by_Glioma_categorical_glioma_progression_overall_",
            i,
            "_significant_results.tsv"
        ))
}


# Combo figures using Maaslin and AncomBC ----

# Taxa table for Maaslin2 output look up
df_taxa_table <- data.frame(tax_table(physeq_SSO_clean_species)) %>%
    rownames_to_column(var = "tax_id") %>%
    dplyr::select(-c(Kingdom:Genus)) %>%
    mutate(Species = gsub("s__", "", Species)) %>%
    mutate(Species = gsub("_", " ", Species))


# Glioma Progression category combine outputs ----

df_combo_res_prog_cat <- data.frame()
for (i in c("Stool", "Saliva", "Oral swab")) {
    df_masslin_out <- read_tsv(paste0(
        "data/processed_data/WGS/maaslin_phyloseq_by_Glioma_categorical_glioma_progression_overall_",
        i,
        "/significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        mutate(feature_copy = feature) %>%
        mutate(feature = gsub("X", "", feature)) %>%
        filter(metadata == "glioma_progression_overall") %>%
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
        "data/processed_data/WGS/ancombc_phyloseq_by_Glioma_categorical_glioma_progression_overall_",
        i,
        "_significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "glioma_progression_overall") %>%
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
        "data/processed_data/WGS/ancombc2_phyloseq_by_Glioma_categorical_glioma_progression_overall_",
        i,
        "_significant_results.tsv"
    )) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "glioma_progression_overall") %>%
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

    df_combo_res_prog_cat <- bind_rows(
        df_combo_res_prog_cat,
        df_masslin_out,
        df_ancombc_out,
        df_ancombc2_out
    )
}

df_combo_res_prog_cat_plt <- df_combo_res_prog_cat %>%
    filter(abs(log2fc) > 2) %>%
    filter(qval < 0.25) %>%
    write_csv(
        "data/processed_data/WGS/df_combo_Glioma_by_progression_categorical_multiple_methods_combined_output_plt.csv"
    )


# Filter for species that are significant with at least 2 methods within each sample type
df_species_filt_combo_res_prog_cat <- df_combo_res_prog_cat_plt %>%
    group_by(sample_type, Species) %>%
    summarize(count = n()) %>%
    filter(count > 1)


# First, check consistency of groupings across methods for each species and sample type
species_grouping_consistency <- df_combo_res_prog_cat_plt %>%
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

df_combo_plot_prog_cat_glioma <- df_combo_res_prog_cat_plt %>%
    inner_join(
        df_species_filt_combo_res_prog_cat,
        by = c("Species", "sample_type")
    ) %>%
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


plot_glioma_progression_heatmap <- ggplot(
    df_combo_plot_prog_cat_glioma,
    aes(y = Species, x = method, color = groupings)
) +
    geom_tile(aes(fill = log2fc), color = "white") +
    scale_fill_gradient2(
        midpoint = 0,
        mid = "white",
        low = "#e8a122",
        high = "#225ea8"
    ) +
    geom_text(
        aes(x = method, y = Species, label = log2fc),
        color = ifelse(
            abs(df_combo_plot_prog_cat_glioma$log2fc) > 7,
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
    labs(title = "Glioma Progression category", x = "Method", y = "Species") +
    labs(
        caption = "Positive log2fc values indicate higher abundance in Y progression group"
    )
ggsave(
    "output/figures/WGS/Glioma_Progression_category_SSO.pdf",
    width = 10,
    height = 14,
    dpi = 300
)


plot_glioma_progression_heatmap_data <- plot_glioma_progression_heatmap$data %>%
    mutate(
        sample_type = factor(
            sample_type,
            levels = c("Saliva", "Oral swab", "Stool")
        )
    ) %>%
    group_by(sample_type) %>%
    distinct(Species) %>%
    arrange(desc(Species), .by_group = T)


plot_glioma_progression_abundance_box <- physeq_SSO_clean_species %>%
    transform_sample_counts(., function(x) (x / sum(x)) * 100) %>%
    fn_clean_tax_table(.) %>%
    ps_filter(tumor_category == "Glioma") %>%
    ps_filter(glioma_progression_overall %in% c("Y", "N")) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    dplyr::select(
        microbiome_sample_type,
        Species,
        glioma_progression_overall,
        Abundance
    ) %>%
    mutate(
        microbiome_sample_type = factor(
            microbiome_sample_type,
            levels = c("Saliva", "Oral swab", "Stool")
        )
    ) %>%
    inner_join(
        plot_glioma_progression_heatmap_data,
        by = c("microbiome_sample_type" = "sample_type", "Species")
    ) %>%
    ggplot(
        .,
        aes(x = Abundance, y = Species, fill = glioma_progression_overall)
    ) +
    geom_boxplot(linewidth = 0.1) +
    geom_point(
        aes(group = glioma_progression_overall),
        position = position_dodge(width = .75)
    ) +
    scale_fill_manual(values = c("N" = "#e8a122", "Y" = "#225ea8")) +
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
    ggh4x::force_panelsizes(rows = c(0.06, 0.26, 0.8)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )


plot_glioma_progression_heatmap +
    plot_glioma_progression_abundance_box +
    plot_layout(
        guides = "collect",
        width = c(0.2, 0.3),
        nrow = 1
    )
ggsave(
    "output/figures/WGS/Glioma_progression_SSO_combo_box.pdf",
    width = 17,
    height = 7,
    dpi = 300
)
