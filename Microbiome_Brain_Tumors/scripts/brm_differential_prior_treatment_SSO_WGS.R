source("src/libraries.R")
source("src/common_functions.R")

# Load Unrarefied data
load("data/processed_data/physeq_SSO_clean.RData")
physeq_SSO_clean

# Metadata
df_metadata <- data.frame(sample_data(physeq_SSO_clean))


physeq_SSO_clean_species <- physeq_SSO_clean %>%
    speedyseq::tax_glom(taxrank = "Species")



# Differential abundance analysis using Maaslin2 ----
# For MET

# PRIOR TREATMENT -----
for (i in c("Stool", "Saliva", "Oral swab")) {
    for (j in c("br_m_prior_treatment")) {
        final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
            speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
            speedyseq::filter_sample_data(tumor_category == "Met")

        df_maaslin_all <- data.frame(otu_table(final_phyloseq_before_rarefaction_species_i) %>% t())

        if (j == "br_m_prior_treatment") {
            df_maaslin_all_meta <- data.frame(sample_data(final_phyloseq_before_rarefaction_species_i)) %>%
                mutate(br_m_prior_treatment = factor(br_m_prior_treatment, levels = c("N", "Y"))) %>%
                filter(br_m_prior_treatment %in% c("N", "Y"))

            # Please refer output/figures/metadata_stats
            maaslin_by_sample <- Maaslin2(
                input_data = df_maaslin_all,
                input_metadata = df_maaslin_all_meta,
                output = paste0("data/processed_data/maaslin_phyloseq_by_MET_with_", j, "_", i),
                analysis_method = "NEGBIN",
                fixed_effects = c("br_m_prior_treatment"),
                plot_scatter = FALSE,
                plot_heatmap = FALSE,
                min_abundance = 100,
                min_prevalence = 0.25,
                cores = 24,
                transform = "NONE",
                reference = c("br_m_prior_treatment,N"),
                normalization = "TMM",
                standardize = TRUE
            )
        }
    }
}




# Differential abundance analysis using AncomBC 1 and 2 ----

for (i in c("Stool", "Saliva", "Oral swab")) {
    for (j in c("br_m_prior_treatment")) {
        if (j == "br_m_prior_treatment") {
            final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
                speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
                speedyseq::filter_sample_data(tumor_category == "Met") %>%
                speedyseq::mutate_sample_data(br_m_prior_treatment = factor(br_m_prior_treatment, levels = c("N", "Y")))



            ancombc_out <- ancombc(
                final_phyloseq_before_rarefaction_species_i,
                tax_level = "Species",
                prv_cut = 0.25,
                formula = "br_m_prior_treatment",
                group = "br_m_prior_treatment",
                struc_zero = TRUE,
                neg_lb = TRUE,
                n_cl = 24
            )

            ancombc_out2 <- ancombc2(
                final_phyloseq_before_rarefaction_species_i,
                tax_level = "Species",
                prv_cut = 0.25,
                fix_formula = "br_m_prior_treatment",
                group = "br_m_prior_treatment",
                struc_zero = TRUE,
                pseudo_sens = TRUE,
                neg_lb = TRUE,
                n_cl = 24
            )

            get_ancombc2_df(ancombc_out2$res, taxa_level = "Species", group_var = "br_m_prior_treatment", group_var_order = c("N", "Y")) %>%
                write_tsv(paste0("data/processed_data/ancombc2_phyloseq_by_MET_with_", j, "_", i, "_significant_results.tsv"))

            get_ancombc_df(ancombc_out$res, taxa_level = "Species", group_var = "br_m_prior_treatment", group_var_order = c("N", "Y")) %>%
                tibble() %>%
                write_tsv(paste0("data/processed_data/ancombc_phyloseq_by_MET_with_", j, "_", i, "_significant_results.tsv"))
        }
    }
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




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# MET br_m_prior_treatment
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_combo_res_MET_prior_treatment <- data.frame()
for (i in c("Stool", "Saliva", "Oral swab")) {
    df_masslin_out <- read_tsv(paste0("data/processed_data/maaslin_phyloseq_by_MET_with_br_m_prior_treatment_", i, "/significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        mutate(feature_copy = feature) %>%
        mutate(feature = gsub("X", "", feature)) %>%
        filter(metadata == "br_m_prior_treatment") %>%
        mutate(log2fc = log2(exp(coef))) %>%
        mutate(groupings = case_when(
            coef > 0 ~ "Y",
            coef < 0 ~ "N"
        )) %>%
        filter(`N.not.0` >= N * 0.25) %>%
        filter(qval < 0.25) %>%
        filter(abs(log2fc) > 2) %>%
        inner_join(df_taxa_table, by = c("feature" = "tax_id")) %>%
        dplyr::select(Species = Species, log2fc, qval, groupings, sample_type) %>%
        mutate(method = "Maaslin2")


    df_ancombc_out <- read_tsv(paste0("data/processed_data/ancombc_phyloseq_by_MET_with_br_m_prior_treatment_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "br_m_prior_treatment") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "Y",
            log2fc < 0 ~ "N"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC")


    df_ancombc2_out <- read_tsv(paste0("data/processed_data/ancombc2_phyloseq_by_MET_with_br_m_prior_treatment_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "br_m_prior_treatment") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "Y",
            log2fc < 0 ~ "N"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC2")

    df_combo_res_MET_prior_treatment <- bind_rows(df_combo_res_MET_prior_treatment, df_masslin_out, df_ancombc_out, df_ancombc2_out)
}

df_combo_plot_MET_prior_treatment <- df_combo_res_MET_prior_treatment %>%
    filter(qval < 0.25) %>%
    mutate(log2fc = round(log2fc, 2)) %>%
    filter(abs(log2fc) > 2) %>%
    write_csv("data/processed_data/df_combo_plot_MET_br_m_prior_treatment_multiple_methods_combined_output_plt.csv")

# Filter for species that are significant with atleast 2 methods within each sample type
df_species_filter_MET_prior_treatment <- df_combo_plot_MET_prior_treatment %>%
    group_by(sample_type, Species) %>%
    summarize(count = n()) %>%
    filter(count > 1)

df_combo_plot_prior_treatment <- df_combo_plot_MET_prior_treatment %>%
    inner_join(df_species_filter_MET_prior_treatment, by = c("sample_type", "Species")) %>%
    mutate(sample_type = factor(sample_type, levels = c("Saliva", "Oral swab", "Stool")))



plot_met_prior_treatment_heatmap <- ggplot(df_combo_plot_prior_treatment, aes(y = Species, x = method, color = groupings)) +
    geom_tile(aes(fill = log2fc), color = "white") +
    scale_fill_gradient2(
        midpoint = 0, limits = c(-11, 11),
        mid = "white", low = "#e8a122", high = "#238443"
    ) +
    geom_text(aes(x = method, y = Species, label = log2fc),
        color = ifelse(abs(df_combo_plot_prior_treatment$log2fc) > 7, "white", "black"),
        fontface = "bold", size = 5
    ) +
    facet_grid(sample_type ~ ., scales = "free_y", space = "free_y") +
    theme_light(base_size = 20) +
    theme(
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    labs(title = "BRM Prior Treatment", x = "Method", y = "Species") +
    labs(caption = "Positive log2fc values indicate higher abundance in Y prior treatment group")
ggsave("output/figures/MET_prior_treatment_SSO.pdf", width = 14, height = 14, dpi = 300)



plot_met_prior_treatment_heatmap_data <- plot_met_prior_treatment_heatmap$data %>%
    mutate(sample_type = factor(sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
    group_by(sample_type) %>%
    distinct(Species) %>%
    arrange(desc(Species), .by_group = T)



plot_met_prior_treatment_abundance_box <- physeq_SSO_clean_species %>%
    transform_sample_counts(., function(x) (x / sum(x)) * 100) %>%
    fn_clean_tax_table(.) %>%
    ps_filter(tumor_category == "Met") %>%
    ps_filter(br_m_prior_treatment %in% c("Y", "N")) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    dplyr::select(microbiome_sample_type, Species, br_m_prior_treatment, Abundance) %>%
    mutate(microbiome_sample_type = factor(microbiome_sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
    inner_join(plot_met_prior_treatment_heatmap_data, by = c("microbiome_sample_type" = "sample_type", "Species")) %>%
    ggplot(., aes(x = Abundance, y = Species, fill = br_m_prior_treatment)) +
    geom_boxplot(linewidth = 0.1) +
    geom_point(aes(group = br_m_prior_treatment), position = position_dodge(width = .75)) +
    scale_fill_manual(values = c("N" = "#e8a122", "Y" = "#238443")) +
    scale_x_log10(labels = scales::comma) +
    facet_grid(microbiome_sample_type ~ ., scales = "free_y") +
    theme_light(base_size = 18) +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    ylab(NULL) +
    ggh4x::force_panelsizes(rows = c(0.85, 0.72, 0.85)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_met_prior_treatment_heatmap + plot_met_prior_treatment_abundance_box +
    plot_layout(guides = "collect", width = c(0.2, 0.3), nrow = 1)
ggsave("output/figures/BRM_MET_prior_treatment_SSO_combo_box.pdf", width = 17, height = 14, dpi = 300)
