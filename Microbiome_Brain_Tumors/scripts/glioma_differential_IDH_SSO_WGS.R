source("src/libraries.R")
source("src/common_functions.R")

# Load Unrarefied data
load("data/processed_data/physeq_SSO_clean.RData")
physeq_SSO_clean

# Metadata
df_metadata <- data.frame(sample_data(physeq_SSO_clean))


physeq_SSO_clean_species <- physeq_SSO_clean %>%
    speedyseq::tax_glom(taxrank = "Species")

# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Glioma IDH mutation controlling for anything significant ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Differential abundance analysis using Maaslin2 ----
# For Glioma
# List for controls is based on the metadata summary wilcox test/fisher test
for (i in c("Stool", "Saliva", "Oral swab")) {
    for (j in c("glioma_idh_mutation", "grade_category")) {
        final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
            speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
            speedyseq::filter_sample_data(tumor_category == "Glioma")

        df_maaslin_all <- data.frame(otu_table(final_phyloseq_before_rarefaction_species_i) %>% t())

        if (j == "glioma_idh_mutation") {
            df_maaslin_all_meta <- data.frame(sample_data(final_phyloseq_before_rarefaction_species_i)) %>%
                mutate(glioma_idh_mutation = factor(glioma_idh_mutation, levels = c("Mutant", "Wild"))) %>%
                filter(glioma_idh_mutation %in% c("Mutant", "Wild")) %>%
                mutate(steroids_at_surgery = case_when(
                    steroids_at_surgery == "Unknown" ~ NA,
                    TRUE ~ steroids_at_surgery
                )) %>%
                mutate(grade_category = case_when(
                    grade_category == "Unknown" ~ NA,
                    TRUE ~ grade_category
                ))

            # Please refer output/figures/metadata_stats
            maaslin_by_sample <- Maaslin2(
                input_data = df_maaslin_all,
                input_metadata = df_maaslin_all_meta,
                output = paste0("data/processed_data/maaslin_phyloseq_by_Glioma_category_with_", j, "_", i),
                analysis_method = "NEGBIN",
                fixed_effects = c(
                    "glioma_idh_mutation",
                    "sex",
                    "age_at_surgery",
                    "steroids_at_surgery",
                    "grade_category"
                ),
                plot_scatter = FALSE,
                plot_heatmap = FALSE,
                min_abundance = 100,
                min_prevalence = 0.25,
                cores = 24,
                transform = "NONE",
                reference = c("glioma_idh_mutation,Wild;
                              sex,F;
                              steroids_at_surgery,No;
                              grade_category,Low"),
                normalization = "TMM",
                standardize = TRUE
            )
        }

        if (j == "grade_category") {
            df_maaslin_all_meta <- data.frame(sample_data(final_phyloseq_before_rarefaction_species_i)) %>%
                mutate(grade_category = factor(grade_category, levels = c("Low", "High"))) %>%
                filter(grade_category %in% c("Low", "High"))

            # Please refer output/figures/metadata_stats
            maaslin_by_sample <- Maaslin2(
                input_data = df_maaslin_all,
                input_metadata = df_maaslin_all_meta,
                output = paste0("data/processed_data/maaslin_phyloseq_by_Glioma_category_with_", j, "_", i),
                analysis_method = "NEGBIN",
                fixed_effects = c("grade_category"),
                plot_scatter = FALSE,
                plot_heatmap = FALSE,
                min_abundance = 100,
                min_prevalence = 0.25,
                cores = 24,
                transform = "NONE",
                reference = c("grade_category,Low"),
                normalization = "TMM",
                standardize = TRUE
            )
        }
    }
}


# Differential abundance analysis using AncomBC 1 and 2 ----

for (i in c("Stool", "Saliva", "Oral swab")) {
    for (j in c("glioma_idh_mutation", "grade_category")) {
        if (j == "glioma_idh_mutation") {
            final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
                speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
                speedyseq::filter_sample_data(tumor_category == "Glioma") %>%
                speedyseq::filter_sample_data(glioma_idh_mutation %in% c("Mutant", "Wild")) %>%
                speedyseq::mutate_sample_data(glioma_idh_mutation = factor(glioma_idh_mutation, levels = c("Wild", "Mutant")))


            ancombc_out <- ancombc(
                final_phyloseq_before_rarefaction_species_i,
                tax_level = "Species",
                prv_cut = 0.25,
                formula = "glioma_idh_mutation + sex + age_at_surgery + steroids_at_surgery + grade_category",
                group = "glioma_idh_mutation",
                struc_zero = TRUE,
                neg_lb = TRUE,
                n_cl = 24
            )

            ancombc_out2 <- ancombc2(
                final_phyloseq_before_rarefaction_species_i,
                tax_level = "Species",
                prv_cut = 0.25,
                fix_formula = "glioma_idh_mutation + sex + age_at_surgery + steroids_at_surgery + grade_category",
                group = "glioma_idh_mutation",
                struc_zero = TRUE,
                pseudo_sens = TRUE,
                neg_lb = TRUE,
                n_cl = 24
            )

            get_ancombc2_df(ancombc_out2$res, taxa_level = "Species", group_var = "glioma_idh_mutation", group_var_order = c("Wild", "Mutant")) %>%
                write_tsv(paste0("data/processed_data/ancombc2_phyloseq_by_Glioma_category_with_", j, "_", i, "_significant_results.tsv"))

            get_ancombc_df(ancombc_out$res, taxa_level = "Species", group_var = "glioma_idh_mutation", group_var_order = c("Wild", "Mutant")) %>%
                tibble() %>%
                write_tsv(paste0("data/processed_data/ancombc_phyloseq_by_Glioma_category_with_", j, "_", i, "_significant_results.tsv"))
        }


        if (j == "grade_category") {
            final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
                speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
                speedyseq::filter_sample_data(tumor_category == "Glioma") %>%
                speedyseq::filter_sample_data(grade_category %in% c("High", "Low")) %>%
                speedyseq::mutate_sample_data(grade_category = factor(grade_category, levels = c("Low", "High")))


            ancombc_out <- ancombc(
                final_phyloseq_before_rarefaction_species_i,
                tax_level = "Species",
                prv_cut = 0.25,
                formula = "grade_category",
                group = "grade_category",
                struc_zero = TRUE,
                neg_lb = TRUE,
                n_cl = 24
            )

            ancombc_out2 <- ancombc2(
                final_phyloseq_before_rarefaction_species_i,
                tax_level = "Species",
                prv_cut = 0.25,
                fix_formula = "grade_category",
                group = "grade_category",
                struc_zero = TRUE,
                pseudo_sens = TRUE,
                neg_lb = TRUE,
                n_cl = 24
            )


            get_ancombc2_df(ancombc_out2$res, taxa_level = "Species", group_var = "grade_category", group_var_order = c("Low", "High")) %>%
                write_tsv(paste0("data/processed_data/ancombc2_phyloseq_by_Glioma_category_with_", j, "_", i, "_significant_results.tsv"))

            get_ancombc_df(ancombc_out$res, taxa_level = "Species", group_var = "grade_category", group_var_order = c("Low", "High")) %>%
                tibble() %>%
                write_tsv(paste0("data/processed_data/ancombc_phyloseq_by_Glioma_category_with_", j, "_", i, "_significant_results.tsv"))
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



df_combo_res_IDH <- data.frame()
for (i in c("Stool", "Saliva", "Oral swab")) {
    df_masslin_out <- read_tsv(paste0("data/processed_data/maaslin_phyloseq_by_Glioma_category_with_glioma_idh_mutation_", i, "/significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        mutate(feature_copy = feature) %>%
        mutate(feature = gsub("X", "", feature)) %>%
        filter(metadata == "glioma_idh_mutation") %>%
        mutate(log2fc = log2(exp(coef))) %>%
        mutate(groupings = case_when(
            coef > 0 ~ "Mutant",
            coef < 0 ~ "Wild"
        )) %>%
        filter(`N.not.0` >= N * 0.25) %>%
        filter(abs(log2fc) > 2) %>%
        inner_join(df_taxa_table, by = c("feature" = "tax_id")) %>%
        dplyr::select(Species = Species, log2fc, qval, groupings, sample_type) %>%
        mutate(method = "Maaslin2")


    df_ancombc_out <- read_tsv(paste0("data/processed_data/ancombc_phyloseq_by_Glioma_category_with_glioma_idh_mutation_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "glioma_idh_mutation") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "Mutant",
            log2fc < 0 ~ "Wild"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC")


    df_ancombc2_out <- read_tsv(paste0("data/processed_data/ancombc2_phyloseq_by_Glioma_category_with_glioma_idh_mutation_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "glioma_idh_mutation") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "Mutant",
            log2fc < 0 ~ "Wild"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC2")

    df_combo_res_IDH <- bind_rows(df_combo_res_IDH, df_masslin_out, df_ancombc_out, df_ancombc2_out)
}

df_combo_plot_IDH <- df_combo_res_IDH %>%
    filter(qval < 0.25) %>%
    mutate(log2fc = round(log2fc, 2)) %>%
    filter(abs(log2fc) > 2) %>%
    mutate(qval = format.pval(qval, digits = 2)) %>%
    write_csv("data/processed_data/df_combo_GBM_by_idh_mutation_multiple_methods_combined_output_plt.csv")

# Filter for species that are significant with atleast 2 methods within each sample type
df_species_filter_idh <- df_combo_plot_IDH %>%
    group_by(sample_type, Species) %>%
    summarize(count = n()) %>%
    filter(count > 1) %>%
    ungroup()

df_combo_plot_idh <- df_combo_plot_IDH %>%
    inner_join(df_species_filter_idh, by = c("sample_type", "Species")) %>%
    mutate(sample_type = factor(sample_type, levels = c("Saliva", "Oral swab", "Stool")))

plt_combo_plot_idh <- ggplot(df_combo_plot_idh, aes(y = Species, x = method, color = groupings)) +
    geom_tile(aes(fill = log2fc), color = "white") +
    scale_fill_gradient2(midpoint = 0, mid = "white", low = "#e8a122", high = "#225ea8") +
    geom_text(aes(x = method, y = Species, label = log2fc),
        color = ifelse(abs(df_combo_plot_idh$log2fc) > 7, "white", "black"),
        fontface = "bold", size = 7
    ) +
    facet_grid(sample_type ~ ., scales = "free_y", space = "free_y") +
    theme_light(base_size = 20) +
    theme(
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    labs(title = "Glioma IDH mutation", x = "Method", y = "Species") +
    labs(caption = "Positive log2fc values indicate higher abundance in Mutant samples")
ggsave("output/figures/Glioma_IDH_mutation_SSO.pdf", width = 15, height = 13, dpi = 300)





box_plt_data <- plt_combo_plot_idh$data %>%
    mutate(sample_type = factor(sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
    group_by(sample_type) %>%
    distinct(Species) %>%
    arrange(desc(Species), .by_group = T)


plot_combo_abundance_box <- physeq_SSO_clean_species %>%
    transform_sample_counts(., function(x) (x / sum(x)) * 100) %>%
    fn_clean_tax_table(.) %>%
    ps_filter(tumor_category == "Glioma") %>%
    ps_filter(glioma_idh_mutation %in% c("Mutant", "Wild")) %>%
    ps_mutate(glioma_idh_mutation = factor(glioma_idh_mutation, levels = c("Wild", "Mutant"))) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    dplyr::select(microbiome_sample_type, Species, glioma_idh_mutation, Abundance) %>%
    mutate(microbiome_sample_type = factor(microbiome_sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
    inner_join(box_plt_data, by = c("microbiome_sample_type" = "sample_type", "Species")) %>%
    ggplot(., aes(x = Abundance, y = Species, fill = glioma_idh_mutation)) +
    geom_boxplot(linewidth = 0.1) +
    geom_point(aes(group = glioma_idh_mutation), position = position_dodge(width = .75)) +
    scale_fill_manual(values = c("Wild" = "#e8a122", "Mutant" = "#225ea8")) +
    scale_x_log10(labels = scales::comma) +
    facet_grid(microbiome_sample_type ~ ., scales = "free_y") +
    theme_light(base_size = 18) +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        # strip.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    ylab(NULL) +
    ggh4x::force_panelsizes(rows = c(0.9, 1.1, 0.9)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("output/figures/Glioma_IDH_mutation_SSO_boxplot.pdf", width = 15, height = 13, dpi = 300)




plt_combo_plot_idh + plot_combo_abundance_box +
    plot_layout(guides = "collect", width = c(0.2, 0.3), nrow = 1)
ggsave("output/figures/Glioma_IDH_mutation_SSO_combo_box.pdf", width = 17, height = 14, dpi = 300)
