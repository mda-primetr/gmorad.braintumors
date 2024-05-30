source("src/libraries.R")
source("src/common_functions.R")

# Load phyloseq data
load("data/processed_data/physeq_SSO_clean.RData")
physeq_SSO_clean

# Metadata
df_metadata <- data.frame(sample_data(physeq_SSO_clean))



physeq_SSO_clean_species <- physeq_SSO_clean %>%
    speedyseq::tax_glom(taxrank = "Species")


df_sample_data_test <- data.frame(sample_data(physeq_SSO_clean)) %>%
    filter(microbiome_sample_type == "Stool") %>%
    dplyr::select(
        sex,
        race,
        steroids_at_surgery,
        tumor_category,
        glioma_time_to_progression,
        glioma_recurrent_status_at_surg,
        grade_category,
        glioma_idh_mutation,
    ) %>%
    mutate(glioma_time_to_progression = as.numeric(glioma_time_to_progression)) %>%
    tbl_summary(by = sex) %>%
    add_p()


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# GLioma by Progression caterogical differential abundance ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Maaslin ---

for (i in c("Stool", "Saliva", "Oral swab")) {
    final_phyloseq_before_rarefaction_species_i <- physeq_SSO_clean_species %>%
        speedyseq::filter_sample_data(microbiome_sample_type == i) %>%
        speedyseq::filter_sample_data(tumor_category == "Glioma") %>%
        speedyseq::filter_sample_data(glioma_progression_overall %in% c("Y", "N")) %>%
        speedyseq::filter_sample_data(glioma_fu_time != "Unknown") %>%
        speedyseq::mutate_sample_data(glioma_fu_time = as.numeric(glioma_fu_time)) %>%
        speedyseq::mutate_sample_data(steroids_at_surgery = factor(steroids_at_surgery, levels = c("N", "Y"))) %>%
        speedyseq::mutate_sample_data(sex = factor(sex, levels = c("F", "M"))) %>%
        speedyseq::mutate_sample_data(grade_category = factor(grade_category, levels = c("Low", "High"))) %>%
        speedyseq::mutate_sample_data(glioma_idh_mutation = factor(glioma_idh_mutation, levels = c("Wild", "Mutant")))


    df_maaslin_all <- data.frame(otu_table(final_phyloseq_before_rarefaction_species_i) %>% t())

    df_maaslin_all_meta <- data.frame(sample_data(final_phyloseq_before_rarefaction_species_i)) %>%
        filter(glioma_progression_overall %in% c("Y", "N"))


    # Please refer output/figures/metadata_stats
    maaslin_by_sample <- Maaslin2(
        input_data = df_maaslin_all,
        input_metadata = df_maaslin_all_meta,
        output = paste0("data/processed_data/maaslin_phyloseq_by_Glioma_categorical_glioma_progression_overall_", i),
        analysis_method = "NEGBIN",
        fixed_effects = c(
            "glioma_progression_overall", "glioma_fu_time",
            "sex", "grade_category", "glioma_idh_mutation", "steroids_at_surgery"
        ),
        plot_scatter = FALSE,
        plot_heatmap = FALSE,
        min_abundance = 100,
        min_prevalence = 0.25,
        cores = 24,
        transform = "NONE",
        reference = c("glioma_progression_overall,N;
                              sex,F;
                              grade_category,Low;
                              glioma_idh_mutation,Wild;
                              steroids_at_surgery,No"),
        normalization = "TMM",
        standardize = TRUE
    )

    # AncomBC version 1----
    ancombc_out <- ancombc(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0.25,
        formula = "glioma_progression_overall + glioma_fu_time + sex + grade_category + glioma_idh_mutation + steroids_at_surgery",
        group = "glioma_progression_overall",
        struc_zero = TRUE,
        neg_lb = TRUE,
        n_cl = 24
    )

    # store AncomBC1 results
    get_ancombc_df(ancombc_out$res, taxa_level = "Species", group_var = "glioma_progression_overall") %>%
        tibble() %>%
        write_tsv(paste0("data/processed_data/ancombc_phyloseq_by_Glioma_categorical_glioma_progression_overall_", i, "_significant_results.tsv"))

    # AncomBC version 2----
    ancombc_out2 <- ancombc2(
        final_phyloseq_before_rarefaction_species_i,
        tax_level = "Species",
        prv_cut = 0.25,
        fix_formula = "glioma_progression_overall + glioma_fu_time + sex + grade_category + glioma_idh_mutation + steroids_at_surgery",
        group = "glioma_progression_overall",
        struc_zero = TRUE,
        pseudo_sens = TRUE,
        neg_lb = TRUE,
        n_cl = 24
    )

    # store AncomBC2 results
    get_ancombc2_df(ancombc_out2$res, taxa_level = "Species", group_var = "glioma_progression_overall") %>%
        write_tsv(paste0("data/processed_data/ancombc2_phyloseq_by_Glioma_categorical_glioma_progression_overall_", i, "_significant_results.tsv"))
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
    df_masslin_out <- read_tsv(paste0("data/processed_data/maaslin_phyloseq_by_Glioma_categorical_glioma_progression_overall_", i, "/significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        mutate(feature_copy = feature) %>%
        mutate(feature = gsub("X", "", feature)) %>%
        filter(metadata == "glioma_progression_overall") %>%
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


    df_ancombc_out <- read_tsv(paste0("data/processed_data/ancombc_phyloseq_by_Glioma_categorical_glioma_progression_overall_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "glioma_progression_overall") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "Y",
            log2fc < 0 ~ "N"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC")


    df_ancombc2_out <- read_tsv(paste0("data/processed_data/ancombc2_phyloseq_by_Glioma_categorical_glioma_progression_overall_", i, "_significant_results.tsv")) %>%
        mutate(sample_type = paste0(i)) %>%
        filter(metadata == "glioma_progression_overall") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "Y",
            log2fc < 0 ~ "N"
        )) %>%
        dplyr::select(Species = feature, log2fc, qval, groupings, sample_type) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        mutate(method = "AncomBC2")

    df_combo_res_prog_cat <- bind_rows(df_combo_res_prog_cat, df_masslin_out, df_ancombc_out, df_ancombc2_out)
}

df_combo_res_prog_cat_plt <- df_combo_res_prog_cat %>%
    filter(abs(log2fc) > 2) %>%
    filter(qval < 0.25) %>%
    write_csv("data/processed_data/df_combo_Glioma_by_progression_categorical_multiple_methods_combined_output_plt.csv")





# Filter for species that are significant with atleast 2 methods within each sample type
df_species_filt_combo_res_prog_cat <- df_combo_res_prog_cat_plt %>%
    group_by(sample_type, Species) %>%
    summarize(count = n()) %>%
    filter(count > 1)

df_combo_plot_prog_cat_glioma <- df_combo_res_prog_cat_plt %>%
    # filter(Species %in% df_species_filt_combo_res_prog_cat$Species) %>%
    inner_join(df_species_filt_combo_res_prog_cat, by = c("Species", "sample_type")) %>%
    mutate(sample_type = factor(sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
    mutate(log2fc = round(log2fc, 2))



plot_glioma_progression_heatmap <- ggplot(df_combo_plot_prog_cat_glioma, aes(y = Species, x = method, color = groupings)) +
    geom_tile(aes(fill = log2fc), color = "white") +
    scale_fill_gradient2(midpoint = 0, mid = "white", low = "#e8a122", high = "#225ea8") +
    geom_text(aes(x = method, y = Species, label = log2fc),
        color = ifelse(abs(df_combo_plot_prog_cat_glioma$log2fc) > 7, "white", "black"),
        fontface = "bold", size = 7
    ) +
    # ggh4x::facet_grid2(sample_type~., scales = "free_y", independent = "y") +
    facet_grid(sample_type ~ ., scales = "free_y", space = "free_y") +
    theme_light(base_size = 20) +
    theme(
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    ) +
    labs(title = "Glioma Progression category", x = "Method", y = "Species") +
    labs(caption = "Positive log2fc values indicate higher abundance in Y progression group")
ggsave("output/figures/Glioma_Progression_category_SSO.pdf", width = 10, height = 14, dpi = 300)


plot_glioma_progression_heatmap_data <- plot_glioma_progression_heatmap$data %>%
    mutate(sample_type = factor(sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
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
    dplyr::select(microbiome_sample_type, Species, glioma_progression_overall, Abundance) %>%
    mutate(microbiome_sample_type = factor(microbiome_sample_type, levels = c("Saliva", "Oral swab", "Stool"))) %>%
    inner_join(plot_glioma_progression_heatmap_data, by = c("microbiome_sample_type" = "sample_type", "Species")) %>%
    ggplot(., aes(x = Abundance, y = Species, fill = glioma_progression_overall)) +
    geom_boxplot(linewidth = 0.1) +
    geom_point(aes(group = glioma_progression_overall), position = position_dodge(width = .75)) +
    scale_fill_manual(values = c("N" = "#e8a122", "Y" = "#225ea8")) +
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
    ggh4x::force_panelsizes(rows = c(0.8, 0.12, 0.9)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())




plot_glioma_progression_heatmap + plot_glioma_progression_abundance_box +
    plot_layout(guides = "collect", width = c(0.2, 0.3), nrow = 1)
ggsave("output/figures/Glioma_progression_SSO_combo_box.pdf", width = 17, height = 14, dpi = 300)





# Check progression outcome for taxa that are significantly different between the groups
# In at atleast 2 methods with KM curves:



# Kaplien Meier plot for PFS ----
res_list <- list()
for (i in unique(plot_glioma_progression_heatmap_data$Species)) {
    for (j in c("Saliva", "Oral swab", "Stool")) {
        tryCatch(
            {
                df_meta_species_glioma <- physeq_SSO_clean_species %>%
                    ps_filter(microbiome_sample_type == j) %>%
                    ps_filter(tumor_category == "Glioma") %>%
                    fn_clean_tax_table(.) %>%
                    transform_sample_counts(., function(x) (x / sum(x)) * 100) %>%
                    psmelt() %>%
                    filter(Species == i) %>%
                    filter(glioma_progression_overall %in% c("Y", "N")) %>%
                    mutate(glioma_progression_overall = case_when(
                        glioma_progression_overall == "Y" ~ 1,
                        glioma_progression_overall == "N" ~ 0
                    )) %>%
                    mutate(glioma_time_to_progression = case_when(
                        is.na(glioma_time_to_progression) ~ as.numeric(glioma_fu_time),
                        TRUE ~ as.numeric(glioma_time_to_progression)
                    )) %>%
                    mutate(med_species = ifelse(Abundance > median(Abundance), "High", "Low"))

                # Get the survival fit ----
                res_list$fit[[paste0(i, "_", j)]] <- survival::survfit(survival::Surv(glioma_time_to_progression, glioma_progression_overall) ~ (med_species), data = df_meta_species_glioma)

                # Get the p value ----
                res_list$pval[[paste0(i, "_", j)]] <- survminer::surv_pvalue(res_list$fit[[paste0(i, "_", j)]])

                # Get the median survival values ----
                res_list$median_survival[[paste0(i, "_", j)]] <- summary(res_list$fit[[paste0(i, "_", j)]])$table %>%
                    data.frame() %>%
                    rownames_to_column(var = "category") %>%
                    dplyr::select(category, median) %>%
                    pivot_wider(names_from = category, values_from = median)



                res_list$plot[[paste0(i, "_", j)]] <- eval(survminer::ggsurvplot(
                    res_list$fit[[paste0(i, "_", j)]],
                    pval = T,
                    conf.int = T,
                    data = df_meta_species_glioma,
                    risk.table = TRUE,
                    fontsize = 8,
                    ggtheme = theme_classic(base_size = 25),
                    title = paste0("PFS: ", i, " in ", j),
                    palette = c("#009E73", "#D55E00")
                ))
            },
            error = function(e) {
                print(paste0(i, " in ", j, " failed"))
            }
        )
    }
}



# Get median and P values
df_surv_median <- enframe(res_list$median_survival) %>%
    unnest(value)

df_surv_pval <- enframe(res_list$pval) %>%
    unnest(value) %>%
    inner_join(df_surv_median, by = c("name" = "name")) %>%
    separate(name, c("species", "microbiome_sample_type"), sep = "_") %>%
    dplyr::select(-pval.txt)






# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Heatmap using the taxa
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Heatmap using the taxa



for (i in c("Stool", "Saliva", "Oral swab")) {
    # Get dataframe used for the ggplot2 plot
    df_combo_species_sample_i <- df_combo_plot_prog_cat_glioma %>%
        filter(sample_type == i)


    physeq_SSO_taxa_clean_species <- fn_clean_tax_table(physeq_SSO_clean_species) %>%
        ps_filter(microbiome_sample_type == i) %>%
        ps_filter(tumor_category == "Glioma")


    physeq_clr <- physeq_SSO_taxa_clean_species %>%
        ps_filter(microbiome_sample_type == i) %>%
        ps_filter(tumor_category == "Glioma") %>%
        microbiome::transform(., transform = "clr") %>%
        subset_taxa(., Species %in% unique(df_combo_species_sample_i$Species)) %>%
        speedyseq::filter_sample_data(glioma_progression_overall %in% c("Y", "N"))



    # Get the number of samples required for each feature to be on the heatmap

    # Taxonomy table
    df_tax <- data.frame(tax_table(physeq_clr)) %>%
        rownames_to_column(var = "tax_id")

    df_mat_complex_heatmap <- data.frame(otu_table(physeq_clr), check.names = FALSE) %>%
        rownames_to_column(var = "tax_id") %>%
        inner_join(df_tax, by = "tax_id") %>%
        inner_join(df_combo_species_sample_i, by = "Species") %>%
        mutate(groupings = case_when(
            log2fc > 0 ~ "Mutant",
            log2fc < 0 ~ "Wild"
        )) %>%
        arrange(log2fc)

    df_mat_complex_heatmap_methods <- df_mat_complex_heatmap %>%
        dplyr::select(Species, log2fc, method) %>%
        pivot_wider(names_from = method, values_from = c(log2fc), values_fill = 0) %>%
        column_to_rownames(var = "Species")


    # Heatmap data
    mat_complex_heatmap <- df_mat_complex_heatmap %>%
        dplyr::select(-c(tax_id, Kingdom:Genus, log2fc:count)) %>%
        distinct(Species, .keep_all = TRUE) %>%
        column_to_rownames(var = "Species")


    # phyloseq metadata
    df_metdata <- data.frame(sample_data(physeq_clr))


    # Column annotation
    ha_column <- HeatmapAnnotation(
        df = data.frame(
            progression = df_metdata$glioma_progression_overall,
            recurrence = df_metdata$glioma_recurrent_status_at_surg,
            glioma_idh_mutation = df_metdata$glioma_idh_mutation,
            batch = df_metdata$shotgun_seq_batch,
            sex = df_metdata$sex,
            grade = df_metdata$grade_category
        ),
        col = list(
            progression = c("Y" = "red", "N" = "blue", "Unknown" = "gray"),
            recurrence = c("Y" = "black", "N" = "gray60", "Unknown" = "gray"),
            glioma_idh_mutation = c("Mutant" = "red", "Wild" = "blue", "Unknown" = "gray"),
            sex = c("F" = "#E41A1C", "M" = "#377EB8", "Unknown" = "gray"),
            batch = c("Batch 1" = "red", "Batch 2" = "blue", "Batch 6" = "orange"),
            grade = c("High" = "orange", Low = "#2aa5a1", "Unknown" = "gray")
        ),
        annotation_name_side = "left"
    )

    # Get Prevalence ----
    df_prev <- get_phyloseq_prev_by_group(physeq_SSO_taxa_clean_species, "glioma_progression_overall", "Species") %>%
        dplyr::select(glioma_progression_overall, Species, prevalence) %>%
        mutate(prevalence = round(prevalence * 100, digits = 2)) %>% # Convert to percentage
        filter(glioma_progression_overall %in% c("Y", "N")) %>%
        pivot_wider(names_from = glioma_progression_overall, values_from = prevalence, values_fill = 0)


    df_prev_plot <- df_prev %>%
        filter(Species %in% c(df_mat_complex_heatmap$Species)) %>%
        mutate(Species = factor(Species)) %>%
        mutate(Species = factor(Species, levels = rownames(mat_complex_heatmap))) %>%
        arrange(Species) %>%
        column_to_rownames(var = "Species")



    # Check differential methods values


    # Check for ANCOMBC values

    if ("AncomBC" %in% names(df_mat_complex_heatmap_methods)) {
        ancombc_data <- df_mat_complex_heatmap_methods$AncomBC
    } else {
        ancombc_data <- rep(0, length(unique(df_mat_complex_heatmap$Species)))
    }

    # Check for ANCOMBC2 values

    if ("AncomBC2" %in% names(df_mat_complex_heatmap_methods)) {
        ancombc2_data <- df_mat_complex_heatmap_methods$AncomBC2
    } else {
        ancombc2_data <- rep(0, length(unique(df_mat_complex_heatmap$Species)))
    }

    # Check for Maaslin2 values

    if ("Maaslin2" %in% names(df_mat_complex_heatmap_methods)) {
        maaslin2_data <- df_mat_complex_heatmap_methods$Maaslin2
    } else {
        maaslin2_data <- rep(0, length(unique(df_mat_complex_heatmap$Species)))
    }



    row_ha2 <- HeatmapAnnotation(
        which = "row",
        AncomBC = anno_simple(
            ancombc_data,
            width = unit(15, "mm"),
            pch = as.character(
                ancombc_data
            ),
            pt_gp = gpar(col = "white"), pt_size = unit(6, "mm"),
            col = colorRamp2(c(-10, 0, 10), c("navy", "white", "firebrick3"))
        ),
        AncomBC2 = anno_simple(
            ancombc2_data,
            width = unit(15, "mm"),
            pch = as.character(ancombc2_data),
            pt_gp = gpar(col = "white"), pt_size = unit(6, "mm"),
            col = colorRamp2(c(-10, 0, 10), c("navy", "white", "firebrick3"))
        ),
        Maaslin2 = anno_simple(
            maaslin2_data,
            width = unit(15, "mm"),
            pch = as.character(maaslin2_data),
            pt_gp = gpar(col = "white"), pt_size = unit(6, "mm"),
            col = colorRamp2(c(-5, 0, 5), c("navy", "white", "firebrick3"))
        ),
        `Prevalence (%)` = anno_barplot(df_prev_plot,
            width = unit(15, "mm"),
            beside = TRUE, attach = TRUE,
            gp = gpar(fill = c("N" = "blue", "Y" = "red"), col = "black"),
            bar_width = 0.3, height = unit(10, "cm")
        ),
        # log2fc = anno_text(
        #       df_mat_complex_heatmap_methods$Maaslin2,
        #       rot = 90,
        #       gp = gpar(fontsize = 10, fontface = "bold", col = "blue")
        # ),


        gap = unit(c(1, 2, 1, 4, 6, 7), "mm"),
        annotation_name_side = "top"
    )



    pdf(file = paste0("output/figures/heatmap/differential_glioma_progression_category_", i, ".pdf"), width = 28, height = 8)
    print(
        Heatmap(mat_complex_heatmap,
            rect_gp = gpar(col = "white", lwd = 1.5),
            col = colorRampPalette(c("#542788", "white", "#b35806"))(50),
            name = "Sample\nAbundances",
            row_names_side = "left",
            cluster_column_slices = T,
            column_gap = unit(5, "mm"),
            cluster_rows = FALSE,
            show_column_names = F,
            top_annotation = ha_column,
            right_annotation = row_ha2,
            column_split = factor(df_metdata$glioma_progression_overall, levels = c("Y", "N")),
            width = ncol(mat_complex_heatmap) * unit(5, "mm"),
            row_names_gp = grid::gpar(fontsize = 20, fontface = "italic")
        )
    )
    dev.off()
}
