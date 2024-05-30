source("src/libraries.R")
source("src/common_functions.R")

load("data/processed_data/physeq_16S_tumor.RData")
load("data/processed_data/physeq_negative_controls.RData")
load("data/processed_data/physeq_samples.RData")




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 1 Apply  proportion cutoffs at the ASV level to the samples
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Load the controls prevalence and proportion data cutoff values ----
df_filt <- read_csv("data/processed_data/otu_cutoffs.csv")


physeq_filt_1 <- physeq_samples %>%
    speedyseq::filter_sample_data(., sample_type %in% c(
        "sample"
    )) %>%
    transform_sample_counts(., function(x) x / sum(x) * 100) %>%
    filter_taxa(., function(x) sum(x) > 0, TRUE) %>%
    filter_taxa(
        .,
        function(x) sum(x) < df_filt$asv_prop_cutoff[[1]], TRUE
    ) %>%
    tax_fix(unknowns = c("g__", "metagenome")) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)

# Export filtered phyloseq object ----
save(physeq_filt_1, file = "data/processed_data/physeq_filt_1.RData")

# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 2 filter by prevalence at the genus level ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Find the prevalence of all the genera in the samples

df_genus_prev <- physeq_filt_1 %>%
    speedyseq::tax_glom(., taxrank = "Genus") %>%
    speedyseq::psmelt() %>%
    dplyr::select(Sample, Genus, Abundance) %>%
    filter(Abundance > 0) %>%
    group_by(Genus) %>%
    summarize(n = n()) %>%
    mutate(prevalence = (n / phyloseq::nsamples(physeq_filt_1)) * 100) %>%
    filter(prevalence > 0) %>%
    filter(prevalence < df_filt$genus_prev_cutoff[[1]])

physeq_filt_2 <- physeq_filt_1 %>%
    subset_taxa(., Genus %in% df_genus_prev$Genus) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)


# Export filtered phyloseq object ----
save(physeq_filt_2, file = "data/processed_data/physeq_filt_2.RData")



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 3 Remove batch unique ASVs from the samples
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_taxa_found_in_more_than_two_batch <- physeq_filt_2 %>%
    psmelt() %>%
    dplyr::select(Genus, x16s_seq_batch, Abundance) %>%
    filter(Abundance > 0) %>%
    group_by(Genus, x16s_seq_batch) %>%
    summarize(n = n()) %>%
    group_by(Genus) %>%
    summarize(n = n()) %>%
    filter(n >= 2)


physeq_filt_3 <- physeq_filt_2 %>%
    subset_taxa(., Genus %in% df_taxa_found_in_more_than_two_batch$Genus) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)

# Export filtered phyloseq object ----
save(physeq_filt_3, file = "data/processed_data/physeq_filt_3.RData")


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 4 Now use scrub
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_asv_filt3 <- psmelt(physeq_filt_3) %>%
    filter(Abundance > 0) %>%
    dplyr::select(Sample, OTU, Abundance) %>%
    distinct(Sample, OTU)

## -------------------------------  Batch 3  -------------------------#
# Get ASVs from the filtered phyloseq object in the previous step
physeq_batch_3 <- physeq_samples %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 3") %>%
    prune_taxa(df_asv_filt3$OTU, .) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_samples(sample_names(.) %in% df_asv_filt3$Sample, .) %>%
    speedyseq::filter_sample_data(!is.na(sample_well))

physeq_neg_ctrl_batch_3 <- physeq_negative_controls %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 3") %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    speedyseq::filter_sample_data(!is.na(sample_well))


# Store list of negative controls ASVs found in Batch 3
list_of_neg_ctrls_batch_3_ASVs <- rownames(tax_table(physeq_neg_ctrl_batch_3))



scrb_out_res_B3 <- fn_run_scrub(physeq_batch_3, physeq_neg_ctrl_batch_3, "Batch 3", sample_names(physeq_filt_3))

physeq_scrubbed_b3 <- fn_update_physeq_with_scrub_results(physeq_batch_3, scrb_out_res_B3) %>%
    speedyseq::filter_sample_data(., sample_type %in% c(
        "sample"
    )) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)


physeq_scrubbed_b3_plus_nc_rem <- prune_taxa(setdiff(
    rownames(tax_table(physeq_scrubbed_b3)),
    list_of_neg_ctrls_batch_3_ASVs
), physeq_scrubbed_b3)

sample_sums(physeq_scrubbed_b3_plus_nc_rem)


## -------------------------------  Batch 5  -------------------------#
physeq_batch_5 <- physeq_samples %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 5") %>%
    prune_taxa(df_asv_filt3$OTU, .) %>%
    prune_samples(sample_names(.) %in% df_asv_filt3$Sample, .) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .)


# Since these are sample controls, we need to explicitly set them as controls in the metadata
physeq_neg_ctrl_batch_5 <- physeq_negative_controls %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 5") %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    speedyseq::mutate_sample_data(is_control = TRUE) %>%
    speedyseq::mutate_sample_data(sample_type_for_well = "sample control") %>%
    speedyseq::filter_sample_data(!is.na(sample_well))


# Store list of sample controls ASVs found in Batch 5 (There are no negative controls in Batch 5)
list_of_sample_control_batch_5_ASVs <- rownames(tax_table(physeq_neg_ctrl_batch_5))


scrb_out_res_B5 <- fn_run_scrub(physeq_batch_5, physeq_neg_ctrl_batch_5, "Batch 5", sample_names(physeq_filt_3))

physeq_scrubbed_b5 <- fn_update_physeq_with_scrub_results(physeq_batch_5, scrb_out_res_B5) %>%
    speedyseq::filter_sample_data(., sample_type %in% c(
        "sample"
    )) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)


## -------------------------------  Batch 4  -------------------------#

# Run SCRUB on Batch 4 ----
physeq_batch_4 <- physeq_samples %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 4") %>%
    speedyseq::filter_sample_data(!is.na(sample_well)) %>%
    prune_taxa(df_asv_filt3$OTU, .) %>%
    prune_samples(sample_names(.) %in% df_asv_filt3$Sample, .) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .)


physeq_neg_ctrl_batch_4 <- physeq_negative_controls %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 4") %>%
    speedyseq::filter_sample_data(!is.na(sample_well)) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .)

# Store list of negative controls found in Batch 4
list_of_negative_control_batch_4_ASVs <- rownames(tax_table(physeq_neg_ctrl_batch_4))



scrb_out_res_B4 <- fn_run_scrub(physeq_batch_4, physeq_neg_ctrl_batch_4, "Batch 4", sample_names(physeq_filt_3))

physeq_scrubbed_b4 <- fn_update_physeq_with_scrub_results(physeq_batch_4, scrb_out_res_B4) %>%
    speedyseq::filter_sample_data(., sample_type %in% c(
        "sample"
    )) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)


physeq_scrubbed_b4_plus_nc_rem <- prune_taxa(
    setdiff(
        rownames(tax_table(physeq_scrubbed_b4)),
        list_of_negative_control_batch_4_ASVs
    ),
    physeq_scrubbed_b4
)


## -------------------------------  Batch 9 -------------------------#
# SCRUB on Batch 9 ----
physeq_batch_9 <- physeq_samples %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 9") %>%
    speedyseq::filter_sample_data(!is.na(sample_well)) %>%
    prune_taxa(df_asv_filt3$OTU, .) %>%
    prune_samples(sample_names(.) %in% df_asv_filt3$Sample, .) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .)


physeq_neg_ctrl_batch_9 <- physeq_negative_controls %>%
    speedyseq::filter_sample_data(x16s_seq_batch == "Batch 9") %>%
    speedyseq::filter_sample_data(!is.na(sample_well)) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .)

# Store list of negative controls found in Batch 9
list_of_negative_control_batch_9_ASVs <- rownames(tax_table(physeq_neg_ctrl_batch_9))



scrb_out_res_B9 <- fn_run_scrub(physeq_batch_9, physeq_neg_ctrl_batch_9, "Batch 9", sample_names(physeq_filt_3))

physeq_scrubbed_b9 <- fn_update_physeq_with_scrub_results(physeq_batch_9, scrb_out_res_B9) %>%
    speedyseq::filter_sample_data(., sample_type %in% c(
        "sample"
    )) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_taxa(setdiff(rownames(tax_table(.)), list_of_negative_control_batch_9_ASVs), .)


physeq_scrubbed_b9_plus_nc_rem <- prune_taxa(
    setdiff(
        rownames(tax_table(physeq_scrubbed_b9)),
        list_of_negative_control_batch_9_ASVs
    ),
    physeq_scrubbed_b9
)


# Save scrubbed phyloseq object ----
physeq_all_scrubbed <- rbind(
    psmelt(physeq_scrubbed_b3),
    psmelt(physeq_scrubbed_b4),
    psmelt(physeq_scrubbed_b5),
    psmelt(physeq_scrubbed_b9)
) %>%
    melt_to_physeq() %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)


# Save the scrubbed phyloseq object ----
save(physeq_all_scrubbed, file = "data/processed_data/physeq_all_scrubbed.RData")

# Save scrubbed and negative control removed phyloseq objects ----
physeq_scrubbed_plus_nc_rem <- rbind(
    psmelt(physeq_scrubbed_b3_plus_nc_rem),
    psmelt(physeq_scrubbed_b4_plus_nc_rem),
    psmelt(physeq_scrubbed_b5), # No negative controls in Batch 5
    psmelt(physeq_scrubbed_b9_plus_nc_rem)
) %>%
    melt_to_physeq() %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)

save(physeq_scrubbed_plus_nc_rem, file = "data/processed_data/physeq_scrubbed_plus_nc_rem.RData")



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Bubble plot of genera after filtering steps
# ────────────────────────────────────────────────────────────────────────────────────────────────────

load("data/processed_data/physeq_all_scrubbed.RData")
physeq_all_scrubbed %>%
    speedyseq::tax_glom(taxrank = "Genus") %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    dplyr::select(Genus, Abundance) %>%
    mutate(Genus = ifelse(Genus == "67-14_unknown", "LKT__Solirubrobacterales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Acidimicrobiia_unknown", "LKT__Acidimicrobiia", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Saccharimonadales_unknown", "LKT__Saccharimonadales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Sphingomonadaceae_unknown", "LKT__Sphingomonadaceae", Genus)) %>%
    group_by(Genus) %>%
    summarize(n = n())


# Get all unique genera ----
all_genera_scrubbed <- psmelt(physeq_all_scrubbed) %>%
    mutate(Genus = ifelse(Genus == "67-14_unknown", "LKT__Solirubrobacterales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Acidimicrobiia_unknown", "LKT__Acidimicrobiia", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Saccharimonadales_unknown", "LKT__Saccharimonadales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Sphingomonadaceae_unknown", "LKT__Sphingomonadaceae", Genus)) %>%
    distinct(Genus)


# Genera that needs to be grayed out ----
genera_order_to_be_grayed <- c(
    "Sphingomonas",
    "Nocardioides",
    "Achromobacter",
    "Aquabacterium",
    "Ferruginibacter",
    "Bacillus",
    "Paenibacillus",
    "Flavobacterium",
    "Microbacterium",
    "Acidimicrobiia_unknown", # This is another name for LKT__Acidimicrobiia
    "LKT__Acidimicrobiia",
    "Devosia",
    "Leuconostoc",
    "Variovorax",
    "Saccharimonadales_unknown", # This is another name for "LKT__Saccharimonadales"
    "LKT__Saccharimonadales",
    "Sphingobacterium",
    "Sphingomonadaceae_unknown", # This is another name for "LKT__Sphingomonadaceae"
    "LKT__Sphingomonadaceae",
    "Kocuria",
    "Rhizobacter",
    "67-14_unknown", # This is another name for LKT__Solirubrobacterales
    "LKT__Solirubrobacterales" # This is from Solirubrobacterales family
)

# Genera that needs to be colored ----
genera_to_be_colored <- setdiff(all_genera_scrubbed$Genus, genera_order_to_be_grayed)

# Get color palette for the genera to be colored ----
color_for_other_genera <- distinct_palette()[1:length(genera_to_be_colored)]

# Set colors for the genera ----
genera_colours <- setNames(
    c(rep("gray", length(genera_order_to_be_grayed)), color_for_other_genera),
    c(genera_order_to_be_grayed, genera_to_be_colored)
)



# Bubble plot  -----
physeq_all_scrubbed %>%
    speedyseq::tax_glom(taxrank = "Genus") %>%
    psmelt() %>%
    filter(!is.na(tumor_category)) %>%
    mutate(tumor_category = case_when(
        tumor_category == "Glioma" ~ "Glioma",
        tumor_category == "Met" ~ "BrM",
        tumor_category == "No tumor" ~ "Normal Brain"
    )) %>%
    mutate(tumor_category = factor(tumor_category, levels = c("Glioma", "BrM", "Normal Brain"))) %>%
    filter(Abundance > 0) %>%
    dplyr::select(Genus, Abundance, dataset, tumor_category) %>%
    mutate(Genus = ifelse(Genus == "67-14_unknown", "LKT__Solirubrobacterales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Acidimicrobiia_unknown", "LKT__Acidimicrobiia", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Saccharimonadales_unknown", "LKT__Saccharimonadales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Sphingomonadaceae_unknown", "LKT__Sphingomonadaceae", Genus)) %>%
    group_by(Genus, dataset, tumor_category) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    group_by(dataset, tumor_category) %>%
    mutate(prop = prop.table(Abundance) * 100) %T>%
    write_csv("data/processed_data/bubble_plot_figure_data_abundance.csv") %>%
    mutate(Genus = factor(Genus, levels = c(genera_order_to_be_grayed, genera_to_be_colored))) %>%
    dplyr::select(Genus, prop, dataset, tumor_category) %>%
    mutate(Genus_cat = ifelse(Genus %in% genera_order_to_be_grayed, "Grayed", "Normal")) %>%
    ggplot(., aes(x = dataset, y = Genus, fill = prop, color = Genus)) +
    geom_point(aes(size = prop), alpha = 0.9) +
    facet_grid(. ~ tumor_category, scales = "free_x", space = "free") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "bold")
    ) +
    scale_color_manual(values = genera_colours) +
    theme(text = element_text(size = 20)) +
    ggtitle("Bubble plot of Genera after scrubbing") +
    xlab("Cohort") +
    ylab("Genera") +
    scale_size_continuous(range = c(3, 12)) +
    guides(color = "none", fill = "none")
ggsave("output/figures/bubble_plot_genera_after_scrubbing.pdf", width = 12, height = 10, dpi = 300)




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Remove genera that are known contaminants or not expected to be present in the samples
# ────────────────────────────────────────────────────────────────────────────────────────────────────

physeq_all_filt <- physeq_scrubbed_plus_nc_rem %>%
    subset_taxa(., !Genus %in% genera_order_to_be_grayed) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)

save(physeq_all_filt, file = "data/processed_data/physeq_all_filt.RData")
