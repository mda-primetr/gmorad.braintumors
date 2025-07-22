source("src/libraries.R")
source("src/common_functions.R")



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Read BacMap for Tumor data before any filteration -----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Pre-filtered phyloseq object ----
physeq_bacmap_tumor_only <- readRDS(
    "data/raw_data/WGS/physeq_bacmap_tumor_prefilter.rds"
)


# Clean taxa tables for spaces or other characters ----
physeq_bacmap_clean_taxa_tumor_prefilter <- fn_clean_tax_table(physeq_bacmap_tumor_only)



# Save the list of taxa after genome completion and coverage filter
physeq_bacmap_clean_taxa_tumor_prefilter %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    distinct(Species, Abundance) %>%
    write_csv(
        "output/tables/WGS/Tumor_Species_abundance_wgs_before_any_filter.csv"
    )



# ────────────────────────────────────────────────────────────────────────────────────────────────────
#-------------------------- Tracking taxa through the filteration steps  -----------------------#
# ────────────────────────────────────────────────────────────────────────────────────────────────────



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. Remove controls -----
# ────────────────────────────────────────────────────────────────────────────────────────────────────



# Calculate the number of species at each step of filteration ----
# Plot top taxa for all sample types -----

# This phyloseq object is already filtered for controls and other contaminants ---
physeq_bacmap_tumor_controls_removed <- readRDS(
    "data/raw_data/WGS/physeq_bacmap_tumor_controls_removed.rds"
)



physeq_bacmap_clean_taxa <- fn_clean_tax_table(physeq_bacmap_tumor_controls_removed)


physeq_bacmap_clean_taxa



Tumor_only_WGS <- physeq_bacmap_clean_taxa %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)



# Get the list of species that are present in the tumor samples
Tumor_only_WGS %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    distinct(Species, Abundance) %>%
    write_csv("output/tables/WGS/Tumor_species_abundance_after_controls_removal_only.csv")



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 2. Remove species not present in the genome completion table
# In other words, only keep species that are found in the genome completion table
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_species_genome_completion <- read_tsv(
    "data/raw_data/WGS/GolnazGBM.tumor.only.genome.completion.04.04.2025.txt"
) %>%
    janitor::clean_names() %>%
    mutate(average_letters_per_reads = total_letters_in_reads / reads)



df_species_genome_completion_filtered <- df_species_genome_completion %>%
    filter(total_letters_in_reads >= 1) %>%
    mutate(name = gsub(" ", ".", name)) %>%
    mutate(name = gsub("\\_", " ", name)) %>%
    dplyr::select(
        name,
        taxa_id,
        total_letters_in_reads,
        max_coverage_in_completed_segments
    )

Tumor_only_WGS_filt1 <- Tumor_only_WGS %>%
    subset_taxa(Species %in% df_species_genome_completion_filtered$name) %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)

# Check species that are not in the genome completion table
setdiff(
    tax_table(Tumor_only_WGS_filt1) %>% data.frame() %>% pull(Species),
    df_species_genome_completion_filtered$name
)

# Save the list of taxa after genome completion and coverage filter
Tumor_only_WGS_filt1 %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    distinct(Species, Abundance) %>%
    write_csv(
        "output/tables/WGS/Tumor_samples_wgs_after_genome_completion_filter.csv"
    )


print("Start removing taxa with coverage > 2")
# ────────────────────────────────────────────────────────────────────────────────────────────────────
#  3. Remove species with coverage > 2
# ────────────────────────────────────────────────────────────────────────────────────────────────────
df_species_303_2x <- df_species_genome_completion_filtered %>%
    filter(max_coverage_in_completed_segments <= 2)

Tumor_only_WGS_filt2 <- Tumor_only_WGS_filt1 %>%
    subset_taxa(Species %in% df_species_303_2x$name) %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)


# Save the list of taxa after genome completion and coverage filter
Tumor_only_WGS_filt2 %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    distinct(Species, Abundance) %>%
    write_csv(
        "output/tables/WGS/Tumor_samples_wgs_after_genome_completion_And_coverage_filter.csv"
    )


# Review the data after filtering
Tumor_only_WGS_filt2 %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    dplyr::select(Species, Abundance, dataset, tumor_category) %>%
    group_by(Species, tumor_category) %>%
    summarize(Abundance = sum(Abundance)) %>%
    pivot_wider(
        names_from = tumor_category,
        values_from = Abundance,
        values_fill = 0
    )


# Get list of species that passed coverage <= 2
list_of_sample_species <- psmelt(Tumor_only_WGS_filt2) %>%
    group_by(Species) %>%
    summarize(Abundance = sum(Abundance))


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 4. Remove median of average letters per read
# ────────────────────────────────────────────────────────────────────────────────────────────────────


# Calculate quartiles for average letters per read
quartiles <- quantile(
    df_species_genome_completion$average_letters_per_reads,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    na.rm = TRUE
)
print(quartiles)

# Comprehensive summary
summary(df_species_genome_completion$average_letters_per_reads)

# Visualize the distribution
# hist(
#     df_species_genome_completion$average_letters_per_reads,
#     main = "Distribution of Average Letters per Read",
#     xlab = "Average Letters per Read",
#     breaks = 30
# )

# # Add vertical lines for quartiles
# abline(
#     v = quartiles,
#     col = c("red", "blue", "green", "blue", "red"),
#     lwd = 2,
#     lty = c(2, 2, 1, 2, 2)
# )
# legend(
#     "topright",
#     legend = c("Min/Max", "Q1/Q3", "Median"),
#     col = c("red", "blue", "green"),
#     lwd = 2,
#     lty = c(2, 2, 1)
# )


# Use median as cutoff
median_cutoff <- quartiles[3]


# Using the median as a cutoff
df_species_less_than_130_letters_average <- df_species_genome_completion %>%
    filter(average_letters_per_reads >= round(median_cutoff)) %>%
    mutate(name = gsub(" ", ".", name)) %>%
    mutate(name = gsub("\\_", " ", name))




Tumor_only_WGS_filt3 <- Tumor_only_WGS_filt2 %>%
    subset_taxa(Species %in% df_species_less_than_130_letters_average$name) %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)


# Save phyloseq object after applying all the filters
physeq_bacmap_tumor_post_filter <- Tumor_only_WGS_filt3
save(physeq_bacmap_tumor_post_filter, file = "data/processed_data/WGS/physeq_bacmap_tumor_post_filter.Rdata")


# Save the list of taxa after genome completion and coverage filter
physeq_bacmap_tumor_post_filter %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    distinct(Species, Abundance) %>%
    write_csv(
        "output/tables/WGS/Tumor_samples_wgs_after_all_filter_including_average_letters_filter.csv"
    )




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Check taxa together on how they are removed by looking at combined tables after filtering steps ---
# ────────────────────────────────────────────────────────────────────────────────────────────────────

rbind(
    read_csv("output/tables/WGS/Tumor_Species_abundance_wgs_before_any_filter.csv") %>%
        group_by(Species) %>%
        summarize(Abundance = sum(Abundance)) %>%
        mutate(filter_step = "Before any filtering"),
    read_csv(
        "output/tables/WGS/Tumor_species_abundance_after_controls_removal_only.csv"
    ) %>%
        group_by(Species) %>%
        summarize(Abundance = sum(Abundance)) %>%
        mutate(filter_step = "After controls removal filter"),
    read_csv(
        "output/tables/WGS/Tumor_samples_wgs_after_genome_completion_filter.csv"
    ) %>%
        group_by(Species) %>%
        summarize(Abundance = sum(Abundance)) %>%
        mutate(filter_step = "After control and genome completion filter"),
    read_csv(
        "output/tables/WGS/Tumor_samples_wgs_after_genome_completion_And_coverage_filter.csv"
    ) %>%
        group_by(Species) %>%
        summarize(Abundance = sum(Abundance)) %>%
        mutate(filter_step = "After control, genome completion, and coverage filter"),
    read_csv(
        "output/tables/WGS/Tumor_samples_wgs_after_all_filter_including_average_letters_filter.csv"
    ) %>%
        group_by(Species) %>%
        summarize(Abundance = sum(Abundance)) %>%
        mutate(filter_step = "After control, genome completion,coverage filter, and median of average letters")
) %>%
    write_csv(paste0(
        "output/tables/WGS/combined_view_Tumor_samples_wgs_after_all_filtering_",
        ".csv"
    ))




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Bubble plot  for WGS
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Genera that needs to be grayed out ----
genera_order_to_be_grayed <- c(
    "Pseudomonas",
    "Sphingobacterium",
    "Jeotgalibaca",
    "Mycolicibacterium",
    "Streptomyces",
    "Sphingomonas"
)



# Get color palette for the genera to be colored ----

# Get color scheme from the 16s plot ----
load("data/processed_data/16S/genera_colors.RData")


genera_colors_16s <- genera_colors %>%
    enframe() %>%
    mutate(value = case_when(
        name %in% genera_order_to_be_grayed ~ "gray",
        TRUE ~ value
    )) %>%
    mutate(seq_tech = "16S") %>%
    mutate(name = ifelse(name == "Prevotella_7", "Prevotella", name))


# Colors for WGS ---
genera_colors_WGS <- rev(distinct_palette(pal = "greenArmytage")[1:length(unique(data.frame(tax_table(Tumor_only_WGS_filt3))$Genus))])


# Create a color assignment dataframe for WGS genera
wgs_colors <- unique(data.frame(tax_table(Tumor_only_WGS_filt3))$Genus) %>%
    enframe() %>%
    rename(Genus = value) %>%
    cbind(genera_colors_WGS) %>%
    left_join(genera_colors_16s, by = c("Genus" = "name")) %>%
    mutate(genera_colors_WGS = case_when(
        !is.na(value) ~ value,
        Genus %in% genera_order_to_be_grayed ~ "gray",
        # Change some colors manually to enhance the contrast
        Genus == "Bifidobacterium" ~ "#c2dbdd",
        Genus == "Lactobacillus" ~ "#761646",
        TRUE ~ genera_colors_WGS
    ))

# Genera of interest and common with 16S data----
genera_of_interest <- c(
    "Fusobacterium",
    "Prevotella",
    "Veillonella",
    "Capnocytophaga",
    "Enterococcus"
)

# Define the manual order for genera
# Put the gray genera at the bottom and the colored ones on top
manual_order <- c(
    # First the genera NOT to be grayed out (those not in genera_order_to_be_grayed)
    genera_of_interest,
    # Add colored genera (not in genera_of_interest and to be grayed out)
    setdiff(
        wgs_colors$Genus[!wgs_colors$Genus %in% c(genera_of_interest)],
        c(genera_of_interest, genera_order_to_be_grayed)
    ),
    # Add gray genera at the bottom
    genera_order_to_be_grayed
)




# Create color palatte for WGS genera ----
wgs_colors_genera <- wgs_colors %>%
    dplyr::select(Genus, genera_colors_WGS) %>%
    deframe()



# Save the list of taxa after genome completion and coverage filter ----
Tumor_only_WGS_filt3 %>%
    speedyseq::tax_glom(taxrank = "Species") %>%
    psmelt() %>%
    filter(Abundance > 0) %T>%
    write_csv(
        "data/processed_data/WGS/bubble_plot_genera_WGS_data_reads_Species.csv"
    ) %>%
    filter(!is.na(tumor_category)) %>%
    write_csv(
        "data/processed_data/WGS/bubble_plot_genera_WGS_after_filtering_data_with_reads_Species.csv"
    )


# Plot the bubble plot for WGS data ----
Tumor_only_WGS_filt3 %>%
    speedyseq::tax_glom(taxrank = "Genus") %>%
    psmelt() %>%
    filter(Abundance > 0) %T>%
    write_csv(
        "output/tables/WGS/bubble_plot_genera_WGS_data_reads.csv"
    ) %>%
    filter(!is.na(tumor_category)) %>%
    mutate(
        tumor_category = case_when(
            tumor_category == "Glioma" ~ "Glioma",
            tumor_category == "Met" ~ "BrM",
            tumor_category == "No tumor" ~ "Normal Brain"
        )
    ) %>%
    mutate(
        tumor_category = factor(
            tumor_category,
            levels = c("Glioma", "BrM", "Normal Brain")
        )
    ) %>%
    filter(Abundance > 0) %T>%
    write_csv(
        "output/tables/WGS/bubble_plot_genera_WGS_after_filtering_data_with_reads.csv"
    ) %>%
    dplyr::select(Genus, Abundance, dataset, tumor_category) %>%
    mutate(
        Genus = ifelse(
            Genus == "67-14_unknown",
            "LKT__Solirubrobacterales",
            Genus
        )
    ) %>%
    mutate(
        Genus = ifelse(
            Genus == "Acidimicrobiia_unknown",
            "LKT__Acidimicrobiia",
            Genus
        )
    ) %>%
    mutate(
        Genus = ifelse(
            Genus == "Saccharimonadales_unknown",
            "LKT__Saccharimonadales",
            Genus
        )
    ) %>%
    mutate(
        Genus = ifelse(
            Genus == "Sphingomonadaceae_unknown",
            "LKT__Sphingomonadaceae",
            Genus
        )
    ) %>%
    group_by(Genus, dataset, tumor_category) %>%
    summarize(Abundance = sum(Abundance)) %>%
    ungroup() %>%
    group_by(dataset, tumor_category) %>%
    mutate(prop = prop.table(Abundance) * 100) %>%
    dplyr::select(Genus, prop, dataset, tumor_category) %>%
    ungroup() %>%
    mutate(
        Genus = factor(
            Genus,
            levels = c(
                rev(manual_order)
            )
        )
    ) %>%
    ggplot(., aes(x = dataset, y = Genus, color = Genus, fill = prop)) +
    geom_point(aes(size = prop), alpha = 0.9) +
    facet_grid(. ~ tumor_category, scales = "free_x", space = "free") +
    theme_bw(base_size = 20) +
    scale_size_continuous(range = c(3, 12)) +
    scale_color_manual(values = wgs_colors_genera) +
    theme(
        strip.background = element_blank(),
        strip.text.y = element_text(size = 24, face = "italic"),
        axis.text.y = element_text(face = "italic"),
        panel.border = element_rect(linewidth = 0.5)
    ) +
    ggtitle("") +
    xlab("Cohort") +
    ylab("Genera") +
    guides(color = "none", fill = "none")
ggsave(
    "output/figures/WGS/bubble_plot_genera_WGS_after_filtering.pdf",
    width = 8,
    height = 10,
    dpi = 300
)
