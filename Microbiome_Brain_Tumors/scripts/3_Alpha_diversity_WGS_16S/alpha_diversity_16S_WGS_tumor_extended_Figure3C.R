source("src/libraries.R")
source("src/common_functions.R")


# Load 16s data after all filters ----
load("data/processed_data/16S/physeq_scrubbed_plus_nc_rem.RData")


# Load WGS tumor data ---
# WGS post filter ----
load(file = "data/processed_data/WGS/physeq_bacmap_tumor_post_filter.Rdata")
physeq_bacmap_tumor_post_filter



# Get taxa that remained after all WGS tumor filters ----
list_of_wgs_remaining_taxa <- read_csv(
    "output/tables/WGS/Tumor_samples_wgs_after_all_filter_including_average_letters_filter.csv"
)




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Plot alpha diversity
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_alpha_WGS <- physeq_bacmap_tumor_post_filter %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    tax_glom(taxrank = "Species") %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    group_by(coded_id, Sample, tumor_category) %>%
    summarize(count_species = n()) %>%
    mutate(sequencing_technology = "WGS")


df_alpha_16S <- physeq_scrubbed_plus_nc_rem %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    group_by(coded_id, Sample, tumor_category) %>%
    summarize(count_species = n()) %>%
    mutate(sequencing_technology = "16S")


# Find samples that are in both datasets ----

df_alpha_WGS_16S_common <- df_alpha_WGS %>%
    filter(coded_id %in% df_alpha_16S$coded_id) %>%
    mutate(sequencing_technology = "WGS")



df_combo_wgs_16s <- rbind(df_alpha_WGS, df_alpha_16S) %>%
    mutate(sequencing_technology = factor(sequencing_technology, levels = c("16S", "WGS")))

df_combo_wgs_16s %>%
    ungroup() %>%
    dplyr::select(tumor_category, count_species, sequencing_technology) %>%
    mutate(count_species = as.numeric(count_species)) %>%
    group_by(sequencing_technology) %>%
    gtsummary::tbl_summary(
        by = sequencing_technology,
        type = c("count_species") ~ "continuous2",
        statistic = all_continuous() ~ "{median} ({min}, {max})" # Explicitly set statistics for numeric variable
    )




# Common samples between WGS and 16S ----
plt_common <- df_combo_wgs_16s %>%
    filter(coded_id %in% df_alpha_WGS_16S_common$coded_id) %>%
    arrange(count_species, coded_id) %>%
    # Create a new ordering for coded_id based on the sum of count_species per coded_id
    mutate(
        coded_id_order = factor(coded_id,
            levels = df_combo_wgs_16s %>%
                filter(coded_id %in% df_alpha_WGS_16S_common$coded_id) %>%
                group_by(coded_id) %>%
                summarize(total_count = sum(count_species)) %>%
                arrange(desc(total_count)) %>%
                pull(coded_id)
        )
    ) %>%
    ggplot(aes(x = coded_id_order, y = count_species, fill = sequencing_technology)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw(base_size = 20) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_blank(), # Removes y-axis text (numbers)
        axis.ticks.y = element_blank(), # Removes y-axis tick marks/dashes
        legend.position = "top",
        legend.title = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 16), expand = c(0, 0)) +
    scale_fill_manual(values = c("16S" = "#5B68B0", "WGS" = "#2FB2AB")) +
    xlab("Patient ID") +
    ylab("")



# Common samples between WGS and 16S ----
plt_separate <- df_combo_wgs_16s %>%
    filter(!coded_id %in% df_alpha_WGS_16S_common$coded_id) %>%
    arrange(count_species, coded_id) %>%
    mutate(coded_id = factor(coded_id, unique(coded_id))) %>%
    ggplot(aes(y = count_species, x = tidytext::reorder_within(coded_id, -count_species, sequencing_technology), fill = sequencing_technology)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw(base_size = 20) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8), # Reduced size for x-axis labels
        legend.position = "top",
        legend.title = element_blank()
    ) +
    facet_grid(~sequencing_technology, scales = "free_x", space = "free") +
    scale_x_reordered() +
    scale_fill_manual(values = c("16S" = "#5B68B0", "WGS" = "#2FB2AB")) +
    xlab("Patient ID") +
    ylab("Richness") +
    theme(strip.text = element_blank()) +
    scale_y_continuous(limits = c(0, 16), expand = c(0, 0))

plt_separate$data %>%
    write_csv("output/figures/alpha_diversity_richness_tumor_samples_separate.csv")
plt_common$data %>%
    write_csv("output/figures/alpha_diversity_richness_tumor_samples_common.csv")


plt_separate + plt_common +
    plot_annotation(
        title = "Alpha diversity (richness) of tumor samples",
        theme = theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
        )
    ) +
    plot_layout(guides = "collect", heights = c(2), widths = c(2, 0.5)) &
    theme(legend.position = "bottom")
ggsave(
    "output/figures/alpha_diversity_richness_tumor_samples.pdf",
    width = 15,
    height = 9,
    dpi = 300
)
