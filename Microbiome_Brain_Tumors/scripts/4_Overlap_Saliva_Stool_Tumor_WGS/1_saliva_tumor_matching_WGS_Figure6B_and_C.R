source("src/libraries.R")
source("src/common_functions.R")



# Read metadata ----
df_metadata <- read_csv("data/metadata/megametadata_combined_16S_WGS_2025_05_01.csv") %>%
    mutate(sample_id = sequencing_id)


df_metadata_saliva_tumor <- df_metadata %>%
    filter(microbiome_sample_type %in% c("Saliva", "Tumor")) %>%
    mutate(sample_id = gsub("\\.", "-", sample_id)) %>%
    dplyr::select(sample_id, coded_id, microbiome_sample_type, tumor_category)


# Read self comparison genome comparison data ----
df_self_matched <- read_tsv("data/raw_data/WGS/SalivaSelfVsNonSelfComparison.04.11.2025.txt")


df_matched <- read_tsv(
    "data/raw_data/WGS/TumorSalivaGenomeFilteredComparison.31Taxa.processedSame.percFiltered.origReconstructions.static.clusters.04.11.2025.txt"
) %>%
    janitor::clean_names() %>%
    dplyr::rename("Saliva" = "other_sample", "Tumor" = "tumor") %>%
    inner_join(df_metadata_saliva_tumor, by = c("Saliva" = "sample_id")) %>%
    inner_join(df_metadata_saliva_tumor, by = c("Tumor" = "sample_id")) %>%
    filter(coded_id.x == coded_id.y) %>%
    filter(taxa_name %in% df_self_matched$taxa_name) %>%
    dplyr::select(
        taxa_name,
        tumor_total_genome_length,
        other_sample_genome_match,
        percent_of_length_matching,
        coded_id = coded_id.x,
        tumor_category = tumor_category.x
    )


df_matched %>%
    mutate(coded_id = substring(coded_id, 1, 10)) %>%
    ggplot(., aes(x = taxa_name, y = percent_of_length_matching)) +
    geom_bar(stat = "identity", width = 0.6) +
    facet_wrap(~coded_id, scales = "free_y") +
    coord_flip() +
    theme_bw(base_size = 17)


# Get color palette for the genera to be colored ----
color_for_species <- distinct_palette()[1:length(unique(df_matched$taxa_name))]
color_for_species[color_for_species == "#1ff8ff"] <- "#85C4Cf"


# Alternative way to plot the data ----
df_matched_order <- df_matched %>%
    mutate(coded_id = substring(coded_id, 1, 4)) %>%
    group_by(tumor_category, coded_id) %>%
    summarize(n_count = n()) %>%
    arrange(tumor_category, desc(n_count))


# Species that passed Tumor filters -----

# Get taxa that remained after all WGS tumor filters ----
list_of_wgs_remaining_taxa <- read_csv(
    "output/tables/WGS/Tumor_samples_wgs_after_all_filter_including_average_letters_filter.csv"
)





df_plt_overlap <- df_matched %>%
    mutate(coded_id = substring(coded_id, 1, 4)) %>%
    mutate(taxa_name = gsub("\\s", ".", taxa_name)) %>%
    filter(taxa_name %in% list_of_wgs_remaining_taxa$Species) %>%
    mutate(taxa_name = gsub("\\.", " ", taxa_name)) %>%
    mutate(coded_id = factor(coded_id, levels = df_matched_order$coded_id)) %T>%
    write_csv(
        "output/figures/WGS/matched_percentage_taxa_saliva_tumor_alternative.csv"
    ) %>%
    mutate(coded_id_tumor = paste(coded_id, tumor_category, sep = "_"))

df_plt_overlap$coded_id


plt_saliva_overlap <- df_plt_overlap %>%
    ggplot(
        .,
        aes(
            y = tidytext::reorder_within(
                coded_id,
                percent_of_length_matching,
                taxa_name
            ),
            x = percent_of_length_matching,
            fill = taxa_name
        )
    ) +
    geom_col(stat = "identity", width = 1) +
    facet_grid(
        coded_id + tumor_category ~ .,
        space = "free",
        scales = "free_y",
        switch = "y"
    ) +
    scale_fill_manual(values = color_for_species) +
    theme_bw(base_size = 20) +
    ylab("Patient ID") +
    tidytext::scale_y_reordered() +
    # legend as one column
    scale_x_continuous(expand = c(0, 0)) +
    theme(strip.background = element_blank()) +
    theme(
        # remove x axis ticks
        strip.text.y.left = element_text(angle = 0),
        axis.text.y = element_blank(), # remove y axis labels
        axis.ticks.y = element_blank() # remove y axis ticks
    ) +
    labs(fill = "Species shared between tumor and saliva") +
    xlab("Proportion of length matching") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )


save(plt_saliva_overlap,
    file = "output/figures/WGS/plt_saliva_overlap.RData"
)
