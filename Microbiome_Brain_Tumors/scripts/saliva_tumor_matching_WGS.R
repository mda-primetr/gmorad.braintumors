# Script to compare matches between saliva and tumor samples based on WGS data
source("src/libraries.R")
source("src/common_functions.R")

load("data/processed_data/physeq_WGS_tumor_saliva.RData")

df_metadata <- sample_data(physeq_WGS_tumor_saliva) %>%
    data.frame() %>%
    mutate(sample_id = gsub("\\.", "-", sample_id)) %>%
    dplyr::select(sample_id, coded_id, microbiome_sample_type, tumor_category)


df_metadata_saliva_tumor <- df_metadata %>%
    filter(microbiome_sample_type %in% c("Saliva", "Tumor")) %>%
    pivot_wider(names_from = microbiome_sample_type, values_from = sample_id) %>%
    janitor::clean_names()



df_matched <- read_tsv("data/raw_data/TumorSalivaGenomeComparison.txt") %>%
    janitor::clean_names() %>%
    dplyr::rename("saliva" = "other_sample") %>%
    inner_join(df_metadata_saliva_tumor, by = c("tumor", "saliva"))


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

# Filter the species based on the manual inspection ----
list_species_to_remove <- c(
    "Streptococcus gordonii",
    "Rothia mucilaginosa",
    "Prevotella nigrescens",
    "Prevotella intermedia",
    "Lautropia mirabilis",
    "Lachnospiraceae bacterium oral taxon 096",
    "Haemophilus parainfluenzae",
    "Streptococcus sp. LPB0220",
    "Segatella copri DSM 18205",
    "Segatella copri",
    "Capnocytophaga sp. FDAARGOS_737"
)


df_plt_overlap <- df_matched %>%
    mutate(coded_id = substring(coded_id, 1, 4)) %>%
    mutate(coded_id = factor(coded_id, levels = df_matched_order$coded_id)) %T>%
    write_csv("output/figures/matched_percentage_taxa_saliva_tumor_alternative.csv") %>%
    filter(!taxa_name %in% list_species_to_remove) %>%
    mutate(coded_id_tumor = paste(coded_id, tumor_category, sep = "_"))

df_plt_overlap$coded_id

df_plt_overlap %>%
    ggplot(., aes(y = tidytext::reorder_within(coded_id, percent_of_length_matching, taxa_name), x = percent_of_length_matching, fill = taxa_name)) +
    geom_col(stat = "identity", width = 1) +
    facet_grid(coded_id + tumor_category ~ .,
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
    theme( # remove x axis ticks
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
        panel.border = element_blank(), panel.background = element_blank()
    )
ggsave("output/figures/matched_percentage_taxa_saliva_tumor_alternative.pdf", width = 16, height = 14, dpi = 300)
