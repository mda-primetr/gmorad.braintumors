source("src/libraries.R")
source("src/common_functions.R")

load("data/processed_data/WGS/physeq_bacmap_tumor_post_filter.Rdata")




Tumor_only_WGS <- physeq_bacmap_tumor_post_filter %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)

# Get count and presence of the species by tumor category
df_presence <- Tumor_only_WGS %>%
    psmelt() %>%
    dplyr::select(Species, Abundance, tumor_category) %>%
    group_by(Species, tumor_category) %>%
    summarize(count_n = sum(Abundance > 0)) %>%
    mutate(present = ifelse(count_n > 0, "Yes", "No")) %>%
    pivot_wider(
        names_from = tumor_category,
        values_from = c(count_n, present),
        names_glue = "{tumor_category}_{.value}"
    )

# Get the dataframe for saliva and tumor match ----
df_species_genome_completion <- read_tsv(
    "data/raw_data/WGS/GolnazGBM.tumor.only.genome.completion.04.04.2025.txt"
) %>%
    mutate(Name = gsub("\\s", ".", Name)) %>%
    mutate(Name = gsub("_", " ", Name))




# Get the list of species that are present in the tumor samples
Tumor_only_WGS %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    group_by(Species) %>%
    summarize(
        total_abundance = sum(Abundance),
        median_abundance = median(Abundance),
        min_abundance = min(Abundance),
        max_abundance = max(Abundance)
    ) %>%
    inner_join(df_presence, by = "Species") %>%
    left_join(df_species_genome_completion, by = c("Species" = "Name")) %>%
    write_csv("output/tables/WGS/Tumor_samples_wgs_abundance_genome_info_table_S4.csv")
