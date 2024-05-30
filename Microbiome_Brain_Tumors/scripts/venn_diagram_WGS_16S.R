source("src/libraries.R")
source("src/common_functions.R")

# ────────────────────────────────────────────────────────────────────────────────────────────────────
#-------------------------- Tracking taxa through the filteration steps  -----------------------#
# ────────────────────────────────────────────────────────────────────────────────────────────────────


# Calculate the number of species at each step of filteration ----
# Plot top taxa for all sample types -----

load("data/processed_data/physeq_WGS_tumor_saliva.RData")


physeq_bacmap_clean_taxa <- fn_clean_tax_table(physeq_WGS_tumor_saliva)


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. Pre-filter : Tumor only samples species count
# ────────────────────────────────────────────────────────────────────────────────────────────────────

Tumor_only_WGS <- physeq_bacmap_clean_taxa %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)

# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 2. Remove species < 303 reads length
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_species_303 <- read_tsv("data/metadata/Golnaz.CP05552.GenomeCompletion.05.18.2024.txt") %>%
    janitor::clean_names() %>%
    filter(total_letters_in_reads >= 303) %>%
    mutate(name = gsub(" ", ".", name)) %>%
    pull(name)

Tumor_only_WGS_filt1 <- Tumor_only_WGS %>%
    subset_taxa(Species %in% df_species_303) %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)


# ────────────────────────────────────────────────────────────────────────────────────────────────────
#  3. Remove species with coverage > 2
# ────────────────────────────────────────────────────────────────────────────────────────────────────
df_species_303_2x <- read_tsv("data/metadata/Golnaz.CP05552.GenomeCompletion.05.18.2024.txt") %>%
    janitor::clean_names() %>%
    filter(total_letters_in_reads >= 303) %>%
    filter(max_coverage_in_completed_segments <= 2) %>%
    mutate(name = gsub(" ", ".", name)) %>%
    pull(name)


Tumor_only_WGS_filt2 <- Tumor_only_WGS_filt1 %>%
    subset_taxa(Species %in% df_species_303_2x) %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# 4. Remove taxa from controls
# ────────────────────────────────────────────────────────────────────────────────────────────────────

physeq_bacmap_taxa_controls <- physeq_bacmap_clean_taxa %>%
    subset_samples(microbiome_sample_type == "Control") %>%
    speedyseq::filter_sample_data(sample_id %in% c(
        "DNA_EXT-NTC_L005_R1_001.fastq.gz",
        "DNA_EXT-PTC_L005_R1_001.fastq.gz"
    )) %>%
    prune_taxa(taxa_sums(.) >= 1, .)



# get list of species in the controls
list_of_control_species <- psmelt(physeq_bacmap_taxa_controls) %>%
    group_by(Species) %>%
    summarize(Abundance = sum(Abundance)) %>%
    pull(Species)

# Get list of species that passed filter >= 303 reads and coverage <= 2
list_of_sample_species <- psmelt(Tumor_only_WGS_filt2) %>%
    group_by(Species) %>%
    summarize(Abundance = sum(Abundance)) %>%
    pull(Species)

# Remove the species that are found in the controls from the list of species that passed the genome completion
list_species_minus_controls <- setdiff(list_of_sample_species, list_of_control_species)


Tumor_only_WGS_filt2_and_controls <- Tumor_only_WGS_filt2 %>%
    subset_taxa(., Species %in% list_species_minus_controls) %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    prune_taxa(taxa_sums(.) >= 1, .)



# count species by tumor category ----
psmelt(Tumor_only_WGS_filt2_and_controls) %>%
    filter(Abundance > 0) %>%
    group_by(tumor_category) %>%
    summarize(count = n_distinct(Species))


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Venn diagram for 16S and WGS data Tumor data ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Get the samples and get taxa with abundance >= 1
df_tumor_WGS_data <- Tumor_only_WGS_filt2_and_controls %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    dplyr::select(Genus, Abundance) %>%
    group_by(Genus) %>%
    summarise(Abundance = sum(Abundance)) %>%
    pull(Genus)


# Get 16S tumor data that is filtered for 5 step filter and manual removal ----
load("data/processed_data/physeq_all_filt.RData")

psmelt(physeq_all_filt) %>%
    filter(Abundance > 0) %>%
    write_csv("data/processed_data/physeq_all_filt_melted.csv")



df_16s_tumor_data <- read_csv("data/processed_data/physeq_all_filt_melted.csv") %>%
    filter(microbiome_sample_type == "Tumor") %>%
    dplyr::select(Genus, Abundance) %>%
    group_by(Genus) %>%
    summarise(Abundance = sum(Abundance)) %>%
    # Assuming that Prevotella_7 is the same as Prevotella
    mutate(Genus = ifelse(Genus == "Prevotella_7", "Prevotella", Genus)) %>%
    pull(Genus)




# Create venn diagram for the tumor data for WGS and 16S ----



# List of items
wgs_16s_combo <- list("16s" = df_16s_tumor_data, "WGS" = df_tumor_WGS_data)



ggVennDiagram::ggVennDiagram(wgs_16s_combo, label_size = 7, force_upset = F) +
    theme_void(base_size = 20) +
    scale_fill_distiller(palette = "RdBu") +
    coord_flip() +
    theme(legend.position = "none") +
    theme(
        axis.title.x = element_text(size = 24)
    ) +
    ylab("")
ggsave("output/figures/venn_diagram_tumor_WGS_16S.pdf", width = 6, height = 4.5)





# Genera common in both 16S and WGS data
common_genera_WGS <- intersect(df_tumor_WGS_data, df_16s_tumor_data) %>%
    cbind() %>%
    data.frame() %>%
    dplyr::rename(Genus = ".") %>%
    write_csv("data/processed_data/common_genera_WGS_16S_with_completion_coverage.csv")

# Genera in WGS but not in 16S data
exclusive_genera_WGS <- setdiff(df_tumor_WGS_data, df_16s_tumor_data) %>%
    cbind() %>%
    data.frame() %>%
    dplyr::rename(Genus = ".") %>%
    write_csv("data/processed_data/Genera_exclusive_WGS_with_completion_coverage.csv")


# Genera in 16S but not in WGS data
setdiff(df_16s_tumor_data, df_tumor_WGS_data) %>%
    cbind() %>%
    data.frame() %>%
    dplyr::rename(Genus = ".") %>%
    write_csv("data/processed_data/Genera_exclusive_16S_with_completion_coverage.csv")


# Get species in the WGS data but not in the 16S data
Tumor_only_WGS_filt2_and_controls %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    filter(Genus %in% exclusive_genera_WGS$Genus) %>%
    distinct(Species) %>%
    arrange(Species) %>%
    write_csv("data/processed_data/Species_exclusive_WGS_with_completion_coverage.csv")


# List of species that have their genera common in both 16S and WGS data
Tumor_only_WGS_filt2_and_controls %>%
    ps_filter(microbiome_sample_type == "Tumor") %>%
    prune_samples(sample_sums(.) >= 1, .) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    filter(Genus %in% common_genera_WGS$Genus) %>%
    distinct(Species) %>%
    write_csv("data/processed_data/Species_Common_WGS_16S_genus_with_completion_coverage.csv")
