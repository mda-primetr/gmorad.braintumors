# Script for alluvial plot to show the filteration steps in 16S data
source("src/libraries.R")
source("src/common_functions.R")

# Load samples only ----

load("data/processed_data/physeq_samples.RData")


# Load after first filter ---
load("data/processed_data/physeq_filt_1.RData")

# Load after second filter ---
load("data/processed_data/physeq_filt_2.RData")

# Load after third filter ---
load("data/processed_data/physeq_filt_3.RData")



# Physeq samples is already in counts format so it can be used directly for alluvial plot

# Create function that returns a list of ASVs that are present in the physeq object > 0

get_asvs <- function(physeq) {
    df_list <- physeq %>%
        psmelt() %>%
        filter(Abundance > 0) %>%
        dplyr::select(OTU, Abundance)
    return(df_list$OTU)
}

# Before any filtering ----
physeq_samples_melt <- physeq_samples %>%
    speedyseq::tax_glom(taxrank = "Genus", .) %>%
    psmelt() %>%
    mutate(filter_step = "Before Filter") %>%
    dplyr::select(Genus, Abundance, filter_step)

# First physeq filter ----
physeq_filt_1_melt <- physeq_samples %>%
    prune_taxa(get_asvs(physeq_filt_1), .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    speedyseq::tax_glom(taxrank = "Genus", .) %>%
    psmelt() %>%
    mutate(filter_step = "Proportion filter") %>%
    dplyr::select(Genus, Abundance, filter_step)

# Second physeq filter ----
physeq_filt_2_melt <- physeq_samples %>%
    prune_taxa(get_asvs(physeq_filt_2), .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    speedyseq::tax_glom(taxrank = "Genus", .) %>%
    psmelt() %>%
    mutate(filter_step = "Prevalence filter") %>%
    dplyr::select(Genus, Abundance, filter_step)



# Third physeq filter ----
physeq_filt_3_melt <- physeq_samples %>%
    prune_taxa(get_asvs(physeq_filt_3), .) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    speedyseq::tax_glom(taxrank = "Genus", .) %>%
    psmelt() %>%
    mutate(filter_step = "unique batch filter") %>%
    dplyr::select(Genus, Abundance, filter_step)



# Scrubbed phyloseq object ----
load("data/processed_data/physeq_all_scrubbed.RData")

# Scrubbed phyloseq object ----
physeq_scrub_melt <- physeq_all_scrubbed %>%
    speedyseq::tax_glom(taxrank = "Genus", .) %>%
    psmelt() %>%
    mutate(filter_step = "SCRuB") %>%
    dplyr::select(Genus, Abundance, filter_step)

# Load scrubbed and negative controls removed phyloseq object ----
load("data/processed_data/physeq_scrubbed_plus_nc_rem.RData")

physeq_scrubbed_plus_nc_rem_melt <- physeq_scrubbed_plus_nc_rem %>%
    speedyseq::tax_glom(taxrank = "Genus", .) %>%
    psmelt() %>%
    mutate(filter_step = "Negative control removal") %>%
    dplyr::select(Genus, Abundance, filter_step)




# Combine all data ----
all_data <- rbind(
    physeq_samples_melt,
    physeq_filt_1_melt,
    physeq_filt_2_melt,
    physeq_filt_3_melt,
    physeq_scrub_melt,
    physeq_scrubbed_plus_nc_rem_melt
)

df_alluvial_filter <- all_data %>%
    dplyr::select("Genus", "Abundance", filter_step) %>%
    mutate(Genus = ifelse(Genus == "67-14_unknown", "LKT__Solirubrobacterales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Acidimicrobiia_unknown", "LKT__Acidimicrobiia", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Saccharimonadales_unknown", "LKT__Saccharimonadales", Genus)) %>%
    mutate(Genus = ifelse(Genus == "Sphingomonadaceae_unknown", "LKT__Sphingomonadaceae", Genus)) %>%
    group_by(Genus, filter_step) %>%
    summarize(sum_no_of_reads = sum(Abundance)) %>%
    group_by(filter_step) %>%
    mutate(prop = sum_no_of_reads / sum(sum_no_of_reads) * 100) %>%
    dplyr::select(site = filter_step, Genus2 = Genus, sum_no_of_reads, prop) %>%
    mutate(Genus2 = factor(Genus2)) %>%
    arrange(prop) %>%
    filter(prop > 0) %>%
    mutate(site = factor(site, levels = c(
        "Before Filter",
        "Proportion filter",
        "Prevalence filter",
        "unique batch filter",
        "SCRuB",
        "Negative control removal"
    )))


myPal <- tax_palette(physeq_all_scrubbed, rank = "Genus", pal = "brewerPlus", n = 41)


ggplot(
    df_alluvial_filter,
    aes(
        x = site, stratum = reorder(Genus2, prop), alluvium = Genus2,
        y = prop,
        fill = Genus2
    )
) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow(width = 1 / 6) +
    geom_stratum(alpha = .5, width = 1 / 4, discern = TRUE) +
    ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(site) == 1, as.character(Genus2), NA)),
        stat = "stratum", size = 3, direction = "x", alpha = 0.6, nudge_x = -0.1, max.overlaps = 3,
    ) +
    ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(site) == 2, as.character(Genus2), NA)),
        stat = "stratum", size = 3, direction = "x", alpha = 0.6, nudge_x = -0.1, max.overlaps = 3,
    ) +
    ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(site) == 3, as.character(Genus2), NA)),
        stat = "stratum", size = 3, direction = "x", alpha = 0.6, nudge_x = -0.1, max.overlaps = 3,
    ) +
    ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(site) == 4, as.character(Genus2), NA)),
        stat = "stratum", size = 3, direction = "x", alpha = 0.6, nudge_x = -0.1, max.overlaps = 3,
    ) +
    ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(site) == 5, as.character(Genus2), NA)),
        stat = "stratum", size = 3, direction = "x", alpha = 0.6, nudge_x = -0.1, max.overlaps = 3,
    ) +
    ggrepel::geom_text_repel(
        aes(label = ifelse(as.numeric(site) == 6, as.character(Genus2), NA)),
        stat = "stratum", size = 3, direction = "x", alpha = 0.6, nudge_x = -0.1, max.overlaps = 3,
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = myPal) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(text = element_text(size = 14)) +
    theme(legend.position = "none") +
    coord_cartesian(clip = "off")
ggsave("output/figures/alluvial_plot_to_show_filter_step.pdf", width = 17, height = 10, dpi = 300)
