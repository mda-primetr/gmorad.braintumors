source("src/libraries.R")
source("src/common_functions.R")



# Combine stool and saliva overlap plots ----
load(file = "output/figures/WGS/plt_saliva_overlap.RData")
load(file = "output/figures/WGS/matched_percentage_taxa_stool_tumor_alternative.RData")



# Save the data for the overlap plot ----
df_combo_saliva_stool_matched_tumor <- rbind(
    plt_saliva_overlap$data %>%
        mutate(sample_type = "Saliva"),
    plt_stool_overlap$data %>%
        mutate(sample_type = "Stool")
) %>%
    write_csv("data/processed_data/WGS/Tumor_Saliva_Stool_overlap_plt_data.csv")


# Get color palette for the genera to be colored ----
color_for_species <- distinct_palette()[1:length(unique(df_combo_saliva_stool_matched_tumor$taxa_name))]
color_for_species[color_for_species == "#1ff8ff"] <- "#85C4Cf"


# Saliva and stool overlap with Tumor plot ----
plt_combo_1 <- df_combo_saliva_stool_matched_tumor %>%
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
        coded_id + tumor_category + sample_type ~ .,
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
    labs(fill = "Species shared between tumor and saliva or stool") +
    xlab("Proportion of length matching") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )
ggsave(
    filename = "output/figures/WGS/matched_percentage_taxa_saliva_stool_tumor_alternative.pdf",
    width = 16,
    height = 14
)



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Plot these species that had overlapped plot above
# Show other species that are not in the overlap plot
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_combo_compare <- df_combo_saliva_stool_matched_tumor %>%
    mutate(taxa_name = gsub("\\s", "\\.", taxa_name)) %>%
    rename(sample_type_site = sample_type)



# Load WGS tumor phyloseq object after all filter ----
load(file = "data/processed_data/WGS/physeq_bacmap_tumor_post_filter.Rdata")


# Get taxa that matched tumor and saliva ----

df_combo_with_taxa_category <- rbind(
    # Not in overlap with saliva or stool
    physeq_bacmap_tumor_post_filter %>%
        tax_glom(taxrank = "Species") %>%
        psmelt() %>%
        filter(Abundance > 0) %>%
        mutate(coded_id = substring(coded_id, 1, 4)) %>%
        anti_join(df_combo_compare, c(
            "coded_id" = "coded_id",
            "Species" = "taxa_name"
        )) %>%
        dplyr::select(
            coded_id,
            Species,
            Abundance,
            tumor_category,
        ) %>%
        group_by(coded_id, Species, tumor_category) %>%
        summarise(Abundance = sum(Abundance)) %>%
        mutate(taxa_status = "Not in overlap with saliva or stool"),
    # Overlap with saliva or stool
    physeq_bacmap_tumor_post_filter %>%
        tax_glom(taxrank = "Species") %>%
        psmelt() %>%
        filter(Abundance > 0) %>%
        mutate(coded_id = substring(coded_id, 1, 4)) %>%
        inner_join(df_combo_compare, c(
            "coded_id" = "coded_id",
            "Species" = "taxa_name"
        )) %>%
        dplyr::select(
            coded_id,
            Species,
            Abundance,
            sample_type_site,
            tumor_category = tumor_category.x
        ) %>%
        group_by(coded_id, Species, tumor_category, sample_type_site) %>%
        summarise(Abundance = sum(Abundance)) %>%
        mutate(taxa_status = case_when(
            sample_type_site == "Stool" ~ "Overlap with stool",
            sample_type_site == "Saliva" ~ "Overlap with saliva",
            TRUE ~ "Not in overlap with saliva or stool"
        ))
) %>%
    group_by(coded_id, tumor_category, taxa_status) %>%
    summarise(richness = n_distinct(Species))



# Extract the order of the coded_id from the first plot
built <- ggplot_build(plt_combo_1)
panels <- unique(built$data[[1]]$PANEL)
coded_id_order <- built$layout$layout$coded_id[order(built$layout$layout$PANEL)]

# Plot the data ----
plt_combo_2 <- df_combo_with_taxa_category %>%
    filter(coded_id %in% df_combo_saliva_stool_matched_tumor$coded_id) %>%
    left_join(df_combo_compare, by = c("coded_id")) %>%
    distinct(coded_id, taxa_status, .keep_all = TRUE) %>%
    mutate(taxa_status = factor(taxa_status, levels = c("Overlap with saliva", "Overlap with stool", "Not in overlap with saliva or stool"))) %>%
    mutate(coded_id = factor(coded_id, levels = coded_id_order)) %>%
    ggplot(., aes(x = coded_id, y = richness, fill = taxa_status)) +
    geom_bar(stat = "identity", width = 0.6) +
    theme_bw(base_size = 20) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
    ) +
    theme(strip.background = element_blank()) +
    scale_fill_manual(values = c(
        "Overlap with stool" = "#2FB2AB",
        "Overlap with saliva" = "#5B68B0",
        "Not in overlap with saliva or stool" = "gray"
    ))
ggsave(
    filename = "output/figures/WGS/matched_percentage_taxa_saliva_stool_tumor_alternative_2.pdf",
    width = 18,
    height = 10
)



plt_combo_2$data %>%
    dplyr::select(coded_id, tumor_category = tumor_category.x, richness, taxa_status) %>%
    write_csv("data/processed_data/WGS/matched_percentage_taxa_saliva_stool_tumor_alternative_2.csv")
