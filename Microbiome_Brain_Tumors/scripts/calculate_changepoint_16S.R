source("src/libraries.R")
source("src/common_functions.R")

load("data/processed_data/physeq_16S_tumor.RData")


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Negative controls only ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

df_controls <- physeq_16S_tumor %>%
    speedyseq::filter_sample_data(sample_type %in% c(
        "Negative control",
        "Water control",
        "Library control",
        "sample control"
    )) %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    filter(x16s_seq_batch %in% c("Batch 3", "Batch 4", "Batch 5", "Batch 9"))


# Save negative controls as a separate file
physeq_negative_controls <- physeq_16S_tumor %>%
    speedyseq::filter_sample_data(sample_type %in% c(
        "Negative control",
        "Water control",
        "Library control",
        "sample control"
    )) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)

save(physeq_negative_controls, file = "data/processed_data/physeq_negative_controls.RData")




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Samples only ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────


df_samples <- physeq_16S_tumor %>%
    speedyseq::filter_sample_data(!sample_type %in% c(
        "Negative control",
        "Water control",
        "Library control",
        "sample control"
    )) %>%
    psmelt() %>%
    filter(Abundance > 0)


df_samples_all <- physeq_16S_tumor %>%
    speedyseq::filter_sample_data(!sample_type %in% c(
        "Negative control",
        "Water control",
        "Library control",
        "sample control"
    )) %>%
    prune_samples(sample_sums(.) > 0, .)


# Save samples as a separate file
physeq_samples <- physeq_16S_tumor %>%
    speedyseq::filter_sample_data(!sample_type %in% c(
        "Negative control",
        "Water control",
        "Library control",
        "sample control"
    )) %>%
    prune_samples(sample_sums(.) > 0, .) %>%
    prune_taxa(taxa_sums(.) > 0, .)


save(physeq_samples, file = "data/processed_data/physeq_samples.RData")

# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Proportion cutoff using controls at ASV level ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Get unique number of samples in the dataframe
df_control_n <- df_controls %>%
    dplyr::select(Sample) %>%
    unique() %>%
    nrow()

# Get prevalence of OTUs in the controls
df_prevalence <- table(df_controls$OTU) %>%
    data.frame() %>%
    mutate(Freq = (Freq / df_control_n[[1]]) * 100) %>%
    rename(OTU = Var1, Prevalence = Freq)

# Get the proportion of OTUs in the controls
df_proportion_prev <- df_controls %>%
    dplyr::select(OTU, Abundance) %>%
    group_by(OTU) %>%
    summarize(Abundance = sum(Abundance)) %>%
    mutate(Percentage = prop.table(Abundance) * 100) %>%
    inner_join(df_prevalence, by = "OTU") %>%
    rename(OTU = OTU, Abundance = Abundance, Percentage = Percentage)




pt_prop_asv <- ggcptplot(sort(df_proportion_prev$Percentage, decreasing = T),
    Q = 1,
    cptline_size = 2,
    change_in = "mean_var", cp_method = "BinSeg"
) +
    labs(title = "Changes in mean variance (BinSeg)") +
    theme_light(base_size = 20) +
    ylab("Abundance (%)") +
    labs(caption = "Proportion cutoff at ASV level") +
    xlab("OTU")

ggsave("output/figures/prop_cutoff_asv.pdf", width = 10, height = 10, dpi = 300)

prop_cutoff <- pt_prop_asv$data %>%
    filter(!is.na(cp_value))

pt_prop_asv_zoomed <- pt_prop_asv +
    ylim(c(0, 0.25)) +
    labs(title = "Changes in mean variance (BinSeg)") +
    theme_light(base_size = 20) +
    ylab("Abundance (%)") +
    labs(caption = "Proportion cutoff at ASV level") +
    xlab("OTU") +
    geom_segment(
        aes(x = -Inf, xend = prop_cutoff$row_number, y = prop_cutoff$cp_value, yend = prop_cutoff$cp_value),
        linetype = "dashed", color = "blue"
    ) +
    geom_text(aes(x = prop_cutoff$row_number, y = prop_cutoff$cp_value),
        label = glue::glue(
            {
                round(prop_cutoff$cp_value, 4)
            },
            "%"
        ), color = "black", nudge_x = 200,
        size = 6
    )
ggsave("output/figures/prop_cutoff_asv_zoomed.pdf", width = 10, height = 10, dpi = 300)




prop_cutoff <- pt_prop_asv$data %>%
    filter(!is.na(cp_value))


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Prevalance cutoffs at Genus level ----
# ────────────────────────────────────────────────────────────────────────────────────────────────────
df_controls_genus <- physeq_16S_tumor %>%
    speedyseq::filter_sample_data(sample_type %in% c(
        "Negative control",
        "Water control",
        "Library control",
        "sample control"
    )) %>%
    speedyseq::tax_glom(taxrank = "Genus") %>%
    psmelt() %>%
    filter(Abundance > 0) %>%
    filter(x16s_seq_batch %in% c("Batch 3", "Batch 4", "Batch 5", "Batch 9")) %>%
    group_by(Genus, Sample) %>%
    summarize(Abundance = sum(Abundance))


# Get prevalence of Genus in the controls
df_prevalence_genus <- table(df_controls_genus$Genus) %>%
    data.frame() %>%
    mutate(Freq = (Freq / df_control_n[[1]]) * 100) %>%
    rename(Genus = Var1, Prevalence = Freq)


# Get the proportion of OTUs in the controls
df_proportion_prev_genus <- df_controls_genus %>%
    dplyr::select(Genus, Abundance) %>%
    group_by(Genus) %>%
    summarize(Abundance = sum(Abundance)) %>%
    mutate(Percentage = prop.table(Abundance) * 100) %>%
    inner_join(df_prevalence_genus, by = "Genus") %>%
    rename(Genus = Genus, Abundance = Abundance, Percentage = Percentage)


pt_prev_genus <- ggcptplot(sort(df_proportion_prev_genus$Prevalence, decreasing = T),
    Q = 1,
    change_in = "mean_var", cp_method = "BinSeg"
)

pt_prev_genus +
    labs(title = "Changes in mean variance (BinSeg)") +
    theme_bw(base_size = 20) +
    labs(caption = "Prevalance cutoff at Genus level") +
    ylab("Prevalence") +
    xlab("OTU") +
    geom_text(
        data = pt_prev_genus$data %>%
            filter(!is.na(cp_value)), aes(label = round(raw_value, digits = 3), x = row_number, y = raw_value),
        size = 6,
        inherit.aes = FALSE, nudge_y = 1.5
    )
ggsave("output/figures/prevalence_cutoff_genus_level.pdf", width = 10, height = 10, dpi = 300)


prev_cutoff <- pt_prev_genus$data %>%
    filter(!is.na(cp_value))

pt_prev_genus +
    xlim(c(0, 225)) +
    labs(title = "Changes in mean variance (BinSeg)") +
    theme_light(base_size = 20) +
    ylab("Abundance (%)") +
    labs(caption = "Prevalance cutoff at Genus level") +
    xlab("OTU") +
    geom_segment(
        aes(x = -Inf, xend = prev_cutoff$row_number, y = prev_cutoff$cp_value, yend = prev_cutoff$cp_value),
        linetype = "dashed", color = "blue"
    ) +
    geom_text(aes(x = prev_cutoff$row_number, y = prev_cutoff$cp_value),
        label = glue::glue(
            {
                round(prev_cutoff$cp_value, 4)
            },
            "%"
        ), color = "black", nudge_x = 20,
        size = 6
    )
ggsave("output/figures/prevalence_cutoff_genus_level_zoomed.pdf", width = 10, height = 10, dpi = 300)







# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Save the results for the prevalence and proportion cutoffs
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Combine the two cutoffs into one file
data.frame(
    genus_prev_cutoff = prev_cutoff[["raw_value"]],
    asv_prop_cutoff = prop_cutoff[["raw_value"]]
) %>%
    write_csv("data/processed_data/otu_cutoffs.csv")
