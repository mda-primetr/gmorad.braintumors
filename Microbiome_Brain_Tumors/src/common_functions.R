library(ggprism)
library(ggrepel)
library(ggtext)
library(gt)
library(gtExtras)
library(microViz)
library(phyloseq)
library(speedyseq)
library(tidyverse)
library(vegan)


#' Function to get data frame for Top N taxa at specified level without any groupings
#'
#' @param physeq phyloseq object
#' @param no_of_taxa Integer to indicate the number of taxa for the grouping
#' @param taxa_level Taxa level to group the dataframe
#'
#' @return dataframe containing top taxa at specified taxa level


get_top_n_taxa_df <- function(physeq, no_of_taxa, taxa_level) {
    df_meta <- data.frame(sample_data(physeq)) %>%
        rownames_to_column(var = "SampleID")

    df_top_n_taxa <- speedyseq::psmelt(physeq) %>%
        group_by(Sample, {{ taxa_level }} := fct_lump(.data[[taxa_level]], n = no_of_taxa, w = Abundance)) %>%
        dplyr::summarize(count = sum(Abundance)) %>%
        mutate(prop = ((count) / sum(count)) * 100) %>%
        inner_join(df_meta, by = c("Sample" = "SampleID")) %>%
        ungroup()

    df_top_n_taxa
}


#' Function to get the top N taxa proportion by grouping
#' @param physeq phyloseq object
#' @param no_of_taxa Integer to indicate the number of taxa for the grouping
#' @param taxa_level Taxa level to group the dataframe
#' @param groups any one grouping variable in the phyloseq metadata
#'
#' @return dataframe containing top N taxa with proportion grouped by the specified variable

get_top_n_taxa_by_group_df <- function(physeq, no_of_taxa, taxa_level, groups) {
    # Glom the physeq to the taxa level
    physeq <- speedyseq::tax_glom(physeq, taxrank = taxa_level)

    # Get the meta data
    df_meta <- data.frame(sample_data(physeq)) %>%
        rownames_to_column(var = "SampleID")

    # Get the top overall taxa without grouping
    df_top_n_taxa <- speedyseq::psmelt(physeq) %>%
        group_by(Sample, {{ taxa_level }} := fct_lump(.data[[taxa_level]], n = no_of_taxa, w = Abundance)) %>%
        dplyr::summarize(count = sum(Abundance)) %>%
        mutate(prop = ((count) / sum(count)) * 100) %>%
        inner_join(df_meta, by = c("Sample" = "SampleID")) %>%
        ungroup()

    # Get the list of taxa to be filtered to be looked into grouping
    df_taxa_unique <- unique(df_top_n_taxa[[taxa_level]])

    # Get the final table with proportions
    df_taxa_by_group <- speedyseq::psmelt(physeq) %>%
        dplyr::select(Sample, {{ taxa_level }}, {{ groups }}, Abundance) %>%
        group_by(Sample, !!sym(groups), !!sym(taxa_level), .drop = FALSE) %>% # Order matters
        dplyr::summarize(cnt = sum(Abundance)) %>%
        mutate(prop = (cnt / sum(cnt)) * 100) %>%
        filter(!!sym(taxa_level) %in% df_taxa_unique)

    df_taxa_by_group
}




# Function to open table in excel
view_xl <- function(.data) {
    if (interactive()) {
        tmp <- tempfile(fileext = ".csv")
        data.table::fwrite(.data, tmp)
        browseURL(tmp, browser = 'open -a "Microsoft Excel"')
    }
}


# convert psmelted data to phyloseq object

melt_to_physeq <- function(melted_physeq) {
    df_taxa <- melted_physeq %>%
        dplyr::select(OTU, any_of(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) %>%
        dplyr::distinct(OTU, .keep_all = TRUE) %>%
        column_to_rownames("OTU") %>%
        as.matrix()

    df_OTU <- melted_physeq %>%
        dplyr::select(OTU, sample_id, Abundance) %>%
        pivot_wider(names_from = sample_id, values_from = Abundance, values_fill = 0) %>%
        column_to_rownames("OTU") %>%
        as.matrix()

    df_meta <- melted_physeq %>%
        dplyr::select(-any_of(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Abundance"))) %>%
        distinct(sample_id, .keep_all = TRUE) %>%
        column_to_rownames("sample_id")

    OTU <- otu_table(df_OTU, taxa_are_rows = TRUE)
    TAX <- tax_table(df_taxa)
    physeq <- phyloseq(OTU, TAX, sample_data(df_meta))

    return(physeq)
}



# Prevalence by sample type
#' Get prevalence of the taxa by sample parameters
#' @param physeq phyloseq object
#' @param sample_param column name or metadata column from the phyloseq sample_data
#' @param taxa_rank Taxa rank on which we want to calculate prevalance
get_phyloseq_prev_by_group <- function(physeq, sample_param, taxa_rank) {
    df_prev_sam <- sample_data(physeq) %>%
        group_by(.data[[sample_param]]) %>%
        dplyr::summarise(sam_sum = n())

    data.frame(speedyseq::psmelt(physeq)) %>%
        filter(Abundance > 0) %>%
        group_by(.data[[sample_param]], .data[[taxa_rank]]) %>%
        dplyr::summarise(counts = sum(Abundance > 0, na.rm = TRUE)) %>%
        inner_join(df_prev_sam, by = c(sample_param)) %>%
        mutate(prevalence = counts / sam_sum) %>%
        ungroup()
}



# Function to run scrub on the phyloseq object -----
# @param phy_samples phyloseq object for the samples
# @param phy_controls phyloseq object for the controls
# @param batch_no batch number for the samples
# @param lst_samples list of samples to be used for the scrub

fn_run_scrub <- function(phy_samples, phy_controls, batch_no, lst_samples) {
    # Get sample phyloseq
    p_sam <- phy_samples %>%
        speedyseq::filter_sample_data(x16s_seq_batch == batch_no) %>%
        prune_samples(sample_names(.) %in% lst_samples, .) %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        prune_samples(sample_sums(.) > 0, .)

    # Get control phyloseq
    p_ctrl <- phy_controls %>%
        speedyseq::filter_sample_data(x16s_seq_batch == batch_no) %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        prune_samples(sample_sums(.) > 0, .)

    # Combine the phyloseq objects for the samples and controls
    p_sam_ctrl <- rbind(psmelt(p_sam), psmelt(p_ctrl)) %>%
        melt_to_physeq() %>%
        speedyseq::filter_sample_data(!is.na(sample_well)) %>%
        prune_taxa(taxa_sums(.) > 0, .) %>%
        prune_samples(sample_sums(.) > 0, .)

    # Create the metadata for the SCRuB
    df_sam_scrb <- sample_data(p_sam_ctrl) %>%
        data.frame() %>%
        mutate(sample_well = gsub(" ", "", sample_well)) %>%
        dplyr::select(is_control, sample_type = sample_type_for_well, sample_well) %>%
        filter(!is.na(sample_well))

    # Create the taxa matrix for the SCRuB
    df_taxa_scrb <- otu_table(p_sam_ctrl) %>%
        t() %>%
        data.frame() %>%
        as.matrix()

    # Run the SCRuB
    scrb_output <- SCRuB::SCRuB(df_taxa_scrb,
        metadata = df_sam_scrb,
        verbose = TRUE
    )

    return(scrb_output)
}


#' Update phyloseq object with scrub results
#'
#' This function updates a phyloseq object by replacing the OTU table with the decontaminated samples from scrub results.
#'
#' @param physeq The original phyloseq object.
#' @param scrub_results The scrub results containing the decontaminated samples.
#'
#' @return The updated phyloseq object with the OTU table replaced by the decontaminated samples.
fn_update_physeq_with_scrub_results <- function(physeq, scrub_results) {
    otu_table(physeq) <- t(as.matrix(scrub_results$decontaminated_samples)) %>% otu_table(taxa_are_rows = T)
    return(physeq)
}



# Clean phyloseq taxa table ----
fn_clean_tax_table <- function(physeq) {
    physeq_clean <- tax_table(physeq) %>%
        data.frame() %>%
        mutate(Kingdom = gsub("k__", "", Kingdom)) %>%
        mutate(Phylum = gsub("p__", "", Phylum)) %>%
        mutate(Class = gsub("c__", "", Class)) %>%
        mutate(Order = gsub("o__", "", Order)) %>%
        mutate(Family = gsub("f__", "", Family)) %>%
        mutate(Genus = gsub("g__", "", Genus)) %>%
        mutate(Species = gsub("s__", "", Species)) %>%
        mutate(Species = gsub("_", " ", Species)) %>%
        as.matrix()

    tax_table(physeq) <- physeq_clean
    physeq
}



# For ancomBC 1
get_ancombc_df <- function(ancombc_out, taxa_level, group_var, group_var_order) {
    # Extract lfc values
    df_ancombc_lfc <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("lfc_")) %>%
        dplyr::select(-c("lfc_intercept"), taxon = ends_with("taxon")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "lfc") %>%
        mutate(metadata = str_remove(metadata, "lfc_"))

    # Extract se values
    df_ancombc_se <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("se_")) %>%
        dplyr::select(-c("se_intercept"), taxon = ends_with("taxon")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "stderr") %>%
        mutate(metadata = str_remove(metadata, "se_"))


    # Extract p-values
    df_ancombc_pval <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("p_")) %>%
        dplyr::select(-c("p_val_intercept"), taxon = ends_with("taxon")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "pval") %>%
        mutate(metadata = str_remove(metadata, "p_val_"))

    # Extract q-values
    df_ancombc_qval <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("q_")) %>%
        dplyr::select(-c("q_val_intercept"), taxon = ends_with("taxon")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "qval") %>%
        mutate(metadata = str_remove(metadata, "q_val_"))



    df_for_heatmap <- df_ancombc_lfc %>%
        inner_join(df_ancombc_se, by = c("taxon", "metadata")) %>%
        inner_join(df_ancombc_pval, by = c("taxon", "metadata")) %>%
        inner_join(df_ancombc_qval, by = c("taxon", "metadata")) %>%
        rowwise() %>%
        mutate(value = toupper(strsplit(metadata, "_(?!.*_)", perl = TRUE)[[1]][2])) %>%
        mutate(metadata = strsplit(metadata, "_(?!.*_)", perl = TRUE)[[1]][1]) %>%
        dplyr::select(feature = taxon, metadata, value, log2fc = lfc, stderr, pval, qval) %>%
        mutate(qval = as.numeric(qval)) %>%
        mutate(pval = as.numeric(pval)) %>%
        mutate(stderr = as.numeric(stderr))
}


# For AncomBC 2 ----
get_ancombc2_df <- function(ancombc_out, taxa_level, group_var, group_var_order) {
    # Extract lfc values
    df_ancombc_lfc <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("lfc_"), taxon) %>%
        dplyr::select(-c("lfc_intercept")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "lfc") %>%
        mutate(metadata = str_remove(metadata, "lfc_"))

    # Extract se values
    df_ancombc_se <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("se_"), taxon) %>%
        dplyr::select(-c("se_intercept")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "stderr") %>%
        mutate(metadata = str_remove(metadata, "se_"))


    # Extract p-values
    df_ancombc_pval <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("p_"), taxon) %>%
        dplyr::select(-c("p_intercept")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "pval") %>%
        mutate(metadata = str_remove(metadata, "p_"))

    # Extract q-values
    df_ancombc_qval <- ancombc_out %>%
        data.frame() %>%
        janitor::clean_names() %>%
        select(starts_with("q_"), taxon) %>%
        dplyr::select(-c("q_intercept")) %>%
        pivot_longer(-c("taxon"), names_to = "metadata", values_to = "qval") %>%
        mutate(metadata = str_remove(metadata, "q_"))



    df_for_heatmap <- df_ancombc_lfc %>%
        inner_join(df_ancombc_se, by = c("taxon", "metadata")) %>%
        inner_join(df_ancombc_pval, by = c("taxon", "metadata")) %>%
        inner_join(df_ancombc_qval, by = c("taxon", "metadata")) %>%
        rowwise() %>%
        mutate(value = toupper(strsplit(metadata, "_(?!.*_)", perl = TRUE)[[1]][2])) %>%
        mutate(metadata = strsplit(metadata, "_(?!.*_)", perl = TRUE)[[1]][1]) %>%
        dplyr::select(feature = taxon, metadata, value, log2fc = lfc, stderr, pval, qval)
}
