source("src/libraries.R")

load(file = "data/processed_data/physeq_SSO_clean.RData")
load(file = "data/processed_data/physeq_16S_tumor.RData")
load(file = "data/processed_data/physeq_WGS_tumor_saliva.RData")



df_sso_clean <- psmelt(physeq_SSO_clean)
df_16S_tumor <- psmelt(physeq_16S_tumor) %>%
    mutate(Kingdom = paste0(""), .after = "data_16s") %>%
    mutate(Species = paste0(""))
df_WGS_tumor_saliva <- psmelt(physeq_WGS_tumor_saliva)


rbind(df_sso_clean, df_16S_tumor, df_WGS_tumor_saliva) %>%
    filter(Abundance > 0) %>%
    write_csv("data/processed_data/all_samples_taxa_with_metadata.csv")
