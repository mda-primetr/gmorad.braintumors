#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
source("src/Functions.R")

# Set reproducibility seed
set.seed(42)

#----------------------------------------------------------------------------------------------------
# Figure 5D – Annotated spatial plot of tumor and neighborhoods
#----------------------------------------------------------------------------------------------------

# Define neighborhood colors
neighborhood_colors <- c(
  "16S+ Tumor" = "#0EFBFF", # bright cyan
  "16S+ Neighborhood" = "#FF00FD", # magenta
  "16S- Neighborhood" = "#218C24",# dark green
  "Excluded" = "lightgrey" # background/other
)

# Upload polygon coordinates with metadata
polygons_annotated<- read_csv("Input_files/SMI/BrM/figure5d_polygons.csv")

# Plot spatial polygons with group-based fill
plot_neighborhood_polygons_by_label(
  polygons_annotated = polygons_annotated,
  fill_colors = neighborhood_colors,
  show_borders = FALSE,  
  crop_left_um = 50,
  shrink_dist = -0.8,
  plot_width = 3,
  plot_height = 2,
  legend_position = "none",
  output_file = "Output_files/SMI/Figures/Figure5D.pdf")


#----------------------------------------------------------------------------------------------------
# Figure 5G1 – Microbe–Host Interaction
#----------------------------------------------------------------------------------------------------

# Upload polygon coordinates with metadata
polygons_annotated<- read_csv("Input_files/SMI/BrM/figure5g1_polygons.csv")

# Spatial plot
plot_auc_neighborhood_polygons(
  polygons_annotated = polygons_annotated,
  auc_col = MicrobeHost_Interaction_AUC,
  crop_left_um = 50,
  shrink_dist = -1.0,
  show_borders = FALSE,
  tumor_outline_thickness = 0.1,
  fill_label = "Microbial-Host",
  legend_position = "bottom",
  plot_width = 3,
  plot_height = 2.5,
  output_file = "Output_files/SMI/Figures/Figure5G1.pdf")

#----------------------------------------------------------------------------------------------------
# Figure 5G2 – Immune Recruitment
#----------------------------------------------------------------------------------------------------

# Upload polygon coordinates with metadata
polygons_annotated<- read_csv("Input_files/SMI/BrM/figure5g2_polygons.csv")

# Spatial plot
plot_auc_neighborhood_polygons(
  polygons_annotated = polygons_annotated,
  auc_col = TumorImmuneRecruitment_AUC,
  crop_left_um = 50,
  shrink_dist = -1.0,
  show_borders = FALSE,
  tumor_outline_thickness = 0.1,
  fill_label = "Immune recruitment",
  legend_position = "bottom",
  plot_width = 3,
  plot_height = 2.5,
  output_file = "Output_files/SMI/Figures/Figure5G2.pdf")

#----------------------------------------------------------------------------------------------------
# Figure 5G3 – Lipid Metabolism
#----------------------------------------------------------------------------------------------------

# Upload polygon coordinates with metadata
polygons_annotated<- read_csv("Input_files/SMI/BrM/figure5g3_polygons.csv")

#Spatial plot
plot_auc_neighborhood_polygons(
  polygons_annotated = polygons_annotated,
  auc_col = TumorLipid_Metabolism_AUC,
  crop_left_um = 50,
  shrink_dist = -1.0,
  show_borders = FALSE,
  tumor_outline_thickness = 0.1,
  fill_label = "Lipid",
  legend_position = "bottom",
  plot_width = 3,
  plot_height = 2.5,
  output_file = "Output_files/SMI/Figures/Figure5G3.pdf")

#----------------------------------------------------------------------------------------------------
# Figure 5G4 – Stress
#----------------------------------------------------------------------------------------------------

# Upload polygon coordinates with metadata
polygons_annotated<- read_csv("Input_files/SMI/BrM/figure5g4_polygons.csv")

# Spatial plot
plot_auc_neighborhood_polygons(
  polygons_annotated = polygons_annotated,
  auc_col = stress_apoptosis_genes_AUC,
  crop_left_um = 50,
  shrink_dist = -1.0,
  show_borders = FALSE,
  tumor_outline_thickness = 0.1,
  fill_label = "Stress",
  legend_position = "bottom",
  plot_width = 3,
  plot_height = 2.5,
  output_file = "Output_files/SMI/Figures/Figure5G4.pdf")

#----------------------------------------------------------------------------------------------------
# Figure 5G5 – Angiogenesis
#----------------------------------------------------------------------------------------------------

# Upload polygon coordinates with metadata
polygons_annotated<- read_csv("Input_files/SMI/BrM/figure5g5_polygons.csv")

#Spatial plot
plot_auc_neighborhood_polygons(
  polygons_annotated = polygons_annotated,
  auc_col = TumorAngiogenesis_AUC,
  crop_left_um = 50,
  shrink_dist = -1.0,
  show_borders = FALSE,
  tumor_outline_thickness = 0.1,
  fill_label = "Angiogenesis",
  legend_position = "bottom",
  plot_width = 3,
  plot_height = 2.5,
  output_file = "Output_files/SMI/Figures/Figure5G5.pdf")
