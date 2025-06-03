#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

#----------------------------------------------------------------------------------------------------
# Load and prepare heatmap matrix
#----------------------------------------------------------------------------------------------------

# Read CSV file with z-scored values for each gene module across conditions
heat_mat <- read_csv("Input_files/SMI/Glioma/figure5B_heatmap_matrix.csv") %>%
  dplyr::select(`16S- Neighborhood`, `16S+ Neighborhood`, `16S+ Tumor`) %>%
  as.matrix()

# Ensure matrix is numeric
mode(heat_mat) = "numeric"

# Define gene set labels (row names for the heatmap)
gene_sets <- c(
  "Glycolipid transport",
  "Macrophage proliferation",
  "Neuron development",
  "Axon development",
  "Glial cell migration",
  "Response to oxidative stress",
  "Regulation of innate immune\nresponses to cytosolic DNA",
  "Microbe-Host interactions",
  "Immune recruitment",
  "Lipid metabolism",
  "Stress response",
  "Angiogenesis",
  "Myeloid inflammation",
  "T cell activation",
  "T cell exhaustion"
)

rownames(heat_mat) <- gene_sets

#----------------------------------------------------------------------------------------------------
# Plot heatmap
#----------------------------------------------------------------------------------------------------

pdf(
  file = "Output_files/SMI/Figures/FigureS5B.pdf",
  width = 3.25,
  height = 6
)

ht <- Heatmap(
  heat_mat,
  # transpose the matrix
  name = "Z-score",
  col = colorRamp2(c(-1.5, 0, 1.5), c("#08306B", "white", "#99000D")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 90,
  row_names_gp = gpar(fontsize = 11),
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 11),
  rect_gp = gpar(col = "black", lwd = 0.5),
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    direction = "vertical",
    legend_height = unit(0.10, "cm") ,
    legend_width = unit(2, "cm")  # ← Make this longer
  )
)
draw(ht, heatmap_legend_side = "right")

dev.off()
