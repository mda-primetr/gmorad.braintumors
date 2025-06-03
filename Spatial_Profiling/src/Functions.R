#----------------------------------------------------------------------------------------------------
# Commonly used functions
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Volcano plot (Gene)
#----------------------------------------------------------------------------------------------------

add_volcano_gene <- function(df, theme_size) {
    df$Color <- "Not Significant"
    df$Color[df$`Pr(>|t|)` < 0.05 & df$log2FC > 0.58] <- "Up-regulated"
    df$Color[df$`Pr(>|t|)` < 0.05 & df$log2FC < -0.58] <- "Down-regulated"
    df$Color <- factor(df$Color,
        levels = c(
            "Not Significant", "Up-regulated",
            "Down-regulated"
        )
    )

    plt <- ggplot(df, aes(
        x = log2FC, y = -log10(`Pr(>|t|)`),
        color = Color, label = Gene
    )) +
        geom_vline(xintercept = c(0.58, -0.58), lty = "dashed") +
        geom_hline(yintercept = -log10(0.05), lty = "dashed") +
        geom_hline(yintercept = -log10(max(df$`Pr(>|t|)`[df$FDR <= 0.05])), lty = "dashed") +
        geom_point(size = 2) +
        labs(
            y = expression("-log"[10] ~ "(p-value)"),
            color = "Significance"
        ) +
        scale_color_manual(values = c(
            `Not Significant` = "#9E9E9E",
            `Up-regulated` = "darkred",
            `Down-regulated` = "darkblue"
        )) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        theme_bw(base_size = {{ theme_size }}) +
        theme(
            legend.position = "none", legend.key.size = unit(0.3, "cm"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
        )

    plt
}

#----------------------------------------------------------------------------------------------------
# Volcano plot (Protein)
#----------------------------------------------------------------------------------------------------

add_volcano_protein <- function(df, theme_size) {
    df$Color <- "Not Significant"
    df$Color[df$FDR < 0.05 & df$log2FC > 0.58] <- "Up-regulated"
    df$Color[df$FDR < 0.05 & df$log2FC < -0.58] <- "Down-regulated"
    df$Color <- factor(df$Color,
        levels = c(
            "Not Significant", "Up-regulated",
            "Down-regulated"
        )
    )

    plt <- ggplot(df, aes(
        x = log2FC, y = -log10(FDR),
        color = Color, label = Protein
    )) +
        geom_vline(xintercept = c(0.58, -0.58), lty = "dashed") +
        geom_hline(yintercept = -log10(0.05), lty = "dashed") +
        # geom_hline(yintercept = -log10(max(df$`Pr(>|t|)`[df$FDR <= 0.05])), lty = "dashed") +
        geom_point(alpha = 0.7, size = 2) +
        labs(
            y = expression("-log"[10] ~ "(FDR)"),
            color = "Significance"
        ) +
        scale_color_manual(values = c(
            `Not Significant` = "#9E9E9E",
            `Up-regulated` = "darkred",
            `Down-regulated` = "darkblue"
        )) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        theme_bw(base_size = {{ theme_size }}) +
        theme(
            legend.position = "none", legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )
    plt
}

#----------------------------------------------------------------------------------------------------
# Signal-to-Noise (Protein - nanoString)
#----------------------------------------------------------------------------------------------------

snrOrder <- function(object, neg.names) {
    if (analyte(object) != "Protein") {
        stop("This function is only meant for protein data")
    }

    if (is.null(neg.names)) {
        neg.names <- iggNames(object)
    }

    raw <- exprs(object)

    # estimate background:
    negfactor <- apply(
        raw[neg.names, , drop = FALSE], 2,
        function(x) {
            pmax(mean(x), 1)
        }
    )

    # calc snr
    snr <- sweep(raw, 2, negfactor, "/")

    igginds <- which(is.element(rownames(snr), neg.names))
    o <- c(igginds, setdiff(order(apply(snr, 1, median)), igginds))

    return(snr[o, , drop = FALSE])
}


#----------------------------------------------------------------------------------------------------
# Plot Spatial Neighborhood Polygons by Category
#----------------------------------------------------------------------------------------------------

plot_neighborhood_polygons_by_label <- function(polygons_annotated,
                                                fill_colors,
                                                crop_left_um = 50,
                                                shrink_dist = -0.8,
                                                border_color = "grey70",
                                                border_size = 0.01,
                                                show_borders = TRUE,
                                                show_bounding_box = TRUE,
                                                legend_position = "bottom",
                                                output_file = NULL,
                                                plot_width = 3,
                                                plot_height = 3) {
  # Flip Y and convert to closed polygons
  polygons_closed <- polygons_annotated %>%
    mutate(y_um_flipped = -y_um) %>%
    group_by(cell_id, neighborhood_intra_label) %>%
    reframe(geometry = {
      x <- x_um
      y <- y_um_flipped
      if (length(x) >= 4 && all(!is.na(x)) && all(!is.na(y))) {
        if (x[1] != x[length(x)] || y[1] != y[length(y)]) {
          x <- c(x, x[1])
          y <- c(y, y[1])
        }
        list(st_polygon(list(cbind(x, y))))
      } else {
        list(st_polygon())
      }
    }) %>%
    st_as_sf()
  
  # Shrink
  polygons_shrunk <- polygons_closed %>%
    st_make_valid() %>%
    mutate(geometry = st_buffer(geometry, dist = shrink_dist)) %>%
    filter(!st_is_empty(geometry))
  
  # Crop
  xmin_cutoff <- st_bbox(polygons_shrunk)["xmin"] + crop_left_um
  polygons_cropped <- polygons_shrunk %>%
    mutate(bbox_xmin = map_dbl(geometry, ~ st_bbox(.x)["xmin"])) %>%
    filter(bbox_xmin > xmin_cutoff) %>%
    dplyr::select(-bbox_xmin)
  
  # Get plot bounds
  bounds <- st_bbox(polygons_cropped)
  
  border_col <- if (show_borders)
    border_color
  else
    NA
  border_sz  <- if (show_borders)
    border_size
  else
    0
  
  # Plot
  p <- ggplot() +
    geom_sf(
      data = polygons_cropped,
      aes(fill = neighborhood_intra_label),
      color = border_col,
      size = border_sz
    ) +
    scale_fill_manual(values = fill_colors, name = "Neighborhood") +
    coord_sf(
      xlim = c(bounds["xmin"], bounds["xmax"]),
      ylim = c(bounds["ymin"], bounds["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = legend_position,
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      legend.direction = if (legend_position == "bottom")
        "horizontal"
      else
        "vertical"
    )
  
  if (show_bounding_box) {
    p <- p + geom_rect(
      aes(
        xmin = bounds["xmin"],
        xmax = bounds["xmax"],
        ymin = bounds["ymin"],
        ymax = bounds["ymax"]
      ),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.5
    )
  }
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = plot_width, height = plot_height)
  }
  
  return(p)
}

#----------------------------------------------------------------------------------------------------
# Plot Spatial Neighborhood Polygons by AUCell score
#----------------------------------------------------------------------------------------------------

plot_auc_neighborhood_polygons <- function(polygons_annotated,
                                           auc_col,
                                           crop_left_um = 50,
                                           shrink_dist = -0.8,
                                           fill_label = "AUC Score",
                                           fill_palette = "magma",
                                           border_color = "grey70",
                                           border_size = 0.01,
                                           grey_color = "lightgrey",
                                           tumor_color = "#00E5FF",
                                           show_borders = TRUE,
                                           tumor_outline_thickness = 0.1,
                                           show_bounding_box = TRUE,
                                           legend_position = "bottom",
                                           use_gradient_for_tumor = FALSE,
                                           output_file = NULL,
                                           plot_width = 3,
                                           plot_height = 3) {
  auc_sym <- ensym(auc_col)
  
  # Attach AUC score using cell_id as the join key
  annotated_with_auc <- polygons_annotated %>%
    mutate(cell_id = paste0("c_1_", fov, "_", cellID)) %>%
    left_join(
      polygons_annotated %>%
        mutate(cell_id = paste0("c_1_", fov, "_", cellID)) %>%
        distinct(cell_id, auc_score = !!auc_sym),
      by = "cell_id"
    )
  
  # Convert to closed polygon geometries and flip y-axis
  polygons_closed <- annotated_with_auc %>%
    mutate(y_um_flipped = -y_um) %>%
    group_by(cell_id, neighborhood_intra_label, auc_score) %>%
    reframe(geometry = {
      x <- x_um
      y <- y_um_flipped
      if (length(x) >= 4 && all(!is.na(x)) && all(!is.na(y))) {
        if (x[1] != x[length(x)] || y[1] != y[length(y)]) {
          x <- c(x, x[1])
          y <- c(y, y[1])
        }
        list(st_polygon(list(cbind(x, y))))
      } else {
        list(st_polygon())
      }
    }) %>%
    st_as_sf()
  
  # Shrink each polygon slightly and ensure valid geometry
  polygons_shrunk <- polygons_closed %>%
    st_make_valid() %>%
    mutate(geometry = st_buffer(geometry, dist = shrink_dist)) %>%
    filter(!st_is_empty(geometry))
  
  # Crop by x_min threshold to exclude far-left polygons (zoom-in)
  xmin_cutoff <- st_bbox(polygons_shrunk)["xmin"] + crop_left_um
  polygons_cropped <- polygons_shrunk %>%
    mutate(bbox_xmin = map_dbl(geometry, ~ st_bbox(.x)["xmin"])) %>%
    filter(bbox_xmin > xmin_cutoff) %>%
    dplyr::select(-bbox_xmin)
  
  # Subset polygons by label
  excluded_cells <- polygons_cropped %>% filter(neighborhood_intra_label == "Excluded")
  scored_cells   <- polygons_cropped %>% filter(neighborhood_intra_label %in% c("16S+ Neighborhood", "16S- Neighborhood"))
  tumor_cells    <- polygons_cropped %>% filter(neighborhood_intra_label == "16S+ Tumor")
  bounds         <- st_bbox(polygons_cropped)
  
  # Border style logic
  border_col     <- if (show_borders)
    border_color
  else
    NA
  border_sz      <- if (show_borders)
    border_size
  else
    0
  tumor_border   <- if (show_borders)
    "grey80"
  else
    NA
  
  # Build plot
  p <- ggplot() +
    geom_sf(
      data = excluded_cells,
      fill = grey_color,
      color = border_col,
      size = border_sz
    ) +
    geom_sf(
      data = scored_cells,
      aes(fill = auc_score),
      color = border_col,
      size = border_sz
    )
  
  # Tumor cells: fixed or gradient fill
  if (use_gradient_for_tumor) {
    p <- p + geom_sf(
      data = tumor_cells,
      aes(fill = auc_score),
      color = tumor_border,
      size = tumor_outline_thickness
    )
  } else {
    p <- p + geom_sf(
      data = tumor_cells,
      fill = tumor_color,
      color = tumor_border,
      size = tumor_outline_thickness
    )
  }
  
  # Optional bounding box
  if (show_bounding_box) {
    p <- p + geom_rect(
      aes(
        xmin = bounds["xmin"],
        xmax = bounds["xmax"],
        ymin = bounds["ymin"],
        ymax = bounds["ymax"]
      ),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.5
    )
  }
  
  # Color scale and formatting
  p <- p +
    scale_fill_viridis_c(
      option = fill_palette,
      name = fill_label,
      na.value = grey_color,
      guide = "colorbar"
    ) +
    coord_sf(
      xlim = c(bounds["xmin"], bounds["xmax"]),
      ylim = c(bounds["ymin"], bounds["ymax"]),
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = legend_position,
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      legend.direction = if (legend_position == "bottom")
        "horizontal"
      else
        "vertical"
    )
  
  # Adjust legend for bottom layout
  if (legend_position == "bottom") {
    p <- p + guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(plot_width * 0.6, "in"),
      barheight = unit(0.2, "in")
    ))
  }
  
  # Save if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = plot_width, height = plot_height)
  }
  
  return(p)
}
