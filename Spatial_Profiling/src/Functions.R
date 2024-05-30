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
# Volcano plot (rotein)
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
