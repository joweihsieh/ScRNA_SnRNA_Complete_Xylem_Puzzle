library(Seurat)
library(magrittr)
library(dplyr)
ori_par <- par(no.readonly = TRUE)


# Define functions =============================================================
#' Output figure
output_png_figure <- function(
    plotting_function,
    # x, y, col, main = "",
    output_figure = FALSE,
    output_path = "temp.png",
    output_without_margin = FALSE,
    ...
) {
    if (output_figure) {
        png(output_path,
            pointsize = 10, res = 300,
            width = 20, height = 15, units = "cm")
    }

    if (output_without_margin) {
        par(mai = c(0, 0, 0, 0))
    } else {
        par(mai = ori_par$mai)
    }

    plotting_function(
        output_figure = output_figure,
        output_path = output_path,
        output_without_margin = output_without_margin,
        ...
    )

    par(mai = ori_par$mai)

    if (output_figure) dev.off()
}

#' Plot seurat UMAP colored on each sample
plotting_function <- function(
    plotting_df,
    output_without_margin,
    ...
) {
    is_not_TenX_Ptr <- (plotting_df$Sample != "TenX_Ptr")
    plotting_df$SS_Color[is_not_TenX_Ptr] <-
        paste0(plotting_df$SS_Color[is_not_TenX_Ptr], "00")
    x <- plotting_df$UMAP.1
    y <- plotting_df$UMAP.2
    col <- plotting_df$SS_Color

    randam_order <- sample(length(x))
    plot(
        x[randam_order], y[randam_order],
        col = col[randam_order], pch = 20, cex = 0.5,
        xlab = "UMAP_1", ylab = "UMAP_2",
        main = ifelse(output_without_margin, "", "main"),
        axes = !output_without_margin, las = 1
    )
}


# Set parameters ===============================================================
#' Get input parameters from command line
input_MS_plotting_csv <- snakemake@input$MS_plotting_csv
output_figure_folder <- snakemake@output$figure_folder

# input_MS_plotting_csv <-
#     "results/Multi_species_analysis/all_plotting_tables_addSS/plotting_PtrPal_seed_42_md_0.3_nn_30.csv"


# Implementation ===============================================================
#' Create output directory
if (!dir.exists(output_figure_folder)) {
    dir.create(output_figure_folder, recursive = TRUE)
}

#' Input MS plotting information
MS_plotting_df <- read.csv(input_MS_plotting_csv)

#' Output figure
output_png_figure(
    plotting_function,
    plotting_df = MS_plotting_df,
    output_figure = TRUE,
    output_path =
        paste0(
            output_figure_folder,
            "/With_TenX_Ptr_SScolor_Clear.png"
        ),
    output_without_margin = TRUE,
    output_without_legend = TRUE
)
