suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(tools)
    library(RColorBrewer)
    ori_par <- par(no.readonly = TRUE)
})


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

#' Plot seurat UMAP colored by correlation
plot_with_correlation <- function(
    plotting_df = SS_MQD_plotting_df,
    cor_colname,
    output_without_margin,
    sorted_order,
    ...
) {
    x <- plotting_df$UMAP.1
    y <- plotting_df$UMAP.2
    cor_vector <- plotting_df[, cor_colname] # "Correlation_PtrFiber1_1"

    max_bound <- quantile(cor_vector, 0.98)
    upper_bound <- max_bound * 0.9
    lower_bound <- max_bound * 0.66
    color_index <-
        ((cor_vector - lower_bound) / (upper_bound - lower_bound)) %>%
        pmax(0) %>%
        pmin(1) %>%
        multiply_by(500) %>%
        round(digits = 0) %>%
        add(1)

    # color_tick <- rev(brewer.pal(11, "Spectral"))
    # color_pool <- colorRampPalette(color_tick[7:11])(501)
    color_tick <- c("#EEF2F9", "#C44233")
    color_pool <- colorRampPalette(color_tick)(501)

    if (sorted_order) {
        plot_order <- order(cor_vector)
    } else {
        plot_order <- sample(length(x))
    }
    plot(x[plot_order], y[plot_order],
        col = color_pool[color_index[plot_order]],
        pch = 20, cex = 0.5,
        xlab = "UMAP_1", ylab = "UMAP_2",
        main = ifelse(
            output_without_margin, "",
            paste0(
                "Min:", round(min(cor_vector), 2), " ",
                "LB:", round(lower_bound, 2), " ",
                "UB:", round(upper_bound, 2), " ",
                "Max:", round(max(cor_vector), 2)
            )
        ),
        axes = !output_without_margin, las = 1
    )

    par(xpd = TRUE)
    current_xlim <- par()$usr[1:2]
    current_ylim <- par()$usr[3:4]
    legend_x_range <- c(sum(current_xlim * c(0.25, 0.75)), current_xlim[2])
    legend_x_vector <-
        seq(
            from = legend_x_range[1],
            to = legend_x_range[2],
            length.out = length(color_pool) + 1
        )
    legend_y_range <-
        c(sum(current_ylim * c(-0.1, 1.1)), sum(current_ylim * c(-0.15, 1.15)))
    for (i in seq_along(color_pool)) {
        lines(
            legend_x_vector[c(i, i + 1)],
            rep(legend_y_range[1], 2),
            col = color_pool[i],
            lwd = 20, lend = "square"
        )
    }
    text(
        legend_x_range[1], legend_y_range[2],
        sprintf("%.2f", lower_bound), cex = 1
    )
    text(
        legend_x_range[2], legend_y_range[2],
        sprintf("%.2f", upper_bound), cex = 1
    )
    par(xpd = FALSE)
}


# Set parameters ===============================================================
#' Get input parameters from command line
input_SS_MQD_plotting_csv <- snakemake@input$SS_MQD_plotting_csv
output_figure_folder <- snakemake@output$figure_folder


# Implementation ===============================================================
#' Create output directory
if (!dir.exists(output_figure_folder)) {
    dir.create(output_figure_folder, recursive = TRUE)
}

#' Input plotting information
SS_MQD_plotting_df <- read.csv(input_SS_MQD_plotting_csv)
is_correlation_column <- grepl("Correlation_", colnames(SS_MQD_plotting_df))

# Output figures
correlation_columns <- colnames(SS_MQD_plotting_df)[is_correlation_column]
for (selected_column in correlation_columns) {
    output_png_figure(
        plotting_function = plot_with_correlation,
        plotting_df = SS_MQD_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, ".png"),
        output_without_margin = FALSE,
        sorted_order = FALSE
    )
    output_png_figure(
        plotting_function = plot_with_correlation,
        plotting_df = SS_MQD_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, "_Clear.png"),
        output_without_margin = TRUE,
        sorted_order = FALSE
    )
    output_png_figure(
        plotting_function = plot_with_correlation,
        plotting_df = SS_MQD_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, "_Sorted.png"),
        output_without_margin = FALSE,
        sorted_order = TRUE
    )
    output_png_figure(
        plotting_function = plot_with_correlation,
        plotting_df = SS_MQD_plotting_df,
        cor_colname = selected_column,
        output_figure = TRUE,
        output_path =
            paste0(output_figure_folder, "/", selected_column, "_Sorted_Clear.png"),
        output_without_margin = TRUE,
        sorted_order = TRUE
    )
}
