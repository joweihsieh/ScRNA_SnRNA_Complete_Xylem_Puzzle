suppressPackageStartupMessages({
  library(dplyr)
  library(slingshot)
  library(magrittr)
  library(Seurat)
  library(reshape2)
  library(RColorBrewer)

})


get_adj <- function(lin, cluster_ids) {
  cluster_ids <- sort(unique(as.character(cluster_ids))) 
  adj <- matrix(0, length(cluster_ids), length(cluster_ids))
  rownames(adj) <- cluster_ids
  colnames(adj) <- cluster_ids
  for (L in lin) {
    L <- as.character(L)  
    for (Ci in seq(length(L) - 1)) {
      adj[L[Ci], L[Ci + 1]] <- 1
      adj[L[Ci + 1], L[Ci]] <- 1
    }
  }
  return(adj)
}


#This function is written in https://menugget.blogspot.com/
#2011/11/define-color-steps-for-colorramppalette.html
color_palette <- function(steps, n_steps_between = NULL, ...) {

  if (is.null(n_steps_between)) n_steps_between <- rep(0, (length(steps) - 1))
  if (length(n_steps_between) != length(steps) - 1) stop
  ("Must have one less n_steps_between value than steps")

  fill_steps <- cumsum(rep(1, length(steps)) + c(0, n_steps_between))
  rgb <- matrix(NA, nrow = 3, ncol = fill_steps[length(fill_steps)])
  rgb[, fill_steps] <- col2rgb(steps)

  for (i in which(n_steps_between > 0)){
    col_start <- rgb[, fill_steps[i]]
    col_end <- rgb[, fill_steps[i + 1]]
    for (j in seq(3)){
      vals <- seq(col_start[j],
                  col_end[j],
                  length.out = n_steps_between[i] +
                    2)[2:(2 + n_steps_between[i] - 1)]
      rgb[j, (fill_steps[i] + 1):(fill_steps[i + 1] - 1)] <- vals
    }
  }

  new_steps <- rgb(rgb[1, ], rgb[2, ], rgb[3, ], maxColorValue = 255)
  pal <- colorRampPalette(new_steps, ...)
  return(pal)
}

cal_cluster_center <- function(plotting_df) {
  rd <- as.matrix(plotting_df[, c("UMAP.1", "UMAP.2")])
  cl <- factor(plotting_df$Cluster)
  x_val <- sapply(levels(cl), function(c) mean(rd[cl == c, 1]))
  y_val <- sapply(levels(cl), function(c) mean(rd[cl == c, 2]))
  cluster_center <- cbind(x_val, y_val)

  return(cluster_center)
}

keep_center_cells <- function(plotting_df, percent_keep) {
  new_df <- data.frame()
  clusters <- levels(factor(plotting_df$Cluster))
  for (cluster in clusters) {
    df_cells <- plotting_df[plotting_df$Cluster == cluster, ]
    center <- c(mean(df_cells$UMAP.1), mean(df_cells$UMAP.1))
    dis_to_center <- apply(df_cells[, c("UMAP.1", "UMAP.2")],
                           1, function(x) euclidean_dist(x, center))
    keep_dis <- quantile(dis_to_center, percent_keep)
    keep_df <- df_cells[dis_to_center < keep_dis, ]
    new_df <- rbind(new_df, keep_df)

  }
  return(new_df)
}

euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))

get_cruve_color <- function(crvs, lineage, plotting_df) {
  cluster_center <- cal_cluster_center(plotting_df)
  lin <- crvs@lineages[[lineage]]
  crv <- crvs@curves[[lineage]]$s
  #col <- sapply(lin,
                #function(x) plotting_df[plotting_df$Cluster == x, "Color"][1])


  col <- sapply(lin, function(x) {
    unique(plotting_df$Color[plotting_df$Cluster == x])[1]
  })
  anchors <- c(1)
  print(lin)
  if (length(lin) > 2) {
    for (i in 2:(length(lin) - 1)) {
      c_center <- cluster_center[lin[i], ]
      dis_to_c <- apply(crv, 1, function(x) euclidean_dist(x, c_center))
      anchors <- append(anchors, which.min(dis_to_c))
    }
  }
  anchors <- append(anchors, nrow(crv))
  steps <- c()

  for (i in 1:(length(anchors) - 1)) {
    steps <- append(steps, anchors[i + 1] - anchors[i])
  }

  pal <- color_palette(col, steps, space = "rgb")
  cols <- pal(nrow(crv))
  color_curve <- list(pos = crv, col = cols)
  return(color_curve)
}

get_cruve_last <- function(crvs, lineage, plotting_df) {
  radius_keep <- 1
  cluster_center <- cal_cluster_center(plotting_df)
  lin <- crvs@lineages[[lineage]]
  crv <- crvs@curves[[lineage]]$s
  last_cluster <- tail(lin, n = 1)
  last_cluster_center <- cluster_center[last_cluster, ]
  col <- plotting_df[plotting_df$Cluster == last_cluster, "Color"][1]
  dis_to_c <- apply(crv, 1, function(x) euclidean_dist(x, last_cluster_center))
  keep <- crv[dis_to_c < radius_keep, ]
  return(list(pos = keep, col = col))
}

plot_umap <- function(
    plotting_df, cluster = 1) {
  plotting_df$Color <- "gray"
  plotting_df[plotting_df$Cluster == cluster, "Color"] <- "black"

  output_without_margin <- FALSE
  png(filename = paste0("UMAP_cla_", cluster, ".png"),
      width = 2800, height = 2000, res = 400)
  plot(
    x = plotting_df$UMAP.1,
    y = plotting_df$UMAP.2,
    #col = plotting_df$Color,
    col = "gray",
    pch = 20,
    cex = 0.2,
    axes = !output_without_margin, las = 1,
    xlim = c(-10, 10),
    ylim = c(-10, 10),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )
  title(main = paste("UMAP cla Cluster ", cluster))
  dev.off()
}

plot_umap_color <- function(
    plotting_df, title = NA) {
  cex <- ifelse(plotting_df$Color == "gray", 0.2, 0.5)
  print(table(plotting_df$Cluster))
  output_without_margin <- FALSE
  png(filename = paste0("UMAP_", title, ".png"),
      width = 2800, height = 2000, res = 400)
  plot(
    x = plotting_df$UMAP.1,
    y = plotting_df$UMAP.2,
    col = plotting_df$Color,
    pch = 20,
    cex = cex,
    axes = !output_without_margin, las = 1,
    xlim = c(-10, 10),
    ylim = c(-10, 10),
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )
  title(main = paste("UMAP ", title))
  dev.off()
}
plot_slingshot <- function(
    slingshot_df,                     
    lins = NULL,
    lins_last = NULL,
    plotting_df = NULL,              
    title = "Test",
    test_mode = FALSE) {

  output_without_margin <- FALSE
  png(filename = paste0("Slingshot_", title, ".png"),
      width = 2800, height = 2000, res = 400)

  ## 
  plot(
    x = plotting_df$UMAP.1,
    y = plotting_df$UMAP.2,
    col = "gray",
    pch = 20,
    cex = 0.2,
    axes = !output_without_margin, las = 1,
    ylab = "UMAP.2",
    xlab = "UMAP.1"
  )

  ## 
  if (length(lins) > 0) {
    rd <- as.matrix(slingshot_df[, c("UMAP.1", "UMAP.2")])
    cl <- slingshot_df$Cluster
    adj <- get_adj(lins, cl)
    sds <- newSlingshotDataSet(rd, cl, lineages = lins, adjacency = adj)
    crvs <- SlingshotDataSet(getCurves(sds, approx_points = 300))

    for (Lineage in names(lins)) {
      color_curve <- get_cruve_color(crvs, Lineage, plotting_df)
      points(color_curve$pos, col = color_curve$col, pch = 20, cex = 4)
    }
  }

  ## 
  if (length(lins_last) > 0) {
    rd <- as.matrix(slingshot_df[, c("UMAP.1", "UMAP.2")])
    cl <- slingshot_df$Cluster
    adj <- get_adj(lins_last, cl)
    sds <- newSlingshotDataSet(rd, cl, lineages = lins_last, adjacency = adj)
    crvs <- SlingshotDataSet(getCurves(sds))

    for (Lineage in names(lins_last)) {
      color_curve <- get_cruve_last(crvs, Lineage, plotting_df)
      points(color_curve$pos, col = color_curve$col, pch = 20, cex = 4)
    }
  }

  ## 
  if (test_mode) {
    cluster_center <- cal_cluster_center(plotting_df)
    for (i in rownames(cluster_center)) {
      x <- cluster_center[i, 1]
      y <- cluster_center[i, 2]
      text(x, y, labels = i, font = 6)
    }
  }

  title(main = title)
  dev.off()
}

project_umap2 <- read.csv("/Users/joweihsieh/Dropbox/YCL/Single_Cell_snRNA_Quanzi/code/hero_shot.csv")


library(RColorBrewer)

#cluster_ids <- as.character(sort(unique(project_umap2$Cluster)))
project_umap2$Cluster <- as.character(project_umap2$Cluster)
project_umap2 <- project_umap2[, c("Species", "Cluster", "Barcode", "umap_1","umap_2","Color")]
colnames(project_umap2) <- c("Species", "Cluster", "Barcode", "UMAP.1",  "UMAP.2", "Color")  

project_umap2$Color <- ifelse(project_umap2$Cluster == "11", "#B2D8FF", project_umap2$Color)
project_umap2$Color <- ifelse(project_umap2$Cluster == "12", "#FFF8D0", project_umap2$Color)




target_clusters <- c("6", "4", "2", "1", "11", "12")

set.seed(42)
slingshot_df <- project_umap2 %>%
  filter(Cluster %in% target_clusters) %>%
  group_by(Cluster) %>%
  slice_sample(n = 445) %>%
  ungroup() %>%
  keep_center_cells(1)

plotting_df <- project_umap2

plot_slingshot(
  slingshot_df = slingshot_df,
  lins = list(Lineage1 = c("6", "4", "2", "1", "11", "12")),
  plotting_df = plotting_df,
  title = "Tung + Quanzi - inferred from even-sampled-2",
  test_mode = FALSE
)
