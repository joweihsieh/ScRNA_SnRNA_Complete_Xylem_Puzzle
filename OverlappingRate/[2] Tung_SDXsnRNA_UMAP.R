setwd("/home/woodydrylab/FileShare/temporary_data_deposit/SC_SnRNA_DATA/results/overlapping")

suppressPackageStartupMessages({
  library(Seurat)
})

runUMAP_getEmbeds <- function(rdsFilePath,
                              dims = 1:30,
                              seed.use = 42,
                              min.dist = 0.3,
                              n.neighbors = 30,
                              metric = "cosine") {
  obj <- readRDS(rdsFilePath)

  obj <- RunUMAP(
    object = obj,
    reduction = "pca",
    dims = dims,
    seed.use = seed.use,
    min.dist = min.dist,
    n.neighbors = n.neighbors,
    umap.method = "uwot",
    metric = metric
  )
  emb <- obj@reductions$umap@cell.embeddings
  return(emb)
}

make_pair_plot_df <- function(umap_emb, sp1, sp2) {
  rn <- rownames(umap_emb)
  idx1 <- grepl(sp1, rn, fixed = TRUE)
  idx2 <- grepl(sp2, rn, fixed = TRUE)

  df1 <- data.frame(UMAP1 = umap_emb[idx1, 1],
                    UMAP2 = umap_emb[idx1, 2],
                    Group = "A")  
  df2 <- data.frame(UMAP1 = umap_emb[idx2, 1],
                    UMAP2 = umap_emb[idx2, 2],
                    Group = "B")  

  rbind(df1, df2)
}

plot_umap_pair <- function(plot_df, prefix,
                           point_cex = 0.25,
                           img_w = 2800, img_h = 2000, img_res = 400) {

  col_black <- "black"
  col_gold  <- "#DAA520"

  # 自動軸範圍，並留一點邊界
  xr <- range(plot_df$UMAP1, na.rm = TRUE)
  yr <- range(plot_df$UMAP2, na.rm = TRUE)
  pad_x <- diff(xr) * 0.05
  pad_y <- diff(yr) * 0.05
  xlim <- c(xr[1] - pad_x, xr[2] + pad_x)
  ylim <- c(yr[1] - pad_y, yr[2] + pad_y)

  png(filename = paste0("UMAP_overlap_", prefix, ".png"),
      width = img_w, height = img_h, res = img_res)

  par(mar = c(4.5, 4.5, 3.5, 1.5))

  with(plot_df[plot_df$Group == "B", ],
       plot(UMAP1, UMAP2,
            col = col_gold, pch = 16, cex = point_cex,
            xlab = "UMAP 1", ylab = "UMAP 2",
            xlim = xlim, ylim = ylim, axes = TRUE))

  with(plot_df[plot_df$Group == "A", ],
       points(UMAP1, UMAP2,
              col = col_black, pch = 16, cex = point_cex))

  title(main = paste0("UMAP overlap  ", prefix))


  dev.off()
}

run_and_plot_pair <- function(rdsFilePath, prefix, sp1, sp2) {
  emb <- runUMAP_getEmbeds(rdsFilePath)
  plot_df <- make_pair_plot_df(emb, sp1, sp2)
  plot_umap_pair(plot_df, prefix)
}

run_and_plot_pair(
  rdsFilePath = "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis/all_data_rds/integration_scRNATung_snRNABio2_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds",
  prefix = "scRNATung_snRNABio2",
  sp1 = "TenX_Ptr_",
  sp2 = "TenX_NucXyBio2"
)

run_and_plot_pair(
  rdsFilePath = "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis/all_data_rds/integration_scRNATung_snRNABio3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds",
  prefix = "scRNATung_snRNABio3",
  sp1 = "TenX_Ptr_",
  sp2 = "TenX_NucXyBio3"
)

run_and_plot_pair(
  rdsFilePath = "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Multi_species_analysis/all_data_rds/integration_snRNA_Bio2_Bio3_an_2000_ft_200_kan_5_seed_42_md_0.3_nn_30_base.rds",
  prefix = "snRNA_Bio2_Bio3",
  sp1 = "TenX_NucXyBio2",
  sp2 = "TenX_NucXyBio3"
)
