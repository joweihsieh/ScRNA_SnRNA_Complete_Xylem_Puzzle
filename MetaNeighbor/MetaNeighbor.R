#/usr/bin/R
#/home/f06b22037/SSD2/utility/miniconda3/envs/snakemake/bin/R

#conda install conda-forge::r-seurat

setwd("/home/woodydrylab/FileShare/temporary_data_deposit/SC_SnRNA_DATA/results/MetaNeighbor/")
suppressPackageStartupMessages({
  library(Seurat)
  library(MetaNeighbor)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(pheatmap)
  library(dplyr)
  library(Matrix)
  library(readr)
  library(tibble)
})

set.seed(42)
outdir <- getwd()
deg_csv <- "/home/woodydrylab/FileShare/temporary_data_deposit/SC_SnRNA_DATA/results/MetaNeighbor/20250520_Input for shinyGO.csv"
logfc_cut <- 1
fdr_cut   <- 0.05
topN_up   <- 150
topN_dn   <- 150
cap_per_cluster <- 150



# ====== 1) ======
bio2_path <- "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/cellranger_reanalysis_TenX_NucXyBio2/outs/filtered_feature_bc_matrix"
bio3_path <- "/home/f06b22037/SSD2/JW/1136project_SingleCell/results/Single_species_analysis/cellranger_reanalysis_TenX_NucXyBio3/outs/filtered_feature_bc_matrix"

Bio2 <- Read10X(bio2_path)
Bio3 <- Read10X(bio3_path)

# ====== 2)  ======
obj2 <- CreateSeuratObject(counts = Bio2, project = 'Quanzi2', min.cells = 3, min.features = 200)
obj3 <- CreateSeuratObject(counts = Bio3, project = 'Quanzi3', min.cells = 3, min.features = 200)

obj2$replicate <- "Quanzi2"
obj3$replicate <- "Quanzi3"


cluster_df0 <- read.csv("hero_shot.csv", header = TRUE, stringsAsFactors = FALSE)
cluster_dfBio2 <- cluster_df0[cluster_df0$Species %in% c("Quanzi2"),]
cluster_dfBio3 <- cluster_df0[cluster_df0$Species %in% c("Quanzi3"),]


# ====== 3) NormalizeData + FindVariableFeatures ======
obj2 <- NormalizeData(obj2, normalization.method = 'LogNormalize', scale.factor = 10000)
obj3 <- NormalizeData(obj3, normalization.method = 'LogNormalize', scale.factor = 10000)

obj2 <- FindVariableFeatures(obj2, selection.method = 'vst', nfeatures = 2000)
obj3 <- FindVariableFeatures(obj3, selection.method = 'vst', nfeatures = 2000)

# ====== 4) cluster label ======

obj2$cluster_late <- cluster_dfBio2$Cluster[match(colnames(obj2), cluster_dfBio2$Barcode)]
obj3$cluster_late <- cluster_dfBio3$Cluster[match(colnames(obj3), cluster_dfBio3$Barcode)]


# ====== 5) Matrix for MetaNeighbor ======

obj2 <- subset(obj2, cells = colnames(obj2)[!is.na(obj2$cluster_late)])
obj3 <- subset(obj3, cells = colnames(obj3)[!is.na(obj3$cluster_late)])

## 
stopifnot(exists("obj2"), exists("obj3"))
stopifnot("cluster_late" %in% colnames(obj2@meta.data),
          "cluster_late" %in% colnames(obj3@meta.data))

##Seurat4 
#expr2 <- GetAssayData(obj2, assay = "RNA", slot = "data")
#expr3 <- GetAssayData(obj3, assay = "RNA", slot = "data")

##Seurat5
expr2 <- GetAssayData(obj2, assay = "RNA", layer = "data")
expr3 <- GetAssayData(obj3, assay = "RNA", layer = "data")


common <- intersect(rownames(expr2), rownames(expr3))

## 
detect_ok <- function(m, genes, p = 0.05) genes[rowMeans(m[genes, , drop = FALSE] > 0) >= p]
keep <- intersect(detect_ok(expr2, common, 0.05),
                  detect_ok(expr3, common, 0.05))

expr2 <- expr2[keep, , drop = FALSE]
expr3 <- expr3[keep, , drop = FALSE]

## 
colnames(expr2) <- paste0("Quanzi2_", colnames(expr2))
colnames(expr3) <- paste0("Quanzi3_", colnames(expr3))
expr_all <- cbind(expr2, expr3)

## 
lab2 <- as.character(obj2$cluster_late[sub("^Quanzi2_", "", colnames(expr2))])
lab3 <- as.character(obj3$cluster_late[sub("^Quanzi3_", "", colnames(expr3))])
cell_type <- c(lab2, lab3)
study_id  <- c(rep("Quanzi2", length(lab2)), rep("Quanzi3", length(lab3)))
names(cell_type) <- colnames(expr_all)
names(study_id)  <- colnames(expr_all)

## 
se_all <- SummarizedExperiment(
  assays  = SimpleList(expr = expr_all), 
  colData = DataFrame(study_id = study_id,
                      cell_type = factor(cell_type),
                      row.names = colnames(expr_all))
)

## 
make_even_se <- function(se, cap = 150) {
  tab <- table(se$study_id, se$cell_type)
  keep <- unlist(lapply(rownames(tab), function(study) {
    idx   <- which(se$study_id == study)
    labs  <- se$cell_type[idx]
    cells <- colnames(se)[idx]
    cl_names <- colnames(tab)[tab[study, ] > 0]
    unlist(lapply(cl_names, function(cl) {
      v <- cells[labs == cl]
      k <- min(length(v), cap)
      if (k > 0) sample(v, k) else character(0)
    }), use.names = FALSE)
  }), use.names = FALSE)
  se[, keep, drop = FALSE]
}
se_even <- make_even_se(se_all, cap = cap_per_cluster)

## HVG for MetaNeighbor
hvg_union <- union(VariableFeatures(obj2), VariableFeatures(obj3))
var_genes_even <- intersect(hvg_union, rownames(se_even))

## 

## 
run_and_save <- function(se, var_genes, tag, outdir = ".", cutoff = 0.9) {
  vg <- intersect(var_genes, rownames(se))
  if (length(vg) < 10) stop(paste0(tag, " Less than 10"))
  au <- MetaNeighborUS(var_genes = vg,
                       dat       = se,
                       study_id  = se$study_id,
                       cell_type = se$cell_type)
  write.csv(au, file.path(outdir, paste0("AUROC_", tag, "_20250910.csv")), quote = FALSE)

  invisible(au)
}

au_even_hvg  <- run_and_save(se_even, var_genes_even,       "Even150_HVG",    outdir)
#se <- se_even
#var_genes <- var_genes_even
#tag<-"Even150_HVG"