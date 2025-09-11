# Jr-Fong Dang
library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)
oriPar = par(no.readonly=T)
setwd("~/local/Diskarray/Rdata")

# Read seurat data start
Tung_xylem.data = Read10X(data.dir = '/Single_species_analysis/cellranger_reanalysis_TenX_Ptr/outs/filtered_feature_bc_matrix')
Quanzi_xylem_2.data = Read10X(data.dir = '/Single_species_analysis/cellranger_reanalysis_TenX_NucXyBio2/outs/filtered_feature_bc_matrix')
Quanzi_xylem_3.data = Read10X(data.dir = '/Single_species_analysis/cellranger_reanalysis_TenX_NucXyBio3/outs/filtered_feature_bc_matrix')
colnames(Tung_xylem.data) <- paste0( "Tung_", colnames(Tung_xylem.data))
colnames(Quanzi_xylem_2.data) <- paste0( "Quanzi2_", colnames(Quanzi_xylem_2.data))
colnames(Quanzi_xylem_3.data) <- paste0( "Quanzi3_", colnames(Quanzi_xylem_3.data))

#Setup the Seurat Object
Tung_xylem = CreateSeuratObject(counts = Tung_xylem.data, project = 'Tung_Xylem', min.cells = 3, min.features = 200)
Quanzi_xylem_2 = CreateSeuratObject(counts = Quanzi_xylem_2.data, project = 'Quanzi_Xylem_2', min.cells = 3, min.features = 200)
Quanzi_xylem_3 = CreateSeuratObject(counts = Quanzi_xylem_3.data, project = 'Quanzi_Xylem_3', min.cells = 3, min.features = 200)

#Normalize the data
Tung_xylem = NormalizeData(Tung_xylem, normalization.method = 'LogNormalize', scale.factor = 10000)
Quanzi_xylem_2 = NormalizeData(Quanzi_xylem_2, normalization.method = 'LogNormalize', scale.factor = 10000)
Quanzi_xylem_3 = NormalizeData(Quanzi_xylem_3, normalization.method = 'LogNormalize', scale.factor = 10000)

#Identify the highly variable features (feature selection)
Tung_xylem = FindVariableFeatures(Tung_xylem, selection.method = 'vst', nfeatures = 2000)
Quanzi_xylem_2 = FindVariableFeatures(Quanzi_xylem_2, selection.method = 'vst', nfeatures = 2000)
Quanzi_xylem_3 = FindVariableFeatures(Quanzi_xylem_3, selection.method = 'vst', nfeatures = 2000)

#Find integration anchors and integrate data
integration_anchors = FindIntegrationAnchors(object.list = list(Tung_xylem,Quanzi_xylem_3,Quanzi_xylem_2),
                                             anchor.features = 2000,
                                             scale = TRUE,
                                             reduction = 'cca',
                                             l2.norm = TRUE,
                                             k.filter = 100,
                                             k.anchor = 5)

Combined_object = IntegrateData(anchorset = integration_anchors)

#Run the standard workflow for visualization and clustering
Combined_object = ScaleData(Combined_object)
Combined_object = RunPCA(Combined_object, npcs = 30)
Combined_object = RunUMAP(Combined_object, reduction = 'pca', dims = 1:30)
Combined_object = FindNeighbors(Combined_object, reduction = 'pca', dims = 1:30, k.param = 3)
Combined_object  = FindClusters(Combined_object, resolution = 0.5)

projection_UMAP = Combined_object@reductions$umap@cell.embeddings %>% as.data.frame
projection_PC30 = Combined_object@reductions$pca@cell.embeddings %>% as.data.frame