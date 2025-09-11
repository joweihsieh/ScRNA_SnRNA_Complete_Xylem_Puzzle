# Jr-Fong Dang
target_list = read.csv("/home/rstudio/local/Diskarray/shinyGO_analysis.csv")

Combined_object@meta.data$kmeans <- NA

idx <- match(target_list$Barcode, rownames(Combined_object@meta.data))
Combined_object@meta.data$kmeans[idx] <- target_list$Cluster

Combined_object@meta.data[["C1_to_others"]] <- ifelse(
  is.na(Combined_object@meta.data[["kmeans"]]), 
  NA, 
  ifelse(Combined_object@meta.data[["kmeans"]] == 1, "C1", "others")
)

DefaultAssay(Combined_object) <- "RNA"
Combined_object <- JoinLayers(Combined_object)

result <- FindMarkers(
  Combined_object,
  group.by = "C1_to_others",
  ident.1 = "C1",
  ident.2 = "others",
  slot = "counts"
)

filtered_result_p <- result %>% filter(p_val_adj < 0.05)
filtered_result_p$trend <- ifelse(filtered_result_p$avg_log2FC > 0, "Up", "Down")
