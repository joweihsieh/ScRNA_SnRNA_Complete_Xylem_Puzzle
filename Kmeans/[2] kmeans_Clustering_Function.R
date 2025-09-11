#Jr-Fong Dang
projection_PC30 = Combined_object@reductions$pca@cell.embeddings %>% as.data.frame

Cluster<-kmeans(
  projection_PC30[, 1:30],
  centers = 19, nstart = 50
)$cluster
kmeans<-cbind(projection_PC30[, 1:30],Cluster)