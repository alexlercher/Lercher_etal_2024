# Functions for subset group cluster annotations
# Brad Rosenberg

# Function for subseting on a major group and re-running PCA
subset_major_group <- function(seurat,
                               group_to_test){
  require(Seurat)
  # Subset to major cell group, re-run PCA, and re-run clustering
  seurat <- subset(seurat,
                         ident = major_cell_groups[[group_to_test]]
  )
  # Remove intermediate cluster columns from metadata
  seurat@meta.data[,
                         grep(
                           colnames(seurat@meta.data),
                           pattern = "^integrated_snn_res*")
  ] <- NULL

  # Run PCA with viral genes and tcr/bcr variabel genes removed from variable features
  DefaultAssay(seurat) <- "integrated"
  features <- seurat@assays$integrated@var.features
  seurat <- RunPCA(
    seurat, npcs = 100, assay = "integrated",
    features = features[!features %in% c(iav_genes, vdj_genes)],
    verbose = FALSE
  )

  # For downstream visualization
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat, assay = "RNA")
  seurat <- ScaleData(seurat, assay = "RNA")
  return(seurat)
}

determine_cluster_res <- function(seurat, ndims, maxres){
  require(Seurat)
  DefaultAssay(seurat) <- "integrated"
  seurat <- FindNeighbors(seurat, dims = 1:ndims, verbose = FALSE) %>%
    RunUMAP(reduction = "pca",
            dims = 1:ndims,
            return.model = TRUE,
            verbose = FALSE) %>%
    FindClusters(algorithm = 3,
                 resolution = seq(from = 0.2,
                                  to = maxres,
                                  by = 0.2),
                 verbose = FALSE
                 )
  return(seurat)
}

# Incoporate minor group annotations and constituent subcluster info
# to Seurat object to facilitate downstream DGE testing and annotation
add_minor_group_metadata <- function (seurat, minor_groups){
  require(reshape2)
  require(Seurat)
  df_minor_groups <- melt(minor_groups)
  colnames(df_minor_groups) <- c("cluster", "minor_group")
  df_minor_groups$cluster <- as.character(df_minor_groups$cluster)
  df_metadata <- seurat@meta.data
  df_metadata$cell_barcode <- rownames(df_metadata)
  df_metadata$subset_clusters <- as.character(df_metadata$subset_clusters)
  df_metadata_toadd <- inner_join(
    x = df_metadata,
    y = df_minor_groups,
    by = c("subset_clusters" = "cluster")
  )
  rownames(df_metadata_toadd) <- df_metadata_toadd$cell_barcode
  seurat <- AddMetaData(seurat, df_metadata_toadd[,c("subset_clusters", "minor_group")])
  seurat$subset_clusters <- factor(seurat$subset_clusters)
  seurat$minor_group <- factor(seurat$minor_group)
  return(seurat)
}
