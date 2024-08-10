# Construct master annotation table for passage back to full Seurat object
buildmasterannotation <- function(seurat = group_subset,
                                  df_masterannotation = NULL
                                  ){
df_metadata <- group_subset@meta.data
df_metadata$major_group <- group_to_test
df_metadata <- df_metadata[, c("major_group", "minor_group")]

if(is.null(df_masterannotation)) {
  df_masterannotation <- df_metadata
  } else {
    df_masterannotation <- rbind(df_masterannotation, df_metadata)
    }

return(df_masterannotation)
}
