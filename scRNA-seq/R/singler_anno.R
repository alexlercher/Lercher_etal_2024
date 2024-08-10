# Function for adding SingleR annotation

singler_anno <- function(SeuratObj,
                         type = "bulk",
                         reference,
                         labels,
                         name){
  # Libraries
  require(SingleCellExperiment)
  require(SingleR)

  if(type == "sce"){
    # Convert to single-cell object SingleR will recognize
    reference <- DietSeurat(reference)
    reference <- as.SingleCellExperiment(reference)
  }

  # SingleR
  predictions <- SingleR(
    test = GetAssayData(object = SeuratObj, slot = "data"),
    ref = reference,
    labels = labels,
    aggr.ref = TRUE
  )
  SeuratObj <- AddMetaData(
    object = SeuratObj,
    metadata = predictions$labels,
    col.name = name)
  return(SeuratObj)
}

