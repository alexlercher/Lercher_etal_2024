# Save formatting for commonly used plots
# Emma DeGrace
# 2022-01-18

#### Discrete UMAP ####
discrete_umap <- function(
  object,
  group.by,
  cols,
  pt.size = NULL,
  label = FALSE,
  repel = FALSE,
  shuffle = FALSE,
  legend.position = "right",
  title = ""
  ) {
  # Required Libraries
  require(ggplot2)
  require(gtools)
  require(Seurat)

  # Make plot
  p <- DimPlot(
    object = object,
    reduction = "umap",
    group.by = group.by,
    pt.size = pt.size,
    cols = cols,
    label = label,
    repel = repel,
    shuffle = shuffle
  ) +
    coord_fixed(ratio = 0.75) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = legend.position
    )
  return(p)
}

#### Feature UMAP ####
feature_umap <- function(
  object,
  features,
  pt.size = NULL,
  label = FALSE,
  order = FALSE,
  title = "",
  feature.name = "",
  legend.position = "right"
) {
  # Required Libraries
  require(ggplot2)
  require(gtools)
  require(Seurat)

  # Make plot
  p <- FeaturePlot(
    object = object,
    reduction = "umap",
    features = features,
    pt.size = pt.size,
    label = label,
    order = order
  ) +
    coord_fixed(ratio = 0.75) +
    scale_color_viridis(option = "E", name = feature.name) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = legend.position
    )
  return(p)
}


#### Marker Dotplot ####
marker_dotplot <- function(
  object,
  features,
  group.by = group.by,
  title = ""
) {
  # Required Libraries
  require(ggplot2)
  require(gtools)
  require(Seurat)

  # Make plot
  p <- DotPlot(object = object, features = features, group.by = group.by, assay = "RNA") +
    coord_fixed() + coord_flip() +
    scale_size(range = c(2, 6)) +
    scale_color_gradientn(colours = brewer_pal(palette = "Blues", direction = 1)(8)) +
    labs(title = title) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 270),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  return(p)
}

#
#### Fill Barplot ####
fill_barplot <- function(
  data,
  x,
  x.limits = NULL,
  y.limits = NULL,
  fill,
  position = "stack",
  cols,
  x.title = "",
  y.title = "",
  title = "",
  caption = "",
  fill.name = "",
  legend.position = "right"
) {

  # Required Libraries
  require(ggplot2)
  require(gtools)
  require(Seurat)

      # Make plot
    p <- ggplot(
      data = data,
      aes(x = {{ x }}, fill = {{ fill }}) # {{}} look for variable within the data
      ) +
      geom_bar(colour = "black", position = position) +
      scale_fill_discrete(type = cols, name = fill.name) +
      scale_x_discrete(expand = c(0, 0), limits = x.limits) +
      scale_y_continuous(expand = c(0, 0), limits = y.limits) +
      labs(title = title, x = x.title, y = y.title, caption = caption) +
      theme(
        legend.position = legend.position,
        rect = element_blank(),
        axis.line = element_line(),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.minor.y = element_line(color = "grey")
      )

    return(p)
  }

#### Value Barplot ####
value_barplot <- function(
  data,
  x,
  y,
  value = FALSE,
  x.limits = NULL,
  y.limits = NULL,
  fill,
  position = "stack",
  cols,
  x.title = "",
  y.title,
  title = "",
  caption = "",
  fill.name = "",
  legend.position = "right"
) {

  if (value == TRUE) {
    # Make plot
    p <- ggplot(
      data = data,
      aes(x = {{ x }}, y = value) # {{}} look for variable within the data
    ) +
      geom_bar(stat = "identity", aes(fill = {{ fill }}), colour = "black", position = position) +
      scale_fill_discrete(type = cols, name = fill.name) +
      scale_x_discrete(expand = c(0, 0), limits = x.limits) +
      scale_y_continuous(expand = c(0, 0), limits = y.limits) +
      labs(title = title, x = x.title, y = y.title, caption = caption) +
      theme(
        legend.position = legend.position,
        rect = element_blank(),
        axis.line = element_line(),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.minor.y = element_line(color = "grey")
      )
  } else {
    # Make plot
    p <- ggplot(
      data = data,
      aes(x = {{ x }}, y = {{ y }}) # {{}} look for variable within the data
    ) +
      geom_bar(stat = "identity", aes(fill = {{ fill }}), colour = "black", position = position) +
      scale_fill_discrete(type = cols, name = fill.name) +
      scale_x_discrete(expand = c(0, 0), limits = x.limits) +
      scale_y_continuous(expand = c(0, 0), limits = y.limits) +
      labs(title = title, x = x.title, y = y.title, caption = caption) +
      theme(
        legend.position = legend.position,
        rect = element_blank(),
        axis.line = element_line(),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.minor.y = element_line(color = "grey")
      )
  }

  return(p)
}

