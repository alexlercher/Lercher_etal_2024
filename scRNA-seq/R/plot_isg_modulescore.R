# Emma DeGrace
# 2023-05-25
# Script for ISG Seurat AddModuleScore

# Always start with session
renv::restore()

# Load required libraries
library(cowplot)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggsci)
library(ggnewscale)
library(introdataviz)
library(pals)
library(RColorBrewer)
library(Seurat)
library(SeuratObject)

# Set file paths
path_fig <- here::here("analysis/figures/")
path_rawdata <- here::here("analysis/data/raw_data/")
path_deriveddata <- here::here("analysis/data/derived_data/")

# Source in R functions
source(here::here("R/ggplot_shortcuts.R"))
source(here::here("R/load_palettes.R"))

# Load annotated Seurat object
balf_obj <- readRDS(file = paste0(path_deriveddata, "sarscov2iavbalf_scRNASeq_AL.RDS"))

# Read gene list to be included in the module
isg_geneset <- read.table(file = paste0(path_deriveddata, "ISG_symbol_ENTREZID.tsv"), header = TRUE)
isg_caption <- isg_geneset[,"SYMBOL"] %>% unique() %in% rownames(balf_obj) %>%
  table() %>% .["TRUE"] %>% paste(
  "NB: In the ISG gene list of", length(unique(isg_geneset[,"SYMBOL"])), "genes,",
  . , "of the genes are included in the GEX matrix input."
  )
isg_geneset <- isg_geneset[isg_geneset[,"SYMBOL"] %in% rownames(balf_obj),]

# Normalize RNA counts
balf_obj <- NormalizeData(
  object = balf_obj,
  assay = "RNA",
  verbose = FALSE
)

# Calculate per-cell ISG module score
balf_obj <- AddModuleScore(
  object =  balf_obj,
  features = list(isg_geneset[,"SYMBOL"]),
  name = "isg_geneset_score",
  assay = "RNA"
)

# Reformat names for plotting
balf_obj$minor_pretty <- dplyr::recode_factor(
  balf_obj$minor_group,
  CD4 = "CD4 T cells",
  CD8 = "CD8 T cells",
  gdT = "gdT cells",
  NKT = "NKT cells",
  B = "B cells",
  macrophage_cluster1 = "macrophage\ncluster1",
  macrophage_cluster2 = "macrophage\ncluster2"
)
balf_obj$group_pretty <- dplyr::recode_factor(
  balf_obj$group,
  PR8_only = "naive/PR8",
  MA10andPR8 = "SARS2/PR8"
)

# Plot module score as heatmap
df <- balf_obj@meta.data[,c("minor_pretty", "group_pretty","isg_geneset_score1")]
df$split <- paste(df$minor_pretty, df$group_pretty, sep = "_")
df <- aggregate(
  x = df$isg_geneset_score1, # Value to summarize
  by = list(df$minor_pretty, df$group_pretty, df$split),
  data = df,
  FUN = median,
  simplify = TRUE
)
p <- ggplot(
  data = df,
  mapping = aes(x = fct_reorder(Group.1, x, .desc = TRUE), y = fct_rev(Group.2))
) +
  geom_tile(aes(fill = x), color = "black") +
  coord_equal() +
  scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu")) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, size = 12),
    axis.text.y = element_text(size = 12, vjust = 0),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid = element_blank(),
    legend.key = element_rect(color = "black")
  ) +
  labs(
    title = "ISG Gene Set Expression Scoring",
    x = "",
    y = "",
    fill = "Median\nISG\ngene set\nscore",
    caption = isg_caption
  ) +
  new_scale_fill() +
  geom_tile(aes(y = 2.7, fill = Group.1), color = "black", height = 0.25) +
  scale_fill_manual(values = pal_minor, guide = "none")
ggsave(
  filename = paste0(path_fig, "ISGModuleScore_Heatmap.pdf"),
  p, width = 12, height = 8, limitsize = FALSE
)
p

# Plot module score as split violin
p <- ggplot(
  data = balf_obj@meta.data,
  aes(
    x = fct_reorder(minor_pretty, isg_geneset_score1, .desc = TRUE),
    y = isg_geneset_score1,
    fill = group_pretty
  )) +
  geom_point(
    aes(fill = group_pretty, color = group_pretty),
    size = 0.5, alpha = .3, show.legend = FALSE,
    position = position_jitterdodge(jitter.width=0.75)
  ) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE, scale = "width") +
  scale_color_manual(values = pal_group) +
  scale_fill_manual(values = pal_group) +
  theme(
    rect = element_blank(),
    axis.line = element_line(),
    panel.grid.major.y = element_line(color = "grey"),
    panel.grid.minor.y = element_line(color = "grey"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Distribution of ISG Gene Set Expression Scores among Cell Types",
    x = "",
    y = "ISG gene set score per cell",
    fill = "",
    caption = isg_caption
  )
ggsave(
  filename = paste0(path_fig, "ISGModuleScore_VlnH.pdf"),
  p, width = 10, height = 6, limitsize = FALSE
)
p
p <- ggplot(
  data = balf_obj@meta.data,
  aes(
    x = fct_reorder(minor_pretty, isg_geneset_score1, .desc = FALSE),
    y = isg_geneset_score1,
    fill = group_pretty
  )) +
  geom_point(
    aes(fill = group_pretty, color = group_pretty),
    size = 0.5, alpha = .3, show.legend = FALSE,
    position = position_jitterdodge(jitter.width=0.75)
  ) +
  introdataviz::geom_split_violin(alpha = .5, trim = FALSE, scale = "width") +
  scale_color_manual(values = pal_group) +
  scale_fill_manual(values = pal_group) +
  coord_flip() +
  theme(
    rect = element_blank(),
    axis.line = element_line(),
    panel.grid.major.x = element_line(color = "grey"),
    panel.grid.minor.x = element_line(color = "grey"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Distribution of ISG Gene Set Expression Scores",
    x = "",
    y = "ISG gene set score per cell",
    fill = "",
    caption = isg_caption
  )
ggsave(
  filename = paste0(path_fig, "ISGModuleScore_VlnV.pdf"),
  p, width = 5, height = 8, limitsize = FALSE
)
p
