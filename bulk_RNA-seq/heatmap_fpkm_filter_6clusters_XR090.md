XR090_heatmap_fpkm_filter_6clusters
================
Alexander Lercher
2025-06-19

This script will do hierarchical clustering (6 clusters total) of DEGs
identified via DESeq2, assign information whether a gene is a known
interferon stimulated gene (ISG) and plot the data as heatmap as well as
save a the list of genes including cluster information as .tsv file.

``` r
#--------------------------------------------------------------------
# LOAD PACKAGES
#--------------------------------------------------------------------
library(readr)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ purrr     1.0.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.3.0
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(dplyr)
library(pheatmap)
```

``` r
#--------------------------------------------------------------------
# DATA IMPORT AND CLEANUP
#--------------------------------------------------------------------
list_of_files <- list.files(path = "input/DESeq2",
                            recursive = TRUE,
                            pattern = ".tsv",
                            full.names = TRUE)

data <- readr::read_tsv(list_of_files, id = NULL)
```

    ## Rows: 211683 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): gene_id, comparison_name
    ## dbl (7): comparison, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
data_fpkm <- read.delim("input/NormSeqData/XR090_FPKM.tsv")

# what are the comparisons made
data_summary <- data %>%
  group_by(comparison) %>%  
  distinct(comparison_name)
```

``` r
#--------------------------------------------------------------------
# FILTER for and identify SIGNIFICANT DEG
#--------------------------------------------------------------------
# define log2FC and adjPval cutoffs
log2FC_cutoff = 1
adjpval_cutoff = 0.05

# subset data by cutoff criteria
significant_DEG <- data %>%
  group_by(comparison) %>%
  filter(abs(log2FoldChange) > log2FC_cutoff & padj < adjpval_cutoff)

# count significant genes per condition
count_DEG <-  significant_DEG %>%
  group_by(comparison, comparison_name) %>%
  count(nrow(comparison)) %>%
  dplyr::rename(DEG_count = n) %>%
  print()
```

    ## # A tibble: 12 × 3
    ## # Groups:   comparison, comparison_name [12]
    ##    comparison comparison_name                     DEG_count
    ##         <dbl> <chr>                                   <int>
    ##  1          1 MA10rec_polyIC_vs_MA10rec_ctrl           1451
    ##  2          2 MA10naive_polyIC_vs_MA10naive_ctrl        456
    ##  3          3 MA10rec_ctrl_vs_MA10naive_ctrl           1777
    ##  4          4 MA10rec_polyIC_vs_MA10naive_polyIC       2673
    ##  5          5 PR8rec_polyIC_vs_PR8rec_ctrl             1741
    ##  6          6 PR8naive_polyIC_vs_PR8naive_ctrl          934
    ##  7          7 PR8rec_ctrl_vs_PR8naive_ctrl             1294
    ##  8          8 PR8rec_polyIC_vs_PR8naive_polyIC         1980
    ##  9          9 MA10naive_ctrl_vs_PR8naive_ctrl           389
    ## 10         10 MA10naive_polyIC_vs_PR8naive_polyIC       383
    ## 11         11 MA10rec_ctrl_vs_PR8rec_ctrl              1623
    ## 12         12 MA10rec_polyIC_vs_PR8rec_polyIC          1247

``` r
# put DEG of individual comparisons into list
list_DEG_comparisons <- list()
for(i in 1:length(count_DEG$comparison)){
  subset_comparison <- subset(significant_DEG, significant_DEG$comparison == i)
  comparison_name <- unique(subset_comparison$comparison_name)
  list_DEG_comparisons[[paste0(i,"_",comparison_name)]] <- subset_comparison
  rm(subset_comparison,comparison_name)
}
```

``` r
#--------------------------------------------------------------------
# CURATE FPKM data
#--------------------------------------------------------------------
# add means for all conditions
data_fpkm$mean_MA10naive_ctrl <-apply(data_fpkm[,2:4],1,mean)
data_fpkm$mean_MA10naive_polyIC <-apply(data_fpkm[,5:7],1,mean)
data_fpkm$mean_MA10rec_ctrl <-apply(data_fpkm[,8:10],1,mean)
data_fpkm$mean_MA10rec_polyIC <-apply(data_fpkm[,11:13],1,mean)
data_fpkm$mean_PR8naive_ctrl <-apply(data_fpkm[,14:16],1,mean)
data_fpkm$mean_PR8naive_polyIC <-apply(data_fpkm[,17:19],1,mean)
data_fpkm$mean_PR8rec_ctrl <-apply(data_fpkm[,20:22],1,mean)
data_fpkm$mean_PR8rec_polyIC <-apply(data_fpkm[,23:25],1,mean)

# rename column gene to gene_id column
data_fpkm <- data_fpkm %>%
  dplyr::rename(gene_id = gene)
```

``` r
#--------------------------------------------------------------------
# LOAD ISG INFO
#--------------------------------------------------------------------
# import ISG list from Schoggins Nature paper, generated by Alex Popa
# problem is that all ISGs are in full caps
ISG_list  <- read.delim("input/ISG_List.csv")

# here i will put the gene names in lower case
# gene names are separated into first letter and rest
# first letter is capitalized, rest is kept lowercase
ISG_list = ISG_list %>%
  mutate(gene_low = tolower(gene),
         first_letter = substr(gene_low,1,1),
         first_letter_cap = toupper(first_letter),
         rest = substr(gene_low,2,20)
  )

# save ISG list as dataframe with the column name "gene"
# add a TRUE column for all ISGs
ISG_list = as.data.frame(paste(ISG_list$first_letter_cap,ISG_list$rest,sep=""))
colnames(ISG_list) = c("gene_id")
ISG_list$ISG = T
```

``` r
#--------------------------------------------------------------------
# DO HEATMAPS only for CONDITIONS of INTEREST
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# SUBSET to MA10 only COMPARISONS
#--------------------------------------------------------------------
# define comparisons of interest
MA10only <- data_summary %>%
  filter(grepl("MA10", comparison_name) & !grepl("PR8", comparison_name)) %>%
  print()
```

    ## # A tibble: 4 × 2
    ## # Groups:   comparison [4]
    ##   comparison comparison_name                   
    ##        <dbl> <chr>                             
    ## 1          1 MA10rec_polyIC_vs_MA10rec_ctrl    
    ## 2          2 MA10naive_polyIC_vs_MA10naive_ctrl
    ## 3          3 MA10rec_ctrl_vs_MA10naive_ctrl    
    ## 4          4 MA10rec_polyIC_vs_MA10naive_polyIC

``` r
# FPKM values for only for MA10 comparisons
MA10only_fpkm <- data_fpkm %>%
  dplyr::select(contains(c("gene","MA10")) & !contains(c("mean")))

# calculate average FPKM level across all samples and define FPKM cutoff
MA10only_fpkm <- MA10only_fpkm %>%
  mutate(average_fpkm = rowSums(MA10only_fpkm[2:length(MA10only_fpkm)]/(length(MA10only_fpkm)-1)))

FPKM_cutoff <- 1

MA10only_fpkm %>%
  ggplot(aes(x=average_fpkm)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
    scale_x_continuous(trans='log10') +
    ggtitle("Average FPKM before cutoff") +
    geom_vline(xintercept = FPKM_cutoff) +
    theme_bw()
```

    ## Warning in scale_x_continuous(trans = "log10"): log-10 transformation
    ## introduced infinite values.

    ## Warning: Removed 168 rows containing non-finite outside the scale range
    ## (`stat_density()`).

![](XR090_heatmap_fpkm_filter_6clusters_files/figure-gfm/Do%20Heatmaps%20only%20for%20conditions%20of%20interest-1.png)<!-- -->

``` r
# filter according to FPKM cutoff
MA10only_fpkm <- MA10only_fpkm %>%
  filter(average_fpkm > FPKM_cutoff)

MA10only_fpkm %>%
  ggplot(aes(x=average_fpkm)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  scale_x_continuous(trans='log10') +
  ggtitle("Average FPKM after cutoff") +
  geom_vline(xintercept = FPKM_cutoff) +
  theme_bw()
```

![](XR090_heatmap_fpkm_filter_6clusters_files/figure-gfm/Do%20Heatmaps%20only%20for%20conditions%20of%20interest-2.png)<!-- -->

``` r
# drop average FPKM column
MA10only_fpkm <- MA10only_fpkm %>%
  dplyr::select(!average_fpkm)

# loop over comparisons of interest
# 1) cluster genes
# 2) plot heatmaps
# 3) add cluster info to gene list

DEG_comparisons_cluster_list <- list()
DEG_comparisons_cluster_DF <- data.frame()

for(i in 1:4){
# select condition of interest
conditionOI <- list_DEG_comparisons[[i]]

# extract DEG from condition of interest
genesOI <- data.frame(gene_id = conditionOI$gene_id)

# get FPKM for DEG of interest
fpkmOI <- left_join(genesOI,MA10only_fpkm)
fpkmOI <- na.omit(fpkmOI)
row.names(fpkmOI) <- fpkmOI$gene_id
fpkmOI$gene_id <- NULL

#--------------------------------------------------------------------
# PLOT HEATMAP with ROW Z SCORES
#--------------------------------------------------------------------
# define formula to calculate row Z score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# calculate row Z score for FPKMs for each row (gene) of interest
fpkmOI_norm <- t(apply(fpkmOI, 1, cal_z_score))

# annotate sample columns
my_sample_col <- data.frame(sample = c(rep("MA10naive_ctrl",3),
                                       rep("MA10naive_polyIC",3),
                                       rep("MA10rec_ctrl",3),
                                       rep("MA10rec_polyIC",3)))
row.names(my_sample_col) <- colnames(fpkmOI)

# calculate dendrogram
hc <- hclust(dist(fpkmOI_norm), method = "complete")
#as.dendrogram(hc) %>%  #no need to print dendrogram
#  plot(horiz = T)

# obtain gene names as per dendrogram
#genes_order <- rev(row.names(fpkmOI_norm)[hc$order]) #no need to print gene names

# add cluster information to heatmap
numberofclusters = 6
my_gene_col <-cutree(hc, k = numberofclusters)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "1",
                                           ifelse(test = my_gene_col == 2, yes = "2",
                                                  ifelse(test = my_gene_col == 3, yes = "3",
                                                         ifelse(test = my_gene_col == 4, yes = "4",
                                                                ifelse(test = my_gene_col == 5, yes = "5", no = "6"
                                                         ))))))

# merge ISG_list heatmap row/gene names to add ISG info
my_gene_col$gene_id <- rownames(my_gene_col)
my_gene_col <- my_gene_col %>%
  left_join(ISG_list) %>%
  replace(is.na(.), FALSE)
my_gene_col$ISG <-as.integer(my_gene_col$ISG)
rownames(my_gene_col) <- my_gene_col$gene_id
my_gene_col$gene_id <- NULL

pheatmap(fpkmOI_norm, annotation_col = my_sample_col, annotation_row = my_gene_col, 
         show_rownames = F, main = paste0(i,":",as.character(MA10only[i,2])),
         #color = colorRampPalette(c("midnightblue", "ivory", "firebrick3"))(50),
         cutree_rows = numberofclusters,
         filename = paste0("output/fpkm_filter/",i,"_",as.character(MA10only[i,2]),"_fpkm_filter_6clusters.pdf"))

# add cluster and ISG information to DEG list
my_gene_col$gene_id <- rownames(my_gene_col) 
conditionOI_cluster <- conditionOI %>%
  left_join(my_gene_col) %>%
  na.omit()
DEG_comparisons_cluster_list[[paste0(i,"_",as.character(MA10only[i,2]))]] <- conditionOI_cluster
DEG_comparisons_cluster_DF <- rbind(DEG_comparisons_cluster_DF,conditionOI_cluster)

write_tsv(conditionOI_cluster, paste0("output/fpkm_filter/",i,"_",as.character(MA10only[i,2]),"_cluster_fpkm_filter_6clusters.tsv"))

}
```

    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`
    ## Joining with `by = join_by(gene_id)`

``` r
#--------------------------------------------------------------------
# CLUSTER GENES SIGNIFICANT IN AT LEAST ONE MA10 CONDITION
#--------------------------------------------------------------------
# get all significant DEG for MA10 comparisons
sig_DEG_MA10 <- significant_DEG %>%
  group_by(comparison_name) %>%
  filter(grepl("MA10",comparison_name) & !grepl("PR8",comparison_name))

# reduce to unique significatn DEG for MA10 comparisons
unique_sig_DEG_MA10 <- data.frame(gene_id = unique(sig_DEG_MA10$gene_id))

# get FPKM for DEG of interest
fpkmOI <- left_join(unique_sig_DEG_MA10,MA10only_fpkm)
```

    ## Joining with `by = join_by(gene_id)`

``` r
fpkmOI <- na.omit(fpkmOI)
row.names(fpkmOI) <- fpkmOI$gene_id
fpkmOI$gene_id <- NULL

# define formula to calculate row Z score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# calculate row Z score for FPKMs for each row (gene) of interest
fpkmOI_norm <- t(apply(fpkmOI, 1, cal_z_score))

# annotate sample columns
my_sample_col <- data.frame(sample = c(rep("MA10naive_ctrl",3),
                                       rep("MA10naive_polyIC",3),
                                       rep("MA10rec_ctrl",3),
                                       rep("MA10rec_polyIC",3)))
row.names(my_sample_col) <- colnames(fpkmOI)

# calculate dendrogram
hc <- hclust(dist(fpkmOI_norm), method = "complete")

# add cluster information to heatmap
numberofclusters = 6
my_gene_col <-cutree(hc, k = numberofclusters)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "1",
                                           ifelse(test = my_gene_col == 2, yes = "2",
                                                  ifelse(test = my_gene_col == 3, yes = "3",
                                                         ifelse(test = my_gene_col == 4, yes = "4",
                                                                ifelse(test = my_gene_col == 5, yes = "5", no = "6"
                                                                ))))))

# merge ISG_list heatmap row/gene names to add ISG info
my_gene_col$gene_id <- rownames(my_gene_col)
my_gene_col <- my_gene_col %>%
  left_join(ISG_list) %>%
  replace(is.na(.), FALSE)
```

    ## Joining with `by = join_by(gene_id)`

``` r
my_gene_col$ISG <-as.integer(my_gene_col$ISG)
rownames(my_gene_col) <- my_gene_col$gene_id
my_gene_col$gene_id <- NULL

pheatmap(fpkmOI_norm, annotation_col = my_sample_col, annotation_row = my_gene_col, 
         show_rownames = F, main = "significant in at least one MA10 comparison",
         #color = colorRampPalette(c("midnightblue", "ivory", "firebrick3"))(50),
         cutree_rows = numberofclusters,
         filename = "output/fpkm_filter/all_MA10_DEG_fpkm_filter_6clusters.pdf")

# add cluster and ISG information to DEG list
my_gene_col$gene_id <- rownames(my_gene_col) 
sig_DEG_MA10_cluster <- sig_DEG_MA10 %>%
  left_join(my_gene_col) %>%
  na.omit()
```

    ## Joining with `by = join_by(gene_id)`

``` r
write_tsv(sig_DEG_MA10_cluster, "output/fpkm_filter/all_MA10_DEG_cluster_fpkm_filter_6clusters.tsv")
```

``` r
# which R packages and versions?
if ("devtools" %in% installed.packages()) devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.0 (2025-04-11)
    ##  os       macOS Sonoma 14.7.3
    ##  system   aarch64, darwin20
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       America/New_York
    ##  date     2025-06-19
    ##  pandoc   3.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
    ##  quarto   1.5.57 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package      * version date (UTC) lib source
    ##  bit            4.6.0   2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64          4.6.0-1 2025-01-16 [1] CRAN (R 4.5.0)
    ##  cachem         1.1.0   2024-05-16 [1] CRAN (R 4.5.0)
    ##  cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
    ##  crayon         1.5.3   2024-06-20 [1] CRAN (R 4.5.0)
    ##  devtools       2.4.5   2022-10-11 [1] CRAN (R 4.5.0)
    ##  digest         0.6.37  2024-08-19 [1] CRAN (R 4.5.0)
    ##  dplyr        * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
    ##  ellipsis       0.3.2   2021-04-29 [1] CRAN (R 4.5.0)
    ##  evaluate       1.0.3   2025-01-10 [1] CRAN (R 4.5.0)
    ##  farver         2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
    ##  forcats      * 1.0.0   2023-01-29 [1] CRAN (R 4.5.0)
    ##  fs             1.6.6   2025-04-12 [1] CRAN (R 4.5.0)
    ##  generics       0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
    ##  ggplot2      * 3.5.2   2025-04-09 [1] CRAN (R 4.5.0)
    ##  glue           1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
    ##  gtable         0.3.6   2024-10-25 [1] CRAN (R 4.5.0)
    ##  hms            1.1.3   2023-03-21 [1] CRAN (R 4.5.0)
    ##  htmltools      0.5.8.1 2024-04-04 [1] CRAN (R 4.5.0)
    ##  htmlwidgets    1.6.4   2023-12-06 [1] CRAN (R 4.5.0)
    ##  httpuv         1.6.16  2025-04-16 [1] CRAN (R 4.5.0)
    ##  knitr          1.50    2025-03-16 [1] CRAN (R 4.5.0)
    ##  labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
    ##  later          1.4.2   2025-04-08 [1] CRAN (R 4.5.0)
    ##  lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
    ##  lubridate    * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr       2.0.3   2022-03-30 [1] CRAN (R 4.5.0)
    ##  memoise        2.0.1   2021-11-26 [1] CRAN (R 4.5.0)
    ##  mime           0.13    2025-03-17 [1] CRAN (R 4.5.0)
    ##  miniUI         0.1.2   2025-04-17 [1] CRAN (R 4.5.0)
    ##  pheatmap     * 1.0.13  2025-06-05 [1] CRAN (R 4.5.0)
    ##  pillar         1.10.2  2025-04-05 [1] CRAN (R 4.5.0)
    ##  pkgbuild       1.4.8   2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload        1.4.0   2024-06-28 [1] CRAN (R 4.5.0)
    ##  profvis        0.4.0   2024-09-20 [1] CRAN (R 4.5.0)
    ##  promises       1.3.3   2025-05-29 [1] CRAN (R 4.5.0)
    ##  purrr        * 1.0.4   2025-02-05 [1] CRAN (R 4.5.0)
    ##  R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
    ##  RColorBrewer   1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
    ##  Rcpp           1.0.14  2025-01-12 [1] CRAN (R 4.5.0)
    ##  readr        * 2.1.5   2024-01-10 [1] CRAN (R 4.5.0)
    ##  remotes        2.5.0   2024-03-17 [1] CRAN (R 4.5.0)
    ##  rlang          1.1.6   2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown      2.29    2024-11-04 [1] CRAN (R 4.5.0)
    ##  rstudioapi     0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
    ##  scales         1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
    ##  sessioninfo    1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
    ##  shiny          1.10.0  2024-12-14 [1] CRAN (R 4.5.0)
    ##  stringi        1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr      * 1.5.1   2023-11-14 [1] CRAN (R 4.5.0)
    ##  tibble       * 3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr        * 1.3.1   2024-01-24 [1] CRAN (R 4.5.0)
    ##  tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange     0.3.0   2024-01-18 [1] CRAN (R 4.5.0)
    ##  tzdb           0.5.0   2025-03-15 [1] CRAN (R 4.5.0)
    ##  urlchecker     1.0.1   2021-11-30 [1] CRAN (R 4.5.0)
    ##  usethis        3.1.0   2024-11-26 [1] CRAN (R 4.5.0)
    ##  utf8           1.2.6   2025-06-08 [1] CRAN (R 4.5.0)
    ##  vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
    ##  vroom          1.6.5   2023-12-05 [1] CRAN (R 4.5.0)
    ##  withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun           0.52    2025-04-02 [1] CRAN (R 4.5.0)
    ##  xtable         1.8-4   2019-04-21 [1] CRAN (R 4.5.0)
    ##  yaml           2.3.10  2024-07-26 [1] CRAN (R 4.5.0)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
