XR090_GO_enrichment_fpkm_filter_6clusters
================
Alexander Lercher
2025-06-19

This script will Perform GO enrichment analyses for clusters of interest
(here cluster 1 and 2) of hierarchically clustered (6 clusters total)
genes identified as DEG via DESeq2.

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
library(clusterProfiler)
```

    ## 
    ## clusterProfiler v4.17.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang,
    ## W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize
    ## multiomics data. Nature Protocols. 2024, 19(11):3292-3320
    ## 
    ## Attaching package: 'clusterProfiler'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
```

    ## Loading required package: GenomicFeatures
    ## Loading required package: BiocGenerics
    ## Loading required package: generics
    ## 
    ## Attaching package: 'generics'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     as.difftime
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     explain
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union
    ## 
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min
    ## 
    ## Loading required package: S4Vectors
    ## Loading required package: stats4
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     rename
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     second, second<-
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand
    ## 
    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     slice
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     %within%
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce
    ## 
    ## Loading required package: GenomeInfoDb
    ## Loading required package: GenomicRanges
    ## Loading required package: AnnotationDbi
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## 
    ## Attaching package: 'AnnotationDbi'
    ## 
    ## The following object is masked from 'package:clusterProfiler':
    ## 
    ##     select
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(org.Mm.eg.db)
```

    ## 

``` r
library(DOSE)  
```

    ## DOSE v4.3.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an
    ## R/Bioconductor package for Disease Ontology Semantic and Enrichment
    ## analysis. Bioinformatics. 2015, 31(4):608-609

``` r
library(enrichplot)
```

    ## enrichplot v1.29.1 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu and Shengqi
    ## Wang. GOSemSim: an R package for measuring semantic similarity among GO
    ## terms and gene products. Bioinformatics. 2010, 26(7):976-978

``` r
library(clusterProfiler)
library(viridis)
```

    ## Loading required package: viridisLite

``` r
#--------------------------------------------------------------------
# DATA IMPORT AND CLEANUP
#--------------------------------------------------------------------
# import significant DEG in at least one MA10 comparison
data_all_MA10 <- read.delim("output/fpkm_filter/all_MA10_DEG_cluster_fpkm_filter_6clusters.tsv")

# group data (genes) by cluster
# remove duplicate genes (some genes significant across mult comparisons)
data_unique <- data_all_MA10 %>%
  arrange(cluster) %>%
  group_by(cluster) %>%
  arrange(log2FoldChange) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  dplyr::rename(SYMBOL = gene_id) %>%
  arrange(cluster)
```

``` r
#--------------------------------------------------------------------
# PRINT TOP 10 MARKER GENES PER CLUSTER
#--------------------------------------------------------------------
# identify top 10 marker genes per cluster
top10_markers_MA10_cluster <- data_unique %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = log2FoldChange) %>%
  arrange(cluster) %>%
  dplyr::select(cluster,SYMBOL) %>%
  print(n=50)
```

    ## # A tibble: 58 × 2
    ## # Groups:   cluster [6]
    ##    cluster SYMBOL       
    ##      <int> <chr>        
    ##  1       1 Vsig10       
    ##  2       1 Cdca5        
    ##  3       1 Sirpb1c      
    ##  4       1 Fads2        
    ##  5       1 Sgsm2        
    ##  6       1 Abcd2        
    ##  7       1 Tmem41a      
    ##  8       1 Sowahc       
    ##  9       1 Ttc38        
    ## 10       1 Cfh          
    ## 11       2 Serpinb2     
    ## 12       2 Iigp1        
    ## 13       2 Mefv         
    ## 14       2 Ifi205       
    ## 15       2 Apol9a       
    ## 16       2 Apol9b       
    ## 17       2 Ms4a4a       
    ## 18       2 Usp18        
    ## 19       2 Trim30b      
    ## 20       2 Cxcl10       
    ## 21       3 Cd93         
    ## 22       3 Slc39a10     
    ## 23       3 Wdr90        
    ## 24       3 Wdhd1        
    ## 25       3 Peli2        
    ## 26       3 Megf9        
    ## 27       3 Htra3        
    ## 28       3 H1f2         
    ## 29       3 Eef2k        
    ## 30       3 Tsc22d3      
    ## 31       4 Selenot      
    ## 32       4 Abcb1a       
    ## 33       4 1600014C10Rik
    ## 34       4 Peli1        
    ## 35       4 Ccnd1        
    ## 36       4 Itpr1        
    ## 37       4 Il15         
    ## 38       4 AW011738     
    ## 39       4 Dcp2         
    ## 40       4 Tnfaip3      
    ## 41       5 Il18         
    ## 42       5 Prr5l        
    ## 43       5 Mtmr7        
    ## 44       5 Cd83         
    ## 45       5 Iqgap2       
    ## 46       5 Fcgr4        
    ## 47       5 Chst15       
    ## 48       5 Cdc42ep5     
    ## 49       6 Nkg7         
    ## 50       6 Lcn2         
    ## # ℹ 8 more rows

``` r
#--------------------------------------------------------------------
# CONVERT GENE SYMBOLS TO ENTREZID
#--------------------------------------------------------------------
# for enrichment analyses, add ENTREZID to SYMBOL
mm <- org.Mm.eg.db
my.symbols <- unique(data_unique$SYMBOL)
gene_list_ENTREZID <- AnnotationDbi::select(mm, 
                                            keys = my.symbols,
                                            columns = c("ENTREZID", "SYMBOL"),
                                            keytype = "SYMBOL")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
colnames(gene_list_ENTREZID) <- c("SYMBOL","ENTREZID")

data_unique <- left_join(data_unique,gene_list_ENTREZID)
```

    ## Joining with `by = join_by(SYMBOL)`

``` r
#--------------------------------------------------------------------
# RUN GO ENRICHMENT for clusters
#--------------------------------------------------------------------
# GO enrichment analyses
# note there might be too few cells for GO enrichment in some clusters
number_of_clusters = list("1","2")
min_GO_pval <- 0.001

GO_plot_list <- list()
for(i in 1:length(number_of_clusters)){
  clusterOI <- number_of_clusters[[i]]
  clusterOI <- as.numeric(clusterOI)
  genes_of_cluster <- data_unique %>%
    filter(cluster == as.character(i)) %>%
    distinct(ENTREZID) %>%
    pull(ENTREZID)
  
  yy <- enrichGO(gene = genes_of_cluster,
                 OrgDb = org.Mm.eg.db,
                 keyType = "ENTREZID",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 readable = TRUE)
  yy <- pairwise_termsim(yy)
  write_tsv(yy@result, paste0("output/fpkm_filter/GO_enrichment_cluster_",clusterOI,"_fpkm_filter_6clusters.tsv"))
  p2 <- emapplot(yy, pie = "count", layout = "nicely") +
    ggtitle(paste0("GO BP Enrichment Map - Cluster ", i)) +
    theme(aspect.ratio = 1,
          plot.title = element_text(size = 10)) +
    scale_color_viridis(option = "plasma", direction = -1)
  GO_plot_list[[paste0("cluster_",clusterOI)]] <- p2
  ggsave(p2, file=paste0("output/fpkm_filter/GO_enrichment_cluster_",clusterOI,"_fpkm_filter_6clusters.pdf"),dpi = 300, units = c("cm"),width = 30, height = 30)
  print(paste0("cluster_",clusterOI," done!"))
  #rm(genes_of_cluster,yy,p2)
}
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## [1] "cluster_1 done!"

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

    ## [1] "cluster_2 done!"

``` r
GO_plot_list
```

    ## $cluster_1

![](XR090_KEGG_GO_enrichment_fpkm_filter_6clusters_files/figure-gfm/Run%20GO%20Enrichment%20for%20clusters-1.png)<!-- -->

    ## 
    ## $cluster_2

![](XR090_KEGG_GO_enrichment_fpkm_filter_6clusters_files/figure-gfm/Run%20GO%20Enrichment%20for%20clusters-2.png)<!-- -->

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
    ##  package                            * version   date (UTC) lib source
    ##  abind                                1.4-8     2024-09-12 [1] CRAN (R 4.5.0)
    ##  AnnotationDbi                      * 1.71.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  ape                                  5.8-1     2024-12-16 [1] CRAN (R 4.5.0)
    ##  aplot                                0.2.5     2025-02-27 [1] CRAN (R 4.5.0)
    ##  Biobase                            * 2.69.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  BiocGenerics                       * 0.55.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  BiocIO                               1.19.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  BiocParallel                         1.43.3    2025-05-29 [1] Bioconductor 3.22 (R 4.5.0)
    ##  Biostrings                           2.77.1    2025-05-18 [1] Bioconductor 3.22 (R 4.5.0)
    ##  bit                                  4.6.0     2025-03-06 [1] CRAN (R 4.5.0)
    ##  bit64                                4.6.0-1   2025-01-16 [1] CRAN (R 4.5.0)
    ##  bitops                               1.0-9     2024-10-03 [1] CRAN (R 4.5.0)
    ##  blob                                 1.2.4     2023-03-17 [1] CRAN (R 4.5.0)
    ##  cachem                               1.1.0     2024-05-16 [1] CRAN (R 4.5.0)
    ##  cli                                  3.6.5     2025-04-23 [1] CRAN (R 4.5.0)
    ##  clusterProfiler                    * 4.17.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  codetools                            0.2-20    2024-03-31 [1] CRAN (R 4.5.0)
    ##  cowplot                              1.1.3     2024-01-22 [1] CRAN (R 4.5.0)
    ##  crayon                               1.5.3     2024-06-20 [1] CRAN (R 4.5.0)
    ##  curl                                 6.3.0     2025-06-06 [1] CRAN (R 4.5.0)
    ##  data.table                           1.17.4    2025-05-26 [1] CRAN (R 4.5.0)
    ##  DBI                                  1.2.3     2024-06-02 [1] CRAN (R 4.5.0)
    ##  DelayedArray                         0.35.1    2025-05-28 [1] Bioconductor 3.22 (R 4.5.0)
    ##  devtools                             2.4.5     2022-10-11 [1] CRAN (R 4.5.0)
    ##  digest                               0.6.37    2024-08-19 [1] CRAN (R 4.5.0)
    ##  DOSE                               * 4.3.0     2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  dplyr                              * 1.1.4     2023-11-17 [1] CRAN (R 4.5.0)
    ##  ellipsis                             0.3.2     2021-04-29 [1] CRAN (R 4.5.0)
    ##  enrichplot                         * 1.29.1    2025-04-28 [1] Bioconductor 3.22 (R 4.5.0)
    ##  evaluate                             1.0.3     2025-01-10 [1] CRAN (R 4.5.0)
    ##  farver                               2.1.2     2024-05-13 [1] CRAN (R 4.5.0)
    ##  fastmap                              1.2.0     2024-05-15 [1] CRAN (R 4.5.0)
    ##  fastmatch                            1.1-6     2024-12-23 [1] CRAN (R 4.5.0)
    ##  fgsea                                1.35.2    2025-06-04 [1] Bioconductor 3.22 (R 4.5.0)
    ##  forcats                            * 1.0.0     2023-01-29 [1] CRAN (R 4.5.0)
    ##  fs                                   1.6.6     2025-04-12 [1] CRAN (R 4.5.0)
    ##  generics                           * 0.1.4     2025-05-09 [1] CRAN (R 4.5.0)
    ##  GenomeInfoDb                       * 1.45.4    2025-05-28 [1] Bioconductor 3.22 (R 4.5.0)
    ##  GenomicAlignments                    1.45.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  GenomicFeatures                    * 1.61.3    2025-06-02 [1] Bioconductor 3.22 (R 4.5.0)
    ##  GenomicRanges                      * 1.61.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  ggfun                                0.1.8     2024-12-03 [1] CRAN (R 4.5.0)
    ##  ggplot2                            * 3.5.2     2025-04-09 [1] CRAN (R 4.5.0)
    ##  ggplotify                            0.1.2     2023-08-09 [1] CRAN (R 4.5.0)
    ##  ggrepel                              0.9.6     2024-09-07 [1] CRAN (R 4.5.0)
    ##  ggtangle                             0.0.6     2024-12-18 [1] CRAN (R 4.5.0)
    ##  ggtree                               3.17.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  glue                                 1.8.0     2024-09-30 [1] CRAN (R 4.5.0)
    ##  GO.db                                3.21.0    2025-06-10 [1] Bioconductor
    ##  GOSemSim                             2.35.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  gridExtra                            2.3       2017-09-09 [1] CRAN (R 4.5.0)
    ##  gridGraphics                         0.5-1     2020-12-13 [1] CRAN (R 4.5.0)
    ##  gson                                 0.1.0     2023-03-07 [1] CRAN (R 4.5.0)
    ##  gtable                               0.3.6     2024-10-25 [1] CRAN (R 4.5.0)
    ##  hms                                  1.1.3     2023-03-21 [1] CRAN (R 4.5.0)
    ##  htmltools                            0.5.8.1   2024-04-04 [1] CRAN (R 4.5.0)
    ##  htmlwidgets                          1.6.4     2023-12-06 [1] CRAN (R 4.5.0)
    ##  httpuv                               1.6.16    2025-04-16 [1] CRAN (R 4.5.0)
    ##  httr                                 1.4.7     2023-08-15 [1] CRAN (R 4.5.0)
    ##  igraph                               2.1.4     2025-01-23 [1] CRAN (R 4.5.0)
    ##  IRanges                            * 2.43.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  jsonlite                             2.0.0     2025-03-27 [1] CRAN (R 4.5.0)
    ##  KEGGREST                             1.49.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  knitr                                1.50      2025-03-16 [1] CRAN (R 4.5.0)
    ##  labeling                             0.4.3     2023-08-29 [1] CRAN (R 4.5.0)
    ##  later                                1.4.2     2025-04-08 [1] CRAN (R 4.5.0)
    ##  lattice                              0.22-7    2025-04-02 [1] CRAN (R 4.5.0)
    ##  lazyeval                             0.2.2     2019-03-15 [1] CRAN (R 4.5.0)
    ##  lifecycle                            1.0.4     2023-11-07 [1] CRAN (R 4.5.0)
    ##  lubridate                          * 1.9.4     2024-12-08 [1] CRAN (R 4.5.0)
    ##  magrittr                             2.0.3     2022-03-30 [1] CRAN (R 4.5.0)
    ##  Matrix                               1.7-3     2025-03-11 [1] CRAN (R 4.5.0)
    ##  MatrixGenerics                       1.21.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  matrixStats                          1.5.0     2025-01-07 [1] CRAN (R 4.5.0)
    ##  memoise                              2.0.1     2021-11-26 [1] CRAN (R 4.5.0)
    ##  mime                                 0.13      2025-03-17 [1] CRAN (R 4.5.0)
    ##  miniUI                               0.1.2     2025-04-17 [1] CRAN (R 4.5.0)
    ##  nlme                                 3.1-168   2025-03-31 [1] CRAN (R 4.5.0)
    ##  org.Mm.eg.db                       * 3.21.0    2025-06-05 [1] Bioconductor
    ##  patchwork                            1.3.0     2024-09-16 [1] CRAN (R 4.5.0)
    ##  pillar                               1.10.2    2025-04-05 [1] CRAN (R 4.5.0)
    ##  pkgbuild                             1.4.8     2025-05-26 [1] CRAN (R 4.5.0)
    ##  pkgconfig                            2.0.3     2019-09-22 [1] CRAN (R 4.5.0)
    ##  pkgload                              1.4.0     2024-06-28 [1] CRAN (R 4.5.0)
    ##  plyr                                 1.8.9     2023-10-02 [1] CRAN (R 4.5.0)
    ##  png                                  0.1-8     2022-11-29 [1] CRAN (R 4.5.0)
    ##  profvis                              0.4.0     2024-09-20 [1] CRAN (R 4.5.0)
    ##  promises                             1.3.3     2025-05-29 [1] CRAN (R 4.5.0)
    ##  purrr                              * 1.0.4     2025-02-05 [1] CRAN (R 4.5.0)
    ##  qvalue                               2.41.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  R.methodsS3                          1.8.2     2022-06-13 [1] CRAN (R 4.5.0)
    ##  R.oo                                 1.27.1    2025-05-02 [1] CRAN (R 4.5.0)
    ##  R.utils                              2.13.0    2025-02-24 [1] CRAN (R 4.5.0)
    ##  R6                                   2.6.1     2025-02-15 [1] CRAN (R 4.5.0)
    ##  ragg                                 1.4.0     2025-04-10 [1] CRAN (R 4.5.0)
    ##  RColorBrewer                         1.1-3     2022-04-03 [1] CRAN (R 4.5.0)
    ##  Rcpp                                 1.0.14    2025-01-12 [1] CRAN (R 4.5.0)
    ##  RCurl                                1.98-1.17 2025-03-22 [1] CRAN (R 4.5.0)
    ##  readr                              * 2.1.5     2024-01-10 [1] CRAN (R 4.5.0)
    ##  remotes                              2.5.0     2024-03-17 [1] CRAN (R 4.5.0)
    ##  reshape2                             1.4.4     2020-04-09 [1] CRAN (R 4.5.0)
    ##  restfulr                             0.0.15    2022-06-16 [1] CRAN (R 4.5.0)
    ##  rjson                                0.2.23    2024-09-16 [1] CRAN (R 4.5.0)
    ##  rlang                                1.1.6     2025-04-11 [1] CRAN (R 4.5.0)
    ##  rmarkdown                            2.29      2024-11-04 [1] CRAN (R 4.5.0)
    ##  Rsamtools                            2.25.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  RSQLite                              2.4.1     2025-06-08 [1] CRAN (R 4.5.0)
    ##  rstudioapi                           0.17.1    2024-10-22 [1] CRAN (R 4.5.0)
    ##  rtracklayer                          1.69.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  S4Arrays                             1.9.1     2025-05-28 [1] Bioconductor 3.22 (R 4.5.0)
    ##  S4Vectors                          * 0.47.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  scales                               1.4.0     2025-04-24 [1] CRAN (R 4.5.0)
    ##  sessioninfo                          1.2.3     2025-02-05 [1] CRAN (R 4.5.0)
    ##  shiny                                1.10.0    2024-12-14 [1] CRAN (R 4.5.0)
    ##  SparseArray                          1.9.0     2025-05-28 [1] Bioconductor 3.22 (R 4.5.0)
    ##  stringi                              1.8.7     2025-03-27 [1] CRAN (R 4.5.0)
    ##  stringr                            * 1.5.1     2023-11-14 [1] CRAN (R 4.5.0)
    ##  SummarizedExperiment                 1.39.0    2025-05-28 [1] Bioconductor 3.22 (R 4.5.0)
    ##  systemfonts                          1.2.3     2025-04-30 [1] CRAN (R 4.5.0)
    ##  textshaping                          1.0.1     2025-05-01 [1] CRAN (R 4.5.0)
    ##  tibble                             * 3.3.0     2025-06-08 [1] CRAN (R 4.5.0)
    ##  tidyr                              * 1.3.1     2024-01-24 [1] CRAN (R 4.5.0)
    ##  tidyselect                           1.2.1     2024-03-11 [1] CRAN (R 4.5.0)
    ##  tidytree                             0.4.6     2023-12-12 [1] CRAN (R 4.5.0)
    ##  tidyverse                          * 2.0.0     2023-02-22 [1] CRAN (R 4.5.0)
    ##  timechange                           0.3.0     2024-01-18 [1] CRAN (R 4.5.0)
    ##  treeio                               1.33.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  TxDb.Mmusculus.UCSC.mm10.knownGene * 3.10.0    2025-06-06 [1] Bioconductor
    ##  tzdb                                 0.5.0     2025-03-15 [1] CRAN (R 4.5.0)
    ##  UCSC.utils                           1.5.0     2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  urlchecker                           1.0.1     2021-11-30 [1] CRAN (R 4.5.0)
    ##  usethis                              3.1.0     2024-11-26 [1] CRAN (R 4.5.0)
    ##  utf8                                 1.2.6     2025-06-08 [1] CRAN (R 4.5.0)
    ##  vctrs                                0.6.5     2023-12-01 [1] CRAN (R 4.5.0)
    ##  viridis                            * 0.6.5     2024-01-29 [1] CRAN (R 4.5.0)
    ##  viridisLite                        * 0.4.2     2023-05-02 [1] CRAN (R 4.5.0)
    ##  vroom                                1.6.5     2023-12-05 [1] CRAN (R 4.5.0)
    ##  withr                                3.0.2     2024-10-28 [1] CRAN (R 4.5.0)
    ##  xfun                                 0.52      2025-04-02 [1] CRAN (R 4.5.0)
    ##  XML                                  3.99-0.18 2025-01-01 [1] CRAN (R 4.5.0)
    ##  xtable                               1.8-4     2019-04-21 [1] CRAN (R 4.5.0)
    ##  XVector                              0.49.0    2025-04-27 [1] Bioconductor 3.22 (R 4.5.0)
    ##  yaml                                 2.3.10    2024-07-26 [1] CRAN (R 4.5.0)
    ##  yulab.utils                          0.2.0     2025-01-29 [1] CRAN (R 4.5.0)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
