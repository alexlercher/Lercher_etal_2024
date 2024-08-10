#### Script to generate mouse gene sets ####
## Emma DeGrace
## 2021-12-17

library(biomaRt)
library(magrittr)
library(dplyr)
library(Seurat)

path_supportingdata <- here::here("analysis/data/raw_data/supporting_data/") # STORE GENE SET HERE

# Stress gene set from O’Flanagan et al. https://doi.org/10.1186/s13059-019-1830-0
oflanagan_list <- read.csv(file = paste0(path_supportingdata, "Campbell_StressCoreGeneSet.csv"))
oflanagan_stress_geneset <- list(
  as.character(
    oflanagan_list$gene_symbol[!is.na(oflanagan_list$gene_symbol)]
  )
)

# Basic function to convert human to mouse gene names
# From: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertHumanGeneList <- function(x) {
  human <- useMart(
    biomart="ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host="uswest.ensembl.org",
    ensemblRedirect = FALSE)
  mouse <- useMart(
    biomart="ENSEMBL_MART_ENSEMBL",
    dataset = "mmusculus_gene_ensembl",
    host="uswest.ensembl.org",
    ensemblRedirect = FALSE)

  genesV2 <- getLDS(attributes = c("hgnc_symbol"),
                    filters = "hgnc_symbol",
                    values = x,
                    mart = human,
                    attributesL = c("mgi_symbol"),
                    martL = mouse,
                    uniqueRows = T
  )

  humanx <- unique(genesV2[, 2])
  return(humanx)
}

# Convert human gene symbols to mouse gene symbols
mm_sgenes <- convertHumanGeneList(cc.genes.updated.2019$s.genes) %>% make.names()
mm_g2mgenes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes) %>% make.names()
oflanagan_stress_geneset <- lapply(oflanagan_stress_geneset, convertHumanGeneList)

# Save as RData
save(mm_sgenes, mm_g2mgenes, file = paste0(path_supportingdata, "mm_cc_genes.RData"))
save(oflanagan_stress_geneset, file = paste0(path_supportingdata, "StressCoreGeneSet_MM.RData"))
