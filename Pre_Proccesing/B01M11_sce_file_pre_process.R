### using filtered results from cellranger
library(scater)
library(DropletUtils)

B01M11_sce <- read10xCounts(getwd(), col.names=T)

## B01M11 is AID F Immunized
## rename "cells"
colnames(B01M11_sce) <- gsub("-1","_B01M11",colnames(B01M11_sce))
## add Condition
colData(B01M11_sce)$Condition <- rep("Immunized", 2081)
## add Sex
colData(B01M11_sce)$Sex <- rep("female", 2081)
## add mouse ID
colData(B01M11_sce)$mouseID <- rep("M11", 2081)

B01M11_sce$cell_id <- colnames(B01M11_sce)

saveRDS(B01M11_sce,'B01M11.rds')
