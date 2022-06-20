### using filtered results from cellranger
library(scater)
library(DropletUtils)

B10M21_sce <- read10xCounts(getwd(), col.names=T)

## B10M21 is AID F Autoimmune
## rename "cells"
colnames(B10M21_sce) <- gsub("-1","_B10M21",colnames(B10M21_sce))
## add Condition
colData(B10M21_sce)$Condition <- rep("Autoimmune", 1713)
## add Sex
colData(B10M21_sce)$Sex <- rep("female", 1713)
## add mouse ID
colData(B10M21_sce)$mouseID <- rep("M21", 1713)

B10M21_sce$cell_id <- colnames(B10M21_sce)

saveRDS(B10M21_sce,'B10M21.rds')
