### using filtered results from cellranger
library(scater)
library(DropletUtils)

B06M13_sce <- read10xCounts(getwd(), col.names=T)

## B06M13 is AID F Autoimmune
## rename "cells"
colnames(B06M13_sce) <- gsub("-1","_B06M13",colnames(B06M13_sce))
## add Condition
colData(B06M13_sce)$Condition <- rep("Autoimmune", 5618)
## add Sex
colData(B06M13_sce)$Sex <- rep("female", 5618)
## add mouse ID
colData(B06M13_sce)$mouseID <- rep("M13", 5618)

B06M13_sce$cell_id <- colnames(B06M13_sce)

saveRDS(B06M13_sce,'B06M13.rds')
