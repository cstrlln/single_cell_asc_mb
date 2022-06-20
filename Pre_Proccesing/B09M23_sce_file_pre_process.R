### using filtered results from cellranger
library(scater)
library(DropletUtils)

B09M23_sce <- read10xCounts(getwd(), col.names=T)

## B09M23 is AID M Immunized
## rename "cells"
colnames(B09M23_sce) <- gsub("-1","_B09M23",colnames(B09M23_sce))
## add Condition
colData(B09M23_sce)$Condition <- rep("Immunized", 583)
## add Sex
colData(B09M23_sce)$Sex <- rep("male", 583)
## add mouse ID
colData(B09M23_sce)$mouseID <- rep("M23", 583)

B09M23_sce$cell_id <- colnames(B09M23_sce)

saveRDS(B09M23_sce,'B09M23.rds')
