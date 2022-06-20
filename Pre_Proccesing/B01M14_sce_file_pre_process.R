### using filtered results from cellranger
library(scater)
library(DropletUtils)

B01M14_sce <- read10xCounts(getwd(), col.names=T)

## B01M14 is AID M Immunized
## rename "cells"
colnames(B01M14_sce) <- gsub("-1","_B01M14",colnames(B01M14_sce))
## add Condition
colData(B01M14_sce)$Condition <- rep("Immunized", 2653)
## add Sex
colData(B01M14_sce)$Sex <- rep("male", 2653)
## add mouse ID
colData(B01M14_sce)$mouseID <- rep("M14", 2653)

B01M14_sce$cell_id <- colnames(B01M14_sce)

saveRDS(B01M14_sce,'B01M14.rds')
