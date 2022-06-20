### using filtered results from cellranger
library(scater)
library(DropletUtils)


B34M34_sce <- read10xCounts(getwd(), col.names=T)

## B34M34 is AID F Autoimmune
## rename "cells"
colnames(B34M34_sce) <- gsub("-1","_B34M34",colnames(B34M34_sce))
## add Condition
colData(B34M34_sce)$Condition <- rep("Autoimmune", 2528)
## add Sex
colData(B34M34_sce)$Sex <- rep("female", 2528)
## add mouse ID
colData(B34M34_sce)$mouseID <- rep("M34", 2528)

B34M34_sce$cell_id <- colnames(B34M34_sce)

saveRDS(B34M34_sce,'B34M34.rds')