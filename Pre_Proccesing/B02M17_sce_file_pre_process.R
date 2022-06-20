### using filtered results from cellranger
library(scater)
library(DropletUtils)

B02M17_sce <- read10xCounts(getwd(), col.names=T)

## B02M17 is AID M Autoimmune
## rename "cells"
colnames(B02M17_sce) <- gsub("-1","_B02M17",colnames(B02M17_sce))
## add Condition
colData(B02M17_sce)$Condition <- rep("Autoimmune", 4878)
## add Sex
colData(B02M17_sce)$Sex <- rep("male", 4878)
## add mouse ID
colData(B02M17_sce)$mouseID <- rep("M17", 4878)

B02M17_sce$cell_id <- colnames(B02M17_sce)

saveRDS(B02M17_sce,'B02M17.rds')

