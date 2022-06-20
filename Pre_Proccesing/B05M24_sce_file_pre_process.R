### using filtered results from cellranger
library(scater)
library(DropletUtils)

B05M24_sce <- read10xCounts(getwd(), col.names=T)

## B05M24 is AID M autoimmune
## rename "cells"
colnames(B05M24_sce) <- gsub("-1","_B05M24",colnames(B05M24_sce))
## add Condition
colData(B05M24_sce)$Condition <- rep("Autoimmune", 2844)
## add Sex
colData(B05M24_sce)$Sex <- rep("male", 2844)
## add mouse ID
colData(B05M24_sce)$mouseID <- rep("M24", 2844)

B05M24_sce$cell_id <- colnames(B05M24_sce)

saveRDS(B05M24_sce,'B01M11.rds')

