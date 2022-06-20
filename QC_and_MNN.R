suppressPackageStartupMessages({
  require(scater)
  require(scran)
  require(batchelor)
  require(mbkmeans)
  require(biomaRt)
  require(DelayedArray)
  require(patchwork) ## for putting ggplots together
  require(BiocParallel)
  require(slingshot)
  require(RColorBrewer)
  require(dittoSeq)
  require(UpSetR)
  require(pheatmap)})
#library(ggplotify)

B01M11 <- readRDS("B01M11.rds")
B01M14 <- readRDS("B01M14.rds")
B02M17 <- readRDS("B02M17.rds")
B05M24 <- readRDS("B05M24.rds")
B06M13 <- readRDS("B06M13.rds")
B09M23 <- readRDS("B09M23.rds")
B10M21 <- readRDS("B10M21.rds")
B34M34 <- readRDS('B34M34.rds')

keep_genes <- Reduce(intersect, list(rownames(B01M11), rownames(B01M14), rownames(B02M17) ,rownames(B05M24), rownames(B06M13), rownames(B09M23), rownames(B10M21), rownames(B34M34)))

B01M11_keep <- B01M11[keep_genes,]
B01M14_keep <- B01M14[keep_genes,]
B02M17_keep <- B02M17[keep_genes,]
B05M24_keep <- B05M24[keep_genes,]
B06M13_keep <- B06M13[keep_genes,]
B09M23_keep <- B09M23[keep_genes,]
B10M21_keep <- B10M21[keep_genes,]
B34M34_keep <- B34M34[keep_genes,]

counts_mix <- cbind(counts(B01M11_keep), counts(B01M14_keep), counts(B02M17_keep), counts(B05M24_keep), counts(B06M13_keep), counts(B09M23_keep), counts(B10M21_keep), counts(B34M34_keep))

sce <- SingleCellExperiment( 
  assays = list(counts = counts_mix),  
  rowData = rowData(B01M11_keep),
  colData = rbind(colData(B01M11_keep), colData(B01M14_keep), colData(B02M17_keep), colData(B05M24_keep), colData(B06M13_keep), colData(B09M23_keep), colData(B10M21_keep), colData(B34M34_keep)))



sce <- getBMFeatureAnnos(sce, ids = rownames(sce), attributes = c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position", "chromosome_name", "percentage_gene_gc_content", "external_gene_name", "gene_biotype"),dataset ="mmusculus_gene_ensembl", host="https://may2021.archive.ensembl.org")



rowData(sce["EYFP",])$gene_biotype <- "protein_coding"
rowData(sce["CREREC",])$gene_biotype <- "protein_coding"

rowData(sce["EYFP",])$mgi_symbol <- "eyfp"
rowData(sce["CREREC",])$mgi_symbol <- "crerec"


sce_original <- sce

table(rowData(sce)$gene_biotype)

ig_list <- c("IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_D_pseudogene", "IG_J_gene", "IG_LV_gene", "IG_pseudogene", "IG_V_gene", "IG_V_pseudogene")

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="https://may2021.archive.ensembl.org")

test_ig <- getBM(filters=c("ensembl_gene_id","biotype"), values = list(rownames(sce), ig_list), attributes=c("mgi_symbol", "gene_biotype", "ensembl_gene_id"), mart= ensembl)

table(rownames(sce)%in% test_ig$ensembl_gene_id) 

ig_genes <- which(rowData(sce)$mgi_symbol %in% test_ig$mgi_symbol)

ribo_genes <- grepl('^Rp[ls].*', rowData(sce)$mgi_symbol)

test_ribo <- getBM(filters=c("ensembl_gene_id","go"), values = list(rownames(sce), "GO:0005840"), attributes=c("mgi_symbol", "gene_biotype", "ensembl_gene_id"), mart= ensembl)

summary(test_ribo) 

sub <- sce[ribo_genes,]
ribo_big <- rownames(sub) # has only the ribo genes from grepl rp[ls]
ribo_big <- rowData(sub)$mgi_symbol
ribo_small <- test_ribo["mgi_symbol"]
ribo_small <- ribo_small[,1] # to have it also as characters

all_ribo <- c(ribo_big, ribo_small)
all_ribo <- unique(all_ribo) # to get only unique values

final_ribo <- rowData(sce)$mgi_symbol %in% all_ribo


## mito genes

mt_genes <- which(rowData(sce)$chromosome_name == "MT")
length(mt_genes)

mito_genes <- grepl('^mt', rowData(sce)$mgi_symbol, ignore.case=T)

table(mito_genes)

tr_factors <- getBM(filters=c("ensembl_gene_id","go"), values = list(rownames(sce), "GO:0003700"), attributes=c("ensembl_gene_id", "mgi_symbol", "gene_biotype"), mart= ensembl)
        
        
tr_factors <- tr_factors["mgi_symbol"]
        
tr_factors <- tr_factors[,1]
        
        
        
## t cells genes
        
tcell_list <- c("TR_C_gene","TR_D_gene" ,"TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene")
        
test_tcell <- getBM(filters=c("ensembl_gene_id","biotype"), values = list(rownames(sce), tcell_list), attributes=c("mgi_symbol", "gene_biotype", "ensembl_gene_id"), mart= ensembl)
        
table(rownames(sce)%in% test_tcell$ensembl_gene_id) 
        
tcell_genes <- rowData(sce)$mgi_symbol %in% test_tcell$mgi_symbol
        
ig_genes <- rowData(sce)$mgi_symbol %in% test_ig$mgi_symbol
        
ctrls <- list(mito = mito_genes, ribo = final_ribo, ig = ig_genes, tcr = tcell_genes)
        
sce_qc <- perCellQCMetrics(sce, subsets=ctrls)
        
colData(sce) <- cbind(colData(sce), sce_qc)
        
  
        
### I would like to make a list of cells that have t cell receptors:
### this code identifies cells which have a count > 1 for tr genes, and subset the cell#.
        
tr_pos_cells <- summary(counts(sce)[rownames(sce)%in%test_tcell$ensembl_gene_id,])[2]$j
        
## now select only single cell# info
        
tr_pos_cells <- unique(tr_pos_cells)
        
## using the cell# we can extract the actual cell barcode name:
        
tr_pos_cells <- colnames(sce[,tr_pos_cells])
        
## removing tr_pos_cells
        
sce_tneg <- sce[,!(colnames(sce)%in%tr_pos_cells)]
        
## right before removing tcells
   
sce <- sce_tneg
        
qc_lib <- sce$sum > 1000
metadata(sce)$qc_libsize <- qc_lib
discard <- qc_lib
filtered <- sce[,discard]

qc_detected <- filtered$detected > 1000
metadata(filtered)$qc_detected <- qc_detected
discard <- qc_detected
filteredb <- filtered[,discard]
high.mito <- filteredb$subsets_mito_percent > 5
metadata(filteredb)$high_mito <- high.mito
discard <- high.mito
filteredc <- filteredb[,!discard]
        
high.ribo <- isOutlier(filteredc$subsets_ribo_percent, nmads=2, type="higher")

high.ribo <- filteredc$subsets_ribo_percent > 40
metadata(filteredc)$high_ribo <- high.ribo
discard <- high.ribo
filteredd <- filteredc[,!discard]

# test_ig already contains list of ig associated genes on the dataset
sce_qc_igless <- filteredd[!(rowData(filteredd)$gene_biotype %in% ig_list),]
## verifying IG genes are really clean
        
igrest <- grepl('^ig[hkl].*', rowData(sce_qc_igless)$mgi_symbol, ignore.case=T)

goodgenes <- c("protein_coding")

sce_qc_igless_safe <- sce_qc_igless
        
sce_qc_igless <- sce_qc_igless[(rowData(sce_qc_igless)$gene_biotype %in% goodgenes),]
        
table(rowData(sce_qc_igless)$gene_biotype)
        
dim(sce_qc_igless)
        
## feature QC
        
per.feat <- perFeatureQCMetrics(sce_qc_igless)
        
ave <- calculateAverage(sce_qc_igless)
filtergenes <- Matrix::rowSums(assay(sce_qc_igless)>1)>5
sce_qc_filt <- sce_qc_igless[filtergenes,]
dim(sce_qc_filt)
sce <- sce_qc_filt

   
set.seed(1234)
qclust<- quickCluster(sce,  BPPARAM = MulticoreParam(8), min.mean=0.1) 
table(qclust)
        
###
        
sce <- computeSumFactors(sce, cluster=qclust) ## with more than 4 it fails .
        
sce <- logNormCounts(sce)
        
tr_factors <- getBM(filters=c("ensembl_gene_id","go"), values = list(rownames(sce), "GO:0003700"), attributes=c("ensembl_gene_id", "mgi_symbol", "gene_biotype"), mart= useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="https://may2021.archive.ensembl.org"))
tr_factors <- tr_factors["ensembl_gene_id"]
tr_factors <- tr_factors[,1]
        
##
design <- model.matrix(~Sex + Condition, colData(sce))
dec <- modelGeneVar(sce, design = design)
hvg <- getTopHVGs(dec, var.threshold=0)
genelist <- union(hvg, tr_factors)
        
### now run PCA
        
set.seed(1234)
sce2 <- runPCA(sce, subset_row = genelist)
        
sce2 <- runUMAP(sce2, dimred = 'PCA')
        
set.seed(1234)
holder <- runUMAP(sce2, dimred = 'PCA', n_dimred = 6,BPPARAM=MulticoreParam(8))
reducedDim(sce2, "UMAP_n6") <- reducedDim(holder, "UMAP")
rm(holder)
        
sce_mnn_simple <- sce2
set.seed(1234)
mnn_simple_out <- batchelor::fastMNN(sce_mnn_simple,
                                             subset.row = genelist,
                                             batch = sce_mnn_simple$mouseID)

reducedDim(sce_mnn_simple, 'MNN') <- reducedDim(mnn_simple_out, 'corrected')

set.seed(1234)
holder <- runUMAP(sce_mnn_simple, dimred = "MNN",BPPARAM=MulticoreParam(8))
reducedDim(sce_mnn_simple, "UMAP_MNN") <- reducedDim(holder, "UMAP")
rm(holder)
        
set.seed(1234)
holder <- runUMAP(sce_mnn_simple, dimred = "MNN", n_dimred=6,BPPARAM=MulticoreParam(8))
reducedDim(sce_mnn_simple, "UMAP_MNN_n6") <- reducedDim(holder, "UMAP")
rm(holder)
        
sce <- sce_mnn_simple

saveRDS(sce,'sce_MNN.rds')
        
