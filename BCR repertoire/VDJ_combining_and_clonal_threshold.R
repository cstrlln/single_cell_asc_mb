library(shazam)

## adding metadata to each VDJ file
## put heave chain info together

B01M11 <- read.table("B01M11_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B01M11$cell_id <- gsub("-1","_B01M11", B01M11$cell_id)
B01M11$MID <- rep('M11', length(B01M11$cell_id))
B01M11$Condition <- rep('Immunized', length(B01M11$cell_id))

B09M23 <- read.table("B09M23_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B09M23$cell_id <- gsub("-1","_B09M23", B09M23$cell_id)
B09M23$MID <- rep('M23', length(B09M23$cell_id))
B09M23$Condition <- rep('Immunized', length(B09M23$cell_id))

B01M14 <- read.table("B01M14_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B01M14$cell_id <- gsub("-1","_B01M14", B01M14$cell_id)
B01M14$MID <- rep('M14', length(B01M14$cell_id))
B01M14$Condition <- rep('Immunized', length(B01M14$cell_id))

##

B05M24 <- read.table("B05M24_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B05M24$cell_id <- gsub("-1","_B05M24", B05M24$cell_id)
B05M24$MID <- rep('M24', length(B05M24$cell_id))
B05M24$Condition <- rep('Autoimmune', length(B05M24$cell_id))

B06M13 <- read.table("B06M13_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B06M13$cell_id <- gsub("-1","_B06M13", B06M13$cell_id)
B06M13$MID <- rep('M13', length(B06M13$cell_id))
B06M13$Condition <- rep('Autoimmune', length(B06M13$cell_id))

B10M21 <- read.table("B10M21_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B10M21$cell_id <- gsub("-1","_B10M21", B10M21$cell_id)
B10M21$MID <- rep('M21', length(B10M21$cell_id))
B10M21$Condition <- rep('Autoimmune', length(B10M21$cell_id))

B02M17 <- read.table("B02M17_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B02M17$cell_id <- gsub("-1","_B02M17", B02M17$cell_id)
B02M17$MID <- rep('M17', length(B02M17$cell_id))
B02M17$Condition <- rep('Autoimmune', length(B02M17$cell_id))

B34M34 <- read.table("B34M34_annotated_db-pass.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")
B34M34$cell_id <- gsub("-1","_B34M34", B34M34$cell_id)
B34M34$MID <- rep('M34', length(B34M34$cell_id))
B34M34$Condition <- rep('Autoimmune', length(B34M34$cell_id))
########

# combining data
all <- rbind(B01M11, B09M23, B01M14, B05M24, B06M13, B10M21, B02M17, B34M34)

all_clean <- all[all$PRODUCTIVE_10X == 'True' ,]

##

write.table(all_clean,"all_clean.tab",sep="\t",row.names=FALSE)

##

db <- read.table("all_clean.tab", header=T, fill=T, stringsAsFactors=F, sep="\t")

###


### take out any duplicated cell out by first looking in the heavy compartment. The duplications
### are coming from cells with multiple contigs.

db_heavy <- db[db$locus=='IGH',]
db_light <- db[db$locus!='IGH',]


duplicated_heavy <- db_heavy$cell_id[duplicated(db_heavy$cell_id)]
duplicated_light <- db_light$cell_id[duplicated(db_light$cell_id)]


db_heavy = db_heavy[order(db_heavy[,'cell_id'],-db_heavy[,'umi_count']),]
db_heavy = db_heavy[!duplicated(db_heavy$cell_id),]

##

db2 <- rbind(db_heavy,db_light)

##
## looking for threshold for clonal assignment

dist_mid <- distToNearest(db2, cellIdColumn="cell_id", locusColumn="locus", VJthenLen=FALSE, onlyHeavy=FALSE, model="ham", normalize="len", fields = 'MID')
output_mid <- findThreshold(dist_mid$dist_nearest, method="density")
threshold_mid <- output_mid@threshold

#threshold_mid
#[1] 0.1691791

## used a distance of 0.1691791 in the following steps with Immcantation in Python
write.table(db2,"all_db_clean.tab",sep="\t",row.names=FALSE)



