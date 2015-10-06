######################################################################
######################################################################
##
## Rebecca Lowdon | 6 Oct 2015
## 


library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)
# args<-commandArgs(TRUE)

#####
## Read in an assemble data into matrix m

files <- paste("fl_in")  #args[1]
fl <- read.table(files,header=FALSE)

d <- matrix( ncol = 5)
colnames(d)<-c("EID","rs","regulatory","transcribed","silent")

for( i in 1:length(fl$V1) ){
#for( i in 1:3 ){
  try(dat <- read.delim2(paste(fl$V1[i]),skip=1,header=FALSE),TRUE)
  try(dat$V12 <- c(i))
  try(dat <- cbind(dat[,12],dat[,4],dat[,9:11]))
  # try(dat <- select(dat, 12,4,9:11))
  try(colnames(dat)<-c("EID","rs","regulatory","transcribed","silent"))
  try(d <- rbind(d, dat),TRUE)
}  

d <- filter(d, !is.na(d[,1]))  ## to remove the row of "NA" created when matrix initialized
m <- melt(d,measure.vars = c(3,4,5))
m <- cbind(m, as.numeric(m[,4]))  ## Make 0/1 a numeric value
colnames(m) <- c("EID","rs","variable","value","VAL")

#######
## Build COUNTS tables and FREQUENCY heatmaps for each RE class

## Regulatory class
r <- filter(m, variable=="regulatory")
rw <- dcast(r, rs ~ EID, value.var = "VAL")

lmat = rbind(c(0,3),c(2,1),c(0,4))
lhei = c(1,4,1.2)
lwid = c(1.5,4)

png("regFreq_allEID_unsorted.png")
heatmap.2(as.matrix(rw[,2:128]), dendrogram = "none", trace="none",
          Rowv = FALSE, Colv = FALSE,
          col=c("lightgray","red"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nRegulatory elements by EID")
dev.off();

png("regFreq_allEID.png")
heatmap.2(as.matrix(rw[,2:128]), dendrogram = "both", trace="none",
          col=c("lightgray","red"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nRegulatory elements by EID")
dev.off();

## Transcribed class
t <- filter(m, variable=="transcribed")
tw <- dcast(t, rs ~ EID, value.var = "VAL")

png("txFreq_allEID_unsorted.png")
heatmap.2(as.matrix(tw[,2:128]), dendrogram = "none", trace="none",
          Rowv = FALSE, Colv = FALSE,
          col=c("lightgray","green"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nTranscribed regions by EID")
dev.off();

png("txFreq_allEID.png")
heatmap.2(as.matrix(tw[,2:128]), dendrogram = "both", trace="none",
          col=c("lightgray","green"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nTranscribed regions by EID")
dev.off();

## Silent class
s <- filter(m, variable=="silent")
sw <- dcast(s, rs ~ EID, value.var = "VAL")

png("silFreq_allEID_unsorted.png")
heatmap.2(as.matrix(sw[,2:128]), dendrogram = "none", trace="none",
          Rowv = FALSE, Colv = FALSE,
          col=c("lightgray","darkblue"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nQuiescent regions by EID")
dev.off();

png("silFreq_allEID.png")
heatmap.2(as.matrix(sw[,2:128]), dendrogram = "both", trace="none",
          col=c("lightgray","darkblue"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nQuiescent regions by EID")
dev.off();

#####
## Import vector of RE space for EIDs and 
## build ENRICHMENT matrix and heatmap for each class

reg_f <- paste("../../chromHMM_REspace/cellType_byClass/allEID/vector/allEID_regul_1-2-3-6-7.vect")
reg_vect <- t(read.delim(reg_f,header=F))
tx_f <- paste("../../chromHMM_REspace/cellType_byClass/allEID/vector/allEID_transcribed_4-5.vect")
tx_vect <- t(read.delim(tx_f,header=F))
s_f <- paste("../../chromHMM_REspace/cellType_byClass/allEID/vector/allEID_quies_repeat.vect")
s_vect <- t(read.delim(s_f,header=F))
vects <- rbind(reg_vect,tx_vect,s_vect)

reg_sum <- colSums(rw[,2:128])
tx_sum <- colSums(tw[,2:128])
s_sum <- colSums(sw[,2:128])
sums <- rbind(reg_sum,tx_sum,s_sum)

mat <- (sums/vects)/( 1016 /3095691400)

lmat = rbind(c(0,3),c(2,1),c(0,4))
lhei = c(1,4,1.2)
lwid = c(1.5,4)

png("enrichm_allEID.png")
heatmap.2(t(mat), cexCol = 1, dendrogram = "both", trace="none",
          lmat=lmat, lhei = lhei, lwid = lwid,
          main="Enrichment of GWAS SNPs in\nRegulatory element classes by EID")
dev.off();
