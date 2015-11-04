######################################################################
######################################################################
##
## Rebecca Lowdon | 6-27 Oct 2015
## Code for plotting frequency and enrichment plots of GWAS cancer SNPs. 
##
######################################################################
######################################################################

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

d <- filter(d, !is.na(d[,1]))       ## to remove the row of "NA" created when matrix initialized
m <- melt(d,measure.vars = c(3,4,5))
m <- cbind(m, as.numeric(m[,4]))    ## Make 0/1 a numeric value
colnames(m) <- c("EID","rs","variable","value","VAL")

#######
## Row colors for RowSideColors arg

rowcols = c(rep("#006666",1),rep("#00FFFF",8),rep("#0099FF",5),rep("#0000FF",9),
            rep("#FF0000",14),rep("#FF00FF",9),rep("#6600FF",4),
            rep("#99FF66",1),rep("#FFCC99",8),rep("#FFFF33",2),
            rep("#FF9900",2),rep("#660066",10),rep("#990000",1),
            rep("#FF6666",5),rep("#99CCFF",5),rep("#33FF33",4),
            rep("#00CC00",12),rep("#009999",11),rep("#000000",16))

#######
## Build COUNTS tables and FREQUENCY heatmaps for each RE class

## Regulatory class: read in, cast data to wide
r <- filter(m, variable=="regulatory")
rw <- dcast(r, rs ~ EID, value.var = "VAL")

## lmat for no row side colors or legend:
# lmat = rbind(c(0,3),c(2,1),c(0,4))
# lhei = c(1,4,1.2)
# lwid = c(1.5,4)
# png("regFreq_allEID_unsorted.png")
# heatmap.2(as.matrix(rw[,2:128]), dendrogram = "none", trace="none",
#           Colv = FALSE, Rowv = FALSE, 
#           col=c("lightgray","red"), lmat=lmat, lhei = lhei, lwid = lwid,
#           main="Counts of GWAS SNPs in\nRegulatory elements by EID")
# dev.off();

## lmat for RowSideColors and legends on top right:
lmat=rbind(c(0,0,4,0), c(0,0,1,0), c(0,3,2,0),c(0,0,5,0))
lhei=c(1.5,0.5,5,1.5)
lwid=c(1,1.5,4,2)

png("regFreq_allEID_sorted.png")
heatmap.2(as.matrix(rw[,2:128]), dendrogram = "both", trace="none",
          ColSideColors = rowcols,
          col=c("lightgray","red"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nRegulatory elements by EID")
dev.off();

## Transcribed class
t <- filter(m, variable=="transcribed")
tw <- dcast(t, rs ~ EID, value.var = "VAL")

# png("txFreq_allEID_unsorted.png")
# heatmap.2(as.matrix(tw[,2:128]), dendrogram = "none", trace="none",
#           Rowv = FALSE, Colv = FALSE,
#           col=c("lightgray","green"), lmat=lmat, lhei = lhei, lwid = lwid,
#           main="Counts of GWAS SNPs in\nTranscribed regions by EID")
# dev.off();

png("txFreq_allEID_sorted.png")
heatmap.2(as.matrix(tw[,2:128]), dendrogram = "both", trace="none",
          ColSideColors = rowcols,
          col=c("lightgray","green"), lmat=lmat, lhei = lhei, lwid = lwid,
          main="Counts of GWAS SNPs in\nTranscribed regions by EID")
dev.off();

## Silent class
s <- filter(m, variable=="silent")
sw <- dcast(s, rs ~ EID, value.var = "VAL")

# png("silFreq_allEID_unsorted.png")
# heatmap.2(as.matrix(sw[,2:128]), dendrogram = "none", trace="none",
#           Rowv = FALSE, Colv = FALSE,
#           col=c("lightgray","darkblue"), lmat=lmat, lhei = lhei, lwid = lwid,
#           main="Counts of GWAS SNPs in\nQuiescent regions by EID")
# dev.off();

png("silFreq_allEID_sorted.png")
heatmap.2(as.matrix(sw[,2:128]), dendrogram = "both", trace="none",
          ColSideColors = rowcols,
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

## lmat without rowsidecols
# lmat = rbind(c(0,3),c(2,1),c(0,4))
# lhei = c(1,4,1.2)
# lwid = c(1.5,4)

col = c("darkgray","white","lightyellow","yellow","gold","orange","red","darkred")

lmat=rbind(c(0,0,0,4,0), c(0,3,1,2,0),c(0,0,0,5,0))
lhei=c(1,5,1)
lwid=c(.5,1,0.3,5,2)

png("enrichm_allEID.png")
heatmap.2(t(mat), cexCol = 1, dendrogram = "both", trace="none",
               lmat=lmat, lhei = lhei, lwid = lwid, 
               col=col, 
               breaks=c(0,0.9,1.1,1.5,1.75,2,2.25,2.75,3),
               RowSideColors = rowcols,
               main="  Enrichment of cancer-related GWAS SNPs\nin Regulatory element classes by EID")
dev.off();

#####
## EID colors legend

leg_text <- c("IMR90", "ESC", "iPSC","ES-derived",
              "Blood & T-cell","HSC & B-cell","Mesenchyme",
              "Myosat","Epithelial","Neuropsh",
              "Thymus","Brain","Adipose",
              "Muscle","Heart","Smooth Muscle",
              "Digestive","Other","ENCODE")   # category labels
leg_col <- c("#006666", "#00FFFF", "#0099FF","#0000FF",
             "#FF0000","#FF00FF","#6600FF",
             "#99FF66","#FFCC99","#FFFF33",
             "#FF9900","#660066","#990000",
             "#FF6666","#99CCFF","#33FF33",
             "#00CC00","#009999","#000000")    # color key

par(lend = 1)              # square line ends for the color legend
legend("topright",         # location of the legend on the heatmap plot
       legend = leg_text,
       col = leg_col,
       lty= 1,             # line style
       lwd = 10            # line width
)

