######################################################################
######################################################################
##
## Rebecca Lowdon | 6-27 Oct 2015
## Code for plotting frequency and enrichment plots of GWAS cancer SNPs. 
## This version uses heatmap.3 and plots col and row side colors.
##
######################################################################
######################################################################

library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)
# args<-commandArgs(TRUE)

##########
## Read in an assemble data into matrix m

files <- paste("fl_noENCODE_in")  #args[1]
fl <- read.table(files,header=FALSE)

d <- matrix( ncol = 6)
colnames(d)<-c("EID","rs","pubmedid","regulatory","transcribed","silent")

for( i in 1:length(fl$V1) ){
  try(dat <- read.delim2(paste(fl$V1[i]),skip=1,header=FALSE),TRUE)
  try(dat$V12 <- c(i))
  try(dat <- cbind(dat[,12],dat[,4],dat[,5],dat[,9:11]))
  try(colnames(dat)<-c("EID","rs","pubmedid","regulatory","transcribed","silent"))
  try(d <- rbind(d, dat),TRUE)
}  

d <- filter(d, !is.na(d[,1]))       ## to remove the row of "NA" created when matrix initialized
m <- melt(d,measure.vars = c(4,5,6))
m <- cbind(m, as.numeric(m[,5]))    ## Make 0/1 a numeric value
colnames(m) <- c("EID","rs","pubmedid","variable","value","VAL")
write.table(m, file="m_111_noENCODE.txt",quote=FALSE,sep = "\t",row.names = FALSE)

##########
## Build COUNTS tables for each RE class

## Regulatory class: read in, cast data to wide
r <- filter(m, variable=="regulatory")
rw <- dcast(r, rs + pubmedid ~ EID, value.var = "VAL")
write.table(rw, file = "rw_111_noENCODE.txt",quote=FALSE,sep = "\t",row.names = FALSE)

## Transcribed class
t <- filter(m, variable=="transcribed")
tw <- dcast(t, rs + pubmedid ~ EID, value.var = "VAL")

## Silent class
s <- filter(m, variable=="silent")
sw <- dcast(s, rs + pubmedid ~ EID, value.var = "VAL")

##########
## Colors for ColSideColors arg for CELL EIDs
## First bulid custom palatte:
pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
          "#A65628","#F781BF","#7FC97F","#BEAED4","#FDC086","#FFFF99",
          "#1B9E77","#F0027F","#BF5B17","#666666","#000000","#33FF33")

## Read in eid_group.txt and assign a color from topo.colors() to each cell group:
cells <- read.delim2(file = "eid_noENCODE_group.txt", header=FALSE)
un_cells <- data.frame(unique(cells$V2))
un_cells$cols <- pal #terrian.colors(length(un_cells[,1]))   ## assign each trait a unique color
colnames(un_cells) <- c("V2","cols")

a <- right_join(cells,un_cells,by="V2")      ## use dplyr join to assign cell colors to EIDs
a <- as.matrix(a[ order(a[,1]), ])
cell_cols <- as.matrix(a[,3])

##########
## colors for RowSideColors arg -- PUBMEDID:
## Read in pubmedid.txt and assign a color from rainbow() to each:
ids <- read.delim2("pubmedids.txt",header=FALSE )
ids$cols <- rainbow(length(ids$V1))

## Create a matrix of colors for pubmedids:
id_cols <- matrix( nrow = length(rw$pubmedid) )

for( i in 1:length(rw$pubmedid)) {        ## iterating through the list of pumedids,
    for( j in 1:length(ids$V2)) {         ## and through the id - color assignments,
      if( rw$pubmedid[i]==ids$V2[j] ){    ## if the current pubmedid == the id in the ref table,
        id_cols[i,1]  <- ids$cols[j]      ## assign that id's color to the id_cols vector
      }
    }
}

## colors for RowSideColors arg -- TRAITS:
## Read in pmid_wPhenotype.txt and assign a color from topo.colors() to each:
traits <- read.delim2(file = "pmid_wPhenotype.txt", header=FALSE)
tcols <- topo.colors(93)                  ## there are 93 unique traits

## Assign each trait a unique color:
n=2
for( i in 2:length(traits$V1)-1 ) {
      j = i+1  
      traits$col[1] <- tcols[1]
      traits$col[2] <- tcols[2]
      if( traits$V2[i] != traits$V2[j]) {
            n=n+1
            traits$col[j] <- tcols[n]
          } else {
            traits$col[j] <- tcols[n] 
          }
}

## Create a matrix of colors for traits:
tr_cols <- matrix( nrow = length(rw$pubmedid) )

for( i in 1:length(rw$pubmedid)) {       ## iterating through the list of pumedids,
  for( j in 1:length(traits$V1)) {       ## and through the traits - color assignments,
    if( rw$pubmedid[i]==traits$V1[j] ){  ## if the current pubmedid == the id in the ref table,
      tr_cols[i,1]  <- traits$col[j]     ## assign that trait's color to the traits_cols vector
    }
  }
}

rcolors <- as.matrix(t(cbind(id_cols,tr_cols)))   ## cbind, transform vector, and make into a matrix.
rownames(rcolors) <- c("id_colors","trait_colors")


##########
## Plotting FREQUENCY heatmpas
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

heatmap.3(as.matrix(rw[,3:113]), dendrogram = "both", trace="none",
          hclustfun = myclust, distfun = mydist,
#          Colv = FALSE, 
#          Rowv = FALSE, 
         ColSideColors = cell_cols,
         ColSideColorsSize = 2, RowSideColorsSize = 2,
         RowSideColors = rcolors,
         col=c("lightgray","red"), #lmat=lmat, lhei = lhei, lwid = lwid,
         main="Counts of GWAS SNPs in\nRegulatory elements by EID")


##########
## Import vector of RE space for EIDs and 
## build ENRICHMENT matrix and heatmap for each class

reg_f <- paste("../../chromHMM_REspace/cellType_byClass/allEID/vector/allEID_regul_1-2-3-6-7.vect")
reg_vect <- t(read.delim(reg_f,header=F))
tx_f <- paste("../../chromHMM_REspace/cellType_byClass/allEID/vector/allEID_transcribed_4-5.vect")
tx_vect <- t(read.delim(tx_f,header=F))
s_f <- paste("../../chromHMM_REspace/cellType_byClass/allEID/vector/allEID_quies_repeat.vect")
s_vect <- t(read.delim(s_f,header=F))
vects <- rbind(reg_vect[1:111],tx_vect[1:111],s_vect[1:111])

reg_sum <- colSums(rw[,3:113])
tx_sum <- colSums(tw[,3:113])
s_sum <- colSums(sw[,3:113])
sums <- rbind(reg_sum,tx_sum,s_sum)

mat <- (sums/vects)/( 1016 /3095691400)

col = c("darkgray","white","lightyellow","yellow","gold","orange","red","darkred")

heatmap.3(mat, dendrogram = "column", trace="none",
          hclustfun = myclust, distfun = mydist,
          #Colv = FALSE, 
          Rowv = FALSE, 
          ColSideColors = cell_cols, ColSideColorsSize = 2, 
          col=col, 
          labRow = c("R","T","S"), labCol = c(""),
          main="Enrichment of GWAS SNPs in each element class by EID")

##########
## EID colors legend

leg_text <- c("IMR90", "ESC", "iPSC","ES-derived",
              "Blood & T-cell","HSC & B-cell","Mesenchyme",
              "Myosat","Epithelial","Neuropsh",
              "Thymus","Brain","Adipose",
              "Muscle","Heart","Smooth Muscle",
              "Digestive","Other")   # category labels
leg_col <- c("#4DAF4A", "#E41A1C", "#984EA3","#377EB8",
             "#F781BF","#A65628","#FF7F00",
             "#7FC97F","#FFFF33","#BEAED4",
             "#33FF33","#F0027F","#FDC086",
             "#000000","#FFFF99","#666666",
             "#BF5B17","#F0027F")    # color key

par(lend = 1)              # square line ends for the color legend
legend("center",         # location of the legend on the heatmap plot
       legend = leg_text,
       col = leg_col,
       lty= 1,             # line style
       lwd = 10            # line width
)

