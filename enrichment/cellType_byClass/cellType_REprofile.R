######################################################################
######################################################################
##
## Rebecca Lowdon | 2 Oct 2015
## 
## This code reads in the SNP annotations calculated from 
## the counts.sh script [INPUT] and assembles them into a matrix 
## with 5 columns: EID (#), rs accession, regulatory, transcribed, and silent.
## The last 3 columns have a value of 0 or 1 depending on which class 
## the SNP falls in the given EID annotation file.
##
## Usage: 
## * interactive: edit line 27 to include file with list of input files
## * command line: /apps/bin/Rscript cellType_REprofile.R <file list>
##
######################################################################
######################################################################

library(dplyr)
library(ggplot2)
library(reshape2)
# args<-commandArgs(TRUE)

#####
## Read in an assemble data into matrix m

files <- paste("fl_in")  #args[1]
fl <- read.table(files,header=FALSE)

d <- matrix( ncol = 5)
colnames(d)<-c("EID","rs","regulatory","transcribed","silent")

#for( i in 1:length(fl$V1) ){
for( i in 1:3 ){
  try(dat <- read.delim2(paste(fl$V1[i]),skip=1,header=FALSE),TRUE)
  try(dat$V12 <- c(i))
  try(dat <- cbind(dat[,12],dat[,4],dat[,9:11]))
  try(colnames(dat)<-c("EID","rs","regulatory","transcribed","silent"))
  try(d <- rbind(d, dat),TRUE)
}  

d <- filter(d, !is.na(d[,1]))  ## to remove the row of "NA" created when matrix initialized
m <- melt(d,measure.vars = c(3,4,5))

#####
## Plot rs x element class, split by EID

g <- ggplot(m, aes(variable, rs)) + geom_tile(aes(fill=value)) + facet_wrap(~EID)
g <- g + theme(axis.text.y = element_blank(), axis.text.x=element_text(angle=30,vjust=0.5)) + 
  scale_x_discrete(name="Regulatory element class") 

png("EID_REprofiles.png")
g
dev.off();

#####
## Plot EID x rs accession, split by class

g <- ggplot(m, aes(EID, rs)) + geom_tile(aes(fill=value)) + facet_wrap(~variable)
g <- g + theme(axis.text.y = element_blank())

png("EID_by_RS_profiles.png")
g
dev.off();
