##########
##
## Programmer: Rebecca
## Date: August 2015
## Purpose: Calculates enrichment of a given SNP set in a cell type x genomic element space.
##
## Usage: /apps/bin/Rscipt calc_enrichm.R <SNP set> <file list of CT genomic elements> 
##
##########

library(dplyr)
args<-commandArgs(TRUE)

inf <- args[1]  #paste("chrGWAS_2015-09-17_canc_hg19.bed") #args[1]
files <- args[2]  #paste("E003_elements")  #args[2]
fl <- read.table(files,header=FALSE)
outf <- paste(inf,"_vs_",files,"_counts",sep="")
#reportf <- paste(inf,"_report",sep="")

#####
## Use command line annotateBed function to get counts

fl_list = paste(paste(paste(fl$V1[1:length(fl$V1)]," ")),collapse = "")
cmd = paste("annotateBed -counts -i",inf,"-files",fl_list,"-names regul transcr silent >",outf,sep=" ")
system(cmd)

#####
## Data are the counts generated in outf above
dat <- tbl_df(read.delim2(paste(outf),skip=1,header=FALSE))
dat <- dat[,1:11]
colnames(dat)<-c("chr","start","end","rs","pubmedID","intergenic",
                 "negLogPval","context","regul","transcr","silent")

#####
## Math:
## Count # of snps in each genome element class with sum()
snps_regul<- sum(dat$regul)
snps_trans <- sum(dat$transcr)
snps_sil <- sum(dat$silent)
tot_snps <- snps_regul + snps_trans + snps_sil

## Need total bp in each genome element class, specific to cell types input
## Currently still hard-coded. Figure out how to make this an arg.
## allESC bp counts
regl_bp <- 1035073400 #119929600
tx_bp <- 3694280400 #451237000
sl_bp <- 20036177400 #2524524800
tot_bp <- 24765531200 #3095691400

## E003 bp counts
# regl_bp <- 119929600
# tx_bp <- 451237000
# sl_bp <- 2524524800
# tot_bp <- 3095691400

## Do the math! Enrichment calclation
regul_enr <- ( (snps_regul / regl_bp) / (tot_snps / tot_bp) )
tx_enr <- ( (snps_trans/tx_bp) / (tot_snps/tot_bp) )
sil_enr <- ( (snps_sil/sl_bp) / (tot_snps/tot_bp) )

#####
## Write report

sink("report.txt")
cat(paste("annotateBed command to calculate counts table:\n\n"))
cat(paste(cmd),sep="\n",append=TRUE)
cat(paste("\n=========================\n\nEnrichments:\n"),append=TRUE)
cat(paste("\nRegulatory SNP enrichment: ",regul_enr),sep="\n",append=TRUE)
cat(paste("Transcribed regions SNP enrichment: ",tx_enr),sep="\n",append=TRUE)
cat(paste("Silent regions SNP enrichment: ",sil_enr),sep="\n",append=TRUE)
sink()
