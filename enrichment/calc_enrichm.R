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

inf <- args[1]
files <- args[2]
fl <- read.table(files,header=FALSE)
outf <- paste(inf,"_counts",sep="")
#reportf <- paste(inf,"_report",sep="")

### Automation:
### How can I get the code to spit out the right # of fl$V1[i] ?
###    for(i in 1:length(fl$V1)) { a=cat(paste("fl$V1[i], ")) }

## Use command line annotateBed function to get counts
cmd = paste("annotateBed -counts -i",inf,"-files",fl$V1[1],fl$V1[2],fl$V1[3],"-names regul transcr silent >",outf,sep=" ")

system(cmd)

## Data are the counts generated in outf above
dat <- tbl_df(read.delim2(paste(outf),header=TRUE))
colnames(dat)<-c("chr","start","end","rs","trait","context","intergenic","negLogPval","regul","transcr","silent")

## Count # of snps in each genome element class with sum()
snps_regul<- sum(dat$regul)
snps_trans <- sum(dat$transcr)
snps_sil <- sum(dat$silent)
tot_snps <- length(dat$chr)

## Need total bp in each genome element class, specific to cell types input
## Currently still hard-coded. Figure out how to make this an arg.
regl_bp <- 19929600
tx_bp <- 451237000
sl_bp <- 2524524800
tot_bp <- 3095691400

## Do the math! Enrichment calclation
regul_enr <- ( (snps_regul / regl_bp) / (tot_snps / tot_bp) )
tx_enr <- ( (snps_trans/tx_bp) / (tot_snps/tot_bp) )
sil_enr <- ( (snps_sil/sl_bp) / (tot_snps/tot_bp) )

sink("report.txt")
cat(paste("annotateBed command to calculate counts table:\n\n"))
cat(paste(cmd),sep="\n",append=TRUE)
cat(paste("\n=========================\n\nEnrichments:\n"),append=TRUE)
cat(paste("\nRegulatory SNP enrichment: ",regul_enr),sep="\n",append=TRUE)
cat(paste("Transcribed regions SNP enrichment: ",tx_enr),sep="\n",append=TRUE)
cat(paste("Silent regions SNP enrichment: ",sil_enr),sep="\n",append=TRUE)
sink()
