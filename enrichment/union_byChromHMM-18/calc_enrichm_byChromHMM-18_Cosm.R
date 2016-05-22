##########
##
## Programmer: Rebecca
## Date: March 2016 -- Modified from calc_enrichm_Cosm.R (orig. August 2015)
## Purpose: Calculates enrichment of a given SNP set in a cell type x genomic element space.
##
## Copied from calc_enrichm_byState_Cosm.R on 14 March 2016
##
## Usage: /apps/bin/Rscipt calc_enrichm_byState_Cosm.R <SNP set> <file list of CT genomic elements> <file with path to CT *unionStates.report> >
##   * Can be made interative: lines 12,14,15,16
##########

setwd("/bar/rlowdon/cancer_annotation/enrichment/union_byChromHMM-18/")

library(dplyr)
library(ggplot2)
args<-commandArgs(TRUE)

inf <- args[1] #paste("liver_keep_L3-WGS_hg19.bed") #args[1]         ## bedfile of SNPs to query
files <- args[2]  #paste("E003_stateRegions")  #args[2]        
stsfl <- args[3] #paste("E003_report") # args[3]

fl <- read.table(files,header=FALSE)

pref <- paste(inf,"_vs_",files,sep="")
dirPath <- pref
outf <- paste(dirPath,"/",pref,"_counts",sep="")
#if( !dir.exists(dirPath) ){
  dir.create(dirPath)
#}

##########
## Make snps table

stsf <- read.table(stsfl, header = FALSE)
sts <- read.table(paste(stsf$V1),header=FALSE,sep = "\t")
colnames(sts)  <- c("state","num","bp")

##########
## Use command line annotateBed function to get counts

fl_list = paste(paste(paste(fl$V1[1:length(fl$V1)]," ")),collapse = "")
cmd = paste("annotateBed -counts -i",inf,"-files",fl_list,"-names 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18>",outf,sep=" ")
system(cmd)

##########
## Data are the counts generated in outf above. Read in 0/1 values for RE classes.

read_in <- c("NULL","NULL","NULL","NULL","NULL","NULL","NULL",
             "integer","integer","integer",
             "integer","integer","integer",
             "integer","integer","integer",
             "integer","integer","integer",
             "integer","integer","integer",
             "integer","integer","integer")
#read_in <- c("NULL","NULL","NULL",
#             "integer","integer","integer",
#             "integer","integer","integer",
#             "integer","integer","integer",
#             "integer","integer","integer",
#             "integer","integer","integer")
dat <- tbl_df(read.delim2(paste(outf),skip=1,header=FALSE,colClasses = read_in))
colnames(dat)<-c("TssA","TssFlnk","TssFlnkU","TssFlnkD","Tx","TxWk","EnhG1",
                 "EnhG2","EnhA1","EnhA2","EnhWk","ZNF_Rpts","Het","TssBiv",
                 "EnhBiv","ReprPC","ReprPCWk","Quies")

##########
## Math:

## Count # of snps in each genome element class with sum()
## and assign to the correct row in sts table
sts$snp_count[sts$state=="1_TssA"] <- sum(dat$TssA)
sts$snp_count[sts$state=="2_TssFlnk"] <- sum(dat$TssFlnk)
sts$snp_count[sts$state=="3_TssFlnkU"] <- sum(dat$TssFlnkU)
sts$snp_count[sts$state=="4_TssFlnkD"] <- sum(dat$TssFlnkD)
sts$snp_count[sts$state=="5_Tx"] <- sum(dat$Tx)
sts$snp_count[sts$state=="6_TxWk"] <- sum(dat$TxWk)
sts$snp_count[sts$state=="7_EnhG1"] <- sum(dat$EnhG1)
sts$snp_count[sts$state=="8_EnhG2"] <- sum(dat$EnhG2)
sts$snp_count[sts$state=="9_EnhA1"] <- sum(dat$EnhA1)
sts$snp_count[sts$state=="10_EnhA2"] <- sum(dat$EnhA2)
sts$snp_count[sts$state=="11_EnhWk"] <- sum(dat$EnhWk)
sts$snp_count[sts$state=="12_ZNF-Rep"] <- sum(dat$ZNF_Rpts)
sts$snp_count[sts$state=="13_Het"] <- sum(dat$Het)
sts$snp_count[sts$state=="14_TssBiv"] <- sum(dat$TssBiv)
sts$snp_count[sts$state=="15_EnhBiv"] <- sum(dat$EnhBiv)
sts$snp_count[sts$state=="16_ReprPc"] <- sum(dat$ReprPC)
sts$snp_count[sts$state=="17_ReprPcWk"] <- sum(dat$ReprPCWk)
sts$snp_count[sts$state=="18_Quies"] <- sum(dat$Quies)
sts$snp_count[sts$state=="Total"] <- sum(sts$snp_count[1:18])  ## should summ to 2807 for 1 sample in hg19 WGS NCV

##########
## Do the math! Enrichment calclation

#dnom <- (sts$snp_count[sts$state=="Total"] / sts$bp[sts$state=="Total"])
dnom <- (sts$snp_count[sts$state=="Total"] / 48500000 )
sts$numer = sts$snp_count/sts$bp
sts$enr = sts$numer/dnom

##########
## Write reports

sink("report_cmd.txt")
cat(paste("annotateBed command to calculate counts table:\n\n"))
cat(paste(cmd),sep="\n",append=TRUE)
sink()

mvcmd = paste( "mv report_cmd.txt ",dirPath,"/.",sep="")
system( mvcmd )

repf <- paste(dirPath,"/",pref,"_REPORT.tsv",sep="")
write.table(sts, file=repf,quote=FALSE,col.names = TRUE,row.names = FALSE,sep="\t")

##########
## Plot

# Reorder plotting of each state by reorder factor levels.
sts$state <- factor( sts$state, levels=paste(rev(unique(sts$state))) )   #c("1_ActTss","2_FlkTss","3_Trx5-3",
#                                     "4_StrTrx","5_WkTrx","6_GenEnh","7_Enh",
#                                     "8_ZfnRep","9_Hetero","10_BivTss",
#                                     "11_FlkBivTss","12_BivEnh","13_RepPly",
#                                     "14_WkRepPly","15_Quies","Total"))

g <- ggplot(sts[1:18,], aes(x=state,y=enr)) + geom_bar(stat="identity") + coord_flip() +
  geom_hline(yintercept = 1, col="pink")
g <- g + ylab("Enrichment") + xlab("") +
 theme(axis.text.x=element_text(size=16,color="black"),
       axis.text.y=element_text(size=16,color="black"),
       axis.title.y=element_text(size=16))

png(paste(dirPath,"/",stsfl,"_enr.png",sep=""),width=550,height = 450)
g
dev.off();
