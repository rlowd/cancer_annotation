#############
##
## Date: 11 Nov 2015 | Rebecca L.
## Purpose: Assemble & plot Observed/Expected values of Cosm NCVs in individ EIDs by State
##
#############

library(dplyr)
library(ggplot2)
library(reshape2)
# args<-commandArgs(TRUE)

## System cmd to create arg[1]:
# ls CosmicNCV_FS_tID_hg19.bed_vs_E*_stateRegions/CosmicNCV_FS_tID_hg19.bed_vs_E*_stateRegions_REPORT.tsv > reports_fl
# ls CosmNCV_2k_hg19.bed_vs_E*_stateRegions/CosmNCV_2k_hg19.bed_vs_E*_stateRegions_REPORT.tsv > reports_fl
# ls liver_L0_hg19.bed_vs_E*_stateRegions/liver_L0_hg19.bed_vs_E*_stateRegions_REPORT.tsv > reports_fl


##########
## Read in an compile data

fl <- paste("reports_fl") #args[1]
rep_fl <- read.table(fl,header=FALSE)

df <- data.frame( matrix( nrow = 1,ncol = 6 ) )
colnames(df) <- c("state","num","bp","snp_count","ex","EID")

cls <- c("factor","integer","numeric","integer","NULL","NULL")

for( i in 1:length(rep_fl$V1) ) {
	try( dat <- read.delim2(paste(rep_fl$V1[i]),colClasses = cls),TRUE )
	
	## Calculate expectation
	dat$ex <- ( dat$snp_count[dat$state=="Total"] )*( dat$bp / dat$bp[dat$state=="Total"] )  

	## annotate EID to set for dat
	pref <- unlist( strsplit(paste(rep_fl$V1[i]),split = "_") )[6]
	dat$EID <- pref
	
	## row bind new EID values to df
	try( df <- rbind( df,dat[1:18,] ),TRUE )
}


df <- filter(df, !is.na(state))
m <- select(df, c(1,4,5,6)) 

##########
## Colors for cell EIDs
## TO DO: use custom palatte for EID colors
# pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33",
#          "#A65628","#F781BF","#7FC97F","#BEAED4","#FDC086","#FFFF99",
#          "#1B9E77","#F0027F","#BF5B17","#666666","#000000","#33FF33")

## Read in EID vs tissue set table and assign EIDs in melted df to tissue set
cells <- read.delim2(file = "../cellType_byClass/eid_noENCODE_group.txt", header=FALSE)
colnames(cells) <- c("EID","set")
a <- right_join(cells,m,by="EID")      ## use dplyr join to assign cell colors to EIDs


#write.table(a,file="liver_keep-ExomSeq_hg19_obsExp.tsv",sep="\t",col.names=TRUE,row.names=FALSE,quote = FALSE)

#########
## For reading in the *tsv file made on server, above
# a <- read.table("liver_L0_obsExp.tsv",header=TRUE,sep="\t")

##########
## Prepare data for plotting
## Remove E060 and E064:
a <- filter(a, !is.na(set))

## Calculate fold obs/exp:
a$fold <- a$snp_count / a$ex

# Control order of facet_wrap plot by reordering the levels for factors 'state' and 'set'
a$state <- factor(a$state, levels=rev(c("1_TssA","2_TssFlnk","3_TssFlnkU","4_TssFlnkD","5_Tx","6_TxWk","7_EnhG1",
                                    "8_EnhG2","9_EnhA1","10_EnhA2","11_EnhWk","12_ZNF-Rep","13_Het","14_TssBiv","15_EnhBiv",
                                    "16_ReprPc","17_ReprPcWk","18_Quies")))

a$set <- factor(a$set, levels=c("ESC","iPSC","ES-deriv","IMR90","Epithelial","HSC & B-cell",
                                "Blood & T-cell","Mesench","Myosat","Neurosph","Adipose","Heart",
                                "Other","Brain","Digestive","Sm. Muscle","Muscle","Thymus"))

## Remove small-sample bias:
a$fold[a$snp_count==1] <- 1.0
a$fold[a$snp_count==2] <- 1.0

#########
## Scatterplots
## 14 states, all EIDs
g <- filter(a, state!="15_Quies") %>%
  ggplot(., aes(x=ex,y=snp_count))+ 
  geom_point(aes(color=state),size=3,alpha=0.8) + geom_abline(intercept=0, slope=1) + ylim(0,200) + xlim(0,200)

# aes(color=set),

g <- filter(a, state=="9_EnhA1") %>%
  ggplot(., aes(x=ex,y=snp_count,color=set))+ 
  geom_point(size=3) + geom_abline(intercept=0, slope=1) #+ ylim(0,60) + xlim(0,60)

g <- g + ylab("# NCVs Observed") + xlab("# NCVs Expected") +
  theme(axis.text.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

g

png("liver_keep-ExomSeq_hg19_ObsExp_10_BivTss.png",width=2500,height=2000,res=200)
g
dev.off();


###########
## Plot fold change
g <- filter(a,state=="9_EnhA1") %>% 
  ggplot(., aes(x=reorder(set,fold),y=fold,group=set)) + 
  geom_boxplot(size=1.5) + geom_hline(y=1,color="red") 

g <- ggplot(a, aes(y=reorder(set,fold),x=fold,group=set)) + 
  geom_boxplot(size=1.5) + geom_hline(y=1,color="red") + coord_flip() + facet_wrap(~state, scales = "free_x")

g <- g + ylab("Observed/Expected") + xlab("") +
  theme(axis.text.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20)) +
  coord_flip()

g

g <- ggplot(a, aes(x=state, y=fold)) + 
  geom_point(position = position_jitter(width = 0.5)) + coord_flip() + 
  geom_hline(yintercept = 1, col="red") + ylab("Observed/Expected") + xlab("") +
  theme(axis.text.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))

png("liver_keep_WGS_ChromHMM-18_fold_noColor.png",width=1500,height=1500,res=200)
g
dev.off();
