setwd("~/cancer_annotation/diffMotifs")

library(dplyr)
library(ggplot2)
library(reshape2)

d <- read.table("1TssA_motifLO_byPosition.txt",header=F)
d[,4:9] <- as.double(cbind(d[,4],d[,5],d[,6],d[,7],d[,8],d[,9]))
colnames(d) <- c("ID","TF","pos","WT_fwd_lo","MUT_fwd_lo","delta_fwd_lo",
                 "WT_rev_lo","MUT_rev_lo","delta_rev_lo","threshold")

# m <- melt(d, id.vars = c("ID","pos","TF","threshold","delta_rev_lo","delta_fwd_lo"))

fwd <- cbind(d[,1:6],d[,10])
fwd$strand <- "fwd"
colnames(fwd) <- c("ID","TF","pos","WT_lo","MUT_lo","delta","threshold","strand")
rev <- cbind(d[,1:3],d[,7:10])
rev$strand <- "rev"
colnames(rev) <- c("ID","TF","pos","WT_lo","MUT_lo","delta","threshold","strand")

allDat <- rbind(fwd,rev)

ggplot(allDat, aes(delta)) + geom_density(fill = "pink") + aes(y = ..count..) +
  ylab("Count") + xlab("Delta") +
  theme(axis.title.x=element_text(size=24), axis.text.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y=element_text(size=24))

ggplot(allDat, aes(y=threshold, x=TF)) + geom_point() + coord_flip() +
  xlab("Threshold for FPR = 0.001") + ylab("TFs") +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14), axis.text.y=element_text(size=14))

## Gain-of-binding plot
gain <- filter(allDat, (delta> 2)) %>% filter(., MUT_lo>2)

ggplot(gain, aes(x=WT_lo,y=MUT_lo, color=strand)) + geom_point() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14), axis.text.y=element_text(size=14))


## Loss-of-binding plot
loss <- filter(allDat, (delta< -2)) %>% filter(., WT_lo>2)

ggplot(loss, aes(x=WT_lo,y=MUT_lo, color=strand)) + geom_point() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme(axis.title.x=element_blank(), axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=14), axis.text.y=element_text(size=14))


loss$lab <- "loss"
gain$lab <- "gain"
combine <- rbind(loss,gain)

ggplot(combine, aes(delta)) + geom_histogram(bins=40) + 
  geom_vline(xintercept = 1.5, col="red") + geom_vline(xintercept = -1.5, col="red") +
  ylab("Count") + xlab("Delta LOR (MUT - WT)") +
  theme(axis.title.x=element_text(size=24), axis.text.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y=element_text(size=24))

ggplot(combine, aes(x=WT_lo,y=MUT_lo, color=lab,size=1.5)) + geom_point() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  ylab("Mutant allele log odds") + xlab("WT allele log odds") +
  theme(axis.title.x=element_text(size=24), axis.text.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y=element_text(size=24),
        legend.text = element_text(size = 24),legend.title = element_text(size = 24))

g

#################
## Onocogenes
fli1 <- filter(combine, TF=="FLI1")
fli1$lab <- "fli1"
myc <- filter(combine,TF=="Myc")
myc$lab <- "myc"
egr1 <- filter(combine,TF=="EGR1")
egr1$lab <- "egr1"
maxtf <- filter(combine,TF=="MAX")
maxtf$lab <- "max"
etv1 <- filter(combine,TF=="ETV1")
etv1$lab <- "etv1"
nfkb2 <- filter(combine,TF=="NFKB2")
nfkb2$lab <- "nfkb2"
mycmax <- filter(combine,TF=="MAX::MYC")
mycmax$lab <- "mycmax"
runx1 <- filter(combine,TF=="RUNX1")
runx1$lab <- "runx1"
ets1 <- filter(combine,TF=="ETS1")
ets1$lab <- "ets1"

## TS
foxp1 <- filter(combine,TF=="FOXP1")
foxp1$lab <- "foxp1"
ebf1 <- filter(combine,TF=="EBF1")  ## translocations in Hodgkin's lymphoma
ebf1$lab <- "ebf1"
nfkb2 <- filter(combine,TF=="NFKB2")
nfkb2$lab <- "nfkb2"
myod1 <- filter(combine,TF=="Myod1")  ## mutated in rhabdomyosarcoma.
myod1$lab <- "myod1"
foxo3 <- filter(combine,TF=="FOXO3")  ## translocations with MLL in acute secondary lymphoma
foxo3$lab <- "foxo3"
foxa1 <- filter(combine,TF=="FOXA1")  ## FOXO3 translocations with MLL in acute secondary lymphoma
foxa1$lab <- "foxa1"
prdm1 <- filter(combine,TF=="PRDM1")  ## inactivated in Burkitt's lymphoma
prdm1$lab <- "prdm1"
hnf4g <- filter(combine,TF=="HNF4G")  ## inactivated in Burkitt's lymphoma
hnf4g$lab <- "hnf4g"


allDat <- rbind(combine, ets1 )#, maxtf,mycmax )#fli1,myc,egr1,maxtf,etv1,nfkb2,mycmax,runx1,foxp1,ebf1,myod1,foxo3,prdm1)
allDat$lab <- factor(allDat$lab, levels=c("gain","ets1","loss"))#,"max","maycmax"))
                                          "fli1","myc","egr1","max","etv1","nfkb2","mycmax",
                                          "runx1","foxp1","ebf1","myod1","foxo3","prdm1"))

ggplot(allDat, aes(x=WT_lo,y=MUT_lo, color=lab,size=1.5)) + geom_point() +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  ylab("Mutant allele log odds") + xlab("WT allele log odds") +
  theme(axis.title.x=element_text(size=24), axis.text.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y=element_text(size=24),
        legend.text = element_text(size = 24),legend.title = element_text(size = 24))

