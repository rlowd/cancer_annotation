############
# Plot pathway p-values
setwd("~/cancer_annotation/counts_per_element/v77")
wgs <- read.delim("liverWGS_v77_curatedPathways.correctedExpectation.UCSC.report",header = T,sep = "\t",skip = 5)
exseq <- read.delim("liverExSeq_v77_curatedPathways.correctedExpectation.UCSC.report",header = T,sep = "\t",skip = 5)

dat <- cbind(wgs$binom_pval,exseq$binom_pval)
colnames(dat) <- c("WGS","ExomeSeq")

paths <- paste(wgs$Pathway)
rownames(dat) <- paths
head(dat)

library(ggplot2)
library(reshape2)
library(dplyr)

m <- melt(dat)
colnames(m) <- c("pathway","SNV_set","value")

head(m)

labs <- c("ERBB signaling","Apoptosis","Transcriptional misregulation in cancer","Cell cycle",
          "JAK/STAT signaling pathway","MAPK signaling pathway","MTOR signaling pathway","Chormating modification",
          "Histone modification","Mismatch repair","WNT singaling pathway","Nucleotide excision repair",
          "TGF-beta signaling pathway","Chromatin remodeling",
          "Nuclear receptor transcription pathway","Histone deacetylase complex","Histone methyltransferase activity",
          "P53 singaling pathway","RAS protein signal transduction","DNA replication","Base excision repair",
          "Hedgehog signaling pathway","Epigenetic regulation of gene expression", "Notch signaling pathway",
          "Histone deacetylase binding")

g <- ggplot(m, aes(x = SNV_set, y=reorder(pathway,-value), fill=-log10(value))) + geom_tile(color="white") + 
  scale_fill_continuous(name = "-log10(P-value)",low = "gray",high = "red") +
  theme_classic(base_size = 20) + ylab("") + xlab("") + scale_y_discrete(labels=rev(labs)) +
  theme(axis.text.x=element_text(angle=90, vjust = 0.5), legend.position="left")

pdf("combined_pathway_binom_pvals.pdf",width = 8.2,height = 7.5)
g
dev.off();
