####################################################################################
#### NEXT: USE txdb object  to annotate distance to nearest genes for all data! ####
####################################################################################

## Load additional packages
#rm(list=ls())
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)
#library(ChIPpeakAnno)

## Read in test file: 
FitCons80_130 <- fread('FitCons_80-130bp_all-annotations_plusK562-tags.txt')

## add sudo count to starts to make annotations work better
FitCons80_130 <- FitCons80_130[ , start := start +1, ]
  head(FitCons80_130)

fitCons_plus_Annos <- makeGRangesFromDataFrame(FitCons80_130, keep.extra.columns = TRUE)
  head(fitCons_plus_Annos) # has all the columns!

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene #shorthand (for convenience)
txdb #hg19; Entrez Gene ID
# seqlevels(txdb) #got a lot of chromosomes here

# get all the genes
Genes_hg19 <- genes(txdb)
length(Genes_hg19) #22304 genes from same chromosomes (still true 4.9.16)

Exons_hg19 <- exons(txdb) # just want exons
length(Exons_hg19) # 267061 (all exons)
Exons_hg19 

## Next, want to expand exons by 10 bp to also exclude any splice junctions
Exons_10plus <- Exons_hg19
start(Exons_hg19) ## was list of list but now just a list
start(Exons_10plus) <- start(Exons_hg19) - 10
end(Exons_hg19)
end(Exons_10plus) <- end(Exons_hg19) + 10
width(Exons_10plus)
head(width(Exons_10plus) - width(Exons_hg19)) # 20! Intervals are expanded!

Introns_hg19 <- intronsByTranscript(txdb)
length(Introns_hg19) # 75333, by transcripts

FitCons80_130 <- makeGRangesFromDataFrame(fitCons_Annos_nearest_genes, keep.extra.columns = TRUE)

## check/get rid of any remaining exon overlapping regions
Exon_overlaps <- GenomicRanges::overlapsAny(FitCons80_130, Exons_hg19, type = "any", ignore.strand = TRUE) # will output full GRanges 
Exon_overlaps[1]
  length(which(Exon_overlaps == TRUE)) # 26612 (FUCKERS)

Exclude_exons <- FitCons80_130[which(Exon_overlaps == FALSE), ]
  length(FitCons80_130) # 840311
  length(Exclude_exons) # 813699
  length(FitCons80_130) - length(Exclude_exons) # 26612 MATCHES
  
FitCons80_130 <- Exclude_exons  ## recode table with exons TOTALLY removed
head(FitCons80_130)

Intron_overlaps <- GenomicRanges::overlapsAny(FitCons80_130, Introns_hg19, type = 'within') ## what's distribution of inside/up-downstream for these guys?
length(which(Intron_overlaps == TRUE)) # 524846 (ALOT)
length(which(Intron_overlaps == TRUE))/length(FitCons80_130)*100 #64.5% of guys overlap introns for 80-130bp!

## Add column for intron overlap for future reference
fitCons_Annos_nearest_genes <- as.data.frame(FitCons80_130)
fitCons_Annos_nearest_genes <- data.table(fitCons_Annos_nearest_genes)
  head(fitCons_Annos_nearest_genes)

## Add field for intronic overlaps (Boolean)
fitCons_Annos_nearest_genes <- fitCons_Annos_nearest_genes[ , "intronicFeature" := 0,  ] # add column to conditionally format
fitCons_Annos_nearest_test <- fitCons_Annos_nearest_genes[1:100, ]
overlap_test <- Intron_overlaps[1:100]
length(which(overlap_test == TRUE)) # 21 should be true out of 100 in test case
fitCons_Annos_nearest_test <- fitCons_Annos_nearest_test[which(overlap_test == TRUE), "intronicFeature" := 1,  ]
  head(fitCons_Annos_nearest_test) # works!! repeat for full dataset

fitCons_Annos_nearest_genes <- fitCons_Annos_nearest_genes[which(Intron_overlaps == TRUE), "intronicFeature" := 1,  ]
  length(which(fitCons_Annos_nearest_genes$intronicFeature == 1)) #524846 (MATCHES; WORKED!)


table(fitCons_Annos_nearest_genes$insideFeature)
### getting included Features, need to exclude remaining exons!
# downstream     inside   upstream 
# 3464162    4543801    3411869 

write.table(fitCons_Annos_nearest_genes, file = "~/Project3_FIles-Geuvadis/FitCons_ALL_nearest-genes_both-Annos_K562-tags_plus-intronicfeatures.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE)

# pie(table(fitCons_Annos_nearest_genes$insideFeature), main = "fitCon blocks to nearest genes", sub = "~1.4 million blocks\nUCSC hg19 genes, Updated")  # roughly matches pie chart for when using 'TSS.human.GRCh37' 
# # save as: "PieChart_fitCon_80-130bp_nearest-genes_UPDATED"