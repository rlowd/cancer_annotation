setwd("~/cancer_annotation/counts_per_element/v77")
library(GenomicFeatures)
library(GenomicRanges)
library(data.table)
library(org.Hs.eg.db)


# get_ensg <- function( g ){
#   ensg <- unlist(strsplit(x=g, split = ".", fixed = T))[1]
#   return( ensg )
# }


# Gencode Basic TxDb

# gencode <- loadDb(file =  "~/genomes/hg19/Gencodev24/gencode_basic_v24.sqlite")
# gencode.promoters <- promoters(transcripts(gencode, filter = "gene_id"), upstream = 2000,downstream = 500 )
# seqlevels(gencode.promoters, force = T) <- seqlevels(gencode.promoters)[1:22]
# gencode.promoters.nr <- unique(gencode.promoters)
# num_promoters_tested <- length(gencode.promoters.nr)
# 
# 
# # Determine background
# 
# ens_promoters <- data.frame(t(data.frame( lapply(X = mcols(gencode.promoters.nr)$tx_name, FUN = get_ensg) )))
# colnames(ens_promoters) <- "ensg"
# ens_prom_gene_names <- data.table(AnnotationDbi::select(x = org.Hs.eg.db,
#                                   keys = as.character(ens_promoters$ensg), columns = c("SYMBOL"), keytype = "ENSEMBLTRANS"))
# ens_prom_gene_names <- dplyr::filter(ens_prom_gene_names, !is.na(SYMBOL))
# ens_prom_gene_names <- dplyr::filter(ens_prom_gene_names, !grepl("LOC", x = SYMBOL))  # Get rid of pseudogenes
# ens_prom_gene_names <- dplyr::filter(ens_prom_gene_names, !grepl("LINC", x = SYMBOL)) 


# Load UCSC Known genes track

# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txdb <- makeTxDbFromUCSC()
txdb <- loadDb(file = "~/genomes/hg19/annotation/ucscKnownGenes_txdb.sqlite")
ucsc.promoters <- promoters(txdb, upstream = 2000, downstream = 500)
seqlevels(ucsc.promoters, force = T) <- seqlevels(ucsc.promoters)[1:22]
ucsc.promoters.nr <- unique(ucsc.promoters)

prom_gene_names <- data.table(AnnotationDbi::select(x = org.Hs.eg.db,
                              keys = mcols(ucsc.promoters.nr)$tx_name, columns = c("SYMBOL"), keytype = "UCSCKG"))
prom_gene_names <- dplyr::filter(prom_gene_names, !is.na(SYMBOL))
prom_gene_names <- dplyr::filter(prom_gene_names, !grepl("LOC", x = SYMBOL))  # Get rid of pseudogenes
prom_gene_names <- dplyr::filter(prom_gene_names, !grepl("LINC", x = SYMBOL)) 
num_promoters_tested <- length(prom_gene_names$UCSCKG)

# Read in  SNVs

df <- read.table(file = "~/cancer_annotation/COSMIC/WGSv77/liver_keep_ExSeq_hg19_v77_99.bed",header=F,sep="\t")
# df <- read.table(file = "~/cancer_annotation/COSMIC/WGSv77/liver_keep_WGS_hg19_v77_975.bed",header=F,sep="\t")
colnames(df) <- c("chr","start","end","cosmicID","sampleID","ExSeq","FATHMM_score")
SNVs <- makeGRangesFromDataFrame(df, keep.extra.columns = T, ignore.strand = T)
seqlevels(SNVs, force = T) <- seqlevels(SNVs)[1:22]

# Select duplicated SNV positions

# dup3x <- makeGRangesFromDataFrame( 
#   read.delim("dup_3x",header = F,sep = "\t",colClasses = c("NULL","factor","integer","integer"),
#              col.names = c("num","chr","start","end")) )
# dup2x <- makeGRangesFromDataFrame( read.delim("dup",header = F,sep = "\t",col.names = c("chr","start","end") ) )
# dup2x <- makeGRangesFromDataFrame( read.delim("dup_WGS_2x",header = F,sep = "\t",
#                                    col.names = c("chr","start","end") ) )
# duplicated_SNVs <- subsetByOverlaps(query = SNVs, subject = dup3x)


# Assign nc-RE SNVs to Gencode ENSEMBL IDs

tmp <- as.data.frame( mergeByOverlaps( SNVs,ucsc.promoters.nr,ignore.strand=T ) )
noncoding_with_gene_promoters_df <- dplyr::select(tmp, -cosmicID,-sampleID,-ExSeq,-FATHMM_score, -tx_name, -tx_id)
# ensg_ids <- data.frame(t(data.frame( lapply(X = noncoding_with_gene_promoters_df$gencode.promoters.nr.tx_name, FUN = get_ensg) )))
# colnames(ensg_ids) <- "ensg"

noncod_gene_names <- data.table(AnnotationDbi::select(x = org.Hs.eg.db,
                                keys = as.character(noncoding_with_gene_promoters_df$ucsc.promoters.nr.tx_name),
                                columns = c("SYMBOL"), keytype = "UCSCKG"))
noncod_gene_names <- dplyr::filter(noncod_gene_names, !is.na(SYMBOL))
noncod_gene_names <- dplyr::filter(noncod_gene_names, !grepl("LOC", x = SYMBOL))
noncod_gene_names <- dplyr::filter(noncod_gene_names, !grepl("LINC", x = SYMBOL))


# Check that all pathway SYMBOLS are in the promoter search list

annot_gene_sets <- read.table("~/cancer_annotation/TCGA_protected/VEP/genesets",header=F)

check_pathway_members <- function( fl,x ){
  tab <- matrix( nrow = length(fl$V1), ncol = 3)
  for( i in 1:length(fl$V1)){
    pref <- read.table( as.character( fl$V1[i] ), nrows = 1)
    name <- pref$V1
    d <- read.table( as.character(fl$V1[i]), skip = 2, col.names="SYMBOL")
    overlap <- data.frame(unique(merge(x = x, y = d, by = "SYMBOL")))
    tab[i,1] <- paste(name)
    tab[i,2] <- length(unique(overlap$SYMBOL))
    tab[i,3] <- length(d$SYMBOL)
  }
  return( tab )
}

pathway_members_hit <- data.frame( check_pathway_members( annot_gene_sets, prom_gene_names ))
colnames(pathway_members_hit) <- c("Pathway","numGenesHit","numGenesTotal")

# repf <- "pathwayGeneMemebersInUcscKnownGenes.tsv"
# if(!file.exists(repf)) { file.create(repf) }
# dat <- pathway_members_hit
# write.table(dat, file = repf, quote = F,row.names = F,col.names = T,sep = "\t")

# Get pathway hits

count_gene_set_overlap <- function( fl,x ){
  tab <- matrix( nrow = length(fl$V1), ncol = 3)
  for( i in 1:length(fl$V1)){
    pref <- read.table( as.character( fl$V1[i] ), nrows = 1)
    name <- pref$V1
    d <- read.table( as.character(fl$V1[i]), skip = 2, col.names="SYMBOL")
    overlap <- data.frame(unique(merge(x = x, y = d, by = "SYMBOL")))
    tab[i,1] <- paste(name)
    tab[i,2] <- length(overlap$SYMBOL)
    tab[i,3] <- length(d$SYMBOL)
  }
  return( tab )
}

noncod_counts_df <- data.frame( count_gene_set_overlap( annot_gene_sets, noncod_gene_names ))
colnames(noncod_counts_df) <- c("Pathway","numNcvs","numGenes")

noncod_counts_df$numGenesTested <- pathway_members_hit$numGenesHit

# Binomial test

# prom_t <- length((ens_prom_gene_names$SYMBOL))  # total number of all promoter TRANSCRIPTS tested
num_promoters_tested <- length(prom_gene_names$UCSCKG)  # total number of unique PROMOTERs in sample (N)
np <- length((noncod_gene_names$UCSCKG))       # total number trials (n)

noncod_counts_df$numNcvs <- as.numeric(paste(noncod_counts_df$numNcvs))
noncod_counts_df$numGenesTested <- as.numeric(paste(noncod_counts_df$numGenesTested))

noncod_counts_df$p_exp <- noncod_counts_df$numGenesTested/num_promoters_tested  ## expected ratio of sucess
# noncod_counts_df$obs <-  noncod_counts_df$numNcvs/np        ## observed ratio of success

for( i in 1:length(noncod_counts_df$Pathway )){
  x <- noncod_counts_df$numNcvs[i]
  n <- np
  p <- noncod_counts_df$p_exp[i]
  res <- binom.test(x=x,n=n,p = p, alternative = "greater")
  noncod_counts_df$binom_pval[i] <- res$p.value
}

###########
# Write report

repf <- "liverExSeq_v77_curatedPathways.correctedExpectation.UCSC.report"
if(!file.exists(repf)) { file.create(repf) }

dat <- num_promoters_tested
mes <- paste("Total number promoters (by tx) tested: N = ",dat,sep = "")
write.table(mes, file = repf, quote = F,row.names = F,col.names = F)

dat <- np
mes <- paste("Total number promoters hit: np = ",dat,sep = "")
write.table(mes, file = repf, quote = F,row.names = F,col.names = F,append = T)

mes <- paste("\nObserved hits by Pathway, with uncorrected p-value:\n",sep = "")
write.table(mes, file = repf, quote = F,row.names = F,col.names = F,append = T)

dat <- noncod_counts_df
write.table(dat, file = repf, quote = F,row.names = F,col.names = T,append = T,sep = "\t")


##########
## Find specific group overlaps
library(dplyr)

d <- read.table( as.character(annot_gene_sets$V1[20]), skip = 2, col.names="SYMBOL")
pref <- read.table( as.character( annot_gene_sets$V1[20] ), nrows = 1)
overlap <- data.frame(unique(merge(x = noncod_gene_names, y = d, by = "SYMBOL")))
head(overlap)
h <- data.frame(table(overlap$SYMBOL))

hits <- data.table(AnnotationDbi::select(x = org.Hs.eg.db,
                                                      keys = as.character(overlap$SYMBOL),
                                                      columns = c("ENTREZID"), keytype = "SYMBOL"))
head(hits)
rank <- group_by(hits, SYMBOL) %>% summarise(., n())
write.table(x = unique(hits$ENTREZID), file = "hits",quote = F,sep = "\t",row.names = F, col.names = F)
