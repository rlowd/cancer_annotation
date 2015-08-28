library(dplyr)

gwas <- read.delim2("gwas_catalog_v1.0-downloaded_2015-08-18.tsv",header=TRUE)

## Clean data
## Create columns for cleaned genomic positions
## If genome position missing, go to SNPS column and take chr:pos information to make
## cleaned genomic position

cleanDat <- function(g) {
    
    g$CLEAN_CHR <- seq(0,0)
    g$CLEAN_POS <- seq(0,0)
    g$CLEAN_rs <- seq(0,0)
    
    for( i in 1:length(g$CHR_POS) ) {
        if( is.na(g$CHR_POS[i]) ) {
            a <- noquote(unlist(strsplit(as.character(g$SNPS[i]),split="(hr[0-9]+:)|(hrX:)")))
            if( !is.na(a[2]) ) {
                a2 <- noquote(unlist(strsplit(as.character(g$SNPS[i]),split=":")))
                b <- noquote(unlist(strsplit(a2,split="r",fixed=TRUE)))
                g$CLEAN_CHR[i]  <- paste("chr",b[2],sep="")
                
                c <- noquote(unlist(strsplit(a2[2],split="-",fixed=TRUE)))
                if( !is.na(c[2]) ){
                    g$CLEAN_POS[i]  <- as.integer(c[1])
                } else {
                    g$CLEAN_POS[i]  <- as.integer(a[2])
                }
                g$CLEAN_rs[i] <- as.character(g$SNPS[i])
            }
            else if( !is.na(a[1]) ) {
                g$CLEAN_rs[i] <- a[1]
            } 
        }
        else {
            g$CLEAN_CHR[i] <- paste("chr",g$CHR_ID[i],sep="")
            g$CLEAN_POS[i] <- g$CHR_POS[i]
            g$CLEAN_rs[i] <- as.character(g$SNPS[i])
        }
    } 
    n <- length(g)
    m <- n-2
    g[,m:n][g[,m:n]==0] <- NA   ## Make missing values for CLEAN columns NA
    return(g)
}

## Function tidyDat tidys the data and returns 4 dataframes:
## The first has all cleaned data, distinct entries.
## g_chr filters the tidy data for entries with genome position information.
## g_rs filteres tidy data for entries with SNP accessions but NO genome position
## g_na is the remainder

tidyDat <- function(g) {
    g_tidy <- mutate(g, CLEAN_END = CLEAN_POS+1) %>%
        select(CLEAN_CHR,CLEAN_POS,CLEAN_END,CLEAN_rs,DISEASE.TRAIT,CONTEXT,INTERGENIC,PVALUE_MLOG) %>%
        distinct(CLEAN_rs)
    g_chr <- filter(g_tidy, CLEAN_CHR>0)
    g_rs <- filter(g_tidy, !is.na(CLEAN_rs), is.na(CLEAN_POS))
    g_na <- filter(g_tidy, is.na(CLEAN_rs), is.na(CLEAN_POS))
    output <- list(g_tidy,g_chr,g_rs,g_na)
    return(output)
}

clean <- cleanDat(gwas)
tidy <- tidyDat(clean)

#tid <- data.frame(tidy[1])
#chr <- data.frame(tidy[2])
#rs <- data.frame(tidy[3])
#na <- data.frame(tidy[4])

write.table(data.frame(tidy[1]),file = "tidyGWAS_08272015.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(data.frame(tidy[2]),file = "chrGWAS_08272015.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(data.frame(tidy[3]),file = "rsGWAS_08272015.bed",sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)


## bedSort chrGWAS.bed

system2("bedSort", args=c("chrGWAS_08272015.bed","chrGWAS_08272015.bed"))
