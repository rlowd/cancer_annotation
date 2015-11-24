################
## 
## Rebecca L | Modified 29 Oct 2015
##
## Code to clean GWAS data
## Usage: /apps/bin/Rscript cleanGwas.R 
##
################

library(dplyr)

classes <- c( "NULL","integer",rep("NULL",5),"factor",rep("NULL",3),"integer","integer",rep("NULL",8),
             "factor",rep("NULL",2),"factor","integer",
             rep("NULL",2),"factor",rep("NULL",5) )

gwas <- read.delim2("gwas_catalog_v1.0-downloaded_2015-08-18.tsv",header=TRUE, colClasses = classes)

#####
## Clean data
## Create columns for cleaned genomic positions
## If genome position missing, go to SNPS column and take chr:pos information to make
## cleaned genomic position

cleanDat <- function(g) {
    
    g$CLEAN_CHR <- seq(0,0)
    g$CLEAN_POS <- seq(0,0)
    g$CLEAN_rs <- seq(0,0)
    
    for( i in 1:length(g$CHR_POS) ) {
        if( is.na(g$CHR_POS[i]) ) {     ## if there is no CHR_POS value, check if theres is a chr:pos value
            a <- noquote(unlist(strsplit(as.character(g$SNPS[i]),split="(hr[0-9]+:)|(hrX:)")))
            if( !is.na(a[2]) ) {
                a2 <- noquote(unlist(strsplit(as.character(g$SNPS[i]),split=":")))  ## split chr# and pos
                b <- noquote(unlist(strsplit(a2,split="r",fixed=TRUE)))             ## sep chr from #
                g$CLEAN_CHR[i]  <- paste("chr",b[2],sep="")                         ## set chrNUM
                
                c <- noquote(unlist(strsplit(a2[2],split="-",fixed=TRUE)))          ## parse and set POS
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
        else {    ## If values exist in CHR_ID, CHR_POS, and CHR_rs fields use them colums
            if( g$CHR_ID[i] == "23" ){    ## check if "23" in CHR_ID field, if yes replace with 'chrX'
                g$CLEAN_CHR[i] <- paste("chrX")
            }
            else {
              g$CLEAN_CHR[i] <- paste("chr",g$CHR_ID[i],sep="")
            }
          
            g$CLEAN_POS[i] <- g$CHR_POS[i]
            g$CLEAN_rs[i] <- as.character(g$SNPS[i])
        }
    } 
    n <- length(g)
    m <- n-2
    g[,m:n][g[,m:n]==0] <- NA   ## Make missing values for CLEAN columns NA
    return(g)
}

#####
## Function tidyDat tidys the data and returns 4 dataframes:
## The first has all cleaned data, distinct entries.
## g_chr filters the tidy data for entries with genome position information.
## g_rs filteres tidy data for entries with SNP accessions but NO genome position
## g_na is the remainder

tidyDat <- function(g) {
    g_tidy <- mutate(g, CLEAN_END = CLEAN_POS+1) %>%
        select(CLEAN_CHR,CLEAN_POS,CLEAN_END,CLEAN_rs,PUBMEDID,INTERGENIC,PVALUE_MLOG,CONTEXT,DISEASE.TRAIT) %>%
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

##### 
## Output files will be date-stamped

date <- Sys.Date()
# format(date, format="%Y-%m-%d")

t <- paste("tidyGWAS_",date,".bed", sep="")
c <- paste("chrGWAS_",date,".bed", sep="")
r <- paste("rsGWAS_",date,".bed", sep="")

write.table(data.frame(tidy[1]),file = t ,sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(data.frame(tidy[2]),file = c ,sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(data.frame(tidy[3]),file = r ,sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)

## bedSort chrGWAS.bed
system2("bedSort", args=c(c,c))
