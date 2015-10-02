library(dplyr)
library(ggplot2)
library(reshape2)
# args<-commandArgs(TRUE)

files <- paste("fl_in")  #args[1]
fl <- read.table(files,header=FALSE)

d <- matrix( ncol = 5)
colnames(d)<-c("EID","rs","regulatory","transcribed","silent")
for( i in 1:length(fl$V1) ){
  try(dat <- read.delim2(paste(fl$V1[i]),skip=1,header=FALSE),TRUE)
  try(dat$V12 <- c(i))
  try(dat <- select(dat, 12,4,9:11))
  try(colnames(dat)<-c("EID","rs","regulatory","transcribed","silent"))
  try(d <- rbind(d, dat),TRUE)
}  

d <- filter(d, !is.na(d[,1]))  ## to remove the row of "NA" created when matrix initialized
m <- melt(d,measure.vars = c(3,4,5))

g <- ggplot(m, aes(variable, rs)) + geom_tile(aes(fill=value)) + facet_wrap(~EID)
g <- g + theme(axis.text.y = element_blank()) + 
  scale_x_discrete(name="Regulatory element class") 

png("EID_REprofiles.png")
g
dev.off();


# library(pheatmap)
# library(RColorBrewer)
# 
# m <- rbind(c(1,2),c(3,4))
# layout(m)
# 
# par(mfrow=c(2,2))
# for( i in 1:length(fl$V1) ){
#   dat <- read.delim2(paste(fl$V1[1]),skip=1,header=FALSE)
#   dat <- select(dat,4,9:11)
#   colnames(dat)<-c("rs","regulatory","transcribed","silent")
#   mat <- cbind(as.numeric(paste(dat[,2])),as.numeric(paste(dat[,3])),as.numeric(paste(dat[,4])))
#   p <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
#              labels_col = colnames(dat[,2:4]), legend = F)
# } 



