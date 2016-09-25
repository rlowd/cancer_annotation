setwd("~/cancer_annotation/COSMIC/WGSv77/")
d <- read.table("pancreas_confSomMut_WGS_SNVs_per_sampleID",header = F, sep = "\t")
  colnames(d) <-  c("numSNVs","sampleID")
  d$sampleID <- factor(d$sampleID)
  glimpse(d)


ggplot(d, aes(x=reorder(sampleID,numSNVs), y=numSNVs)) + geom_bar(stat="identity",fill="black",color="black") +
  theme(axis.text.x=element_blank()) + scale_y_log10() + 
  theme_bw(base_size = 20)


quants <- data.frame(quantile(d$numSNVs, c(seq(0,1,by = 0.005)))) 
colnames(quants) <- "values"
quants$quantile <- rownames(quants)
glimpse(quants)


ggplot(quants, aes(x=reorder(quantile,values),y=values)) + geom_point() + 
  theme(axis.text.x=element_blank()) + theme_bw(base_size = 20)


labs <- names(quantile(d$numSNVs,c(0.97, 0.975,0.98, 0.985, 0.99)))
zoom <- data.frame(quantile(d$numSNVs,c(0.97, 0.975,0.98, 0.985, 0.99)))
colnames(zoom) <- "values"
zoom$quantile <- labs
zoom


slope <- function( y,x ){
  slopes <- data.frame(ncol(2))
  delta_y <- y[2] - y[1]
  m <- delta_y / 0.25  
  slopes[1,1] <- m
  slopes[1,2] <- paste(x[1],x[2],sep = " : ")
  for( i in 3:length(y)){
    delta_y <- y[i] - y[i-1]
    m <- delta_y / 0.25
    slopes[i-1,1] <- m
    slopes[i-1,2] <- paste(x[i-1],x[i],sep = " : ")
  }
  return(slopes)
}

s <- slope(quants$values,quants$quantile)
thresh <- max(s$V1[1:199])
thresh_val <- thresh*0.5
# ggplot(zoom, aes(x=reorder(quantile,values),y=values)) + geom_point() #+ geom_hline(yintercept = 367)


d_sorted <- d[rev(order(d$numSNVs)),]
head(d_sorted)


png("pancreas_quantile_ExSeq.png",height = 200,width = 600)
ggplot(quants, aes(x=reorder(quantile,values),y=values)) + geom_point() + 
  geom_vline(xintercept = 197.0, col="red") + xlab("quantile") + ylab("NumSNVs") +
  theme_classic(base_size = 20) + theme(axis.text.x=element_blank()) 
dev.off();

png("pancreas_samples_by_numSNVs_ExSeq.png",height = 400,width = 600)
ggplot(d, aes(x=reorder(sampleID,numSNVs), y=numSNVs)) + 
  geom_bar(stat="identity",fill="black",color="black") + 
  geom_vline(xintercept = 554, col="red") +
  theme_classic(base_size = 20) + theme(axis.text.x=element_blank()) + scale_y_log10()
dev.off();
