#!/bin/bash


## What are the sample IDs and their metadata with the NCVs of interest?

echo 20:44401510  > ncvs
echo 12:114401777 >> ncvs

for p in `cat ncvs`;
do
	echo $p
	zcat CosmicWGS_NCV.tsv.gz | grep "$p" - >> ncvs_samples
done

for s in `cut -f 2 ncvs_samples`;
do
	echo $s
	zcat CosmicWGS_SamplesExport.tsv.gz | grep "$s" - | awk '{print $1"\t"$3"\t"$7}' - >> ncvs_samples_metadat
done


## Do these samples have expression information?

#zcat CosmicCompleteGeneExpression.tsv.gz | cut -f 1 | sort - | uniq - > CCGE_samples
#
#for s in `cut -f 1 ncvs_samples_metadat`;
#do
#	echo $s
#	grep "$s" CCGE_samples >> ncvs_CCGE
#done
#
### Get expression data of interest
#
#for s in `cut -f 1 ncvs_samples_metadat`;
#do
#	echo $s
#	zcat CosmicCompleteGeneExpression.tsv.gz | grep "$s" - | awk '{if($3=="FGF5") print $0}' >> ncvs_target_expression
#done