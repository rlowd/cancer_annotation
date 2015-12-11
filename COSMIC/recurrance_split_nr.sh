#!/bin/bash

#######
## Parse filtered NCV set into bedfiles for each recurrence value
## Sort and count unique genome positions in the format chr:start-end
## awk positions for each level of recurrence into separate files
## bedSort, then cat files together to get a single non-redundant list.

#cut -f 1-2 CosmNCV_2k_hg19.bed | sort - | uniq -c - | awk '{print $1"\t"$2"\t"$3}' - | sort -nr -k1 - > d
#
##for i in 1 2 3 4 5 8 10 12 13; 
#for i in 16; 
#do
#	awk -v var="$i" 'OFS="\t" { if($1==var) print $2,$3,$3+1 }' d > ${i}xnr
#	bedSort ${i}xnr CosmNCV_v74_2k_${i}xnr.bed
#done

#rm d
#rm *xnr

#cat CosmNCV_v74_2k_1xnr.bed CosmNCV_v74_2k_2xnr.bed CosmNCV_v74_2k_3xnr.bed \
#CosmNCV_v74_2k_4xnr.bed CosmNCV_v74_2k_5xnr.bed CosmNCV_v74_2k_8xnr.bed CosmNCV_v74_2k_10xnr.bed \
#CosmNCV_v74_2k_12xnr.bed CosmNCV_v74_2k_16xnr.bed > CosmNCV_v74_2k_allxnr.bed


#######
## Analysis of recurrantally mutated positions.
## What sample types are recurrent NCVs found in?

## First grab list of sampleIDs with mutated position of interest
## Then grab sample metadata.

#for p in `cat 2k_3-16xnr.bed`; 
#for p in `cat 16xnr`;
#do
#	echo $p
#	awk '{print $1":"$2"-"$3"\t"$5}' CosmNCV_2k_hg19.bed | grep $p - | cut -f 2 -> ${p}_ids
#	
#	for s in `cat ${p}_ids`;
#	do
#		echo $s
#		zcat v75/CosmicSample.tsv.gz | grep $s - > ${p}_${s}_met
#	done
#done

## Concatenate *metadat files for each variant. e.g.:

#cat chr14\:23887231-23887232_*metadat > chr14_cardiac_miRNA
#cat chr17\:43587569-43587570_*met > chr17_AluSp
#cat chr17\:43587576-43587577*met > chr17_4x_AluSp

for i in `cut -f 4 tert_vars`; 
do 
	echo $i
	grep "$i" CosmicNCV.tsv >> tert_vars_hg38
done