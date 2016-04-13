#!/bin/bash

## Data pre-processing for differential motif analysis
## First get NCVs in states of interest (for keep and removed NCVs) and liftOver to hg38.
## Then form grep-able input, first awk command.
## Next grep the CosmicWGS_NCV.tsv.gz to get the NCVs of interest and their WT and mutated alleles. 
## Print to outfile that will be used in python script
## for creating reference and mutated FASTA sequences.

ln -s ~/cancer_annotation/enrichment/union_byChromHMM-18/allEIDs_14_TssBiv.bed_uniq .
ln -s ~/cancer_annotation/enrichment/union_byChromHMM-18/allEIDs_1_TssA.bed_uniq .

f="allEIDs_1_TssA.bed_uniq"

liftOver $f ~/genomes/hg19/hg19ToHg38.over.chain.gz ${f}_hg38 unm

awk 'BEGIN { FS="r";} {print $2}' ${f}_hg38 | awk 'BEGIN { FS="\t"; } { print $1":"$2"-"$2 }' - > ${f}_hg38-list

for p in `cat ${f}_hg38-list`;
do
	zcat ~/cancer_annotation/COSMIC/WGSv75/CosmicWGS_NCV.tsv.gz | grep "$p" - | awk 'BEGIN {FS="\t";} { print $6"\t"$8"\t"$9 }' - >> ${f}_hg38.alleles
done

sort ${f}_hg38.alleles | uniq - > ${f}_hg38.uniqAlleles

python mutateFasta.py ${f}_hg38.uniqAlleles


