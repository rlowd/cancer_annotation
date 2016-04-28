#!/bin/bash

## Data pre-processing for differential motif analysis
## First get NCVs in states of interest (for keep and removed NCVs) and liftOver to hg38.
## Then form grep-able input, first awk command.
## Next grep the CosmicWGS_NCV.tsv.gz to get the NCVs of interest and their WT and mutated alleles. 
## Print to outfile that will be used in python script
## for creating reference and mutated FASTA sequences.

#ln -s ~/cancer_annotation/enrichment/union_byChromHMM-18/allEIDs_14_TssBiv.bed_uniq .
#ln -s ~/cancer_annotation/enrichment/union_byChromHMM-18/allEIDs_1_TssA.bed_uniq .
#ln -s ~/cancer_annotation/enrichment/union_byChromHMM-18/RM78-EIDs_ExSeq* .
#
#f="RM78-EIDs_ExSeq_1_TssA.bed_uniq"


## For loop to iterate over each file -- get alleles

#ls RM78-EIDs_ExSeq* > rmlist

#for f in `cat rmlist.4`;
#do
#	echo $f
#	liftOver $f ~/genomes/hg19/hg19ToHg38.over.chain.gz ${f}_hg38 unm
#	
#	awk 'BEGIN { FS="r";} {print $2}' ${f}_hg38 | awk 'BEGIN { FS="\t"; } { print $1":"$2"-"$2 }' - > ${f}_hg38-list
#	
#	for p in `cat ${f}_hg38-list`;
#	do
#		zcat ~/cancer_annotation/COSMIC/WGSv75/CosmicWGS_NCV.tsv.gz | grep "$p" - | awk 'BEGIN {FS="\t";} { print $6"\t"$8"\t"$9 }' - >> ${f}_hg38.alleles
#	done
#	
#	sort ${f}_hg38.alleles | uniq - > ${f}_hg38.uniqAlleles
#	
#	suff=".bed_uniq"
#	newdir=${f%$suff}
#	mkdir $newdir
#	mv ${f}_hg38.uniqAlleles ${newdir}/.
#	
#	## Cleanup files
#	rm ${f}_hg38
#	rm ${f}_hg38-list
#	rm ${f}_hg38.alleles
#	mv $f ${newdir}/.
#done

## Loop to run mutateFasta.py

for f in `cat rmlist.4`;
do
	suff=".bed_uniq"
	newdir=${f%$suff}
	
	echo "running mutateFasta.py script on "${newdir}/${f}_hg38.uniqAlleles

	python mutateFasta.py ${newdir}/${f}_hg38.uniqAlleles
	
	## Cleanup files from python script
	rm unm
	rm out
	rm inbed
	rm alleles
	mv fastaWithAlleles ${newdir}/.
done