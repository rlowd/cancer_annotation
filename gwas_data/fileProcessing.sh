#!/bin/bash

###########
##
##  File proccesing for generating chrGWAS_<date>_canc.bed file
##  Bash commands will take cleanGwas.R output (arg $1) and filter for cancer-associated
##  SNPs and report duplicate variants.
##  Then litfOver tool will convert to hg19 coordinates.
##
##  Dev: Aug 2015
##  R Lowdon
##  http://wang.wustl.edu/mediawiki/index.php/RebeccaLowdon_August_2015#Filter_GWAS_SNPs_for_cancer_phenotypes
##
###########

## 1. Grep to select cancer-associated traits:

d=$(date +%Y%m%d)
out=${1%.bed}_canc.bed
hg19=${out%.bed}_hg19.bed

#echo "grep cancer-traits to"
echo $out

grep -i 'cancer' $1 > $out
grep -i 'leukemia' $1 >> $out
grep -i 'lymphoma' $1 >> $out
grep -i 'sarcoma' $1 >> $out
grep -i 'carcinoma' $1 >> $out

## 2. Sort and use uniq command to find duplicate entries:

sort $out > tmp
uniq -c -d tmp dupl.$d
unlink tmp

## 3. Use liftOver tool to convert to hg38

liftOver -bedPlus=5 $out hg38ToHg19.over.chain.gz $hg19 unmapped.$d
echo "liftOver -bedPlus=5 $out hg38ToHg19.over.chain.gz $hg19 unmapped.$d"

## Then remove fields 9+
mv $hg19 tmp
cut -f 1-8 tmp > $hg19
unlink tmp