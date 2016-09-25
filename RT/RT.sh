#!/bin/bash

## June 6, 2016
## Same strategy as naive analysis.
## Using partitioned RT bins

#tail -n+17 RT_CyT49_Liver_D16_All.txt | cut -f 3-6 - > liver_D16_R1.bed
#tail -n+17 RT_CyT49_Liver_D16_All-2.txt | cut -f 3-6 - > liver_D16_R2.bed
#tail -n+17 RT_CyT49_Definitive\ Endoderm_All.txt | cut -f 3-6 - > DE.bed


## Instead of using only probe intervals, use python script to partition into Mb-size regions
#python RT_partitions.py RT_CyT49_Liver_D16_All.txt > liver_D16_all.partitions.txt
#python RT_partitions.py RT_CyT49_Liver_D16_All-2.txt > liver_D16_all-2.partitions.txt


for f in `ls *hg38.bed`;
do
#	intersectBed -c -a liver_D16_all.partitions.txt -b $f > ${f}_vs_liver_R1
#	intersectBed -wo -a $f -b liver_D16_all.partitions.txt> ${f}_vs_liver_R1
	intersectBed -c -a liver_D16_all-2.partitions.txt -b $f > ${f}_vs_liver_R2
#	intersectBed -wo -a $f -b DE.bed > ${f}_vs_DE
done

#for f in `ls *_vs_*`;
#do
#	cut -f 11 $f | sort -nr - | awk '{print NR"\t"$1}' -  > ${f}.valList
#done
