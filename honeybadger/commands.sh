#!/bin/bash

awk 'BEGIN {FS="-|:"; OFS="\t";} {print $1,$2,$3}' prom-DNase_regions.txt > prom-DNase_regions.bed

intersectBed -wao -a prom-DNase_regions.bed -b /home/comp/twlab/rlowdon/cancer_annotation/COSMIC/WGSv75/liver_keep_L3-WGS_hg19.bed > livWgsNcvs_prom-DNase.bed

awk 'BEGIN {OFS="\t";} {print $1":"$2"-"$3}' livWgsNcvs_prom-DNase.bed > livWgsNcvs_prom-DNase.txt


awk 'BEGIN {FS="-|:"; OFS="\t";} {print $1,$2,$3}' enh-DNase_regions.txt > enh-DNase_regions.bed

intersectBed -wa -a enh-DNase_regions.bed -b /home/comp/twlab/rlowdon/cancer_annotation/COSMIC/WGSv75/liver_keep_L3-WGS_hg19.bed > livWgsNcvs_enh-DNase.bed

awk 'BEGIN {OFS="\t";} {print $1":"$2"-"$3}' livWgsNcvs_enh-DNase.bed > livWgsNcvs_enh-DNase.txt