#!/bin/bash

fileDir="../../chromHMM_REspace/cellType_byClass/allEID/bedfile/"
snpFile="CosmicNCV_FS_tID_hg19.bed"
#d=$(date +%Y-%m-%d)
outDir=allEID_vs_CosmNCV_FS_tID_2015-11-06

mkdir $outDir

for i in `seq 3 127`; do
	if (( "$i" < 10 )); then
		pref=E00$i
		echo $pref
	elif (( "$i" > 9 )) && (( "$i" <100 )); then
		pref=E0$i
		echo $pref
	elif (( "$i" > 99 )); then
		pref=E$i
		echo $pref
	fi
		
	regf=$fileDir${pref}_15_coreMarks_dense_regul_1-2-3-6-7.bed
	txf=$fileDir${pref}_15_coreMarks_dense_transcribed_4-5.bed
	quiesf=$fileDir${pref}_15_coreMarks_dense_quies_repeat.bed
	dest=${outDir}/${snpFile%.bed}_vs_$pref

	annotateBed -names regul transcr silent -i $snpFile -files $regf $txf $quiesf > $dest
done

