#!/bin/bash

for i in `ls *vcf.gz`;
do
	`gunzip $i`
	`bgzip ${i%.gz}`
	`tabix -p vcf $i`
done
