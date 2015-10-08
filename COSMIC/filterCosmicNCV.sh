awk  'BEGIN{FS="\t";OFS="\t";} NR>1 {
n = split( $6, t, ":" )
for ( j = 0; j++ < n; )
	split( t[j],s,"-" )
	for ( k = 0; k++ < 2; )
		start = s[k]
		end = start+1
		if((index($7,"somatic") !=0) == 1)
		print "chr"t[1],start,end,$3,$17
}' CosmicNCV.tsv > tmp

bedSort tmp sorted

liftOver -bedPlus=5 sorted hg38ToHg19.over.chain.gz CosmicNCV_hg19.bed unmapped

unlink tmp
unlink sorted