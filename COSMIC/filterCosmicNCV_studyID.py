#!/user/bin/python

################
## Programmer : Rebecca
## Date: October 2015
## Purpose: Filter and pre-process CosmicNCV.tsv. Can filter by study ID to isolate
## tumor origin-specific variants; also filters for only confirmed somatic variants.
##
## Outputs include:
##   * <output pref>.tsv with only variants from the specifiec tissue, 
##     liftOver to hg19 and bedSorted.
##
## Usage: python RE_space.py <CosmicNCV.tsv> <output prefix> 
################


from itertools import izip
import sys,gzip,os

if len(sys.argv) != 3:
    print '''usage: {0} <CosmicNCV.tsv> <output prefix> \n\n'''.format(sys.argv[0])

cosmic = sys.argv[1]
out = sys.argv[2]+".tsv"
#idfile = sys.argv[3] <study ID list>

def parse_data():

    with open( cosmic,"r" ) as inf, open( out,"w" ) as outf:

#        elist = e.readlines()
#        for i in elist:
#            a = i.split("\t")
#            j = 0
#            for 
#        elist = [19,251,351,419]
        
        next( inf )
        for line in inf:

            l = line.split( "\t" )
            if (int(l[16])== 19 or int(l[16])== 251 or int(l[16]) == 351 or int(l[16]) == 419) and l[6]=="Confirmed somatic variant":
                a = l[5].split( ":" )
                c = "chr"+a[0]
                b = a[1].split( "-" )
                st = b[0]
                end = b[1]
                outf.write( c+"\t"+st+"\t"+end+"\t"+l[2]+"\t"+l[16]+"\n" )

    inf.close()
    outf.close()

if __name__=="__main__":
    parse_data()