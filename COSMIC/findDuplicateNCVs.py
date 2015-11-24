#!/user/bin/python

################################
################################
## Programmer : Rebecca
## Date: 18 Nov 2015
## Purpose: Filter and pre-process CosmicNCV.tsv v75
##   * Filter to retain only confirmed somatic variants (col 9)
##	 * Retain only variants from samples that have been WG resequenced
##   * Retain only variants NOT previously detected ("n" in field 10)
##   * Keep sampleID field in output (col 2)
##
## Outputs include:
##   * <output pref>.tsv with only filtered variants, 
##     liftOver to hg19 and bedSorted.
##
## Usage: python RE_space.py <CosmicNCV.tsv> <output prefix> 
################################
################################


import sys,gzip,os

if len(sys.argv) != 3:
    print '''usage: {0} <CosmicNCV.tsv> <output prefix> \n\n'''.format(sys.argv[0])

cosmic = sys.argv[1]
out = sys.argv[2]+"_gt2.bed"


def get_gt2():

    inf = open( cosmic,"r" )
    vars = inf.readlines()
    
    (c0,s0) = ("","")
    
    with open( "tmp","w" ) as tmpf:
       
        for i in range( 0, len(vars)-1 ):

            line = vars[i]
            l = line.split( "\t" )
            c = l[0]
            s = l[1]
                
            nline = vars[i+1]
            nl = nline.split( "\t" )
            nc = nl[0]
            ns = nl[1] 
            
            if ( nc == c and ns == s):
                tmpf.write( line )

    inf.close()
    tmpf.close()


def get_unique():

#    vars = open( "tmp","r" ).readlines()    ## if using both functions
    inf = open( cosmic,"r" )               ## if only want to remove duplicates,
    vars = inf.readlines()                 ## not find variants >=2x in db
    uni = set(vars)
    
    with open( "tmp2","w" ) as tmpf:
        for line in uni:
            tmpf.write( line )
    
    tmpf.close()
    
    cmd = "bedSort tmp2 "+out
    os.system( cmd )
    

if __name__=="__main__":
    get_gt2()
    get_unique()
