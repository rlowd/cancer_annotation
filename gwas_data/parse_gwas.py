#!/usr/bin/python

#################
##
## Programmer: Rebecca
## Date: August 2015
## Purpose: Filter GWAS catalog for cancer-related SNPs
##
#################

import sys

#print sys.argv[0]
#print sys.argv[1]
#print len(sys.argv)

if len(sys.argv) != 3:
    print '''usage: {0} <input GWAS db> <output bedfile> \n\n'''.format(sys.argv[0])
    
def start(X):
    try:
        print 'opening file :',X
        infile = open(X,"r").readlines()
        print 'Total ',len(infile),' lines.'
        return infile
    except IOError,message:
        print >> sys.stderr, "cannot open file",message
        sys.exit(1)
        

def parse_gwas():
    start( sys.argv[1] )
    gw = open( sys.argv[1],"r" ).readlines()
#    l = gw[0].split('\t')
#    print l[7]
#    print '\n'
    
    with open( sys.argv[2],"w" ) as outf:
        for line in gw:
            l = line.split( '\t' )
            if "cancer" in l[7] and l[12] != '':
                st = int(l[12])
                end = st+1
                outf.write( l[11]+'\t'+str(st)+'\t'+str(end)+"\t"+l[7]+'\n' )
            if "cancer" in l[7] and l[12] == '':
                print line 
    
    outf.close()


if __name__=="__main__":
    parse_gwas()