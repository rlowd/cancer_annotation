#!/usr/bin/python

################
## Programmer : Rebecca
## Date: Sept 2015 -- modified from <union_REclass_space.py>
## Purpose: Create output files of "regulatory element space" *for each cell type* 
##          from chromHMM for all roadmap data. "RE classes" are regulatory, transcribed, and silent.
##
## Input: text file with list of chromHMM file locations to parse
## Outputs include:
##   * report file with # of each chromatin state element and total # of bp for that element in hg19.
##   * bedfile for each RE class for each cell type (3x # CTs input) with non overlapping genomic positions in hg19. 
##
## Usage: python RE_space.py <file list> <output prefix>
################


from itertools import izip
import sys,gzip
#import array


if len(sys.argv) != 2:
    print '''usage: {0} <file list> <output prefix>\n\n'''.format(sys.argv[0])

def start(X):
    try:
        print 'opening file :',X
        infile = open(X,"r").readlines()
        print 'Total ',len(infile),' lines.'
        return infile
    except IOError,message:
        print >> sys.stderr, "cannot open file",message
        sys.exit(1)

def read_data():

    vec_pref = sys.argv[2]

## create vector files that will have a list/vector of total bp for each class
    regul_vec = vec_pref+"_regul_1-2-3-6-7.vect"
    tx_vec = vec_pref+"_transcribed_4-5.vect"
    quies_vec = vec_pref+"_quies_repeat.vect"

    print "create output files\n"
    print "outfile prefix: "+vec_pref+"\n"
    stat = pref+"_RE_stats.report"
    rep = pref+"_roadmap_data.report"
    regul = pref+"_regul_1-2-3-6-7.bed"
    tx = pref+"_transcribed_4-5.bed"
    quies = pref+"_quies_repeat.bed"

    Ef = start(sys.argv[1])
    for f in Ef: 
        dat = gzip.open( f.rstrip("\n"),"r" ).readlines()
        repf.write( f+" contains "+str(len(dat))+" lines\n" )
        print "opening "+f+" ...\n"

        pref = f

        regul_num = 0
        regul_bp = 0
        tx_num = 0
        tx_bp = 0
        quies_num = 0
        quies_bp = 0
        repeat_num = 0
        repeat_bp = 0

        print "create output files\n"
        print "outfile prefix: "+pref+"\n"
        regul = pref+"_regul_1-2-3-6-7.bed"
        tx = pref+"_transcribed_4-5.bed"
        quies = pref+"_quies_repeat.bed"

        with open( regul,"w" ) as regulf, open( tx,"w" ) as txf, open( quies,"w" ) as quiesf, open( rep,"w" ) as repf:

            for line in dat:
                l = line.split("\t")
                if int(l[3])== 3 or int(l[3])== 6 or int(l[3])==7 or int(l[3])==1 or int(l[3])==2:
                    regul_num +=1
                    regul_bp += int(l[2]) - int(l[1])
                    regulf.write( line )
                if int(l[3])==4 or int(l[3])==5:
                    tx_num +=1
                    tx_bp += int(l[2]) - int(l[1])
                    txf.write( line )
                if int(l[3])==15 or int(l[3])==9 or int(l[3])==12 or int(l[3])==13 or int(l[3])==14 or int(l[3])==10 or int(l[3])==11 or int(l[3])==8:
                    quies_num +=1
                    quies_bp += int(l[2]) - int(l[1])
                    quiesf.write( line )

        regulf.close()
        txf.close()
        quiesf.close()

        with open( regul_vec,"w" ) as rvf:
            rvf.write(str(regul_num)+"\n")
        rvf.close()
        with open( tx_vec,"w" ) as tvf:
            tvf.write(str(tx_num)+"\n")
        tvf.close()
        with open( quies_vec,"w" ) as qvf:
            qvf.write(str(quies_num)+"\n")
        qvf.close()

if __name__=="__main__":
    read_data()
