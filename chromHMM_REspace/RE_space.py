#!/usr/bin/python

################
## Programmer : Rebecca
## Date: August 2015
## Purpose: Create output files of "regulatory element space" from chromHMM for all roadmap data.
## Outputs include:
##   * report file with # of each chromatin state element and total # of bp for that element in hg19.
##   * bedfile for each chromatin state with non overlapping genomic positions in hg19. 
##
## Usage: python RE_space.py
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

    regul_num = 0
    regul_bp = 0
    tx_num = 0
    tx_bp = 0
    quies_num = 0
    quies_bp = 0
    repeat_num = 0
    repeat_bp = 0
    
    pref = sys.argv[2]

#    num = array.array('i',(0,)*15)
#    bp = array.array('i',(0,)*15)

    print "create output files\n"
    print "outfile prefix: "+pref+"\n"
    stat = pref+"_RE_stats.report"
    rep = pref+"_roadmap_data.report"
    regul = pref+"_regul_1-2-3-6-7.bed"
    tx = pref+"_transcribed_4-5.bed"
    quies = pref+"_quies_repeat.bed"
      
    Ef = start(sys.argv[1])
       
    with open( regul,"w" ) as regulf, open( tx,"w" ) as txf, open( quies,"w" ) as quiesf, open( rep,"w" ) as repf:

        for f in Ef:
            dat = gzip.open( f.rstrip("\n"),"r" ).readlines()
            repf.write( f+" contains "+str(len(dat))+" lines\n" )
            print "opening "+f+" ...\n"
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
    repf.close()

    with open( stat,"w" ) as out:
        out.write( "Total number of cell types (files): "+str(len(Ef)) )
        out.write( "\n\nRegulatory regions\nNumber of sites: "+str(regul_num)+"\nTotal bp: "+str(regul_bp) )
        out.write( "\n\nTranscribed regions\nNumber of sites: "+str(tx_num)+"\nTotal bp: "+str(tx_bp) )
        out.write( "\n\nQuiescent regions\nNumber of sites: "+str(quies_num)+"\nTotal bp: "+str(quies_bp) )
        total_genome_bp = regul_bp + tx_bp + quies_bp
        out.write( "\n\nTotal genome bp: "+str(total_genome_bp) ) 
          

if __name__=="__main__":
    read_data()
