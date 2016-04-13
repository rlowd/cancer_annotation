#!/usr/bin/python

################################################################
################################################################
##
##  Programmer : Rebecca
##  Date: Sept 2015 -- modified from <union_REclass_space.py>
##
##  Purpose: Create output files of "regulatory element space" *for each cell type* 
##           from chromHMM for all roadmap data. "RE classes" are regulatory, transcribed, and silent.
##
##  Input: 
##    * text file with list of chromHMM file locations to parse
##    * prefix for vector file outputs
##    * directory location for ouputs
##
##  Outputs include:
##    * /outdir/vector/ contains 3 files: pref+regul, tx, or quies 
##      contaiing a vector of bp length for that element in each cell type quiered.
##    * /outdir/bedfile/ contains 1 bedfile for each RE class for each cell type (3x # CTs input) 
##      of non overlapping genomic positions in hg19. 
##
##  Usage: python RE_space.py <file list> <output prefix> <output dir>
## 
################################################################
################################################################

from itertools import izip
import sys,gzip,os

if len(sys.argv) != 3:
    print '''usage: {0} <file list> <output prefix> <bedfiles output dir>\n\n'''.format(sys.argv[0])

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

    dirPath = sys.argv[3]
    vec_dir = dirPath+"/vector/"
    vec_pref = vec_dir+sys.argv[2]
    
    if not os.path.exists(vec_dir):
        os.makedirs(vec_dir)

    print "creating vector files at: "+vec_pref+"\n"   
    regul_vec = vec_pref+"_regul_1-2-3-6-7.vect"
    tx_vec = vec_pref+"_transcribed_4-5.vect"
    quies_vec = vec_pref+"_quies_repeat.vect"

    Ef = start(sys.argv[1])
    for f in Ef: 
        dat = gzip.open( f.rstrip("\n"),"r" ).readlines()
        print "opening "+f+" ...\n"
               
        cellName = f.split("/")[9].split(".")[0]   ## was [2] -- [9] for running w/ canAnn.sh
        bed_dir = dirPath+"/bedfile/"
        pref = bed_dir+cellName
        
        if not os.path.exists(bed_dir):
            os.makedirs(bed_dir)

        print "creating bedfiles at: "+pref+"\n"
        regul = pref+"_regul_1-2-3-6-7.bed"
        tx = pref+"_transcribed_4-5.bed"
        quies = pref+"_quies_repeat.bed"
        
        regul_num = 0
        regul_bp = 0
        tx_num = 0
        tx_bp = 0
        quies_num = 0
        quies_bp = 0
        repeat_num = 0
        repeat_bp = 0

        with open( regul,"w" ) as regulf, open( tx,"w" ) as txf, open( quies,"w" ) as quiesf:

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

        with open( regul_vec,"a" ) as rvf:
            rvf.write(str(regul_bp)+"\n")
        rvf.close()
        with open( tx_vec,"a" ) as tvf:
            tvf.write(str(tx_bp)+"\n")
        tvf.close()
        with open( quies_vec,"a" ) as qvf:
            qvf.write(str(quies_bp)+"\n")
        qvf.close()

if __name__=="__main__":
    read_data()
