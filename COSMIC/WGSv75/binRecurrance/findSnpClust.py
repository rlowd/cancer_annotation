#!/usr/bin/python

################################
################################
##
## Programmer : Rebecca
## Date: 10-16 Dec 2015
## Purpose: Count number of NCVs (from filtered set) in a given binsize
##    * Binsize is sliding -- e.g. every 500bp bin starting with each NCV
##
## Input args:
##    1. SNP bed file
##    2. bin size
##    3. output file name
##
## Outputs:
##    1. bedfile with binsize windows and # of SNPs per bin in col4
##    2. *_<binsize>_stat.tsv file with the value of SNPs per bin and count of
##       how many bins have x number of SNPs.
##
## Usage: python findSnpClust.py <SNP file> <binsize> <outfile name>
##
################################
################################


import sys,os
import numpy as np


if len(sys.argv) != 4:
    print '''usage: {0} <SNP file> <binsize> <outfile name>  \n\n'''.format(sys.argv[0])
    

## This function reads in the input file. It first calls the get_chr_sizes() function
## which will set the max bp for a bin end point based on the appropriate chr.
## The infile data is read in with readlines(). This is so the get_bin() function
## can iterate over the data and remember which line it is on.
## Then for each line, this function gets the chr and start position and
## passes this information to the get_bin() function.

def split_bins():
    
    snpfile = sys.argv[1]
    binsize = sys.argv[2]
    outfile = sys.argv[3]+"_bin"+str(binsize)+".bed"
    
    with open( snpfile,"r" ) as inf, open( outfile,"w" ) as outf:
        
        sizes = get_chr_sizes()
        d = inf.readlines()
        count = 0        
        
        for i in range(0,len(d)):
            
            l = d[i].split("\t")
            ch,st = l[0],l[1]
            
            chr_max = sizes[ch]
            count,bin_end = get_bin( i,ch,st,chr_max,d )
                       
            outf.write( ch+"\t"+st+"\t"+str(bin_end)+"\t"+str(count)+"\n" )
            
    inf.close()
    outf.close()


## This function takes 4 inputs from split_bins() --
## i = line number, ch = chromosome value for line i
## st = start position for the NCV in line i
## chr_max = the last chr bp position for the chr in line i
## d = the entire dataset read in with readlines() above.
## This function uses this information to count the number of SNPs in the bin
## given the SNP position does not exceed the binsize or the chr end point.
## Returns the # of SNPs and the bin end position.

def get_bin( i,ch,st,chr_max,d ):

    snp_count = 0
    binsize = int(sys.argv[2])
    bin_end = int(st) + binsize
    
    if len(d)-binsize > i:   ## This line is to check if the line is close to the end of the file
        for n in range(0,binsize-1):     ## Theoretically, only the next # of lines = binsize can possibly have SNPs that should included in this bin.
            cur_st = int(d[(i+n)].split("\t")[1])
            cur_ch = d[(i+n)].split("\t")[0]
        
            if bin_end > (chr_max - (binsize-1)):   ## check that the end of the bin isn't exceeding the chr end position
                bin_end = chr_max                   ## if TRUE, set chr end position as bin_end
        
            if cur_ch == ch and cur_st < (int(st)+binsize):
                snp_count+=1

            if cur_ch != ch:
                continue
    
    else:                   ## If close to the end of the file,
        z = len(d)-i        ## Then only iterate over the remaining lines (z)
        for n in range(0,z):
            cur_st = int(d[(i+n)].split("\t")[1])
            cur_ch = d[(i+n)].split("\t")[0]
        
            if bin_end > (chr_max - (binsize-1)):
                bin_end = chr_max
        
            if  cur_st < (int(st)+binsize) and cur_ch == ch:
                snp_count+=1
            
            if cur_ch != ch:
                continue

    return snp_count,bin_end


## Get end chr position for current chr

def get_chr_sizes():

    f = 'chr_sizes'
    dt={}
    sizes = open( f,"r" ).readlines()
    a = np.genfromtxt(f,delimiter="\t",dtype=None)
    
    for c in range(0,len(a)):
        dt[a[c][0]] = a[c][1]

    return dt


## Count number of bins with x number of SNPs

def get_statistics():

    binsize = sys.argv[2]
    statf = sys.argv[3]+"_bin"+str(binsize)+".bed"
    outf = statf+"_"+str(binsize)+"_stat.tsv"
    
    with open( outf,"w" ) as out:
        out.write( "Bin size\tNum\n" )
    out.close()
    
    string1 = "cut -f 4 "
    string2 = "| sort - | uniq -c - | awk \'{print $2\"\\t\"$1}\' - |sort -n -k1 - >> "
    
    cmd = string1+statf+string2+outf
    os.system( cmd )


if __name__=="__main__":
    split_bins()
    get_statistics()


