#!/usr/bin/python


import sys,os
import numpy as np


if len(sys.argv) != 2:
    print '''usage: {0} <SNP file>  \n\n'''.format(sys.argv[0])
    

def split_bins():
    
    snpfile = sys.argv[1]
    
    with open( snpfile,"r" ) as inf:
        
        sizes = get_chr_sizes()
        d = inf.readlines()
        count = 0        
        
        for i in range(0,len(d)):
            
            l = d[i].split("\t")
            ch,st = l[0],l[1]
            print "i: "+str(i)+"\tst: "+st
            
            chr_max = sizes[ch]
            count = get_bin( i,ch,st,chr_max,d )
                       
            print "returned count: "+str(count)+"\n"
            
    inf.close()


def get_bin( i,ch,st,chr_max,d ):
    snp_count = 0
    
    
    for n in range(1,6):
        cur_st = int(d[(i+n)].split("\t")[1])
        cur_ch = d[(i+n)].split("\t")[0]
        bin_end = int(st) + 200
        
        if bin_end > (chr_max - 199):
            bin_end = chr_max
        
        print "cur_chr: "+cur_ch+"\tcur_st: "+str(cur_st)+"\tbin end: "+str(bin_end)
        
        print cur_st < bin_end
        
        if  cur_st < int(st)+200 and cur_ch == ch:
            snp_count+=1
            print "snp count: "+str(snp_count)
            
    return snp_count


def get_chr_sizes():
    f = '/home/comp/twlab/rlowdon/genomes/hg19/chr_sizes'
    dt={}
    sizes = open( f,"r" ).readlines()
    a = np.genfromtxt(f,delimiter="\t",dtype=None)
    
    for c in range(0,len(a)):
        dt[a[c][0]] = a[c][1]

    return dt
    
if __name__=="__main__":
    split_bins()