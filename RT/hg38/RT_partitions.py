#!/usr/bin/python


#########
## 
## June 3, 2016
## Purpose: To partition the replication timing data into larger bins
## Usage: python RT_partitions.py <infile>
##
## Input: file is RT data downloded from FSU
##
##########

import sys
import numpy as np

if len(sys.argv) != 2:
    print '''usage: {0} <input filename> \n\n'''.format(sys.argv[0])

inf = sys.argv[1]

with open( inf,"r" ) as f:

    for _ in range(16):
        f.next()

    prevLine = f.next()
    cur_sign = 0
    nextl = 1
    
    for line in f:
        
        pl = prevLine.split("\t")
        cl = line.split("\t")
        
        ( p_c,p_g,p_ch,p_st,p_end,p_delta ) = pl
        ( c_c,c_g,c_ch,c_st,c_end,c_delta ) = cl
        

        if( np.sign(float(p_delta.rstrip())) == np.sign(float(c_delta.rstrip())) ):
            cur_sign = np.sign(float(c_delta.rstrip()))
            bin_end = c_end
            
            nextl = 0    ## stay on current interval
                
        elif( np.sign(float(p_delta.rstrip())) != np.sign(float(c_delta.rstrip())) ):
            print "%s\t%s\t%s\t%s" % ( p_ch,p_st,bin_end,cur_sign)

            nextl = 1  # signal to move on to next set of lines in file.
        
        
        
        # Given if should continue to next set of lines or not,
        # assign prevLine as specified and reset count if needed
        if( nextl == 1 ):
            prevLine = line
            cur_sign = np.sign(float(c_delta.rstrip()))
            bin_end = c_end
            
        elif( nextl == 0 ):
            prevLine = prevLine
            
            

        