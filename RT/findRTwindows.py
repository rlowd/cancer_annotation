#!/usr/bin/python

import sys,os

if len(sys.argv) != 2:
    print '''usage: {0} <site value>  \n\n'''.format(sys.argv[0])
    

def combine():
    with open( sys.argv[1] ) as inf:
        for line in inf:
            if line.startswith( "CHR" ):
                l = line.split("\t")
                start,end,val = l[3],l[4],l[5]
                print str(l[3])+"\t"+str(l[5])
                
                



if __name__=="__main__":
    combine()