#!/usr/bin/python

################################
################################
## Programmer : Rebecca
## Date: 15 Dec 2015
## Purpose: Filter and pre-process CosmicNCV.tsv v75
##   * Based on filterCosmicNCV_v74.py script
##
## Outputs include:
##   * CosmicWGS_<feature term>_sampleFeat.tsv - table of metadata for each NCV in input file (CosmicWGS_NCV.tsv.gz)
##   * CosmicWGS_<feature term>_sampleFeat.tsv-NCVsPerSample_byFeature - count of NCVs per sampleID in the feature-filtered NCVs
##   * CosmicWGS_<feature term>_sampleFeat.tsv-numSamples_byFeature - count of number of samples for each unique set of site/histology features
##
## Usage: python annotateSamplesWithFeatureData.py <site/hitsology metadata filter term> 
################################
################################

import sys,gzip,os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.plot([1,2,3])


if len(sys.argv) != 2:
    print '''usage: {0} <site/hitsology metadata filter term>  \n\n'''.format(sys.argv[0])


def get_sampleIDs():

## This function opens the Cosmic sample features file
## and creates a dict with the sample IDs as keys, and
## primary and seconday site and histology metadata as values.
## The dict is then passed to the get_tissue() function to filter NCVs
## according to desired metadata field.

    s = 'CosmicWGS_SamplesExport.tsv.gz'   ## v75 NCVs
    with gzip.open( s,"r" ) as inf:
        
        dt={}
        next( inf )                       ## skip header line
        for line in inf:

            l = line.split( "\t" )
            t = (l[2],l[3],l[6],l[7])
            dt[l[0]] = t
    
    inf.close()
    select_tissue( dt )


def select_tissue( dict ):

## This function takes as an input the sampleID-metadata dict
## and scans it for the tissue metadata desired.
## !! Note !! modify the -if- statement below to change which metadata field is scanned.
## Then the filtered sampleIDs + sample feature values are passed to the set ts_set 
## (and written to an outfile - for QC only).
## The list of sampleIDs and the sampleID-features dict is passed to the next function
## which will use the sampleIDs to select NCVs.
    
    ts = sys.argv[1]
    ts_set = set()
    ts_dict={}
    
    with open( 'tissue_sampleIDs','w' ) as tsf:
        for k,v in dict.iteritems():
            s1,s2,h1,h2 = v
#            if s1==ts:n                    ## turn off this line to run without selecting feature
#                ts_set.add( k )
            tsf.write( k+"\t"+s1+"\t"+s2+"\t"+h1+"\t"+h2+"\n" )
            ts_dict[ k ] = v
    
    tsf.close()
    
    associateTissue( ts_dict )
    

def associateTissue( ts_dict ):

## This function takes as input the tissue set IDs and the NCVs to be filtered.
## After parsing and cleaning, the NCV entry may be written to the tmp outfile
## if it passes filter requirements.
## !! Note !! to change filter requirements, modify the 3rd -if- statement below.
## Filtered NCVs are then sorted and mapped to hg19.

    cosmic = 'CosmicWGS_NCV.tsv.gz'
    out = "CosmicWGS_NCV_"+sys.argv[1]+"_sampleFeat.tsv"
    with gzip.open( cosmic,"r" ) as inf, open( out,"w" ) as outf:
        next( inf )                       ## skip header line
        for line in inf:
            l = line.split( "\t" )
            smplID = l[1]
            if smplID in ts_dict.keys(): 
                s1,s2,h1,h2 = ts_dict[ smplID ] 
                outf.write( smplID+"\t"+s1+"\t"+s2+"\t"+h1+"\t"+h2+"\n" )
    inf.close()
    outf.close()
    quantify( outf )
    
    
def quantify( f ):
    inf = os.path.basename(f.name)
    outf1 = inf+"-NCVsPerSample_byFeature"
    outf2 = inf+"-numSamples_byFeature"
    cmd1 = "sort -k1 -n "+inf+" | uniq -c - | awk 'OFS=\"\t\" {print $1,$2,$3,$4,$5,$6}' - > "+outf1
    os.system( cmd1 )
#    cmd2 = "cut -f 3-6 "+outf1+" | sort - | uniq -c - | awk 'OFS=\"\t\" {print $1,$2,$3,$4,$5}' - > "+outf2
    cmd2 = "cut -f 3,5-6 "+outf1+" | sort - | uniq -c - | awk 'OFS=\"\t\" {print $1,$2,$3,$4}' - > "+outf2
    os.system( cmd2 )


if __name__=="__main__":
    get_sampleIDs()
