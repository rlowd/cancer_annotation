#!/usr/bin/python

################################
################################
## Programmer : Rebecca
## Date: 14 June, 2016 (copied from filterNCVs.py orig date = 10-16 Dec 2015)
## Purpose: Filter and pre-process CosmicNCV.tsv v75
##   * Improves on filterCosmicNCV_v74.py script
##   1 First get sampleIDs from the Samples feature file.
##   2 Select SampleIDs for feature of interest.
##   3 Isolate CosmicNCVs by sample feature.
##   3.1 Also filter for desired level:
##       - remove population variants
##       - remove ExomeSeq OR WG-reseq'd vars
##       - remove non-confirmed somatic vars
##   4 Identify and remove hypermutated variants. Hypermuated threshold = top 2.5%tile
##
## Output files (3):
##   * <filter_term_pref_hg19>.tsv with only filtered variants and bedSorted.
##   * <filter_term_pref_hg19>_keepSampleIDs - list of sampleIDs below NCV threshold.
##   * <filter_term_pref_hg19>_removeSampleIDs - list of sampleIDs above NCV threshold.
##
## Usage: python filterNCV.py <site/hitsology metadata filter term> 
################################
################################

import sys,gzip,os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
plt.plot([1,2,3])


if len(sys.argv) != 2:
    print '''usage: {0} <site value>  \n\n'''.format(sys.argv[0])


## This function opens the Cosmic sample features file
## and creates a dict with the sample IDs as keys, and
## primary and seconday site and histology metadata as values.
## The dict is then passed to the get_tissue() function to filter NCVs
## according to desired metadata field.

def get_sampleIDs():
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



## This function takes as an input the sampleID-metadata dict
## and scans it for the tissue metadata desired.
## !! Note !! modify the -if- statement below to change which metadata field is scanned.
## Then the filtered sampleIDs are assed to the set ts_set 
## (and written to an outfile) and this set is passed to the next function
## which will use the sampleIDs to select NCVs.

def select_tissue( dict ):    
    ts = sys.argv[1]
    ts_set = set()
    with open( 'tissue_sampleIDs','w' ) as tsf:
        for k,v in dict.iteritems():
            s1,s2,h1,h2 = v
#            print s1
            if s1==ts:
                ts_set.add( k )
                tsf.write( k+"\t"+s1+"\t"+s2+"\t"+h1+"\t"+h2+"\n" )
    tsf.close()
#    return ts_set
    
    filterByTissue( ts_set )
    


## This function takes as input the tissue set IDs and the NCVs to be filtered.
## After parsing and cleaning, the NCV entry may be written to the tmp outfile
## if it passes filter requirements.
## !! Note !! to change filter requirements, modify the 3rd -if- statement below.
## Filtered NCVs are then sorted and mapped to hg19.

def filterByTissue( ts_set ):
    cosmic = 'CosmicWGS_NCV.tsv.gz'
#    cosmic = 'ncv'
    filename = "tmp"
    with gzip.open( cosmic,"r" ) as inf, open( filename,"w" ) as tmpf:
        
        next( inf )                       ## skip header line
        for line in inf:
            
            l = line.split( "\t" )
            
            if ( l[10] != "" ):                         # some variants have missing FATHMM score values.
                fs = float(l[10])                               ## if not missing, make the value a float
            if ( l[10] == "" or float(l[10]) < float(0.0) ):    ## if FATHMM value missing or < 0.0, 
                fs = float(0.0)                                 ## assign it to be 0.0.
                
            if ( l[9] == "" ):
                pvar = "y"
            else:
                pvar = l[9]
                
            ## Edit line 115 for each filtering level --------
#            if l[1] in ts_set and l[9]=="n" and l[14]=="y" and l[6]=="Confirmed somatic variant":  ## for WGS seq vars
            if l[1] in ts_set and l[9]=="n" and l[15]=="y" and l[6]=="Confirmed somatic variant":  ## for ExomeSeq vars

                a = l[5].split( ":" )
                c = "chr"+a[0]
                b = a[1].split( "-" )
                st = b[0]
                end = int(st)+1
                tmpf.write( c+"\t"+st+"\t"+str(end)+"\t"+l[2]+"\t"+l[1]+"\t"+pvar+"\t"+str(fs)+"\n" )
                
    inf.close()
    tmpf.close()
    
    remove95( filename )
    

def remove95( filename ):
    
    sdt={}
    with open( filename,"r" ) as inf:
        for line in inf:              ## for loop to count # of times each sampleID is in the NCV file
            l = line.split("\t")
            idv = l[4]
            if idv in sdt:
                sdt[idv] += 1
            else:
                sdt[idv] = 1
        
#    values = sdt.values() 
#    values.sort()
#    q95 = np.percentile( values,95 )
#    
#    print "95%tile value: "
#    print q95
#        
#    filt={}
#    for k,v in sdt.iteritems():
#        if v <= q95:
#            filt[k] = v
#        
#    fvalues = filt.values()
        
#    plt.hist(fvalues)
#    plt.savefig('NCVs_per_sample_hist.png')
    
    find_hypermut( sdt,filename )


def find_hypermut( dt,filename ):
    
    values = dt.values()
    values.sort()
    cutoff = np.percentile( values, 98.5 ) #97.5 )
    
    print "\ncutoff: 99.0%tile\n"
    print cutoff
    
    keep={}
    remove={}
    
    for k,v in dt.iteritems():
        if v < cutoff:
            keep[ k ] = v
        else:
            remove[ k ] = v
               
    kvals = keep.values()
    rvals = remove.values()
    
    kvals.sort()
    rvals.sort()
    
    print "\nkeep vals min, max:\n"
    print str(kvals[0])+"\t"+str(kvals[-1])
    print "\nremove vals min, max:\n"
    print str(rvals[0])+"\t"+str(rvals[-1])
    
    kkeys = keep.keys()
    rkeys = remove.keys()
    
    kkeys.sort()
    rkeys.sort()
    
    kkfn = sys.argv[1]+"_ExSeq_keepSampleIDs_v77_985"
    rkfn = sys.argv[1]+"_ExSeq_removeSampleIDs_v77_985"
    
    with open( kkfn,"w" ) as kkf:
        for i in range(0,len(kkeys)):
            kkf.write( str(kkeys[i])+"\n" )
    kkf.close()
    
    with open( rkfn,"w" ) as rkf:
        for i in range(0,len(rkeys)):
            rkf.write( str(rkeys[i])+"\n" )
    rkf.close()

    remove_hypermut( kkeys,rkeys,filename )
     
#    values = keep.values()
    
#    plt.hist(values)
#    plt.savefig('NCVs_per_sample_KeepHist.png')
   
#    plt.hist(values)
#    plt.savefig('NCVs_per_sample_RemoveHist.png')
    
    



def remove_hypermut( keep,remove,filename ):
    keepname = 'keep'
    removename = 'remove'
    with open( filename,"r" ) as inf, open( keepname,'w' ) as keepf, open( removename,'w' ) as removef:
        sdt={}
        for line in inf:        
            l = line.split( "\t" )
            smplID = l[4]
            if smplID in keep:
                keepf.write( line )
            if smplID in remove:
                removef.write( line )
                     
    cleanup( keepname )
    cleanup( removename )



def cleanup( filename ):
    out = sys.argv[1]+"_"+filename+"_ExSeq_hg19_v77_985.bed"
    cmd0 = "bedSort "+filename+" "+out
    os.system( cmd0 )
    
    cleanup1 = "rm tmp"
    cleanup2 = "rm "+filename
    
#    os.system( cleanup1 )
#    os.system( cleanup2 )
    

if __name__=="__main__":
    get_sampleIDs()
