#!/usr/bin/python


import sys,gzip,os

if len(sys.argv) != 2:
    print '''usage: {0} <site value>  \n\n'''.format(sys.argv[0])


def parse_data():

    s = 'CosmicWGS_SamplesExport.tsv.gz'   ## v75 NCVs
    with gzip.open( s,"r" ) as inf:
        
        dt={}
        next( inf )                       ## skip header line
        for line in inf:

            l = line.split( "\t" )
            t = (l[2],l[3],l[6],l[7])
            dt[l[0]] = t
    
    inf.close()
    get_tissue( dt )

def get_tissue( dict ):
    
    ts = sys.argv[1]
    ts_set = set()
    
    with open( 'tissue_sampleIDs','w' ) as tsf:
        for k,v in dict.iteritems():
            s1,s2,h1,h2 = v
            if h2==ts:
                ts_set.add( k )
                tsf.write( k+"\t"+s1+"\t"+s2+"\t"+h1+"\t"+h2+"\n" )
    
    tsf.close()
    
    filter_NCVs( ts_set )
    

def filter_NCVs( ts_set ):

    cosmic = 'CosmicWGS_NCV.tsv.gz'
    with gzip.open( cosmic,"r" ) as inf, open( "tmp","w" ) as tmpf:
       
        next( inf )                       ## skip header line
        for line in inf:

            l = line.split( "\t" )
            
            if ( l[10] != "" ):                         # some variants have missing FATHMM score values.
                fs = float(l[10])                               ## if not missing, make the value a float
            if ( l[10] == "" or float(l[10]) < float(0.0) ):    ## if FATHMM value missing or < 0.0, 
                fs = float(0.0)                                 ## assign it to be 0.0.
            
            ## Edit line 59 for each filtering level --------
            if l[1] in ts_set and l[9]=="n" and l[14]=="y": # and l[6]=="Confirmed somatic variant"
                a = l[5].split( ":" )
                c = "chr"+a[0]
                b = a[1].split( "-" )
                st = b[0]
                end = int(st)+1
                tmpf.write( c+"\t"+st+"\t"+str(end)+"\t"+l[2]+"\t"+l[1]+"\t"+str(fs)+"\n" )

    inf.close()
    tmpf.close()
    
    
    out = sys.argv[1]+"_L2_hg19.bed"
    
    cmd = "liftOver -bedPlus=5 tmp hg38ToHg19.over.chain.gz hg19 unmapped"
    os.system( cmd )
    cmd2 = "bedSort hg19 "+out
    os.system( cmd2 )



if __name__=="__main__":
    parse_data()