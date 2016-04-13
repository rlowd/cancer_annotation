#!/usr/bin/python


################################
################################
## Programmer : Rebecca
## Date: 20 March 2016
## Purpose: 
##
##
## Usage: python filterNCV.py <uniqAlleles file> 
##
################################
################################


import sys,gzip,os
import numpy as np
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import motifs


## The uniqAlleles file is formatted as:
## col1 = position (chr:start-end)
## col2 = WT allele
## col3 = MUT allele
##
## NOTE: getNCVallele.sh is a wrapper script that creates this file from a 3 column bedfile
## using grep on the gzipped NCV db to get the allele information.


if len(sys.argv) != 2:
    print '''usage: {0} <uniqAlleles file>  \n\n'''.format(sys.argv[0])


## This function takes the input and parses it to get the chr,start,end data,
## plus the WT and MUT alleles. IDs for each unique position-allele set is formatted as follows:
## <position|WT|MUT>. This information is written to the "alleles" file and
## will become the ID for the FASTA sequence (alleles variable) in the next function.
## Then write the position +/-20bp as a 1 line bedfile. Use this bedfile and
## system.os() to call fastaFromBed commandline tool. 
## Finally, the ID and its FASTA sequence are pasted together appended ot the <fasta> tmp file.


def make_fasta():
    uniqfile = sys.argv[1]
    with open( uniqfile,"r" ) as inf:
        for line in inf:
            
            ## Get NCV position information...
            l = line.split("\t")
            pos = l[0]
            ch = l[0].split(":")[0]
            bases = l[0].split(":")[1]
            st,end = ( bases.split("-")[0],bases.split("-")[1] )
            
            ## Get alleles information and print to alleles file
            wt,mut = l[1],l[2].rstrip("\n")
            alleles = pos+"|"+wt+"|"+mut
            
            with open( "alleles","w" ) as allelesf:
                allelesf.write( alleles+"\n" )
            allelesf.close()
            
            ## Write NCV position in bed format, +/-20bp
            ## Then write to "inbed" file, which is used for fastaFromBed
            bedline = "chr"+ch+"\t"+str(int(st)-20)+"\t"+str(int(end)+20)
            
            with open( "inbed","w" ) as bedf:
                bedf.write( bedline+"\n" )
            bedf.close()
            
            ## fastaFromBed should be hg38!!!
            cmd1 = "fastaFromBed -fi /bar/genomes/hg38/hg38.fa -bed inbed -fo out"
            os.system( cmd1 )
            
            ## Put the alleles together with its fasta seq
            ## "fasta" file is what will be used to make mutated seqs
            ## At this point, inspect the "fasta" file
            ## to ensure that the fasts seq makes sense with the given
            ## WT allele, at lest in most cases
            cmd2 = "paste out alleles >> fasta"
            os.system( cmd2 )
     
    inf.close()


## Take the <fasta> file and for each header line, split on tabs and
## save the second element (as "newl"). newl = the position|WT|MUT allele
## tag created in make_fasta(). This tag becomes the new header.
## How you get to this point doesn't really matter -- but having the 
## allele information you want to mutate to in the header line is the important
## part for the make_seq() function to work.


def reformat():
    infile = "fasta"
#    with open( infile,"r" ) as inf, open( "new","w" ) as newf:
    with open( infile,"r" ) as inf, open( "fastaWithAlleles","w" ) as newf:
        for line in inf:
            if line.startswith( ">" ):
                l = line.split("\t")
                newl = l[1].rstrip("\n")
                newf.write( ">"+newl+"\n" )
            else:
                newf.write( line )
    inf.close()
    
    cmd3 = "rm fasta"
    os.system( cmd3 )


## make_seqs() is the meat of this script. This function uses Biopython to create an
## index dict of the fastaWithAlleles file where keys are the header lines and 
## values are the sequence.
## Use the method mutable_seq() to change a seq_record to a data type that can be mutated.
## (To preserve integrity, you can't mutate most seq items in Biopython. You have to
## make it a specific mutable object.)
## Once the mutable_seq is made, use it to create the WT and MUT alleles, as needed.
## I also count how many times the WT allele matches the ref and/or the MUT, for my records.
## Finally, the mutated seq (WT or MUT) are associated with thier ID and added to either
## the WT or MUT dictionaries. These dicts are passed to the next function, find_motifs().


def make_seqs():    
    ## First make a dict of the seq_records
    index_dict = SeqIO.index("fastaWithAlleles","fasta",alphabet=IUPAC.unambiguous_dna)

    agree = 0
    disagree = 0
    exclusive = 0
    
    WTdict = {}
    MUTdict = {}
    
    for k in index_dict.keys():
               
        cur_seq_record = index_dict[k]
        mutable_seq = index_dict[k].seq.tomutable()
        
        WT = cur_seq_record.id.split("|")[1]
        MUT = cur_seq_record.id.split("|")[2]
        ref = mutable_seq[20]
        
        ## Change seq to mutable_seq and get WT and MUT alleles
        ## Check to see if patient WT allele matches reference allele        
        
        ## How often does the ref == WT?
        ## If TRUE, print ref allele to patientWT.fa
        ## and print MUT alleles to patientMUT.fa
        if( ref.upper() == WT ):
            agree += 1
            
            mutable_seq[20]=MUT
            new_seq = mutable_seq.toseq()
            new_seq.id = cur_seq_record.id
            MUTdict[new_seq.id] = new_seq
            
            WTdict[cur_seq_record.id] = cur_seq_record.seq
              
        ## How often does the REF == MUT?
        if( ref.upper() == MUT ):
            disagree += 1
                  
        # How often are all three alleles different?
        if( ref.upper() != MUT and ref.upper() !=WT ):
            exclusive += 1
            
            mutable_seq[20]=WT
            new_seq = mutable_seq.toseq()
            new_seq.id = cur_seq_record.id
            WTdict[new_seq.id] = new_seq
                       
            mutable_seq[20]=MUT
            new_seq = mutable_seq.toseq()
            new_seq.id = cur_seq_record.id
            MUTdict[new_seq.id] = new_seq


    print "agree: "+str(agree)+"\tdisagree: "+str(disagree)+"\texclusive: "+str(exclusive)
    
    find_motifs( WTdict,MUTdict )


## find_motifs() will actually do the motif scanning.
## First, build_motif_db() method is called (see below) which creates a 
## dict where the key = motif name and value = an array where the first element
## is the PSSM for that motif and the second element is the score threshold.
## For each key in motif_dict, save the PSSM, threshold, and PSSM for the reverse complement.
## Then for each FASTA in the WT dict (which is paried with a FASTA of the SAME KEY in the
## MUT dict), use the calculate() method from Biopython to find the PSSM score for 
## _each start position_ for that sequence on both alleles, on both strands.
## Then determine the delta value.
## Based on these values, write values to the outfile depending on criteria expressed
## in the if statement on line 223.


def find_motifs( WT,MUT ):

    ## First call funciton to build a dict where 
    ## key = motif name and value = motif pssm
    motif_dict = build_motif_db()
    
    for mk in motif_dict.keys():
        PSSM = motif_dict[mk][0]
        threshold = motif_dict[mk][1]
        RPSSM = PSSM.reverse_complement()
    
        for k in WT.keys():
            
            WT_fwd_lo = PSSM.calculate( WT[k] )
            WT_rev_lo = RPSSM.calculate( WT[k] )

            MUT_fwd_lo = PSSM.calculate( MUT[k] )
            MUT_rev_lo = RPSSM.calculate( MUT[k] )
            
            delta_fwd_lo = MUT_fwd_lo - WT_fwd_lo
            delta_rev_lo = MUT_rev_lo - WT_rev_lo
                       
            with open( "test-out.txt","a" ) as outf:
                for i in range(0,len(WT_fwd_lo)):
                    if( (WT_fwd_lo[i] > threshold or MUT_fwd_lo[i] > threshold or WT_rev_lo[i] > threshold or MUT_rev_lo[i] > threshold) and (delta_fwd_lo[i] != 0.0 or delta_rev_lo[i] != 0.0) ):
                        outf.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (k,mk,i,WT_fwd_lo[i], MUT_fwd_lo[i],delta_fwd_lo[i], WT_rev_lo[i], MUT_rev_lo[i],delta_rev_lo[i],threshold) )
            outf.close()
            

## Function to bulid the motif_dict needed in the find_motifs() function.
## <jaspar_curated.pfm> file is the JASPAR vertbrates motif file that I have
## curated to include only motifs of interst (currently 99).


def build_motif_db():
    motif_dict= {}
    background = {'A': 0.3, 'C': 0.2, 'T': 0.3, 'G': 0.2}
    jf = open("jaspar_curated.pfm")
    for m in motifs.parse(jf,"jaspar"):
        pwm = m.counts.normalize(pseudocounts={'A': 0.6, 'C': 0.4, 'T': 0.6, 'G': 0.4})
        pssm = pwm.log_odds(background)
        distribution = pssm.distribution(background=background, precision = 10**3)
#        threshold = distribution.threshold_patser()
        threshold = distribution.threshold_fpr(0.001)
        motif_dict[m.name] = (pssm,threshold)
    return motif_dict


if __name__=="__main__":
    make_fasta()
    reformat()
    make_seqs()
    
    