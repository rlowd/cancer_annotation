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


if len(sys.argv) != 2:
    print '''usage: {0} <uniqAlleles file>  \n\n'''.format(sys.argv[0])


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


def make_seqs():    
    ## First make a dict of the seq_records
    index_dict = SeqIO.index("fastaWithAlleles","fasta",alphabet=IUPAC.unambiguous_dna)
#    index_dict = SeqIO.index("new","fasta",alphabet=IUPAC.unambiguous_dna)

    agree = 0
    disagree = 0
    exclusive = 0
    
    WTdict = {}
    MUTdict = {}
    
    for k in index_dict.keys():
        
#        print k+"\t"+index_dict[k]+"\n"
        
        cur_seq_record = index_dict[k]
        mutable_seq = index_dict[k].seq.tomutable()
        
        WT = cur_seq_record.id.split("|")[1]
        MUT = cur_seq_record.id.split("|")[2]
        ref = mutable_seq[20]
        
        ## Change seq to mutable_seq and get WT and MUT alleles
        ## Check to see if patient WT allele matches reference allele        
#        print WT+"\t"+MUT+"\t"+ref
        
        ## How often does the ref == WT?
        ## If TRUE, print ref allele to patientWT.fa
        ## and print MUT alleles to patientMUT.fa
#        print "ref equals WT? "+str( ref.upper() == WT )
        if( ref.upper() == WT ):
            agree += 1
            
            mutable_seq[20]=MUT
            new_seq = mutable_seq.toseq()
            new_seq.id = cur_seq_record.id
            MUTdict[new_seq.id] = new_seq
            
            WTdict[cur_seq_record.id] = cur_seq_record.seq
            
#            print "ref == WT+\n"
#            print "original: "+cur_seq_record.id+"\t"+cur_seq_record.seq
#            print "new MUT: "+new_seq.id+"\t"+MUTdict[new_seq.id]+"\n"
        
        ## How often does the REF == MUT?
#        print "ref equals MUT? "+str( ref.upper() == MUT )
        if( ref.upper() == MUT ):
            disagree += 1
            
#            print "disagree\n"
        
        # How often are all three alleles different?
#        print "all three different? "+str( ref.upper() != MUT and ref.upper() !=WT )
        if( ref.upper() != MUT and ref.upper() !=WT ):
            exclusive += 1
            
            mutable_seq[20]=WT
            new_seq = mutable_seq.toseq()
            new_seq.id = cur_seq_record.id
            WTdict[new_seq.id] = new_seq
            
#            print "ref != WT != MUT+\n"
#            print "original: "+cur_seq_record.id+"\t"+cur_seq_record.seq
#            print "new WT: "+new_seq.id+"\t"+WTdict[new_seq.id]+"\n"
            
            mutable_seq[20]=MUT
            new_seq = mutable_seq.toseq()
            new_seq.id = cur_seq_record.id
            MUTdict[new_seq.id] = new_seq

#            print "new MUT: "+new_seq.id+"\t"+MUTdict[new_seq.id]+"\n"

#        return WTdict,MUTdict
    print "agree: "+str(agree)+"\tdisagree: "+str(disagree)+"\texclusive: "+str(exclusive)
#    print "Num entries in MUT dict: "+str(len(MUTdict))+"\n"
#    print "Num enteries in WT dict: "+str(len(WTdict))+"\n"
    
    
    find_motifs( WTdict,MUTdict )


def find_motifs( WT,MUT ):

    ## First call funciton to build a dict where 
    ## key = motif name and value = motif pssm
    motif_dict = build_motif_db()
    
    for mk in motif_dict.keys():
#        print mk+"\n"+str(motif_dict[mk])+"\n"
        PSSM = motif_dict[mk][0]
        threshold = motif_dict[mk][1]
        RPSSM = PSSM.reverse_complement()
    
        for k in WT.keys():
#            print k+"\tWT: "+WT[k]+"\tMUT: "+MUT[k]
            
            WT_fwd_lo = PSSM.calculate( WT[k] )
            WT_rev_lo = RPSSM.calculate( WT[k] )
            
#            print "\nWT_fwd_lo:\n"
#            print WT_fwd_lo
#            print "\nWT_rev_lo:\n"
#            print WT_rev_lo
        
            MUT_fwd_lo = PSSM.calculate( MUT[k] )
            MUT_rev_lo = RPSSM.calculate( MUT[k] )
            
#            print "\nMUT_fwd_lo:\n"
#            print MUT_fwd_lo
#            print "\nMUT_rev_lo:\n"
#            print MUT_rev_lo
            
            delta_fwd_lo = MUT_fwd_lo - WT_fwd_lo
            delta_rev_lo = MUT_rev_lo - WT_rev_lo
            
#            print "\ndelta_fwd_lo:\n"
#            print delta_fwd_lo
#            print "\ndelta_rev_lo:\n"
#            print delta_rev_lo
            
            with open( "motifLO_byPosition.txt","a" ) as outf:
                for i in range(0,len(WT_fwd_lo)):
                    if( (WT_fwd_lo[i] > threshold or MUT_fwd_lo[i] > threshold or WT_rev_lo[i] > threshold or MUT_rev_lo[i] > threshold) and (delta_fwd_lo[i] != 0.0 or delta_rev_lo[i] != 0.0) ):
                        outf.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (k,mk,i,WT_fwd_lo[i], MUT_fwd_lo[i],delta_fwd_lo[i], WT_rev_lo[i], MUT_rev_lo[i],delta_rev_lo[i],threshold) )
            outf.close()
            


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
    
    