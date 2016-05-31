#!/usr/bin/python


#########
## 
## May 24, 2016
## Purpose: To count the # of NCVs in a particular regulatory element
## Usage: python analyze_element_counts.py <infile>
##
## Input: file is the output of command:
##    intersectBed -wo -a autosomesHg19_enhPriority_merged_enh.bed -b liver_keep_L3-WGS_hg19.bed > merge_enh
##
##########

import sys

if len(sys.argv) != 2:
    print '''usage: {0} <input filename> \n\n'''.format(sys.argv[0])

inf = sys.argv[1]


# For line in file, compare current and previous line.
# Each line/element starts with a mutation count of 1.
# If the next line has the same or contiguous element, count +=1.
# If the next line has the same Cosmic sampleID, check if it is 
# an independent mutation or a dinucleotide mutation.


prevLine = ""          ## save previous line: http://stackoverflow.com/questions/17373118/read-previous-line-in-a-file-python


with open( inf,"r" ) as f:

    chr_pos = ""
    st_pos = ""
    end_pos = ""
    
    for line in f:
        
        # Variable to control if should continue searching in the current RE or 
        # move to next line
        nextl = 1 
        
        # If not the first line,
        # split both current and prevLine on tab
        if( len(prevLine.split("\t")) > 1 ):
            pl = prevLine.split("\t")
            cl = line.split("\t")
            
            
            # Check if current line is in the same reg element as previous line:
            # Either same start or the adjacent elements.
            # If yes, this element gets a count of +=1
            if( pl[1]==cl[1] or pl[2]==cl[1] ):
            
                chr_pos = cl[0]
                st_pos = pl[1]
                end_pos = cl[2]
                count +=1
                
                
                # If in the same element, check if the mutations are in the same patient.
                # If yes, check if they are indpendent mutations (not subsequent nts) -> count +=0
                # or if it is likely dinucleotide mutation -> count -=1
                # If none are true, script continues to next set of lines.
                if( cl[8] == pl[8] ):
                
                    if( cl[5]!=pl[5] and cl[5]!=pl[6] ):
                        continue
                        
                    elif( cl[5] == pl[6] ):
                        count -=1
                        
                nextl = 0  # signal to stay on the current element.
                        
            # If the prevLine and current line are different elements,
            # print information from the previous element saved in variables,
            # trigger nextl to continue to new set of input lines.            
            elif( pl[1]!=cl[1] and pl[2]!=cl[1] ):
                
                print "%s\t%s\t%s\t%s" % ( chr_pos,st_pos,end_pos,count)
                
                chr_pos = cl[0]
                st_pos = cl[1]
                end_pos = cl[2]
                count = 1
                nextl = 1  # signal to move on to next set of lines in file.
        
        
        # Given if should continue to next set of lines or not,
        # assign prevLine as specified and reset count if needed
        if( nextl == 1 ):
            prevLine = line
            count = 1
            
        elif( nextl == 0 ):
            prevLine = prevLine
            count = count

        