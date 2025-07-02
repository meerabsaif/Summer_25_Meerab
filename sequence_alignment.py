#July 2nd, 2025
#Pairwise Sequence Alignment 
"""
1. In this file pairwise aligment is performed using Biopython. 
2. Created a Python Script using sys module (to pass arguments from command line).
3. Throughout the script I implemented modular functions to perform each task:
- Aligning sequences 
- Calculating Similarity %
- Calculating gap frequency 
- Finding conserved regions in sequences 
4. Error Handling was done in each function.
5. Comments are added to explain the code. 
"""

import sys                    #importing sys module to give argument from command line 
from Bio import AlignIO       #importing Align10, to read sequence files 
from Bio import pairwise2     #importing pairwise2, for seq alignment 
from Bio import SeqIO

import seq

def align_seq(filepath, filepath2, match=2, mismatch= -2,gap_open = -2,gap_extended= -2):    #defining fn to align sequences, also setting scoring parameters
    try: 
        seq1 = SeqIO.read(filepath, "fasta").seq
        seq2 = SeqIO.read(filepath2, "fasta").seq
        alignment = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open ,gap_extended)      #performs alignment with provided scores 
        best_alignment = alignment [0]                #selecting best alignment from results, it is present at 0 index

        print("the alignment score is this:", best_alignment.score)
        print("the alignment seq 1 is this:", best_alignment.seqA)
        print("the alignment seq 2 is this:", best_alignment.seqB)
        print("the start of alignment is this:", best_alignment.start)    #point from where the alignment starts 
        print("the end of alignment is this:", best_alignment.end)        #point till alignment is done 

        return best_alignment 
    
    except Exception as error:
        print(f"Error in aligenment: {error}")
        sys.exit(1)

#defining function to calculate similarity b/w aligned sequences 
def similarity (alignment):
    try: 
        seq1= alignment.seqA                 #extracting first aligned seq and storing it in variable seq1
        seq2 = alignment.seqB
        start = alignment.start              #alignment start position
        end = alignment.end                  #alignment end position
        aligned1 = seq1[start:end]           #extracting aligned portion of seq1
        aligned2 = seq2[start:end]           #extracting aligned portion of seq2

        matches = 0                                #creating counter to store base matches 
        for i in range(len(aligned1)):               
            if aligned1[i] == aligned2[i] and aligned1 != '-' :      #ignore gaps (hyphen)
                matches += 1                                     #incrementing matches counter when base match

        print ('Base matches are:', matches)
        length = end - start     #subtracting start from alignment from end , to calculate length of aligned region
        
        #calculating similarity percentage
        similarity =( matches/length ) * 100  if length > 0 else 0               #if length >0 to avoid dividing from zero 
        print("The similarity of alignment is:", similarity, '%')
        return 

    except Exception as error:
        print(f"Error in calculating similarity percentage: {error}")
        sys.exit(1)

#Defining a fn to calculate gap frequency in the aligned sequences
def gap_frequency(alignment): 
    try: 
        seq1= alignment.seqA 
        seq2 = alignment.seqB
        gapsA = seq1.count('-')           #counting gaps in seq1
        gapsB = seq2.count('-')           #we use - to represent gaps
        print("checking", gapsA)

        gap_frequency_1 = gapsA / len(seq1) if len(seq1) > 0 else 0      ##calculating frequency for seq1
        gap_frequency_2 = gapsB / len(seq2) if len(seq2) > 0 else 0 

        print("gap frequency of sequence 1:", gap_frequency_1)
        print("gap frequency of sequence 2:", gap_frequency_2)

        return 

    except Exception as error:
        print(f"Error in calculating gap frequency: {error}")
        sys.exit(1)


def finding_conserved_regions(alignment):
    try: 
        seq1= alignment.seqA            
        seq2 = alignment.seqB
        start = alignment.start          #start position of alignment
        end = alignment.end              #end position of aligment 
        aligned1 = seq1[start:end]
        aligned2 = seq2[start:end]

        conserved_regions = []           #initializng a list to store conserved regions 
        current_region = ''              #a string to build and store conserved regions 
        pointer = None                   #creating a variable called pointer, to track start of a conserved region, setiing it to None initially 

        for i in range(len(aligned1)):
            sequence1= aligned1[i]       #getting base at position i 
            sequence2= aligned2[i]
            if sequence1 == sequence2 and sequence1 != '-' and sequence2 != '-' :     #base matching and no gaps check
                if pointer is None :         #from start 
                    pointer = i 
                current_region += sequence1         #adding matching base to current region 

            else:                                    #if bases done match 
                if len(current_region) >= 20:        #if region >= 20
                    conserved_regions.append(current_region)       #adiing regions to list 
                    # print("Conserved region found:", current_region) 
                    # print("Length of Conserved region:", len(current_region))
                
                current_region = ''
                pointer = None 

        
        if len(current_region) >= 20:           #checking last region
            if current_region not in conserved_regions:
                conserved_regions.append(current_region)

                print("Conserved region found:", current_region) 
                print("Length of Conserved region:",len(current_region))

        if conserved_regions:
            print("Total conserved regions ≥20 bp:",len(conserved_regions))
        else:
            print("No conserved regions of 20 bp or longer found")
        
        return conserved_regions

    except Exception as error:
        print(f"Error in identifying conserved regions: {error}")
        sys.exit(1)



if __name__ == '__main__' : 
    if len(sys.argv) != 3:          #only 2 command line arguments to be passed 
        print("To Use: python alignment.py COX1_gene_homosapiens.fasta COX1_gene_zebrafish.fasta")     #to use
        sys.exit(1)

    filepath= sys.argv[1]
    filepath2= sys.argv [2]


    alignment_results= align_seq(filepath, filepath2)  #aligning sequences and assiging to alignment_results variablr, to call our functions later 
    similarity(alignment_results)
    gap_frequency(alignment_results)
    finding_conserved_regions(alignment_results)

