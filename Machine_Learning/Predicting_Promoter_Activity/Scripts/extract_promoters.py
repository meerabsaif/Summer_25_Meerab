#extract_promoters.py
"""
This script processes genomic data to extract promoter sequences.
They are regions of DNA that are present upstream of genes.
They regulate gene expression. 
I have used a GTF (Gene Transfer Format) file for gene annotation.
For genome sequence i have used a FASTA file. 
"""
from Bio import SeqIO       #to handle biological sequences
import pandas as pd         #for data manipulation 

def parse_gtf(gtf_file):         #defining a fn to parse gtf file and extract promote regions
    promoters = []               #initializing an empty list to store promoter regions data later
    with open(gtf_file) as file:        #opening the file
        for line in file:               #iterating over each line in file
            if line.startswith("#"):     #to skip comment/ metadata lines
                continue
            fields = line.strip().split('\t')       #Removing whitespaces and using tab as a sperater
            #promoters are defind to TSS 
            if fields[2] == "transcript":    #we are filtering for the feature type (transcript), index 2: we have 3rd feild which is transcript
                chrom = fields[0]                #extracting chrom name
                start = int(fields[3])            #start position of transcript
                end = int(fields[4])              #end position of transcript
                strand = fields[6]                #strand which indicated direction of transcription
                info = fields[8]                  #attribute field which contains metadata(gene name)

                #Finding gene_name from last column (attribute)
                gene_name = ""             #initializing empty string to store gene name later
                for item in info.split(';'):        #splitting attribute to a list of key :value pairs on basis of ;
                    if "gene_name" in item:        #checking if attribute has gene_name
                        gene_name = item.strip().split('"')[1]     #splits based on "", and remove whitespace, and extract second element which is gene_name
                        break             #when found, not looking further
                if not gene_name:         #if no gene name found, skip
                    continue              #skipping to next line

                #Defining the promoter region, located upstream of the transcription start site 
                if strand == '+':        #For genes on positive strand, upstream is before the TSS
                    promoter_start = max(1, start - 1000)    #Setting it to start 1000 bp before transcript start, making sure it doesn’t go below position 1
                    promoter_end = start            #Setting promoter end at the transcript start.
                else:              #For genes on negative strand, upstream is after the TSS
                    promoter_start = end              #Setting promoter start at transcript end
                    promoter_end = end + 1000        #extending 1000 bp downstream 

                promoters.append((chrom, promoter_start, promoter_end, strand, gene_name))     #storing promoter data as a tuple
    return promoters         #returning list of promoter region

#this fn loads genome sequence from a FASTA file into a dictionary, keys :chromosome names , values :sequences.
def get_genome_dict(fasta_path):

    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))     #parsing fasta file, and converting parsed records into dictionaries
    return genome          #returning dict of sequences

#extracts DNA sequences for promoter regions from the genome
def extract_promoter_sequences(promoter_regions, genome):      #giving list of promoter regions and genome dictionary as argument
    sequences = []                        #initializing an empty list to store gene names and promoter seuqneces later
    if not genome:                    #if dict is empty
        print("Genome dictionary is empty")
        return []               #returning empty list 

#Retrieving the chromosome sequence
    chrom_key = list(genome.keys())[0]      #get first chrom ID from dict
    chrom_seq = genome[chrom_key].seq       #extract seq of that ID

    for chrom, start, end, strand, gene in promoter_regions:       #Iterating over each promoter to extract its sequence
        promoter_dna_seq = chrom_seq[start-1:end]        #Slicing chromosome seq to get the promoter region
        if strand == '-':       #checking if gene is on neg strand 
            promoter_dna_seq = promoter_dna_seq.reverse_complement()         #using reverse complement to reverse the seq, and replace each base w its compliment

        sequences.append((gene, str(promoter_dna_seq).upper()))        #storing gene name and promoter sequence in upper case, and converting it into a string
    return sequences       #returns list of gene names and promoter sequences 