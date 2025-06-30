#Day 6 | June 30th, 2025
#Task: Convert a Genbank file into fasta file. 


from Bio import SeqIO             #Importing SeqIO from Biopython, it handles reading and writing sequence data in different formats 
from Bio.Seq import Seq

#Reading the GenBank file
seq = SeqIO.read("gen.gbk", "genbank")     #reading a single seq from genbank file, storing it into variable seq 
#                                          #syntax ( 'file name', 'file type') file type tells the format file needs to be parsed in. 

#Writing it to a Fasta file
with open("converted.fasta", "w") as f:    #we write the sequence to fasta file, and also mention the name of the file we want to generate. 'w' write mode 
    SeqIO.write(seq, f, "fasta")           #seq contains the sequence that was read from the GenBank file.
#                                          #f is where the sequence will be written.
#                                          #"fasta" Tells SeqIO to set the output format as Fasta file. (includes header (>,sequence identifier) and sequence.
    print("Genbank file is converted to FASTA")      #if this executes successfully, this message will be printed 

