import sys 

#Writing a code to read fasta file involves understanding how a fasta file is read and involves coding for each operation accordingly 
#A file contains multiple lines, a code is read line by line and a sequence may span on more than one lines.     
#Each sequence first starts with a header > , then comes the ID and then the sequence. 
""" 
Ive created a function which opens a fasta file provided as argument and reads it line by line. 
It performs the following operations: 
1. cleans the file, removing any space or new lines usually present in the fasta file. I have used strip function. 
2. It checks if the file starts with > , if it does it means its a new sequence. 
3. Theres another condition in the same block, if key ID consists of any value it means its not the first sequence. 
4. A sequence was already present, which we save to the dictionary. 
5. To extract the sequence ID, i have used indexing to remove >. seq ID is tracked using the variable key_id. 
6. The value_seq is reset to now store the upcoming new sequence. 
7. If line doesnt starts with > , it means its the same sequence, so in else the line is added to the seq value and stored in the sequence value. 
8. After loop ends, we check with if statement that key ID is not null, if its not then we store the last sequence . 
9. The dictionary is returned at end of function. 
"""

def read_fasta (name):

    seq_Dict = {}              #initializing an empty dictionary 
    key_id = None              #to store our ID's , keep track of the sequence we are dealing with 
    values_seq = ""            #to store sequence, mainly using it to add the multiple line sequences. 

    with open (name, 'r') as x:
        for line in x :
            line=line.strip()         #removing any spaces or new lines 

            if line.startswith ('>') :

                if key_id is not None :
                    seq_Dict [key_id]= values_seq

                key_id = line [1: ]    
                values_seq = ""            #resetting the values, to avoid data mixing. 

            else :
                values_seq = values_seq + line 
 
        if key_id is not None : 
            seq_Dict [key_id]= values_seq
                                                         
    return seq_Dict


def analysis (a):
    dnaseq = a 
    len_dna = len(dnaseq)
    gc_count = (dnaseq.count('G') + dnaseq.count('C'))
    gc_content = gc_count / len_dna *100 if len_dna > 0 else 0 

    return gc_content 

def valid_seq (a):

    valid_sequence = {'A', 'T', 'C', 'G'}
    a = a.upper ()
    for x in a: 
         if x not in valid_sequence:
             return False 
    return True 

        
if __name__ == "__main__":
    if len(sys.argv) != 2 : 
        print ("To use this script: python Script_to_FASTA_file_processing.py <fasta file name>")
        sys.exit (1)

    name = sys.argv [1]
    seq_dict = read_fasta (name)


#Identifying unique nucleotides                  #combined sequences in a string and converted into a set.
    total_nucleotides = set()                    #initializing an empty set  
    for seq_values in seq_dict.values():                     #reading values from the dictionary 
        total_nucleotides.update(seq_values)                  #adds the current sequence to the set of total nucleotides 
    print("Unique Nucleotides:", total_nucleotides)
    
#preparing results
    analysis_results = []                 #creating list, keeping it empty 
    for id, seq in seq_dict.items():       #storing results from the dictionary we created above 
        length = len(seq)                   #assigning values by calling functions 
        gc_content = analysis(seq)
        validity = valid_seq(seq)
        analysis_results.append([id, length, gc_content, validity])       #storing in form of a list 
    
#Writing results to a CSV file
    with open("analysis_results.csv", "w") as f:        # w refers to write mode, also overwrites already existing file 
        f.write("ID,Length,GC_Content,Validity\n")      #writing a heading 
        for result in analysis_results:                 #using loop to write down results for each sequence
            result = id, length, gc_content, validity   #storing in individual variables  
            f.write(f"{id},{length},{gc_content},{validity}\n")      #using f string 
    print("Results saved to 'analysis_results.csv'")
    
    
