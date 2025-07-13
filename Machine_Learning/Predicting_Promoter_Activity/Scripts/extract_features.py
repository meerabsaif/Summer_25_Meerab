"""
This script will compute features from promoter sequences
(GC content, AT content, presence of a TATA box, dinucleotide repeat scores and k-mer frequencies)
they are inputs for the machine learning model.
"""

import re                 #finding dinucleotide repeats
import numpy as np      #For numerical computations
from collections import Counter       #For counting k-mers in sequences.
import math 

#to Count frequency of all k-mers, kmer are the subsequences of length k in a seq
def kmer_counts(seq, k=3):            #defining a fn that takes sequence and k-mer length (3) as an input
    return Counter([seq[i:i+k] for i in range(len(seq)-k+1)])     #list of all k-mers in creayed: by sliding window of size k over the sequence
#Counter counts frequency of each k-mer

#calculating GC content
def gc_content(seq):
    g = seq.count('G')     #counting occurences of G
    c = seq.count('C')     #same as G
    return (g + c) / len(seq) if len(seq) > 0 else 0      #returning gc content after calculating, returns 0 if seq is empty to avoid division by 0 

##calculating GC content
def at_content(seq):
    a = seq.count('A')
    t = seq.count('T')      ##counting occurences
    return (a + t) / len(seq) if len(seq) > 0 else 0     ##returning at content after calculating

#checking is seq has tata box
#tata box is a common promoter element 
def tata_box(seq):        #passing seq as argument 
    return 'TATAAA' in seq       #if the string is present it will return true

#Counting no. of dinucleotide repeat motifs in seq
def dinucleotide_repeat_sc(seq):
    count = 0       #initializing counter for repeat motif.
    for dinuc in ['AT', 'TA', 'CG', 'GC', 'AG', 'GA', 'CT', 'TC']:   #iterating over a list of pairs
        pattern = f"({dinuc}){{3,}}"  #creating pattern for 3 or more repeats of the dinuc
        if re.search(pattern, seq):         #checking if pattern exists
            count += 1              #if it does, incrementing the counter 
    return count          #returning total no. of dinucleotide repeats found

#critical core promoter element. It's a binding site for common transcription factors
def has_gc_box(seq):      #passing seq as argument 
    #Checking for the presence of the GC Box motif (GGGCGG)
    return 'GGGCGG' in seq        #if the string is present it will return true


"""
Checking for a simplified initiator element.
The Inr is at end of the promoter (at the Transcription Start Site)
We will check the last 10 bases of our sequence.
It plays a key role in positioning the machinery that reads the gene
"""
def has_initiator_element(seq):
    tss_region = seq[-10:]
    return 'CA' in tss_region or 'CG' in tss_region            #checking for a pattern like 'CA' or 'CG' right at the very end

"""
Shannon Entropy: Powerful feature from information theory.
A low complexity (low entropy) sequence is very repetitive (e.g., ATATATATAT).
A high complexity (high entropy) sequence is more random and information-rich (e.g., AGTCGTTACGGA).
Promoter regions often have a characteristic complexity, and this can be a very useful feature for the model.

"""
#Calculates the Shannon entropy as a measure of sequence complexity
def shannon_entropy(seq):
    if not seq:
        return 0
    
#Getting counts of each base (A, C, G, T)
    base_counts = Counter(seq)
    seq_length = len(seq)
    
    entropy = 0.0
#Calculating entropy using the formula: H = -sum(p_i * log2(p_i))
    for base in base_counts:
        probability = base_counts[base] / seq_length
        if probability > 0:
            entropy -= probability * math.log2(probability)
            
    return entropy   #returning entropy

#defining a fn to extract features from a list 
def extract_features(sequences):    #containing a seq list which contains gene names and promoter sequences 
    featurelist = []             #initializing an empty list
    for gene, seq in sequences:        #iterating over each pair
        features = {                  #creating dict for each seq
            "gene": gene,                #gene name 
            "gc_content": gc_content(seq),
            "at_content": at_content(seq),
            "tata_box": int(tata_box(seq)),        #converting boolean to int, type casting 
            "dinuc_repeat_count ": dinucleotide_repeat_sc(seq),
            "has_gc_box": int(has_gc_box(seq)),    #converting boolean to int, type casting 
            "has_initiator": int(has_initiator_element(seq)),     #converting boolean to int, type casting 
            "shannon_entropy": shannon_entropy(seq),
        }

#adding the freq of each 3 mer to dict
        kmer_frequency = kmer_counts(seq, 3)      #getting count of each 3 mer in seq
        for kmer in kmer_frequency:              #iterating over each 3 mer in counter
            features[f'kmer_{kmer}'] = kmer_frequency[kmer]    #adding key with count
        featurelist.append(features)        #storing featire dict to list 
    return featurelist          #returning list of feature dictionaries 

