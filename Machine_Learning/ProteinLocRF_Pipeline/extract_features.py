#This file contains functions to extract physicochemical features from protein sequences (in a FASTA file).
#And then assigning synthetic labels for classification. 

#importing libraries
import sys 
import pandas as pd 
from Bio import SeqIO          #will be used to parse files
from Bio.SeqUtils.ProtParam import ProteinAnalysis     #to compute physiochemical properties for protein sequences 
import numpy as np               #for performing numerical operations 
import logging                   #used to log information about the script's execution for debugging and tracking


def protein_features_analysis(fasta_file):          #defining a function that takes fasta file as an argument. 
    try:                                            #to handle potential errors
        all_features = [ ]                          #creatin an empty list to store fetures 
        #file = SeqIO.parse(fasta_file, "fasta")     #(pass parameter, write format type )      #perse: ghanghalna in urdu

        for record in SeqIO.parse(fasta_file, "fasta"):      #iterates over every seq in file
            #record contains id and sequence
            seq_id = record.id                               #extracting ids and storing in a variable 
            sequence = record.seq                            #extracting sequences and storing in a variable
            
            #adding a line before protein analysis           
            #sequence= str(record.seq).replace("-","").replace("","")    #did this for cleaning

            protein = ProteinAnalysis(str(sequence))                          #applying protein analysis on our sequence also converting seq into a string

            features = {                                    #making a dict to store physiochemical features for the protein
                #"id" : record.id,
                "molecular_weight": protein.molecular_weight(),                  #retrieving molecular weight from the variable protein
                "isoelectric point": protein.isoelectric_point(),                #isoelectric point, ph at wich protein has no net charge
                "aromacity": protein.aromaticity(),                              #calculating fraction of aromatic amino acids
                "instability_index": protein.instability_index(),                #Estimating protein stability (lower values : more stable proteins
                "gravy":protein.gravy(),                                         #GRAVY: Grand average of hydropathicity
            }
            all_features.append(features)      #append stores at last of list, adding feature dictionary to all_features (empty list) we created in the begining of fn

            #break                                       #only one value needed, break fn. 
        #df= pd.DataFrame([features])                      #this will only store features of last protein. so we make a list above 

        df = pd.DataFrame(all_features)                #creating a df of the list, row : point at each protein and columns are the features 
        df.to_csv('proteinfeatures.csv', index=False)           #saving the df to a csv file, excluding index
        logging.info(df)                          #logging dataframe content to the log file 
        return df                                  #returning the dataframe
    except FileNotFoundError:                      #error handling
        print("file not found")
        sys.exit(1)
    sys.exit(1)

#we have to predict in which part of cell these proteins work, so we will assign labels 
#this fn assigns random synthetic labels to each protein in teh df, to classify later

def synthetic_labels(df):          #defining a fn that takes df as input (giving features)
    try:
        labels = ['cytoplasm', 'nucleus', 'mitochondria']      #defining a list of labels
        df['labels'] = np.random.choice(labels, size=len(df))  #adding a new column names labels to the dataframe
#                                                              #randomly selecting a label from the column for each row in the df
#size=len(df) : creating an array of length equal to no. of rows in our df

        return df                       #this would return a modified df which has a new column label
    except Exception as e:              #error handling
        print("no dataset found", e)
        sys.exit(1)