# main.py
from extract_promoters import parse_gtf, get_genome_dict, extract_promoter_sequences
from extract_features import extract_features
from model import train_model
import pandas as pd
import numpy as np

#defining File Paths
gtf_file = "data/annotations.gtf"
fasta_file = "data/genome.fa"
expr_file = "data/expression.tsv"

#loading expression data 
def load_expression_data(expression_file):
    """
    This function correctly reads 'expression.tsv' file.
    It groups all the data for the same gene and calculates average expression (TPM) across all tissues.
    """
    #Reading the tab-separated file (.tsv)
    df = pd.read_csv(expression_file, sep='\t')

    #Renaming 'Gene name' to 'gene' to match the other data
    df = df.rename(columns={"Gene name": "gene"})

    #Grouping by gene column and calculating the mean of TPM column
    #gives one average expression value for each unique gene
    gene_expression = df.groupby('gene')['TPM'].mean().reset_index()

    #Applying a log transform to expression values. 
    #common step in biology to make data easier for model to learn 
    gene_expression["expression"] = np.log1p(gene_expression["TPM"])

    #Returning only the columns we need
    return gene_expression[["gene", "expression"]]

#Pipeline

#Extracting Promoter Sequences from genome files
print("Parsing files:")
promoters = parse_gtf(gtf_file)
genome = get_genome_dict(fasta_file)
sequences = extract_promoter_sequences(promoters, genome)
print(f"Extracted {len(sequences)} promoter sequences.")

#Extracting Features from the DNA sequences
print("Extracting features:")
features = extract_features(sequences)
features_df = pd.DataFrame(features)
print("Features extraction Complete.")

#Saving complete set of extracted features to a CSV file for review
features_df.to_csv("extracted_features.csv", index=False)
print("Extracted features saved to extracted_features.csv")

#Loading Expression Data and Merging it w seq features
print("Loading and merging expression data:")
expression_df = load_expression_data(expr_file)

merged_df = pd.merge(features_df, expression_df, on="gene", how="inner") #Merging two dataframes based on gene column, 'how="inner"' means only keeping genes in both files.
print(f"{len(merged_df)} genes with both sequence and expression data.")

#Train the Machine Learning Model
if not merged_df.empty:
    print("\nTraining the model:")
    train_model(merged_df)
    print("\nSaved to 'predicted_gene_activity.csv'.")
else:
    print("\nno genes were found in both the annotation and expression files")