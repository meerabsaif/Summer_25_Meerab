#Project 1 
#1st july 2025
import sys 
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.Restriction import EcoRI 
from Bio.SeqUtils import gc_fraction 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px


def read_fasta (filepath):
    try: 
        seq = list(SeqIO.parse(filepath, 'fasta'))
        for i in seq :
            seq_id = i.id
            seq_sequences = i.seq  
            print("The ID is:", seq_id )
    except FileNotFoundError:
        print("FILE NOT FOUND")
    return seq

#validating the sequence 

def validate_sequence (seq):
    try:
        is_valid=all(base in "ATGC" for base in seq)
        return is_valid
    except Exception as error:
        print(f"inavlid sequence: {error}")
        sys.exit(1) 
 
#sequence_analysis

def gc_content(seq):
    try:
        return gc_fraction(seq) * 100
    except Exception as error:
        print(f"Error calculating GC content: {error}")
        sys.exit(1)

def seq_transcribe(seq):
    try:
        rna = seq.transcribe()
        return rna 
    except Exception as error:
        print(f"Error transcribing sequence: {error}")
        sys.exit(1)   
    
def seq_translate(seq):
    try:
        protein = seq.translate()
        return protein
    except Exception as error:
        print(f"Error translating sequence: {error}")
    sys.exit(1)
    
# def codon_usage(seq):
#     try:
#         #sequence length divisible by 3
#         extra_bases = len(seq) % 3  
#         if extra_bases != 0:
#             seq = seq[:-extra_bases]  

#         #Count each codon 
#         codon_counts = {}
#         for i in range(0, len(seq), 3):  
#             codon = seq[i:i+3]
#             if codon in codon_counts:
#                 codon_counts[codon] += 1
#             else:
#                 codon_counts[codon] = 1

#         # Calculate percentages
#         total = len(seq) // 3  
#         codon_percentages = {}
#         for codon in codon_counts:
#             codon_percentages[codon] = codon_counts[codon] / total

#         return codon_percentages
#     except Exception as e:
#         print(f"Error calculating codon usage: {e}")
#         return {}  

def codon_usage(seq):
    try:
        #Convert to string and clean sequence
        seq = str(seq).upper()
        seq = ''.join(c for c in seq if c in 'ATGC')
        
        #checking seq length
        if len(seq) < 3:
            return {}
        
        extra_bases = len(seq) % 3  
        if extra_bases != 0:
            seq = seq[:-extra_bases]  

        codon_counts = {}
        for i in range(0, len(seq), 3):  
            codon = seq[i:i+3]
            if len(codon) == 3:  #Only counting complete codons
                if codon in codon_counts:
                    codon_counts[codon] += 1
                else:
                    codon_counts[codon] = 1
        #percentages
        total = len(seq) // 3  
        codon_percentages = {}
        for codon in codon_counts:
            codon_percentages[codon] = codon_counts[codon] / total

        return codon_percentages
    except Exception as e:
        print(f"Error calculating codon usage: {e}")
        return {} 

def restriction_sites(seq):    
    try:
        cut_sites = EcoRI.search (seq)   
        return len(cut_sites)
    except Exception as error:
        print(f"Error finding restriction sites: {error}")
    sys.exit(1)
        
def analyze_sequences(sequences):
    results = []
    for i in sequences:
        seq_sequences = i.seq
        if validate_sequence(seq_sequences):
            result = {
                'Sequence_ID': i.id,
                'Length': len(seq_sequences),
                'GC_Content': gc_content(seq_sequences),
                'EcoRI_Sites': (restriction_sites(seq_sequences)),  # Count of EcoRI sites
                'RNA': str(seq_transcribe(seq_sequences))[:30],
                'Protein': str(seq_translate(seq_sequences))[:50],
                'Codon_Usage': codon_usage(seq_sequences)
            }
        results.append(result)
    return results


def save_to_csv(results):
    try:
        df = pd.DataFrame(results)
        df.to_csv('sequence_analysis.csv')
        print("Results saved to sequence_analysis.csv")
        return df
    except Exception as e:
        print(f"Error saving to CSV: {e}")   
    sys.exit(1)    

def visualize_data(df):
    try: 
        #Matplotlib Bar plot
        plt.figure(figsize=(10, 6))
        plt.bar(df["Sequence_ID"],df["GC_Content"])
        plt.xlabel("Sequence ID")
        plt.ylabel('GC Content')
        plt.title('GC content per sequence')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig("gc_content_bar_plot.png", dpi=300)
        plt.close()
        print("Saved as gc_content_bar_plot.png")

        #Matplotlib Histogram
        plt.figure(figsize=(8, 6))
        plt.hist(df["Length"],bins=10)
        plt.xlabel("Sequence Length")
        plt.ylabel('Frequency')
        plt.title('Sequence Length distribution')
        plt.tight_layout()
        plt.savefig("sequence_length_histogram.png", dpi=300)
        plt.close()
        print("Saved: sequence_length_histogram.png")

        if 'Codon_Usage' in df.columns:
            valid_codon_data = [usage for usage in df['Codon_Usage'] if usage and isinstance(usage, dict)]
    
            if valid_codon_data:
                codon_usage_df = pd.DataFrame(valid_codon_data)
                codon_usage_df = codon_usage_df.fillna(0)
        
                plt.figure(figsize=(15, 10))
                sns.heatmap(codon_usage_df.T, 
                   annot=True, 
                   cmap="Blues", 
                   fmt='.3f',yticklabels=False)
                plt.title('Codon Usage Frequencies')
                plt.tight_layout()
                plt.savefig("codon_usage_heatmap.png", dpi=300, bbox_inches='tight')
                plt.close()
                print("Saved: codon_usage_heatmap.png")
            else:
                print("No valid codon usage data found")

        # Plotly scatterplot
        fig = px.scatter(df, x="GC_Content", y="EcoRI_Sites", 
                        hover_data=['Sequence_ID'],
                        title="GC Content vs EcoRI Sites")
        fig.write_image("scatterplot.png")
        print("Saved: scatterplot.png")
            
    except Exception as e:
        print(f"Error in visualization: {e}")
        
        # #Plotly scatterplot
        # fig= px.scatter(df,x="GC_content",y="EcoRI_sites")
        # plt.savefig("scatterplot.png")
        # #fig.show()
        # plt.close()

def main():
    if len(sys.argv) != 2:
        print("To Use: python sac_genome_analysis.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    sequences = read_fasta(fasta_file)   
    results = analyze_sequences(sequences)

    df = pd.DataFrame(results)
    
    df = save_to_csv(results)
    if df is not None: 
        visualize_data(df)
        print("Analysis is completed and visualizations are saved.")


if __name__ == "__main__":
    main()







