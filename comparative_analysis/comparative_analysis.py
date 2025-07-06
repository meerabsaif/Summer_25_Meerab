#Comparative Analysis of SNP Distribution and Gene Expression 
#Organism: Drosophila melanogaster 

import sys                    #importing sys module to give argument from command line 
from Bio import AlignIO        #importing Align10, to read sequence files 
from Bio import SeqIO
from Bio.Seq import Seq 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#def loading_data ():

#parsing gff file: using pandas to read and parse gff3 file
#We will add specific parameters to handle the gff file correctly. (The format is entirely different from a fasta or gbk file)
#It consists of one line per feature, each line containing 9 columns of data. 
#The gff file has values in tabular form, so we sperate info on the basis of tabs (sep='\t')
#Setting the data types for columns with exon coordinates (dtype = int) to avoid float values if any
#The file has lines starting with #, so we comment # as the information is irrelevant (comment = '')
#This file has no headers, so we need to assign names to refer to data and extract information later (names= '')

gff_column = ['chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'] #check ensembl to assign field names (for clarity)
gff_df = pd.read_csv('D_melanogaster_BDGP6.gff3', sep='\t',names = gff_column ,comment= '#' , dtype = {'start':int, 'end': int})  #passed the arguments, reason is discussed above 

#filtering the original dataset for exon features, assigining to a new dataset after filtering  
exon_df = gff_df[gff_df['type'] == 'exon'].copy()       # .copy(). was not used initially. I used it later when i was grouping exons with respective gene id from attribute column. 
#                                                        #pandas was facing an issue determining if exon_df was linked to my original gff_df, or if its a copy of it.  

#extracting gene id from attributes column to link it with exons. (will be utilized when mapping SNP data to exons)
#first we look at how the attribute column is formatted, in my file it looks like this:
#Parent=transcript:FBtr0077815;Name=FBtr0077815-E3;constitutive=1;ensembl_end_phase=1;ensembl_phase=2;exon_id=FBtr0077815-E3;rank=3
#we need to seperate the values, we will use (;) to seperate, and extract id denoted by (parent=) 

def extract_gene_id (attributes):                     #giving attributes as an argument 
    if attributes is not None:                        #checking if attributes in not empty 
        attribute_values = attributes.split(';')      #attributes are separated by (;) so we split them by to obtain seperated list values
        for i in attribute_values:                    #using a loop to read each seperated pair from attributes line 
            if 'Parent=' in i:                        #if it matches 'parents=' key 
                attribute = i.split('=')              #split by = , as id value is present after = 
                if len(attribute) > 1:                #confirming if length of splitted value is atleast 2 parts  
                    id_value = attribute [1]          #used index to return second element of the value, output showed: transcript:FBtr03089
                    gene_id = id_value.split(':')     #to remove 'transcript' we split again on basis of (:) 
                    if len(gene_id) == 2:             #ensuring length 
                      return gene_id[1]               #returning value present at 1, which is our gene id 

    return False                                      #if attribute is empty it returns false 


#adding extracted gene_id column from attribute string to my exon_df dataset by applying the fn 'extract_gene_id'
exon_df['gene_id'] = exon_df['attributes'].apply(extract_gene_id)       #this allowed applying the function after i used .copy() with my modified exon_df dataset
print("Df after adding gene_id column:") 
print(exon_df.head(5))                                                  #printing exon_df, to confirm if gene_id column is added or not. 

# def calculating_exon_length (start, end ):
#     if start is not None and end is not None : 
#         exon_length = end - start + 1 
#         return exon_length 
#     return None 

#Calculating exon length to calculate SNP density later. 
exon_df['exon_length'] = exon_df['end'] - exon_df['start'] + 1       #calculating length directly and adding as column to exon_df 

print("Df after calculating exon length and displaying it as a column :" )
print(exon_df.head(5))

#cleaning dataset = Extracting only essential columns for clarity 
exon_data = exon_df[['chromosome', 'start', 'end', 'gene_id', 'exon_length']]
print("The clean exon dataset:")
print(exon_data.head(5))

#Total exons found in dataset
print("Total exons found:", len(exon_data))

#Saving the exon coordinates clean dataset to a csv file. (For submitting to github and using in SNP mapping)
exon_data.to_csv('exon_coordinates.csv', index=False)
print("The exon data of all genes is saved to exon_coordinates.csv")


# Task 4: Parsing SNP data
#reading a VCF file containing SNP data and converting it into CSV format
#VCF file contains metadata lines, goal 1 is to remove them, they contain file information only 
#commenting (#)and seperating by ('\t')
#selecting relevant columns to prevent from later cleaning 
#also ensuring datatype of each column 
#removing null values

metadata_line_count = 0                     #to count metadata lines (starting with ##) and skip them later 

try:
    with open('snpsites.vcf', "r") as f:             #opening vcf file in read mode 
        for line in f:                               #reading each line in file opened as f 
            if line.startswith("##"):                #if any line starts with ## (showing that it is metadata)
                metadata_line_count += 1             #add if to the counter, if condition above is true
            elif line.startswith("#CHROM"):          #metadata lines end when header starts, so when header line starting with (#CHROM) is found, we stop the metadata line count
                print("VCF header:")                 #to print header 
                print(line.strip())                  #to print header line without any space
                break                                #to print only one line
    print("Metadata lines found:", metadata_line_count)               #to print no. of lines found 

except Exception as error:                                            #error handling
    print("Error reading VCF file:", error)                           #if theres an error reading the file, it prints error 
    exit()


#Reading vcf using pandas 
snp_df = pd.read_csv(
    'snpsites.vcf',
    sep='\t',                                                  #sep='\t': using tabs as separators
    skiprows=metadata_line_count,                              #skipping metadata lines we counted above
    usecols=["#CHROM", "POS", "ID", "REF", "ALT"],             #only loading the columns we need
    dtype={"#CHROM": str, "ID": str, "REF": str, "ALT": str}   #specifying datatypes of each column
)


#cleaning the position data column by finding and removing null values 
try:
    snp_df["POS"] = pd.to_numeric(snp_df["POS"], errors="coerce")           #Converting pos column to numeric, coerce converts invalid values to Nan
    if snp_df["POS"].isna().any():                                          #checking if any Nan value found, isna() would return true or false 
        print("No. of rows with invalid POS values:" , snp_df['POS'].isna().sum())          #sum is used to count Nan values
        snp_df = snp_df.dropna(subset=["POS"]).copy()                                       #.dropna() used to remove rows with nan values| .copy() used to avoid settingwithcopywarnings
        print("Invalid POS values removed")
    snp_df["POS"] = snp_df["POS"].astype(int)                               #after removing nan values we convert pos dtype from numeric to int 
except Exception as e:
    print("Error converting POS to integer: ", e )
    exit()

#Renaming the columns for clarity
snp_df = snp_df.rename(columns={
    "#CHROM": "chromosome",
    "POS": "position",
    "ID": "snp_id",
    "REF": "ref",
    "ALT": "alt"
})

#filtering the dataset and only keeping selected columns
snp_df = snp_df[["chromosome", "position", "snp_id", "ref", "alt"]].copy()         #again used .copy so panda reads it as a copy of the original df

#checking for null values in chromosome column and removing 
if snp_df["chromosome"].isna().any():                                                           #checking for null values
    print("No. of missing values found in chrom column:", snp_df['chromosome'].isna().sum())    #counting null values using sum()
    snp_df = snp_df.dropna(subset="chromosome").copy()                                          #removing rows with missing values

#displaying to check data
print(snp_df.head())                                                        #displaying the first five rows of cleaned dataset
print("Unique chromosomes:", snp_df["chromosome"].unique())                 #displaying all chromosomes found
print("Total SNPs loaded:", len(snp_df))                                    #displaying total SNPs 


#writing output to a csv file 
snp_df.to_csv('snp_data.csv', index=False)
print("The snp data is saved to snp_data.csv")


#Task 5. Mapping SNPs to Exons
#Identifying which SNPs from VCF file fall within the exonic coordinates extracted from the GFF file.  
#We will check if each SNP’s position lies within the start and end of any exon on the same chromosome. 

#printing few rows of both datasets
print(exon_data.head(30))
print(snp_df.head())

from math import ceil      #To calculate expected chunks (added this bcs system was not processing snp mapping)

#Load the exon and SNP DataFrames
exon_file = "exon_coordinates.csv"
snp_file = "snp_data.csv"

#Loading with specified data types
try:
     exon_data = pd.read_csv(            #reading csv file as dataframe
         exon_file,
         dtype={"chromosome": str, "start": int, "end": int, "gene_id": str, "exon_length": int}
     )
     snp_df = pd.read_csv(
         snp_file,                                      #also setting column dtypes below
         dtype={"chromosome": str, "position": int, "snp_id": str, "ref": str, "alt": str}
     )
     print("Exon and SNP DataFrames loaded")
except Exception as e:                                #error handling 
     print(f"Error loading CSV files: {e}")
     exit()

#Filter to chromosomes X and 4 to reduce runtime due to large SNP
#Due to certain computational limitations and time shortage i am restricting further analysis to chrom X and 4 only 

valid_chromosomes = ['X', '4']
exon_data = exon_data[exon_data["chromosome"].isin(valid_chromosomes)].copy()     #filters exon_data to only rows where chromosome is X or 4
snp_df = snp_df[snp_df["chromosome"].isin(valid_chromosomes)].copy()              #Filters snp_df to only rows w chrom X or 4
if exon_data.empty or snp_df.empty:                                 #if datasets become empty after filtering
     print("No exons or SNPs found for chromosomes X, 4")
     exit()

#Validating chromosomes 
print("Total exons after filtering:", len(exon_data))
print("Total SNPs after filtering:", len(snp_df))
print("Exon chromosomes:", exon_data["chromosome"].unique())
print("SNP chromosomes:", snp_df["chromosome"].unique())
for chrom in valid_chromosomes:                             #adding this code, to estimate total chunks to be printed for each chrom
    chrom_snps = snp_df[snp_df["chromosome"] == chrom]
    print(f"Total SNPs for {chrom}: {len(chrom_snps)} (Expected chunks: {ceil(len(chrom_snps) / 500)})")


# Step 6: Map SNPs to exons with chunk size 500
chunk_size = 500  # Specified chunk size
snp_mappings = []
for c in valid_chromosomes:
    print(f"Processing chromosome {c}...")
    chrom_snps = snp_df[snp_df["chromosome"] == c]
    chrom_exons = exon_data[exon_data["chromosome"] == c]
    
    # Skip if no SNPs or exons
    if chrom_snps.empty or chrom_exons.empty:
        print(f"No SNPs or exons found for chromosome: {c}")
        continue
    
    # Process SNPs in chunks of 500
    total_chunks = ceil(len(chrom_snps) / chunk_size)
    for start in range(0, len(chrom_snps), chunk_size):
        chunk_num = start // chunk_size + 1
        print(f"Processing chunk {chunk_num}/{total_chunks} for {c}")
        chunk = chrom_snps.iloc[start:start+chunk_size]
        
        # Merge chunk with exons
        merged = pd.merge(
            chunk,
            chrom_exons,
            on="chromosome",
            how="inner",
            suffixes=("_snp", "_exon")
        )
        
        # Filter for SNPs within exon ranges
        chunk_mappings = merged[
            (merged["position"] >= merged["start"]) & 
            (merged["position"] <= merged["end"])
        ][["snp_id", "chromosome", "position", "gene_id"]]
        
        snp_mappings.append(chunk_mappings)

# Step 7: Create DataFrame from mappings
snp_exon_map_df = pd.concat(snp_mappings) if snp_mappings else pd.DataFrame()

# Step 8: Validate results
if snp_exon_map_df.empty:
    print("No SNPs mapped to exons.")
else:
    print("Total SNP-exon mappings:", len(snp_exon_map_df))
    print("First 5 mappings:")
    print(snp_exon_map_df.head())

# Check for missing values
if not snp_exon_map_df.empty and snp_exon_map_df.isna().any().any():
    print("Warning: Missing values in mappings. Removing them.")
    snp_exon_map_df = snp_exon_map_df.dropna()

# Step 9: Save mappings
output_file = "snp_exon_mappings.csv"
snp_exon_map_df.to_csv(output_file, index=False)
print(f"SNP-exon mappings saved to {output_file}")
#I focused on X and 4 for time efficiency and biological relevance, used Pandas with chunking.

#Calculating SNP Density for Chromosome X and 4
#Snp density is the number of SNPs mapped to a gene’s exons divided by the total exonic length
#Loading SNP mapping file
snp_mapping_file = "snp_exon_mappings.csv"

try:                                #loading csv file into 2 dataframes for applying pandas data manipulation
    exon_data = pd.read_csv(        #ensuring the type they are uploaded as 
    'exon_coordinates.csv',
    dtype={"chromosome": str, "start": int, "end": int, "gene_id": str, "exon_length": int}
    )
    snp_exon_map_df = pd.read_csv(snp_mapping_file, dtype={"snp_id": str, "chromosome": str, "position": int, "gene_id": str})
except Exception as e:
    print("Couldnt load files", e)
    exit()

#Filtering to chromosomes X and 4
snp_exon_map_df = snp_exon_map_df[snp_exon_map_df["chromosome"].isin(valid_chromosomes)]      #Keeping only rows where chromosome is X or 4
exon_data = exon_data[exon_data["chromosome"].isin(valid_chromosomes)]

#Summing exon lengths for each gene
exon_length_df = exon_data.groupby("gene_id")["exon_length"].sum().reset_index()    #using gene_id to group exon data and summing their lengths, using reset index to make it a df
exon_length_df.columns = ["gene_id", "total_exon_length"]                         #renaming columns

#Counting SNPs for each gene
snp_count_df = snp_exon_map_df.groupby("gene_id")["snp_id"].count().reset_index()   #grouping by gene_id and counting snp count per gene
snp_count_df.columns = ["gene_id", "snp_count"]             #renaming column


#Merging exon lengths and SNP counts
snp_density_df = pd.merge(exon_length_df, snp_count_df, on="gene_id", how="left")       #using gene_id to merge exon length with snp count, how="left" means we are including all genes

#Setting snp_count to 0 for genes with no SNPs
snp_density_df["snp_count"] = snp_density_df["snp_count"].fillna(0).astype(int)        #we replace missing snp count value with zero, making it integer

#Calculating SNP density: (snp_count / total_exon_length)
snp_density_df["snp_density"] = snp_density_df["snp_count"] / snp_density_df["total_exon_length"]


#Printing results
print("Total genes:", len(snp_density_df))
print(snp_density_df.head)

#Saving to CSV
snp_density_df.to_csv("snp_density.csv", index=False)      #not to include index
print("Saved to snp_density.csv")                          #prints after saving file as csv

# Task 7: Correlating SNP density with gene expression for chromosomes X and 4
snp_density_file = "snp_density.csv"           #Setting the file name 
snp_density_df = pd.read_csv(snp_density_file)   #Loading SNP data into a DataFrame
expression_file = "expression.tsv"                

# Inspecting file structure to read properly 
with open(expression_file, 'r') as f:           #Openning expression.tsv in read mode
    first_10_lines = [f.readline().strip() for _ in range(10)]         #Reading first 10 lines and removing extra spaces

print("First 10 lines of the file:")        
for i, line in enumerate(first_10_lines):     #Loop through the 10 lines with their index
    print(f"Line {i}: {line}")               #Printing each line with its line no.

# Finding the header line 
header_line = None                    #Creating a variable to store the header line number
with open(expression_file, 'r') as f:         #Open file again in read mode
    for i, line in enumerate(f):           #Loop through each line with its index
        if not line.startswith('#'):            #Checking if the line does not start with #
            header_line = i                 #Saving line no. where header starts
            print(f"\nFound header at line {i}: {line.strip()}")      #Printing the header n line number 
            break              #stopping loop

#Loading the file 
if header_line is not None:       #Checking if header line was found
    #Load without specifying columns first to see what's available
    expression_df = pd.read_csv(expression_file, sep='\t', skiprows=header_line)        #Load TSV file, skipping rows up to header
    print("Available columns:")  
    for i, col in enumerate(expression_df.columns):          #Looping through column names with index
        print(f"  {i}: '{col}'")                             #each column name with its index
    print(f"\nDataframe shape: {expression_df.shape}")       #number of rows and columns in DataFrame
    print("First few rows:")  
    print(expression_df.head())  
else:                                                        #if no header line found
    print("couldnt read header")                             #if header not found
    exit() 

#Finding the correct columns
gene_id_col = None            #Createing variable for gene ID column, assigning None
expression_col = None          #variable for expression column, start with None

#Gene ID column 
for col in expression_df.columns:           #Looping through all column names
    col_lower = col.lower()                 #Converting column name to lowercase for comparison
    if any(term in col_lower for term in ['gene_id', 'gene', 'id', 'geneid']):  #Checking if column name has gene-related words
        gene_id_col = col                                           #Saving the matching column name as gene ID
        print("gene ID column:", col)  
        break  

#Expression column 
for col in expression_df.columns: 
    col_lower = col.lower()  
    if any(term in col_lower for term in ['expression', 'value', 'level', 'count', 'fpkm', 'tpm']):  # Check if column name has expression-related words
        expression_col = col                                     #Save the matching column name as expression
        print("expression column:", col)  
        break  

#If columns not found, use the first two columns
if gene_id_col is None:                        #if no gene ID column was found
    gene_id_col = expression_df.columns[0]    #Using the first column as gene ID
    print(f"Using first column as gene ID: '{gene_id_col}'")  

if expression_col is None:                           #if no expression column was found
    expression_col = expression_df.columns[1]        #Using the second column as expression
    print(f"Using second column as expression: '{expression_col}'")  

#Creating clean expression dataframe
print("Creating clean expression dataframe")  

#Keeping columns we need
clean_expression_df = expression_df[[gene_id_col, expression_col]].copy()  

#Renaming columns to standard names
clean_expression_df = clean_expression_df.rename(columns={  
    gene_id_col: 'gene_id',             #gene ID column name to gene_id
    expression_col: 'expression_level'  #expression column name to expression_level
})

#Cleaning the data
clean_expression_df['gene_id'] = clean_expression_df['gene_id'].astype(str)         #Converting to string type
clean_expression_df['expression_level'] = pd.to_numeric(clean_expression_df['expression_level'], errors='coerce')  #Converting expression_level to numbers, inavlid set as Nan

#Removing rows with missing values
clean_expression_df = clean_expression_df.dropna()  

print("Clean expression data:")      
print(clean_expression_df.head())  
print(f"Shape after cleaning: {clean_expression_df.shape}")  #number of rows and columns after cleaning

#Merging with SNP density data
merged_df = pd.merge(snp_density_df, clean_expression_df, on='gene_id', how='left')  #Merging SNP density and expression data, keep all genes from SNP density

#Filling missing expression values with 0
merged_df['expression_level'] = merged_df['expression_level'].fillna(0)  #Setting missing expression_level values to 0

#Keeping only necessary columns
final_df = merged_df[['gene_id', 'snp_density', 'expression_level']].copy()   

print("Final merged data:")  
print(final_df.head()) 
print("Total genes:", len(final_df))  #total number of genes in final DataFrame

#Saving the results
final_df.to_csv("snp_expression.csv", index=False)         #Saving without index
print("Data saved to snp_expression.csv")  

#Summary
print("Summary:")  
print("Total genes in final dataset:", len(final_df))  #total number of genes in final dataset
print("Genes with expression data:", len(final_df[final_df['expression_level'] > 0]))  #number of genes with expression_level greater than 0
print("Genes without expression data:", len(final_df[final_df['expression_level'] == 0]))  #number of genes with expression_level equal to 0

#Statistics
print("Expression statistics:")  
print("Mean expression level: ", final_df['expression_level'].mean())  
print("Median expression level: ", final_df['expression_level'].median())  
print("Max expression level: ", final_df['expression_level'].max())  
print("Min expression level: ", final_df['expression_level'].min())  

#Many genes have expression_level = 0 due to a mismatch between transcript IDs (FBtr...) in snp_density.csv and gene IDs (FBgn...) in expression.tsv

#Visualizing SNP density and gene expression for chromosomes X and 4
#Using scatter plot and histograms to show the relationship
#Loading data
snp_expression_file = "snp_expression.csv" 
snp_expression_df = pd.read_csv(snp_expression_file)  # Loading data into a DataFrame

#Running Pearson correlation test
import scipy.stats as stats  #Importing scipy.stats for Pearson correlation test
print("Running Pearson correlation test:")  
corr_coeff, p_value = stats.pearsonr(snp_expression_df["snp_density"], snp_expression_df["expression_level"])  # Calculating Pearson correlation and p-value
print("Pearson correlation coefficient:", corr_coeff)            #displaying correlation coefficient
print("P-value:", p_value)                                   #Printing the p-value


#Creating a scatter plot of SNP density vs. expression level
plt.figure(figsize=(8, 6))                 #Setting the plot size to 8x6 inches
sns.scatterplot(x="snp_density", y="expression_level", data=snp_expression_df)               #SNP density on x-axis and expression on y-axis
plt.title("SNP Density vs. Gene Expression (Chromosomes X and 4)")                  #title
plt.xlabel("SNP Density (SNPs/bp)")                                      # Labeling x-axis
plt.ylabel("Expression Level (FPKM)")                        #Labeling y-axis
plt.savefig("snp_expression_scatter.png")               #saving the scatter plot as a PNG file
plt.close()  
print("Scatter plot saved to snp_expression_scatter.png")  

#Creating a histogram for SNP density
plt.figure(figsize=(8, 6))                         #Setting the plot size
sns.histplot(snp_expression_df["snp_density"], bins=30)       #Creating histogram for SNP density with 30 bins
plt.title("Histogram of SNP Density (Chromosomes X and 4)")     #Adding title to the plot
plt.xlabel("SNP Density (SNPs/bp)")                               #Labeling x-axis
plt.ylabel("Number of Genes")                                   #Labeling y-axis
plt.savefig("snp_density_histogram.png")                   #saving as a PNG file
plt.close()  
print("SNP density histogram saved to snp_density_histogram.png")  

#Creating a histogram for expression level
plt.figure(figsize=(8, 6))                          #Setting the plot size 
sns.histplot(snp_expression_df[snp_expression_df["expression_level"] > 0]["expression_level"], bins=30)  
plt.title("Histogram of Gene Expression (Chromosomes X and 4)")  #title 
plt.xlabel("Expression Level (FPKM)")  #x-axis
plt.ylabel("Number of Genes")  #y-axis
plt.savefig("expression_histogram.png")         #saving as PNG file
plt.close()  #
print("Expression histogram saved to expression_histogram.png")  
