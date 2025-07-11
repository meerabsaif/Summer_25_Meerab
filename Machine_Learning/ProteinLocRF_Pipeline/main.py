#thie file is the entry point for pipeline

#whatevr code we will run, it will run through the main file 
#like we import modules, we can also import functions from other files

import logging            #for tracking execution
from extract_features import protein_features_analysis, synthetic_labels     #importing functions from file
from model import model   #importing model fn

#configueing logging system to record information about our pipelines execution
logging.basicConfig(level=logging.INFO,               #setting up logging
                    filename= "outputlogging.log", 
                    format= '%(asctime)s - %(levelname)s - %(message)s' ),
"""
level=logging.INFO: Logs messages at info level or higher
filename=" ": Writes logs to an output file
format='%(asctime)s - %(levelname)s - %(message)s': Specifing the log format : timestamp, log level, and message.
"""

def main():                         #defining main fn
    fasta_file = "proteins.fasta"                  #importing a fasta file, this is only name of file, we loaded file in protein_feature_analysis fn from file extract features file
    #protein_features_analysis(fasta_file)         #calling fn, passing argument 
    logging.info("Features extracted successfully. Dataframe:")      #logging a message to indicate feature extraction to occur next
    dataframe = protein_features_analysis(fasta_file)                 #calling fn, passing file as argument, it processes file and stores features in the dataframe
    print(dataframe)                                          #printing dataframe

    logging.info("assigning labels")               #logging the message that label assignment is about to start 
    labeled_datset= synthetic_labels(dataframe)    #calling fn, passing dataframe. it adds random labels to df and returns modified df to a new variable "labeled dataframe"
    #model(labeled_datset)

    accuracy, classification_rep = model(labeled_datset)      #calling the model and passing labeled df to it, accuracy and classifiction report generated
    print(accuracy, classification_rep)

if __name__ == "__main__":        #chechking if script is the main program
    main()                        #calling main function to execute the pipeline. 