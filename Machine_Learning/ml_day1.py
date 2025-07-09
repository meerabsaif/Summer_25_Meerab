import pandas as pd             #from to import a specific module
import numpy as np              #import, to import complete library
import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight, gc_fraction, MeltingTemp as mt 
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split        
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score        #we can import more metrices 
import math 

data = pd.read_csv("dna_sequences.csv")
print(data.head())                         #doing feature extrction of data   
def extract_features(seq):                 #to find features of each sequence we need to use loop, use loop in main fn not in fn 
                                           #make a fn for one sequence, then call that function in main in form of a loop
    seq = Seq(seq.upper())
    A = seq.count('A')
    T = seq.count('T')    
    G = seq.count('G') 
    C = seq.count('C')         

    N = len(seq)

    gc_content= gc_fraction(seq)
    mol_w = molecular_weight(seq, 'DNA')
    melting_temp = mt.Tm_Wallace(seq)

    #gc skew 
    gc_skew = (G-C) / (G+C) if (G+C) != 0 else 0

    #at skew 
    at_skew = (A-T) / (A+T) if (A+T) != 0 else 0

    #shanon entropy 
    probs = {'A': A/N, "T": T/N, 'G' : G/N, 'C':C/N  }
    shannon_entropy= -sum(p * math.log2(p) for p in probs.values() if p > 0)


    return pd.Series( {
        "length": N,
        "gc_content": gc_content, 
        "mol_weight": mol_w,
        "gc_skew": gc_skew,
        "at_skew": at_skew,
        "shannon_enotropy": shannon_entropy 

    })

def main(): 
    features = data['sequence'].apply(extract_features)     #data loads all sequences, and applies the fn on all sequences 
    print(features.head())       
    # print(features["melting_temp"])
    print("\n")
#we only have features, now we will add labels 

    features = pd.concat([data['label'], features], axis=1)     #axis 1 means vertcal, axis is same as dimension 
    print(features.head())
#now our dataset has features and labels both 


#splitting x and y 
#writing x as capital, y as small

#step 1 , splitting into x and y 
    X= features.drop('label', axis=1)               #dropping labels 
    y= features['label']                           #loading labels 


#step 2 , train test split 
#input data set used for model training, split into test and training data 
#test size 0/2 means 20 percent , 20 percent data not used during training 
#random state, how our data is randomly assigned in train test split 
#v1, v2 = 1,2 

#train_test_split returns 4 values, in the same sequence

    X_train, X_test, y_train, y_test = train_test_split(X,y, test_size=0.2, random_state=42)

    #step3 , loading the model and fit the model 
    #assigning our model to a variable 
    model = RandomForestClassifier(n_estimators=100,random_state=42)
    model.fit(X_train, y_train)             #fitting model, x train= features and y train= labels of training dataset


#step 4 predict the model 
    y_pred= model.predict(X_test)            #using model for prediction, we give only features (x_test), it will predict y

#step 5, evaluating the model 
    print("Confusion Matrix")
    print(confusion_matrix(y_test, y_pred))           #actual labels, y_test. predicted labels: y_pred
    print("\nClassification report:")
    print(classification_report(y_test, y_pred))
    print("\nAccuracy: ", accuracy_score(y_test, y_pred))

if __name__ == "__main__":
    main()



# features : x 
#actual features x_train, test features x_test 


#how to check which seuences it took in training and which sequences it took for test, and also check what preditction was made against which sewuence