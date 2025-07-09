#Project 1 | Machine learning 
#Predicting Drug-like Ligands Using Molecular Descriptors and Atom Encodings
#Predicting whether a molecule is drug-like or non-drug-like using molecular descriptors (physicochemical properties), Atom Encodings and Random Forest Classifier (Machine learning Model)

#Import libraries
import pandas as pd         #to handle data
import numpy as np          #for mathematical operations
from rdkit import Chem      #for reading .sdf files and chemical properties
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors as rdm        #for computing molecular descriptors
from sklearn.ensemble import RandomForestClassifier   #importing model
from sklearn.model_selection import train_test_split,cross_val_score  #to split data
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, mean_absolute_error,mean_squared_error   #evaluation metrics
import os                   #for handling files and folders
import matplotlib.pyplot as plt      #to create plots for visualizing the output
import seaborn as sns 

#Loading .sdf files and creating a dataset with labels
def load_sdf_files(druglike, nondruglike):            #loading .sdf files from two folders
    data = []                                         #creating an empty list to store data later
    
    #Checking if (druglike) folder exists
    if not os.path.exists(druglike):             #print message if folder doesnt exists
        print("Folder:", druglike, "not found")
        return pd.DataFrame(data)                #returning empty dataframe
     
    #Reading all .sdf files present in the folder and labelling 
    for file in os.listdir(druglike):                    #looping through all files in the folder
        if file.endswith(".sdf"):                        #only considering files ending with .sdf which points that file is in sdf format      
            file_path = os.path.join(druglike, file)     #building complete path of file
            compound_name = file.replace(".sdf", "")     #getting molecule name from filename, also removing .sdf 
            data.append({"file_path": file_path, "Molecule": compound_name,"label": 1})       #label 1 for drug like, adding dictionary to data list
    
    #Checking if nondruglike folder exists
    if not os.path.exists(nondruglike):
        print("Folder", nondruglike, "not found")
        return pd.DataFrame(data)
    
    #Reading non-drug-like .sdf files and labelling
    for file in os.listdir(nondruglike):
        if file.endswith(".sdf"):
            file_path = os.path.join(nondruglike, file)
            compound_name = file.replace(".sdf", "")  
            data.append({"file_path": file_path, "Molecule": compound_name,"label": 0})     #label 0 for non-drug like
    
    #Converting list of dictionaries above to a dataframe
    return pd.DataFrame(data)

#Extracting molecular descriptors from .sdf files
#molecular descriptors: physiochemical properties of a compound 
def extract_descriptors(file_path):                   #defining a fn, giving file_path as argument 
    #Reading the .sdf files
    mol = Chem.MolFromMolFile(file_path)              #loading molecule from file using rdkit
    if mol is None:         #check if molecule is valid
        return None         #to skip invalid files
    #Calculating descriptors and atom encodings 
    mol_weight = Descriptors.MolWt(mol)                #Molecular weight
    log_p = Descriptors.MolLogP(mol)                   #Hydrophobicity
    h_donors = Descriptors.NumHDonors(mol)             #Hydrogen bond donors
    h_acceptors = Descriptors.NumHAcceptors(mol)       #Hydrogen bond acceptors
    tpsa = Descriptors.TPSA(mol)                       #Topological polar surface area
    num_aromatic = len(mol.GetSubstructMatches(Chem.MolFromSmarts("c")))  #counting aromatic carbons
    num_rings = mol.GetRingInfo().NumRings()           #counting no. of rings
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)  #Rotatable bonds
    fsp3 = rdm.CalcFractionCSP3(mol)                   #Fraction of sp3 carbons
    heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [6, 1])  #Non-carbon, non-hydrogen atoms
    bertz_ct = Descriptors.BertzCT(mol)                    #Molecular complexity
    # syn_access = sascorer.calculateScore(mol) if 'sascorer' in globals() else 0  #Synthetic accessibility
    heavy_atoms = mol.GetNumHeavyAtoms()               #Non-hydrogen atoms


    #Calculating Lipinski violations and checking if molecule violates it 
    lipinski_violations = 0                 #counter for calculating no. of rule violations
    if mol_weight > 500:
        lipinski_violations += 1
    if log_p > 5:
        lipinski_violations += 1
    if h_donors > 5:
        lipinski_violations += 1
    if h_acceptors > 10:
        lipinski_violations += 1

    #Returning descriptors as a series
    #later will be converted in a dataframe
    #column names: values 
    return pd.Series({
        "Molecule": "",   
        "Aromatic": num_aromatic,
        "Rings": num_rings,
        "LogP": log_p * 100,  
        "HBD": h_donors * 10,  
        "HBA": h_acceptors * 10, 
        "TPSA": tpsa * tpsa,  
        "Mol_Weight": mol_weight * 10,  
        "RB": rotatable_bonds * 10, 
        "Lipinski_Violations": lipinski_violations,
        "Fsp3": fsp3,
        "Heteroatoms": heteroatoms,
        "BertzCT": bertz_ct,
        # "SynAccess": syn_access,
        "HeavyAtoms": heavy_atoms
    })


#alternative to try next time: "rdMolDescriptors.CalcDescriptors" gives all descriptors at once.

#defining main function
def main():
    #Loading data
    print("Loading the files:")
    data = load_sdf_files("drug_like_compounds", "non_drug_like_compounds")  #loading both foalders 
    if data.empty:                               #if foalders are empty print next message
        print("No data is loaded.")
        return
    print("Head of dataset:")
    print(data.head(5)) 
    print("Total Molecules:", len(data))         #gives total no. of molecules
    
    #Extracting descriptors for each compound
    print("Extracting molecular descriptors:")
    features = data["file_path"].apply(extract_descriptors)     #applying extract descriptors function to every line, storing it in a variable called features
    #after this step data has been converted into a dataframe, data is a dataframe

    #Removing rows with none or zero values 
    features = features.dropna()                  #removing rows with null values 
    data = data.loc[features.index]          #features.index gives row index of valid molecules
    #data.loc[features.index] : selecting rows from data whose indexes are in features.index
    #data = data.loc[features.index] : overwriting data to align with features
    #why aligning with features? ans: later we’ll add labels back to features, only works if both have matching rows and indexes (features["Drug"] = data["label"])
    #this above code prevented from mismatched rows, which would have caused errors

    #Adding molecule names and labels to features
    features["Molecule"] = data["Molecule"]       #adding new molecule column to features, to keep track of of the information of each mol
    features["Drug"] = data["label"]              #adding a drug column to the dataset. Value comes from label column from data
    #label: it is our target variable. 
    
    print("Dataset with features and labels:")
    print(features.head(10))                     #displaying dataset to see if our df is proper for applying model and also analyzing it 
    print("Class distribution:")                 #printing how many drug like and non drug like molecules are present
    print(features["Drug"].value_counts())
    
    #preparing data for machine learning 
    #Splitting the dataset into features (X) which are input for model and labels (y) which are target values that our model will predict
    X = features.drop(["Molecule", "Drug"], axis=1)  #Features (X) , the convention is to use X as capital
    #removing these two columns (axis=1 means column, 0 means rows)
    #now x only contains numerical features 
    #molecule name was not useful for predictions
    #it is important to remove target variable from the input

    y = features["Drug"]               #Labels (y) , we use y in lowercase
    #selecting drug column from features, it contails our labels 1 and 0 
    #labeling showing that it is supervised learning

    #Performing 5-fold cross-validation   
    ##to ensure model generalizes to unseen molecules and not just memorizes the datset               
    #it tests models performance, avoids overfitting 
    #more explanation in the notebook
    print("Performing 5-fold cross-validation:")
    model = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
    #we applied the machine learning model to our model
    #why used in notebook, parameters also explained in the notebook

    cv_scores = cross_val_score(model, X, y, cv=5, scoring="accuracy")           #Running 5-fold cross-validation
    """"
    #Spliting dataset into 5 equal parts ("folds").
    For each fold: Using 4 folds for training, 1 for testing.
    Training model and measuring accuracy(metrics) on test fold. Used accuracy as metric for propotion of correct predictions.
    Repeating process 5 times, each fold is used as test data once.
    Returning accuracy for each iteration.
    """
    print("Cross-validation accuracies:", cv_scores)     #printing individual score
    print("Mean CV accuracy:", cv_scores.mean())         #taking average of all scores 
    #eg: average, model predicts correctly ~87.8% time on unseen data

    #Split into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    #X_train, features for training, X_test	features for testing
    #y_train, labels for training, y_test abels for testing (evaluation)
    #random_state=42  (split is reproducible (same split every run) explore more later
   
    #Save test indices and file paths for tracking row in original dataset corresponding to test set
    test_indices = X_test.index          #getting original row index of test set from x_test
    test_files = data.loc[test_indices, ["file_path", "Molecule"]]    #looking at molecule and filepath, and saving to test_files, later will use for prediction w molecule names and file paths
    
    #Training the model
    print("Training the model:")
    #model = RandomForestClassifier(n_estimators=100, random_state=42)      #defined earlier in cross validation, so commented here
    model.fit(X_train, y_train)                    #Building multiple decision trees using X_train (input features) and y_train (labels)
    #Learning patterns to predict molecule druglikeness later 
    
    #Predicting 
    print("Making predictions:")
    y_pred = model.predict(X_test)          #making prediction features for testing which were x_test
    #later we will compare with actual values : y_test

    #Saving predictions with file paths to a dataframe 
    #we will use this df for further analysis and saving
    results = pd.DataFrame({
        "file_path": test_files["file_path"],       #actual file path for test sample
        "Molecule": test_files["Molecule"],         #molecule name for test sample
        "actual_label": y_test,                     #true labels
        "predicted_label": y_pred                   #lables model predicted now
    })

    results.to_csv("predictions.csv", index=False)     #we are saving to csv 
    print("Predictions saved to predictions.csv")
    
    #Evaluation metrics
    print("Calculating evaluation metrics:")       
    print("Confusion Matrix:")                          #table showing TP, TN, FP, FN 
    cm = confusion_matrix(y_test, y_pred)
    print(cm) 
    print("Classification Report:")                     #gives detailed metrics
    print(classification_report(y_test, y_pred))
    print("Accuracy:", accuracy_score(y_test, y_pred))
    print("Mean Absolute Error:", mean_absolute_error(y_test, y_pred))
    #Measures average absolute difference between predicted labels and true labels.
    #counts how many were wrong divided by total

    print("Mean Squared Error:", mean_squared_error(y_test, y_pred))
    #Similar to MAE but squares errors 
    #In classification, not super meaningful as errors are 0 or 1, but still computed
    
    #Visualizing confusion matrix
    plt.figure(figsize=(6, 4))         #creating figure
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=["Non-Drug", "Drug"], yticklabels=["Non-Drug", "Drug"])
    #annot=True: Show numbers on cells, fmt="d": Format no as int, xticklabels or yticklabels: Labels for axes.

    plt.title("Confusion Matrix")       #adding title
    plt.xlabel("Predicted")             #labeling x axis
    plt.ylabel("Actual")                #labeling y axis
    plt.savefig("confusion_matrix.png")         #saving as image
    plt.close()                                 #closing to prevent from mix w next plot
    print("Confusion matrix is saved to confusion_matrix.png")

    #Visualizing feature importance
    feature_names = X.columns                 #x.columns are names of our mol descriptors
    importances = model.feature_importances_      #model.feature_importances_ is Random Forest built-in importance scores for each feature
    plt.figure(figsize=(10, 12))                  #setting fig size
    plt.bar(feature_names, importances, color="blue")       
    plt.title("Feature Importance")
    plt.xlabel("Features")
    plt.ylabel("Importance")
    plt.xticks(rotation=60)            #rotating label tags
    plt.savefig("feature_importance.png")
    plt.close()
    print("Feature importance is saved to feature_importance.png")

    #Save the complete dataset
    features = features[["Molecule", "Drug", "Mol_Weight", "LogP", "HBD", "HBA", "TPSA", "RB", 
                     "Lipinski_Violations", "Aromatic", "Rings", "Fsp3", "Heteroatoms", 
                     "BertzCT", "HeavyAtoms"]]    #mentioning columns we want to save from features
    features.to_csv("finaldataset.csv", index=False)              #saved to csv
    print("Complete dataset saved to finaldataset.csv")

    
#Running the program
if __name__ == "__main__":
    main()