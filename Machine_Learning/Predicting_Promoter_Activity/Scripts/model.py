#importing libraries
from sklearn.ensemble import RandomForestClassifier   #for the main model
from sklearn.model_selection import train_test_split    #tool to split data
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix   #for checking how good the model is
import pandas as pd     #for handling tables of data (dataframes)
import matplotlib.pyplot as plt    #for makin plots
import seaborn as sns     #also for plots


#function to determine if gene expression is High or Low.
def get_activity_label(expression_value): #a function to sort genes
    if expression_value >= 1:          #if the expression number is 1 or bigger
        return 'High'       #we call it 'High'
    else:            #otherwise
        return 'Low'     #it's 'Low'


#defining a function to train the model and fills missing values in the feature df
def train_model(feature_df):     #this is the main function, it takes our feature table
    #replacing nan w 0 and modifying df in place
    feature_df.fillna(0, inplace=True)        #if any cells are empty, just put a 0 there
    
    #Using function above to create the label column
    #it runs the fn on every row in the 'expression' column
    feature_df['label'] = feature_df['expression'].apply(get_activity_label)

    #X has all the features for training, we drop columns we dont need for prediction
    X = feature_df.drop(columns=['gene', 'expression', 'label'])
    y = feature_df['label']    #y is the answer column ('High' or 'Low')
    genes = feature_df['gene']    

    #Splitting data
    #this splits everything into a training set and a testing set
    #test_size=0.2 means 20% of data is for testing, 80% for training
    # andom_state=42 means it splits the same way every time we run it
    X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
        X, y, genes, test_size=0.2, random_state=42
    )

    #Train classification model
    #we're making a model called Random Forest, it uses 100 'decision trees'
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)        #here's where the model actually learns from the training data

    #Predicting
    y_pred = model.predict(X_test)      #now we ask the trained model to make guesses on the test data

    #Evaluating
    acc = accuracy_score(y_test, y_pred)    #checking how many it got right (the accuracy)
    print(f"Accuracy: {acc:.4f}") 
    
    #Generateing the classification report 
    class_report = classification_report(y_test, y_pred)      #makes a detailed report
    print("\nClassification Report:\n", class_report)     #displaying the report on the screen
    
    #Save the report to a text file
    with open("classification_report.txt", "w") as report_file:      #openning a new file to write in
        report_file.write("Classification Report\n")   #put a title 
        report_file.write(class_report)      #write the whole report string to a file
    print("Classification report saved to classification_report.txt") 


    #Combining predictions with genes
    test_features = X_test.copy()     #making a copy of the test features table
    test_features["gene"] = genes_test.values #adding the gene names back
    test_features["actual_label"] = y_test.values   #adding the actual answers column
    test_features["predicted_label"] = y_pred     #adding a column for what our model predicted 

    #Reordering columns
    #this just makes the table look nicer, with the important columns first
    cols = ["gene", "actual_label", "predicted_label"] + [
        col for col in test_features.columns if col not in ["gene", "actual_label", "predicted_label"]
    ]
    test_features = test_features[cols]    #applying the new column order

    #Saving predictions
    #saves the final table to a csv file, so we can look at it in excel
    test_features.to_csv("predicted_gene_activity.csv", index=False)
    print("Predicted labels saved to predicted_gene_activity.csv") #

    return model       #returning the trained model 