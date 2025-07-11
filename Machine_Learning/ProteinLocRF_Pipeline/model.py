# #from sklearn pipeline import pipeline 
# from sklearn.svm import SVR 

#here we will implement the machine learning model, and evaluate using random forest classifier

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler           #StandardScaler: Standardizes features to have zero mean and unit variance, this improves model performance
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report

def model(df):        #defining a function model that takes df with features and labels
    #step 1: seperating features
    X = df.drop(columns=['labels'])      #removing label column from x
    y = df['labels']      #extracting label column as the target variable


    #train test split data
    #splitting data into 80 : 20 ratio fro evaluating models performance on unseen data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    #training (X_train, y_train) and testing (X_test, y_test) sets. test_size=0.2 (20% data for testing) random_state=42 (reproducibility)

#pipepline code can be inserted here and we can move direct to prediction

    #preprocessing step: standardization of features to ensure they are contributing equally
    scaler = StandardScaler()       #initializing a standard scalar object
    
    #First transform
    X_train_scaled_np = scaler.fit_transform(X_train)     #this fits the scalar to the traning data and transform x_train to a scaled numpy array
    X_test_scaled_np = scaler.transform(X_test)           #uses training data's mean and std to the test data and ensure consistency
    
    #converts the scaled numpy arrays back dataframe, preserving column names
    X_train_scaled = pd.DataFrame(X_train_scaled_np, columns=X.columns)   #creates dataframe from scaled training data, we used original feature names
    X_test_scaled = pd.DataFrame(X_test_scaled_np, columns=X.columns)

    #initializing and training random forest classifier on scaled training data
    model = RandomForestClassifier(n_estimators=100, random_state=42)      #making 100 different trees
    model.fit(X_train_scaled, y_train)  #training the model on scaled training features and labels 

    #above is the pipeline code i was referring to

    # p = Pipeline({
    #     ('scaler', StandardScaler()),
    #     ('svr', SVR(kernel='linear')),
    #     ('rf', RandomForestClassifier(n_estimators=100, random_state=42))

    # })

    # p.fit_transform(X_train, y_train)

    #instead of complete block of code. After data splitting till here we can use pipline, short and simple

    #using trained model to predict labels for the test set
    y_pred = model.predict(X_test_scaled)      #this generates prediction for x test scaled and resturns an array of predicted laebls

    #evaluating models performance 
    accuracy = accuracy_score(y_test, y_pred)     #propotio of correct predictions
    classification_rep = classification_report(y_test, y_pred, zero_division=1)     #generates a report

    #zero divison =1 ,  we set the ouput to 1, for the mertic when class has zero instances. This avoids division by zero errors
    
    return accuracy, classification_rep            #returns accuracy and classification report
