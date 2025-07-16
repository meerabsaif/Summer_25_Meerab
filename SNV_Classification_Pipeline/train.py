import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Subset
from snv_dataset import SNVDataset
from model import SNVClassifier
from sklearn.model_selection import KFold
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def train_model(features_df, labels, batch_size=32, lr=1e-3, epochs=5, n_splits=5):
    input_dim = features_df.shape[1]
    logging.info(f"DEBUG: input_dim = {input_dim}, type = {type(input_dim)}")
    
    try:
        #Initializing k-fold cross-validation: 
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        best_model = None
        best_loss = float('inf')
        best_params = {"input_dim": input_dim, "hidden_dim": 128, "lr": lr}
        
        for fold, (train_idx, val_idx) in enumerate(kf.split(features_df)):
            logging.info(f"Starting fold {fold+1}/{n_splits}")
            
            #Creating train and validation datasets
            train_dataset = SNVDataset(features_df.iloc[train_idx], labels[train_idx])
            val_dataset = SNVDataset(features_df.iloc[val_idx], labels[val_idx])
            
            train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
            val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
            
            #Initializing model
            model = SNVClassifier(input_dim).to('cpu')  # Explicitly use CPU for simplicity
            optimizer = torch.optim.Adam(model.parameters(), lr=lr)
            criterion = nn.BCEWithLogitsLoss()
            
            #Training loop
            for epoch in range(epochs):
                model.train()
                total_train_loss = 0.0
                for i, (X, y) in enumerate(train_loader):
                    try:
                        optimizer.zero_grad()
                        preds = model(X)
                        loss = criterion(preds.squeeze(), y)
                        loss.backward()
                        optimizer.step()
                        total_train_loss += loss.item()
                        if i % 100 == 0:
                            logging.info(f"Fold {fold+1}, Epoch {epoch+1}, Batch {i}, Loss: {loss.item():.4f}")
                    except Exception as e:
                        logging.error(f"Error in training batch {i}: {e}")
                        raise
                
                #Validation
                model.eval()
                total_val_loss = 0.0
                with torch.no_grad():
                    for X, y in val_loader:
                        preds = model(X)
                        loss = criterion(preds.squeeze(), y)
                        total_val_loss += loss.item()
                
                avg_val_loss = total_val_loss / len(val_loader)
                logging.info(f"Fold {fold+1}, Epoch {epoch+1}/{epochs}, Train Loss: {total_train_loss/len(train_loader):.4f}, Val Loss: {avg_val_loss:.4f}")
                
                #saving the best model
                if avg_val_loss < best_loss:
                    best_loss = avg_val_loss
                    best_model = model
                    best_params['fold'] = fold + 1
                    best_params['val_loss'] = avg_val_loss
            
        logging.info(f"Best model from fold {best_params['fold']} with validation loss: {best_params['val_loss']:.4f}")
        return best_model, val_loader, best_params
    
    except Exception as e:
        logging.error(f"Error in train_model: {e}")
        raise
