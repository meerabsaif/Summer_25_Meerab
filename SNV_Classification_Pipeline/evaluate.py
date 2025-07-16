
import torch
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, confusion_matrix, classification_report
import matplotlib.pyplot as plt
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def evaluate_model(model, test_loader):
    try:
        model.eval()
        predictions = []
        actuals = []
        probabilities = []
        
        with torch.no_grad():
            for features, labels in test_loader:
                outputs = model(features).squeeze()
                probs = torch.sigmoid(outputs)
                preds = (probs >= 0.5).float()
                
                predictions.extend(preds.tolist())
                actuals.extend(labels.tolist())
                probabilities.extend(probs.tolist())
        
        accuracy = accuracy_score(actuals, predictions)
        precision = precision_score(actuals, predictions, zero_division=0)
        recall = recall_score(actuals, predictions, zero_division=0)
        f1 = f1_score(actuals, predictions, zero_division=0)
        roc_auc = roc_auc_score(actuals, probabilities)
        cm = confusion_matrix(actuals, predictions)
        fpr, tpr, _ = roc_curve(actuals, probabilities)
        
        # Generate classification report
        class_report = classification_report(actuals, predictions, target_names=['Benign', 'Pathogenic'], output_dict=True, zero_division=0)
        logging.info(f"Classification Report:\n{classification_report(actuals, predictions, target_names=['Benign', 'Pathogenic'], zero_division=0)}")
        
        # Create and save classification report as an image
        report_df = pd.DataFrame(class_report).transpose()
        plt.figure(figsize=(8, 4))
        plt.table(cellText=report_df.round(4).values,
                  rowLabels=report_df.index,
                  colLabels=report_df.columns,
                  cellLoc='center',
                  loc='center')
        plt.axis('off')
        plt.title('Classification Report')
        plt.savefig('classification_report.png', bbox_inches='tight', dpi=300)
        plt.close()
        logging.info("Classification report saved as classification_report.png")
        
        logging.info(f"Test Accuracy: {accuracy:.4f}, Precision: {precision:.4f}, Recall: {recall:.4f}, F1: {f1:.4f}, ROC-AUC: {roc_auc:.4f}")
        
        return {'accuracy': accuracy, 'precision': precision, 'recall': recall, 'f1': f1, 'roc_auc': roc_auc, 'cm': cm, 'fpr': fpr, 'tpr': tpr}
    
    except Exception as e:
        logging.error(f"Error in evaluate_model: {e}")
        raise
