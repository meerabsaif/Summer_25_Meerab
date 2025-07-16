import matplotlib.pyplot as plt
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def plot_roc_curve(fpr, tpr, roc_auc, filename='roc_curve.png'):
    try:
        plt.figure()
        plt.plot(fpr, tpr, label=f'ROC Curve (area = {roc_auc:.2f})')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic')
        plt.legend(loc="lower right")
        plt.savefig(filename)
        plt.close()
        logging.info(f"ROC curve saved as {filename}")
    except Exception as e:
        logging.error(f"Error plotting ROC curve: {e}")
        raise

def plot_confusion_matrix(cm, filename='confusion_matrix.png'):
    try:
        plt.figure()
        sns.heatmap(cm, annot=True, fmt='d')
        plt.title('Confusion Matrix')
        plt.ylabel('Actual Label')
        plt.xlabel('Predicted Label')
        plt.savefig(filename)
        plt.close()
        logging.info(f"Confusion matrix saved as {filename}")
    except Exception as e:
        logging.error(f"Error plotting confusion matrix: {e}")
        raise


def plot_feature_importance(ranked_features, top_n=10, filename='feature_importance.png'):
    try:
        features, importances = zip(*ranked_features[:top_n])
        
        plt.figure(figsize=(10, 6))
        plt.barh(range(top_n), importances, align='center')
        plt.yticks(range(top_n), features)
        plt.xlabel('Absolute Pearson Correlation')
        plt.title(f'Top {top_n} Feature Importances')
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        logging.info(f"Feature importance plot saved as {filename}")
    
    except Exception as e:
        logging.error(f"Error plotting feature importance: {e}")
        raise