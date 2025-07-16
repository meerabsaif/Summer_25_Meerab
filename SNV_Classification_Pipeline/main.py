import logging
import pandas as pd
from snv_feature_extractor import SNVFeatureExtractor
from train import train_model
from evaluate import evaluate_model
from visualize import plot_roc_curve, plot_confusion_matrix

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)')

def prepare_data(csv_path, sample_size=None):
    try:
        df = pd.read_csv(csv_path, low_memory=False)
        if sample_size:
            df = df.sample(n=sample_size, random_state=42)
        logging.info(f"Loaded {len(df)} variants from {csv_path}")
        extractor = SNVFeatureExtractor()
        features_df = extractor.process_dataset(df)
        labels = df['label'].astype(float).values
        logging.info(f"Features shape: {features_df.shape}, Labels shape: {labels.shape}")
        logging.info(f"Any NaN in features: {features_df.isna().any().any()}")
        logging.info(f"Feature columns: {features_df.columns.tolist()}")
        return features_df, labels
    except Exception as e:
        logging.error(f"Error preparing data: {e}")
        raise

if __name__ == "__main__":
    try:
        csv_path = 'data/clinvar_processed.csv'
        # Use a smaller sample for testing to avoid memory issues
        features_df, labels = prepare_data(csv_path, sample_size=10000)
        model, test_loader, best_params = train_model(features_df, labels, batch_size=64, lr=1e-3, epochs=5, n_splits=3)
        metrics = evaluate_model(model, test_loader)
        plot_roc_curve(metrics['fpr'], metrics['tpr'], metrics['roc_auc'])
        plot_confusion_matrix(metrics['cm'])
        logging.info(f"Best hyperparameters: {best_params}")
    except Exception as e:
        logging.error(f"Error in main pipeline: {e}")
        raise