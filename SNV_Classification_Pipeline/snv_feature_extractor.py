import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class SNVFeatureExtractor:
    def extract_features(self, df):
        try:
            features = pd.DataFrame(index=df.index)
            features['pos'] = pd.to_numeric(df['POS'], errors='coerce').fillna(0)
            
            # One-hot encode REF and ALT alleles
            for base in ['A', 'T', 'C', 'G']:
                features[f'ref_{base}'] = (df['REF'] == base).astype(int)
                features[f'alt_{base}'] = (df['ALT'] == base).astype(int)

            # Basic variant type features
            features['is_snv'] = ((df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)).astype(int)
            
            # Check for NaN values
            if features.isna().any().any():
                logging.warning("NaN values detected in features, filling with 0")
                features = features.fillna(0)
                
            logging.info(f"Extracted features: {features.columns.tolist()}")
            return features
        except Exception as e:
            logging.error(f"Error extracting features: {e}")
            raise

    def process_dataset(self, df):
        return self.extract_features(df)