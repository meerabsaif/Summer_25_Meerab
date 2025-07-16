import torch
from torch.utils.data import Dataset
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class SNVDataset(Dataset):
    def __init__(self, features_df, labels):
        try:
            # Convert to torch tensors
            self.features = torch.tensor(features_df.values, dtype=torch.float32)
            self.labels = torch.tensor(labels, dtype=torch.float32)
            logging.info(f"SNVDataset created with {len(self.features)} samples")
        except Exception as e:
            logging.error(f"Error creating SNVDataset: {e}")
            raise

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        return self.features[idx], self.labels[idx]