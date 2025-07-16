import torch
import torch.nn as nn
import torch.nn.init as init
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class SNVClassifier(nn.Module):
    def __init__(self, input_dim):
        super(SNVClassifier, self).__init__()
        try:
            self.fc1 = nn.Linear(input_dim, 128)
            self.bn1 = nn.BatchNorm1d(128)
            self.fc2 = nn.Linear(128, 64)
            self.bn2 = nn.BatchNorm1d(64)
            self.fc3 = nn.Linear(64, 1)
            
            init.kaiming_uniform_(self.fc1.weight, nonlinearity='relu')
            init.kaiming_uniform_(self.fc2.weight, nonlinearity='relu')
            logging.info("SNVClassifier initialized")
        except Exception as e:
            logging.error(f"Error initializing SNVClassifier: {e}")
            raise

    def forward(self, x):
        try:
            x = self.bn1(nn.functional.relu(self.fc1(x)))
            x = self.bn2(nn.functional.relu(self.fc2(x)))
            x = self.fc3(x)
            return x
        except Exception as e:
            logging.error(f"Error in SNVClassifier forward pass: {e}")
            raise