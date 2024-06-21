import torch
import torch.nn as nn
import torchvision.transforms as transforms
import torch.nn.functional as F
import torch.optim as optim

# % Define the autoencoder model (DSA has a symmetric layout)        
class SparseAutoencoder(nn.Module):
    def __init__(self, n_channels):
        super(SparseAutoencoder, self).__init__()
    
        # encoder
        self.enc1 = nn.Linear(in_features=n_channels, out_features=15, bias=True)
        self.enc2 = nn.Linear(in_features=15, out_features=10, bias=True)
        self.enc3 = nn.Linear(in_features=10, out_features=3, bias=True)        
    
        # decoder 
        self.dec1 = nn.Linear(in_features=3, out_features=10, bias=True)
        self.dec2 = nn.Linear(in_features=10, out_features=15, bias=True)
        self.dec3 = nn.Linear(in_features=15, out_features=n_channels, bias=True)        
    
    def forward(self, x):
        # encoding
        x = F.sigmoid(self.enc1(x))
        x = F.sigmoid(self.enc2(x))
        x = F.sigmoid(self.enc3(x))        
    
        # decoding
        x = F.sigmoid(self.dec1(x))
        x = F.sigmoid(self.dec2(x))
        x = F.sigmoid(self.dec3(x))
                
        return x