
from torch.utils.data import Dataset, DataLoader
from PIL import Image
import numpy as np

class MyDataset(Dataset): #convert numpy to tensor
    def __init__(self, data, targets, transform=None):
        self.data = data
        # self.targets = torch.LongTensor(targets)
        self.targets = targets
        self.transform = transform
        
    def __getitem__(self, index):
        x = self.data[index]
        y = self.targets[index]   
        
        if self.transform:            
            x = Image.fromarray(self.data[index].astype(np.uint8))
            x = self.transform(x)            
        return x, y
    
    def __len__(self):
        return len(self.data)