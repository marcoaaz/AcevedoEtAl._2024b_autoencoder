
import torch

def get_device(): # get the computation device
    if torch.cuda.is_available():
        device = 'cuda:0' 
        #device = torch.cuda.get_device_name(0)
        print(f"using {device}")
    else:
        device = 'cpu'
        print(f"using {device}")
    return device       