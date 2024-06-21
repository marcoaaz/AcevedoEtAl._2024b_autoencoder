import numpy as np

#for saving the dimensionally reduced image
def rescale(arr, outPct):          
    min_temp = outPct
    max_temp = 100 - outPct
    
    #enhance colouring
    low_p = np.percentile(arr, min_temp) #default = 0.03 (synchrotron), 3 (EDX)
    high_p = np.percentile(arr, max_temp) #default = 99.97
    # low_p = arr.min()
    # high_p = arr.max()
    
    print(f"low p: {low_p}, high p {high_p}")
    new_arr = ((arr - low_p) * (1/(high_p - low_p) * 255)).astype('uint8')
    return new_arr