# -*- coding: utf-8 -*-
"""
Predict RGB image (embedded space) with Deep sparse autoencoder. 
This script loads the DSA model and save an image mosaic

This is the update of the second half (predicting DSA) of the second script for 
dimensionality reduction in Acevedo Zamora et al. (2024). 
The method followed Thomas et al. (2016, 2017, 2022) papers.
The implementation uses CPU 'multi-thread' with low GPU (if available) consumption
The new 'DSA_image_training_v2.py' script must be ran beforehand.  

Cite as: https://doi.org/10.1016/j.chemgeo.2024.121997

Written by: Marco Andres, ACEVEDO ZAMORA
Created on Tue Apr 12 11:13:17 2022
Published (first version): 14-Nov-23
Updated (second version): 29-May-24

Followed sources:

#Implementation note 1:
#If pyvips import is taken outside '__main__':
#(process:69204): GLib-GIO-WARNING **: 15:47:44.631: 
#Unexpectedly, UWP app `Microsoft.OutlookForWindows_1.2023.719.200_x64__8wekyb3d8bbwe' 
#(AUMId `Microsoft.OutlookForWindows_8wekyb3d8bbwe!Microsoft.OutlookforWindows') 
#supports 1 extensions but has no verbs

"""
#!/usr/bin/python3

#Dependencies

import os
import argparse
import time
import glob
import re
import pickle

import numpy as np
import h5py

import matplotlib
import matplotlib.pyplot as plt 
matplotlib.style.use('ggplot')

import torch
import torch.nn as nn
import torchvision.transforms as transforms
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader

from PIL import Image
from tqdm import tqdm

#Developed functions (see folder)
from MyDataset import MyDataset
from SparseAutoencoder import SparseAutoencoder
from rescale import rescale
from get_device import get_device

#region
if __name__ == '__main__':    #lock to only 1 subprocess (Windows) 
    import pyvips #see Note 1    

    #region
    ## User input
    save_recoloured = 1 #save image tiles
    save_stitched = 1 #save mosaic
    percentOut_dsaImage = 5 #percentile out in the output (improves DSA image colour contrast); default=1
    calc_depth = 2**16-1 #for precision 8-bit=255; 16-bit=65535

    # fileDest = f"../{destDir1}/{descrip_pred1}_{imageName_pred1}_sessionVariables.pckl"    
    workingDir = "D:\\Chris_Collaboration\\2024_Ioan_Purdys_reward\\XFM\\153874_P2\\tiff\\Au\\recoloured_trial2_pctOut0.1"

    #input folder from training script
    tag_combo_input = '20Jun24_test1' #folder with trained DSA model (from 'DSA_images_training_v3.py')
    metadataFile_pckl = os.path.join(workingDir, 'dsa_results\\'+  tag_combo_input + 
                                '\\sessionVariables.pckl')
    parentDir = os.path.dirname(metadataFile_pckl)

    #output folder name
    outputFolder2 = os.path.join(parentDir, 'dsa_tiles')
    fileName_mosaic = 'DSA_rgb_pctOutput' + str(percentOut_dsaImage) + '.tif' #stitched image

    #For loading tiles (choose linear or log-transformed pyramid)
    imageDir = os.path.join(workingDir, 'recoloured_pyramid_files\\0')

    #endregion
    
    #region
    ## Helper functions
    def make_dir(destDir):
        image_dir = destDir
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)  

    #Define colour map (falsecolour images)
    img_indexes = pyvips.Image.identity()
    lut = img_indexes.falsecolour() #using standard heatmap
    #256x1 uchar, 3 bands, srgb, pngload

    ## Script
    os.chdir(workingDir)
    print("Current working directory: {0}".format(os.getcwd()))     
    make_dir(outputFolder2) #saving dir   

    device = get_device()

    # Prediction functions
    def embedded_space(images, model_children):
        values = images    
        for i in range(3):
            values = model_children[i](values)
            values = F.sigmoid(values) # sigmoid because we need probability distributions            
        return values
    
    def predict_space(model, dataloader):

        model_children = list(model.children()) # get the layers       

        print('Predicting')
        model.eval()   #inferencing mode
        counter = 0
        with torch.no_grad():
            pixel_batch = []
            for i, data in tqdm(enumerate(dataloader), total=int(len(all_set)/dataloader.batch_size)):
                counter += 1
                img, _ = data
                img = img.to(device)
                img = img.view(img.size(0), -1)
                
                output = embedded_space(img, model_children)
                pixel_batch.append(output)
            # result = torch.cat(pixel_batch, dim=0)
        return pixel_batch           

    #Load image processing and PyTorch model metadata (from 'DSA_images_training_v2.py')    
    f = open(metadataFile_pckl, 'rb')
    obj = pickle.load(f)
    f.close()    
    [workingDir, metadataFile_mat, outputFolder, 
     outPct, LEARNING_RATE, RHO,
     epoch_default, alpha_reg, betha_reg, fraction,
     test_ratio, n_workers, BATCH_SIZE, n_workers_pred,
     BATCH_SIZE_pred, EPOCHS, ALPHA, BETA, ADD_SPARSITY] = obj        
    
    #Overwrite parameters
    BATCH_SIZE_pred = 8192

    # Understanding dataset (tileset)    
    f = h5py.File(metadataFile_mat, 'r') #info about montage and selected stack layers
    metadata_struct = f.get('s')    
    imageHeight = metadata_struct['imageHeight'][()][0]
    imageWidth = metadata_struct['imageWidth'][()][0]
    idx_in = metadata_struct['idx_in'][()].flatten() #same order as 'descriptiveStats.csv'
    n_layers_input = np.sum(idx_in, axis= 0)    
    idx_in_bool = idx_in.astype(bool)  

    #scan tileset
    pattern = re.compile(r".*\\(\d+)_(\d+)\.tif") #Windows
    #pattern = re.compile(r".*/(\d+)_(\d+)\.tif") #Linux
    max_x = 0
    max_y = 0
    for filename in glob.glob(f"{imageDir}/*_*.tif"):
        match = pattern.match(filename)
        if match:        
            max_x = max(max_x, int(match.group(1)))
            max_y = max(max_y, int(match.group(2)))
    
    print(list(metadata_struct.keys()))
    print(f"mosaic of WxHxC= {imageWidth}x{imageHeight}x{n_layers_input} px")
    print(f"mosaic of {max_x + 1}x{max_y + 1} tiles")    

    #Load model
    modelPATH = f"{outputFolder}/model_epochs{EPOCHS}.tar"
    #modelPATH = 'E:\\paper 3_datasets\\data_DMurphy\\91712-81R5w-glass-quadrant cpx\\tiff\\15-Nov-23\\DSA_test1_model_5.tar'    
    #Option: If wanting a manual input (transfer model from other dataset). It must be *.tar not *.pth    

    #Initialize
    n_channels = n_layers_input
    checkpoint = torch.load(modelPATH) #dictionary containing objects
    model_pred = SparseAutoencoder(n_channels=n_channels).to(device)
    optimizer_pred = optim.Adam(model_pred.parameters(), lr=LEARNING_RATE)    
     
    model_pred.load_state_dict(checkpoint['model_state_dict'])
    optimizer_pred.load_state_dict(checkpoint['optimizer_state_dict'])
    epoch = checkpoint['epoch']
    loss = checkpoint['loss']         
    
    model_pred.eval() #set dropout and batch normalisation to evaluation mode
    # model_pred = deepcopy(model_pred.state_dict()) #prevents retraining    

    #Load tiles and resample
    tiles_across = max_x + 1
    tiles_down = max_y + 1
    n_tiles = tiles_down*tiles_across
    
    # k = 0
    # #tiles = []
    # for y in range(0, tiles_down):
    #     for x in range(0, tiles_across):
    #         k += 1

    #         fileName_expression = f"{x}_{y}.tif"
    #         fileName_input = os.path.join(imageDir, fileName_expression)
    #         fileName_output = os.path.join(outputFolder2, fileName_expression)
    #         fileName_parse = f"{fileName_input}"
    #         print(f"Processing tile {k}/{n_tiles}: {fileName_expression}")   

    #         image_temp = pyvips.Image.new_from_file(fileName_parse)
    #         image_np = image_temp.numpy()
    #         dim_shape = image_np.shape  
    #         image_xyz = np.reshape(image_np, (-1, dim_shape[2])) #39 channels
    #         n_rows = image_xyz.shape[0]                            

    #         #Full image tile
    #         data_full = image_xyz[:, idx_in_bool]
    #         #n_rows_sel = data_full.shape[0]
    #         #n_channels = data_full.shape[1]     

    #         #Predict embedding space (data = target data)
    #         #Data transformations 
    #         loader = transforms.Compose([ 
    #             transforms.ToTensor(),
    #             transforms.Normalize((0.5,), (0.5,)), #normalisation for faster training
    #             ])
    #         all_set = MyDataset(data_full, data_full, transform= loader)

    #         alldata_loader = DataLoader(all_set, batch_size= BATCH_SIZE_pred,
    #                                     shuffle= False, num_workers= n_workers_pred)   

    #         pixel_batch = predict_space(model_pred, alldata_loader)    
    #         #print(pixel_batch[0].shape)
    #         result = torch.cat(pixel_batch, dim=0)    
    #         outputs = result.view(dim_shape[0], dim_shape[1], 3).cpu().numpy() #double

    #         #Saving recoloured images (for retrospective feedback): time-consuming        
    #         if save_recoloured == 1:
    #             image_output = pyvips.Image.new_from_array(outputs)                            
    #             image_output.write_to_file(fileName_output) 

    #         #tiles.append(outputs)

    #endregion

    #region
    
    k = 0
    tiles = []
    for y in range(0, tiles_down):
        for x in range(0, tiles_across):
            k += 1

            fileName_expression = f"{x}_{y}.tif"
            fileName_input = os.path.join(outputFolder2, fileName_expression)            
            fileName_parse = f"{fileName_input}"
            print(f"Processing tile {k}/{n_tiles}: {fileName_expression}")   

            image_temp = pyvips.Image.new_from_file(fileName_parse)            

            tiles.append(image_temp)

    #stitching; DSA values from [0-1] with 0.5 mean
    image_stitched = pyvips.Image.arrayjoin(tiles, across=tiles_across)    
    n_bands = image_stitched.bands
    
    #crop background borders (optional)
    left, top, width, height = image_stitched.find_trim(threshold=0.001, background=[0])
    image_cropped = image_stitched.crop(left, top, width, height) #modify accordingly    
        
    #Descr. statistics (follows 'tilingAndStacking_v3.py')
    out = pyvips.Image.stats(image_cropped)
    out1 = out.numpy()    
    #print(out1)

    r, g, b = image_cropped.bandsplit()    
    channel_list_in = []        
    channel_list_in.append(r)
    channel_list_in.append(g) 
    channel_list_in.append(b) 
    channel_list_out = []

    for i in range(n_bands):
        print(i)

        channel_temp = channel_list_in[i]
        
        #Positive image
        statistic_vals = out1[i+1, :]        
        min_val_orig = statistic_vals[0]
        image_positive = channel_temp - min_val_orig 

        #Finding percentiles
        min_val = image_positive.min()
        max_val = image_positive.max()
        
        image_rs1 = (image_positive - min_val) * (calc_depth / (max_val - min_val))         
        image_rs2 = (image_rs1).cast("uint") #'uchar' is for 8-bit        
        th_low = image_rs2.percent(percentOut_dsaImage)
        th_high = image_rs2.percent(100 - percentOut_dsaImage)        
        
        #Finding P threshold 
        th_low_input = th_low*((max_val - min_val)/calc_depth) + min_val 
        th_high_input = th_high*((max_val - min_val)/calc_depth) + min_val                    

        #Capping
        image_rs3 = (image_positive - th_low_input) * (255 / (th_high_input - th_low_input)) 
        image_rs3 = (image_rs3 > 255).ifthenelse(255, image_rs3) #true, false
        image_rs3 = (image_rs3 < 0).ifthenelse(0, image_rs3)

        image_rs3 = image_rs3.cast("uchar") #uint8        

        channel_list_out.append(image_rs3)

    #RGB
    image_rescaled = channel_list_out[0].bandjoin(channel_list_out[1:])    
    
    #Medicine (fixing pyvips vertical flip)
    #image_flipped = image_rescaled.flipver()

    #Saving recoloured mosaic      
    if save_stitched == 1:        
        file_output = os.path.join(parentDir, fileName_mosaic)
                
        image_rescaled.write_to_file(file_output) 
        #image_rescaled.write_to_file(file_output, compression="lzw", tile=True, 
         #              tile_width=512, tile_height=512,  
          #             pyramid=True, subifd=True, bigtiff=True) 

    #endregion
    