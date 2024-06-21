# -*- coding: utf-8 -*-
"""
Deep sparse auto-encoder. 

This is the update of the first half (training DSA) of the second script for 
dimensionality reduction in Acevedo Zamora et al. (2024).
The method followed Thomas et al. (2016, 2017, 2022) papers.
The implementation uses CPU 'multi-thread' with low GPU (if available) consumption
The new 'tilingAndStacking_v3.py' and 'wsi_dimPCA_v1.m' scripts 
must be ran beforehand.

Cite as: https://doi.org/10.1016/j.chemgeo.2024.121997

Written by: Marco Andres, ACEVEDO ZAMORA
Publication (code update): 14-Nov-23
Update: 27-May-24

#Implementation note 1:
#If pyvips import is taken outside '__main__':
#(process:69204): GLib-GIO-WARNING **: 15:47:44.631: 
#Unexpectedly, UWP app `Microsoft.OutlookForWindows_1.2023.719.200_x64__8wekyb3d8bbwe' 
#(AUMId `Microsoft.OutlookForWindows_8wekyb3d8bbwe!Microsoft.OutlookforWindows') 
#supports 1 extensions but has no verbs

"""
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

from tqdm import tqdm

#Developed functions (see folder)
from MyDataset import MyDataset
from SparseAutoencoder import SparseAutoencoder
from get_device import get_device

#region
if __name__ == '__main__':    #lock to only 1 subprocess (Windows) 
    import pyvips #see Note 1

    #region 
    ## User Input
    #use \\ for Windows
    
    workingDir = "D:\\Chris_Collaboration\\2024_Ioan_Purdys_reward\\XFM\\153874_P2\\tiff\\Au\\recoloured_trial2_pctOut0.1"
    
    #input folder from MatLab script (within PCA folder)
    tag_combo = '13-Jun-24_trial1' 
    metadataFile_mat = os.path.join(workingDir, 'pca_results\\'+  tag_combo + 
                                '\\montage_pcaInfo.mat')
    
    #output folder name
    tag_combo_output = '20Jun24_test1' #tag (date, trial)    
    outputFolder = os.path.join(workingDir, 'dsa_results\\'+  tag_combo_output)
    
    #For loading tiles (choose linear or log-transformed pyramid)
    imageDir = os.path.join(workingDir, 'recoloured_pyramid_files\\0')

    #Autoencoder settings (see 'parametrisation_guide.txt')       
    #Model parameters
    LEARNING_RATE = 10e-4 #5e-2 (might diverge!, decrease if necessary)              
    RHO = 0.5 #KL new mean; default 0.5 (larger reduces train loss >> validation loss)
    epoch_default = 10 #default= 5, can get better Loss
    alpha_reg = 0.0001 #default= 0.001
    betha_reg = 10 #default= 10

    #Training parameters:
    fileSize_th = 12 #default<= 4GB (overload RAM risk); fraction of WSI for training
    fraction = 1 #default= 0.05; fraction of data for training
    test_ratio = 0.2 #test/training populations  
    n_workers = 8    
    BATCH_SIZE = 8192 #default=1024, 32 with CPU, 256 with GPU, predict uses=8192 (recommended for training)          
    seed_value = 0 #seed for subsetting arrays in training/test sets            

    #Prediction parameters:    
    n_workers_pred = n_workers #default=0; ideal=8; OverflowError: cannot serialize a bytes object larger than 4 GB
    BATCH_SIZE_pred = 8192 
    #default =8192; larger is more consistent and takes less time, e.g.: 2048 or 8192    

    #imaging settings
    outPct = 0.03 #default= 0.03; increases DSA image colour contrast

    #endregion
    
    #region
    ## Helper functions     

    #Adding to binary to path (avoid conflict w/ other programs)
    vipsbin = r'C:\\Users\\n10832084\\AppData\\Local\\vips-dev-8.12\bin'
    os.environ['PATH'] = vipsbin + ';' + os.environ['PATH']
    print("vips version: " + str(pyvips.version(0))+"."+str(pyvips.version(1))+"."+str(pyvips.version(2)))

    def make_dir(destDir):
        image_dir = destDir
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)     
        
    # Cost regularization functions    
    #Kullback-Leibler penalty
    def kl_divergence(rho, rho_hat):
        rho_hat = torch.mean(F.sigmoid(rho_hat), 1) # sigmoid because we need the probability distributions
        rho = torch.tensor([rho] * len(rho_hat)).to(device)
        output_loss = torch.sum(rho*torch.log(rho/rho_hat) + (1-rho) * torch.log((1-rho)/(1-rho_hat)))
        return output_loss
    
    def sparse_loss2(rho, images, model_children):
        values = images
        loss = 0
        for i in range(len(model_children)):
            values = model_children[i](values)
            loss += kl_divergence(rho, values)
        return loss
    
    #Tikhonov L2
    def sparse_loss1(model):    
        l2_norm = sum(torch.linalg.norm(p, 2) for p in model.parameters())

        #To avoid when learning rate is too high:
        #UserWarning: During the SVD execution, batches 0 failed to converge. 
        #A more accurate method will be used to calculate the SVD as a fallback. 
        #C:\actions-runner\_work\pytorch\pytorch\builder\windows\pytorch\aten\src\ATen\native\cuda\linalg\BatchLinearAlgebraLib.cpp:836.
        return l2_norm
    
    # Training function
    def fit(model, dataloader, epoch):
        
        model_children = list(model.children()) # get the layers       

        print('Training')
        model.train()
        
        counter = 0
        running_loss = 0.0        
        for i, data in tqdm(enumerate(dataloader), total=int(len(trainset)/dataloader.batch_size)):
            counter += 1
            img, _ = data
            img = img.to(device)
            img = img.view(img.size(0), -1)
            
            optimizer.zero_grad()
            outputs = model(img)
            
            mse_loss = criterion(outputs, img)
            if ADD_SPARSITY == 'yes':
                l2_loss = sparse_loss1(model)
                kl_loss = sparse_loss2(RHO, img, model_children)
                
                # add the sparsity penalty
                loss = mse_loss + ALPHA*l2_loss + BETA*kl_loss
            else:
                loss = mse_loss
                
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
        
        epoch_loss = running_loss / counter
        print(f"Train Loss: {epoch_loss:.3f}")
                
        return epoch_loss
    
    # Validation function
    def validate(model, dataloader, epoch):
        print('Validating')
        model.eval()
        
        running_loss = 0.0
        counter = 0
        with torch.no_grad():
            for i, data in tqdm(enumerate(dataloader), total=int(len(testset)/dataloader.batch_size)):
                counter += 1
                img, _ = data
                img = img.to(device)
                img = img.view(img.size(0), -1)
                
                outputs = model(img)
                loss = criterion(outputs, img)
                running_loss += loss.item()
                
        epoch_loss = running_loss / counter
        print(f"Val Loss: {epoch_loss:.3f}")          
        
        return epoch_loss
    
    #endregion    

    #region   
    #%% Script
    np.random.seed(seed_value)  
    torch.manual_seed(seed_value)  

    #Command-line argument parsers 
    ap = argparse.ArgumentParser()
    ap.add_argument('-e', '--epochs', type=int, default= epoch_default,
        help='number of epochs to train our network for')
    ap.add_argument('-l1', '--reg_param1', type=float, default= alpha_reg, 
        help='regularization parameter `lambda`')
    ap.add_argument('-l2', '--reg_param2', type=float, default= betha_reg, 
        help='regularization parameter `lambda`')
    ap.add_argument('-sc', '--add_sparse', type=str, default='yes', 
        help='whether to add sparsity contraint or not') 
    args = vars(ap.parse_args())
    
    EPOCHS = args['epochs']
    ALPHA = args['reg_param1'] 
    BETA = args['reg_param2'] 
    ADD_SPARSITY = args['add_sparse'] #if 'yes', the cost function is regularized    
    print(f"Add sparsity regularization: {ADD_SPARSITY}")   
    
    #Directories
    os.chdir(workingDir)              
    print("Current working directory: {0}".format(os.getcwd()))    
    make_dir(outputFolder) #saving dir   

    fileDest = f"{outputFolder}/sessionVariables.pckl"    
    modelPATH = f"{outputFolder}/model_epochs{EPOCHS}.tar"
    lossPlotPATH = f"{outputFolder}/loss_plot.png"
        
    #Save script metadata
    obj = [workingDir, metadataFile_mat, outputFolder,
           outPct, LEARNING_RATE, RHO, 
           epoch_default, alpha_reg, betha_reg, fraction, 
           test_ratio, n_workers, BATCH_SIZE, n_workers_pred, 
           BATCH_SIZE_pred, EPOCHS, ALPHA, BETA, ADD_SPARSITY]
    
    f = open(fileDest, 'wb')
    pickle.dump(obj, f)
    f.close()            
    
    #endregion

    #region    
    # Loading data    
    f = h5py.File(metadataFile_mat, 'r') #info about montage and selected stack layers
    metadata_struct = f.get('s')    
    imageHeight = metadata_struct['imageHeight'][()][0]
    imageWidth = metadata_struct['imageWidth'][()][0]
    idx_in = metadata_struct['idx_in'][()].flatten() #same order as 'descriptiveStats.csv'
    n_layers_input = np.sum(idx_in, axis= 0)  
    idx_in_bool = idx_in.astype(bool) #layers chosen in 'wsi_dimPCA_v2.m' 

    #print(list(metadata_struct.keys()))
    print(f"mosaic of WxHxC= {imageWidth}x{imageHeight}x{n_layers_input} px")

    # scan tileset    
    pattern = re.compile(r".*\\(\d+)_(\d+)\.tif") #Windows
    #pattern = re.compile(r".*/(\d+)_(\d+)\.tif") #Linux

    max_x = 0
    max_y = 0
    for filename in glob.glob(f"{imageDir}/*_*.tif"):
        match = pattern.match(filename)
        if match:        
            max_x = max(max_x, int(match.group(1)))
            max_y = max(max_y, int(match.group(2)))

    print(f"mosaic of {max_x + 1}x{max_y + 1} tiles")

    #Predict file size
    bitDepth = (2**8)
    fileSize = (10**(-9))*(imageHeight*imageWidth*n_layers_input*bitDepth/8) #GB    

    if fileSize > fileSize_th:
        fraction = fileSize_th/fileSize    
    else:
        fraction = 1
    print(f"sampling fraction= {fraction}")    

    #Load tiles and resample
    tiles_across = max_x + 1
    tiles_down = max_y + 1
    tile_arrays = []
    for y in range(0, tiles_down):
        for x in range(0, tiles_across):

            fileName_expression = f"{x}_{y}.tif"
            fileName = os.path.join(imageDir, fileName_expression)
            fileName_parse = f"{fileName}"
            image_temp = pyvips.Image.new_from_file(fileName_parse)

            image_np = image_temp.numpy()
            dim_shape = image_np.shape  
            image_xyz = np.reshape(image_np, (-1, dim_shape[2])) #39 channels
            n_rows = image_xyz.shape[0]                

            to_row = int(n_rows*(fraction))
            indices = np.random.permutation(n_rows) #scrambling
            sampling_idx = indices[:to_row]        
                        
            image_sampled = image_xyz[sampling_idx, :][:, idx_in_bool] #selected data (data4 = target)    
            #print(image_sampled.shape)
            tile_arrays.append(image_sampled)

    #Full dataset
    data4 = np.concatenate([i for i in tile_arrays])
    n_rows_sel = data4.shape[0]
    n_channels = data4.shape[1]    
    print(n_channels)

    #endregion

    #region
    #Dataset splitting
    to_row = int(n_rows_sel*(1 - test_ratio))
    indices = np.random.permutation(n_rows_sel)
    training_idx, test_idx = indices[:to_row], indices[to_row:]
    training_x, test_x = data4[training_idx,:], data4[test_idx,:]
    training_y, test_y = data4[training_idx, :], data4[test_idx, :]
    
    #Data transformations 
    loader = transforms.Compose([ 
        transforms.ToTensor(),
        transforms.Normalize((0.5,), (0.5,)), #normalisation for faster training
        ])
    trainset = MyDataset(training_x, training_y, transform= loader) #convert to tensor (uint8 required by PIL)
    testset = MyDataset(test_x, test_y, transform= loader)
    
    #Arranging in batches     
    trainloader = DataLoader(trainset, batch_size=BATCH_SIZE,
                             shuffle=False, num_workers= n_workers)
     
    testloader = DataLoader(testset, batch_size=BATCH_SIZE, 
                            shuffle=False, num_workers= n_workers)          
    
    device = get_device()
    model = SparseAutoencoder(n_channels=n_channels).to(device) # print(model)        
    criterion = nn.MSELoss() #loss function= mean squared error
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE) #optimizer    
    model_children = list(model.children()) # get the layers       
    
    # %% Train and test            
   
    train_loss = []
    val_loss = []    
    start = time.time()
    for epoch in range(EPOCHS):
        print(f"Epoch {epoch+1} of {EPOCHS}")
        train_epoch_loss = fit(model, trainloader, epoch)
        val_epoch_loss = validate(model, testloader, epoch)
        
        train_loss.append(train_epoch_loss)
        val_loss.append(val_epoch_loss)
    
    end = time.time()
    print(f"DSA training and validation took {(end-start)/60:.3} min")
    
    #Save loss plot
    plt.figure(figsize=(10, 7))
    plt.plot(train_loss, color='orange', label='training loss')
    plt.plot(val_loss, color='red', label='validation loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(lossPlotPATH)
    plt.show()

    #Save the model    
    torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': train_loss
            }, modelPATH)   
    
    #endregion
    
    