# -*- coding: utf-8 -*-
"""
Deep sparse auto-encoder following Thomas et al. (2016, 2017, 2022) papers.
Written by: Marco Andres, ACEVEDO ZAMORA
Place: QUT
Last update: 14-Nov-2023

"""
import os
import argparse
import time
import numpy as np

from scipy.io import savemat
import h5py
import pickle

import matplotlib
import matplotlib.pyplot as plt 
matplotlib.style.use('ggplot')

import torch
import torch.nn as nn
import torchvision.transforms as transforms
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
# from torchvision.utils import save_image

from PIL import Image
from tqdm import tqdm


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
    
    
if __name__ == '__main__':    #lock to only 1 subprocess (Windows) 
    
    #% Helper functions     
    def make_dir(destDir):
        image_dir = '../' + destDir + '/images'
        if not os.path.exists(image_dir):
            os.makedirs(image_dir)

    
    def get_device(): # get the computation device
        if torch.cuda.is_available():
            device = 'cuda:0' 
            print(f"using {device}")
        else:
            device = 'cpu'
            print(f"using {device}")
        return device    
        
    # % Define the autoencoder model (DSA has a symmetric layout)        
    class SparseAutoencoder(nn.Module):
        def __init__(self):
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
        
    # Cost regularization functions    
    #Kullback-Leibler penalty
    def kl_divergence(rho, rho_hat):
        rho_hat = torch.mean(F.sigmoid(rho_hat), 1) # sigmoid because we need the probability distributions
        rho = torch.tensor([rho] * len(rho_hat)).to(device)
        output_loss = torch.sum(rho*torch.log(rho/rho_hat) + (1-rho) * torch.log((1-rho)/(1-rho_hat)))
        return output_loss
    
    def sparse_loss2(rho, images):
        values = images
        loss = 0
        for i in range(len(model_children)):
            values = model_children[i](values)
            loss += kl_divergence(rho, values)
        return loss
    
    #Tikhonov L2
    def sparse_loss1(model):    
        l2_norm = sum(torch.linalg.norm(p, 2) for p in model.parameters())
        return l2_norm
    
    # Training function
    def fit(model, dataloader, epoch):
        print('Training')
        model.train()
        
        running_loss = 0.0
        counter = 0
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
                kl_loss = sparse_loss2(RHO, img)
                
                # add the sparsity penalty
                loss = mse_loss + ALPHA*l2_loss + BETA*kl_loss
            else:
                loss = mse_loss
                
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
        
        epoch_loss = running_loss / counter
        print(f"Train Loss: {epoch_loss:.3f}")
        
        # # save the reconstructed images 
        # save_decoded_image(outputs.cpu().data, f"../outputs/images/train{epoch}.png")
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
        
        # # save the reconstructed images 
        # outputs = outputs.view(outputs.size(0), 1, 28, 28).cpu().data
        # save_image(outputs, f"../outputs/images/reconstruction{epoch}.png")
        return epoch_loss
    
    # Prediction functions
    def embedded_space(images):
        values = images    
        for i in range(3):
            values = model_children[i](values)
            values = F.sigmoid(values) # sigmoid because we need probability distributions            
        return values
    
    def predict_space(model, dataloader):
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
                
                output = embedded_space(img)
                pixel_batch.append(output)
            # result = torch.cat(pixel_batch, dim=0)
        return pixel_batch
    
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
    
    # for saving the reconstructed images
    # def save_decoded_image(img, name):
    #     img = img.view(img.size(0), 1, 28, 28)
    #     save_image(img, name)        
    
    # %% User Input
    
    # Change the current working directory (use \\ for Windows machine)      
    # workingDir = 'C:\\Users\\n10832084\\OneDrive - Queensland University of Technology\\Desktop\\E3_christoph\\E3_66039_XFM scanning\\synchrotron_stack'     
    # matlabFile = 'stack_21febList_0.1_denoised.mat'
    # destDir = '1f_5x5MedF_2'
    # imgName_pred = 'MedF_1'
    
    workingDir = 'E:\\paper 3_datasets\\harzburgite_synchrotron_christoph\\tiff_HiRes\\synchrotron_stacks2'     
    matlabFile = 'stack_temp_ext_notMedF.mat'
    destDir = '17-Nov-23' #destination directory (date)    
    descrip_pred = 'DSA' 
    imgName_pred = 'test3_balz-energySettings' #trial number (some parameters)
    outPct = 0.03 #default= 0.03
    
    #Autoencoder settings        
    #Model parameters
    LEARNING_RATE = 5e-2              
    RHO = 0.5 #KL new mean; default 0.5 (larger reduces train loss >> validation loss)
    epoch_default = 6 #default= 5, can get better Loss
    alpha_reg = 0.0001 #default= 0.001
    betha_reg = 10 #default= 10
    
    #Training parameters:
    fraction = 1 #default= 0.05; fraction of data for training
    test_ratio = 0.2 #test/training populations  
    n_workers = 8    
    BATCH_SIZE = 8192 #default=1024, 32 with CPU, 256 with GPU, predict uses=8192          
            
    #Prediction parameters:    
    n_workers_pred = n_workers #default=0; ideal=8; OverflowError: cannot serialize a bytes object larger than 4 GB
    BATCH_SIZE_pred = 8192 
    #default =8192; larger is more consistent and takes less time, e.g.: 2048 or 8192    
    
    #Notes:        
    #LEARNING_RATE default= 1e-4 (might be dumb!, useful for SEM-EDX); 
    #1e-3 (if important stuff is lost); 
    #1e-2 (might diverge!, useful for Synchrotron)
    
    # ALPHA = 0.001 
    #Tikhonov L2 reg., 0.01 to default 0.0001; 
    #at smaller alpha, training approches validation
    #weight-decay term prevents overfitting

    # BETA = 10 
    #KL reg., 0.001; default 100; works 0.1; 10 for PIXL
    #measure how different two distributions are (mean of activations vs. rho)
    #compensate large number of hidden units against sparsity (nodes > than input)        
    
    # %% Script
    
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
    
    os.chdir(workingDir)
    make_dir(destDir) #saving dir             
    print("Current working directory: {0}".format(os.getcwd()))    
    
    #Save script metadata
    obj = [workingDir, matlabFile, destDir, descrip_pred, 
           imgName_pred, outPct, LEARNING_RATE, RHO, 
           epoch_default, alpha_reg, betha_reg, fraction, 
           test_ratio, n_workers, BATCH_SIZE, n_workers_pred, 
           BATCH_SIZE_pred, EPOCHS, ALPHA, BETA, ADD_SPARSITY]
    
    fileDest = f"../{destDir}/{descrip_pred}_{imgName_pred}_sessionVariables.pckl"    
    f = open(fileDest, 'wb')
    pickle.dump(obj, f)
    f.close()            
    
    # %% Load metadata
    
    #manual (optional)
    destDir1 = '17-Nov-23'
    descrip_pred1 = 'DSA'
    imageName_pred1 = 'test2'
    fileDest = f"../{destDir1}/{descrip_pred1}_{imageName_pred1}_sessionVariables.pckl"
    f = open(fileDest, 'rb')
    obj = pickle.load(f)
    f.close()    
    [workingDir, matlabFile, destDir, descrip_pred, 
           imgName_pred, outPct, LEARNING_RATE, RHO, 
           epoch_default, alpha_reg, betha_reg, fraction, 
           test_ratio, n_workers, BATCH_SIZE, n_workers_pred, 
           BATCH_SIZE_pred, EPOCHS, ALPHA, BETA, ADD_SPARSITY] = obj            
    
    # %%
    # Loading data    
    f = h5py.File(matlabFile, 'r') #'all_modes2.mat'
    data = f.get('stack_temp') #'img_stack_multiplex'
    data = np.array(data) # For converting to a NumPy array    
    data2 = np.transpose(data) #matlab to python convention
    
    # x_tl = 587 #optional (to avoid fringes in WSI)
    # y_tl = 71
    # bb_w = 5095
    # bb_h = 8921
    # x_br = x_tl + bb_w
    # y_br = y_tl + bb_h
    # data2 = data2[y_tl:y_br, x_tl:x_br, :]
    
    dim_shape = data2.shape  
    data3 = np.reshape(data2, (-1, dim_shape[2])) #for image input; 39 channels
    n_rows = data3.shape[0]
    data_full = data3
    targets_full = data_full
    
    #Resampling dataset (economize RAM)
    seed_value = 0
    torch.manual_seed(seed_value)
    np.random.seed(seed_value)  
    
    to_row = int(n_rows*(fraction))
    indices = np.random.permutation(n_rows) #scrambling
    sampling_idx = indices[:to_row]
    data4 = data3[sampling_idx, :] #selected data
    n_rows = data4.shape[0]
    n_channels = data4.shape[1]
    targets = data4
    
    #Dataset splitting
    to_row = int(n_rows*(1 - test_ratio))
    indices = np.random.permutation(n_rows)
    training_idx, test_idx = indices[:to_row], indices[to_row:]
    training_x, test_x = data4[training_idx,:], data4[test_idx,:]
    training_y, test_y = targets[training_idx, :], targets[test_idx, :]
    
    #Data transformations 
    loader = transforms.Compose([ 
        transforms.ToTensor(),
        transforms.Normalize((0.5,), (0.5,)), #normalisation for faster training
        ])
    trainset = MyDataset(training_x, training_y, transform= loader) #convert to tensor
    testset = MyDataset(test_x, test_y, transform= loader)
    
    #Arranging in batches     
    trainloader = DataLoader(
         trainset, 
         batch_size=BATCH_SIZE,
         shuffle=False, #True not mandatory
         num_workers= n_workers
     )
     
    testloader = DataLoader(
         testset, 
         batch_size=BATCH_SIZE, 
         shuffle=False,
         num_workers= n_workers
     )          
    
    device = get_device()
    model = SparseAutoencoder().to(device) # print(model)        
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
    print(f"{(end-start)/60:.3} minutes")
    
    modelPATH = f"../{destDir}/{descrip_pred}_{imgName_pred}_model_{EPOCHS}.tar"
    lossPlotPATH = f"../{destDir}/{descrip_pred}_{imgName_pred}_loss.png"
    
    # save the model    
    torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': train_loss
            }, modelPATH)
    
    #loss plots
    plt.figure(figsize=(10, 7))
    plt.plot(train_loss, color='orange', label='training loss')
    plt.plot(val_loss, color='red', label='validation loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    plt.savefig(lossPlotPATH)
    plt.show()
    
    # %%Predict: Loading model and save embedded space image
    #Note: this section can run independently from script           
    # del modelPATH
    # del checkpoint
    # del model_pred
    
    #Info pre-loaded
    dim_shape = data2.shape            
    modelPATH = f"../{destDir}/{descrip_pred}_{imgName_pred}_model_{EPOCHS}.tar" #must be *.tar not *.pth     
        
    #Input model path manually: Transfer model from other dataset
    # modelPATH = 'E:\\paper 3_datasets\\data_DMurphy\\91712-81R5w-glass-quadrant cpx\\tiff\\15-Nov-23\\DSA_test1_model_5.tar'
    
    #Initialize
    checkpoint = torch.load(modelPATH) #dictionary containing objects
    model_pred = SparseAutoencoder().to(device)
    optimizer_pred = optim.Adam(model_pred.parameters(), lr=LEARNING_RATE)    
    model_children = list(model_pred.children()) # get the layers       
    #Load    
    model_pred.load_state_dict(checkpoint['model_state_dict'])
    optimizer_pred.load_state_dict(checkpoint['optimizer_state_dict'])
    epoch = checkpoint['epoch']
    loss = checkpoint['loss']    
    
    model_pred.eval() #set dropout and batch normalisation to evaluation mode
    # model_pred = deepcopy(model_pred.state_dict()) #prevents retraining
    
    #Info:         
    print("Model's state_dict:")
    for param_tensor in model_pred.state_dict():
        print(param_tensor, "\t", model_pred.state_dict()[param_tensor].size())
        
    print("Optimizer's state_dict:")
    for var_name in optimizer_pred.state_dict():
        print(var_name, "\t", optimizer_pred.state_dict()[var_name])
        
    #Predict embedding space   
    all_set = MyDataset(data_full, targets_full, transform= loader)
    alldata_loader = DataLoader(
        all_set, 
        batch_size= BATCH_SIZE_pred, 
        shuffle= False,
        num_workers= n_workers_pred
    )         
        
    # Run model forward
    pixel_batch = predict_space(model_pred, alldata_loader)    
    print(pixel_batch[0].shape)
    result = torch.cat(pixel_batch, dim=0)    
    
    # %Generate image    
    outputs = result.view(dim_shape[0], dim_shape[1], 3).cpu().numpy() #matrix              
    R = rescale(outputs[:, :, 0], outPct)
    G = rescale(outputs[:, :, 1], outPct)
    B = rescale(outputs[:, :, 2], outPct)
    img_rgb = np.stack((R, G, B), axis=2) #RGB image   
    
    # %save the reconstructed images         
    fileName = f"../{destDir}/{descrip_pred}_{imgName_pred}_mtx.mat"
    savemat(fileName, dict(x= outputs[:, :, 0], y= outputs[:, :, 1], z= outputs[:, :, 2]))

    img_rgb = Image.fromarray(img_rgb)    
    img_rgb.save(f"../{destDir}/{descrip_pred}_{imgName_pred}.tif")
    # img_rgb.save(f"E:\\Alienware_March 22\\current work\\data_DMurphy\\91690-81r5w-glass-overview\\tiff\\output_ext_0.5dataset\\{descrip}_trainedOriginal.tif")   
   
    plt.figure(figsize=(10, 7))
    plt.imshow(img_rgb, interpolation="nearest")
   
    # %% User Feedback
   
    n_bins = 40 #'auto'
    alphaFace = 0.3 #0.7
    alphaVal = 1 #0.75, grid
    rwidthVal = .85 #0.85
    
    #Channel histograms     
    n, bins, patches = plt.hist(x= R.flatten(), bins= n_bins, alpha= alphaFace, rwidth= rwidthVal, facecolor= (1, 0, 0, 0.5), edgecolor='red', linewidth=1.2)
    plt.grid(axis='y', alpha= alphaVal, )    
    
    n, bins, patches = plt.hist(x= G.flatten(), bins= n_bins, alpha= alphaFace, rwidth= rwidthVal, facecolor= (0, 1, 0, 0.5), edgecolor='green', linewidth=1.2)
    plt.grid(axis='y', alpha=alphaVal)    
    
    n, bins, patches = plt.hist(x= B.flatten(), bins= n_bins, alpha= alphaFace, rwidth= rwidthVal, facecolor= (0, 0, 1, 0.5), edgecolor='blue', linewidth=1.2)
    plt.grid(axis='y', alpha= alphaVal)
    
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Channel histograms')
    
# %% K-means for mini-phase map

# import pandas as pd
# from sklearn.cluster import KMeans

# X = img_rgb
# kmeans = KMeans(n_clusters=4, init='k-means++', max_iter=300, n_init=10, random_state=0)
# kmeans.fit(X)
# y_kmeans = kmeans.predict(X)
# x_centre = kmeans.cluster_centers_[:, 0]
# y_centre = kmeans.cluster_centers_[:, 1]
# z_centre = kmeans.cluster_centers_[:, 2]
# centres = np.stack((x_centre, y_centre, z_centre), axis=1)

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# ax.scatter(img_rgb[:, 0], img_rgb[:, 1], img_rgb[:, 2], marker='.', s=3, c= y_kmeans, cmap='viridis')
# ax.scatter(x_centre, y_centre, z_centre, s=100, c='red')
# ax.set_xlabel('R')
# ax.set_ylabel('G')
# ax.set_zlabel('B')
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot()
# ax.scatter(img_rgb[:, 0], img_rgb[:, 1], marker='.', s=3, c= y_kmeans, cmap='viridis')
# ax.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], s=100, c='red')
# ax.set_xlabel('R')
# ax.set_ylabel('G')
# plt.show()

# np.savetxt('kmeans_pred2.csv', (y_kmeans), fmt='%1d', delimiter=',')
# np.savetxt('embedded_space2.csv', img_rgb, fmt='%3d', delimiter=',')
# np.savetxt('centres2.csv', centres, fmt='%e', delimiter=',')    