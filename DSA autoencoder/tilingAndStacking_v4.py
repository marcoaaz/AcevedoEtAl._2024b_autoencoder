# -*- coding: utf-8 -*-
"""
Tiling and stacking chemical layers of very large (WSI) X-ray maps

This is the update of the first half of the first script for dimensionality reduction 
in Acevedo Zamora et al. (2024). 
The update follows 'stitch_stretch_batch.py' (Acevedo Zamora and Kamber, 2023) in
implementing pyvips 

Cite as: https://doi.org/10.1016/j.chemgeo.2024.121997

Written by: Marco Andres, ACEVEDO ZAMORA
Created on Tue Apr 12 11:13:17 2022
Published (first version): 14-Nov-23
Updated (second version for whole-slide imaging): 19-Jun-24

Followed sources:
https://libvips.github.io/pyvips/vimage.html
https://www.libvips.org/API/current/Examples.html
https://github.com/libvips/libvips/issues/2600   
https://forum.image.sc/t/reading-regions-of-tif-files-with-more-than-3-channels/93299/6
https://forum.image.sc/t/pyvips-2-2-is-out-with-improved-numpy-and-pil-integration/66664/3
https://www.libvips.org/API/current/libvips-arithmetic.html#vips-stats

"""
#!/usr/bin/python3

#Dependencies   
import os
import sys
import pyvips
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import shutil

## User input

imageFolder = "D:\\Chris_Collaboration\\2024_Ioan_Purdys_reward\\XFM\\153874_P2\\tiff\\hybrid"
sourceFolder = os.path.split(imageFolder)

fileList = glob.glob(f"{imageFolder}/153874_*-*.tiff") #perfect match needed with 'pattern'
pattern = re.compile(r".+\\153874_(.+)-(.+)\.tiff") #for Windows \\
#pattern = re.compile(r".*/(.*)_(\d+)_(\d+)\.tiff") #for Linux

trial_name = 'trial1'
pctOut = 0.1 #percentile out in the input (for colour contrast)
filterSize = 5 #for smoothness
calc_depth = 2**16-1 #for precision 8-bit=255; 16-bit=65535
save_recoloured = 0 #quite time-consuming for each layer
save_stack_transformed = 1 #linear/log (recoloured) pyramid
save_stack = 1 #original pyramid (takes longer since it is float)


## Functions

def percentilesAndCapping(img_input, pctOut):
    
    #Finding percentiles
    min_val = img_input.min()
    max_val = img_input.max()
    image_rs1 = (img_input - min_val) * (calc_depth / (max_val - min_val))         
    image_rs2 = (image_rs1).cast("uint") #'uchar' is for 8-bit        
    th_low = image_rs2.percent(pctOut)
    th_high = image_rs2.percent(100 - pctOut)        
    #Finding P threshold in img_input
    th_low_input = th_low*((max_val - min_val)/calc_depth) + min_val 
    th_high_input = th_high*((max_val - min_val)/calc_depth) + min_val      

    #prevent division by zero (constant image)
    if th_low_input == th_high_input:
        th_high_input += 1            

    #Capping
    image_rs3 = (img_input - th_low_input) * (255 / (th_high_input - th_low_input)) 
    image_rs3 = (image_rs3 > 255).ifthenelse(255, image_rs3) #true, false
    image_rs3 = (image_rs3 < 0).ifthenelse(0, image_rs3)

    return image_rs3

def remove(path):
    """ param <path> could either be relative or absolute. """
    if os.path.isfile(path) or os.path.islink(path):
        os.remove(path)  # remove the file
    elif os.path.isdir(path):
        shutil.rmtree(path)  # remove dir and all contains
    else:
        raise ValueError("file {} is not a file or dir.".format(path))
    
## Script    

n_files = len(fileList)
print(f"There are {n_files} original images")    
#print(np.transpose(np.array(fileList)))

#Define colour map (falsecolour images)
img_indexes = pyvips.Image.identity()
lut = img_indexes.falsecolour() #using standard heatmap
#256x1 uchar, 3 bands, srgb, pngload

## Script
destDir = os.path.join(imageFolder, 
                       'recoloured_' + trial_name + '_pctOut' + str(pctOut))

if not os.path.exists(destDir):
    try:
        os.mkdir(destDir)
    except OSError as e:
        print(f"An error has occurred: {e}")
        raise

# scan image set (perfect match needed)
path_list = []
experiment_list = []
element_list = []
out2 = [] #for csv
stack_layers = [] #linear data
stack_layers_log = [] #natural log() data
k = 0
for filename in fileList:
#for filename in fileList[0:3]: #trial

    match = pattern.match(filename)

    if match:                
        k += 1

        #sample = match.group(1)
        #element = match.group(2) #default

        #element = match.group(1) #most common

        experiment = match.group(1) #using hybrid folders
        element = match.group(2) 

        print(f"Processing experiment {experiment} match {k}: {element}")         

        path_list.append(filename)
        experiment_list.append(experiment)
        element_list.append(element)
        
        image = pyvips.Image.new_from_file(filename) 
        #optional: path, page=i, access="sequential" 
        
        #Descr. statistics
        out = pyvips.Image.stats(image)
        out1 = out.numpy()
        statistic_vals = out1[0, :]
        out2.append(statistic_vals)        

        #Transforming input
        min_val_orig = statistic_vals[0]
        image_positive = image - min_val_orig + 1
        image_log = image_positive.log()                

        #Finding percentiles, capping and stretching image histogram top and bottom
        image_rs = percentilesAndCapping(image_positive, pctOut)
        image_log_rs = percentilesAndCapping(image_log, pctOut)

        #Median filter
        image_med = image_rs.median(filterSize) #linear
        image_log_med = image_log_rs.median(filterSize) #log
        image_med = image_med.cast("uchar") #uint8        
        image_log_med = image_log_med.cast("uchar") #uint8        

        #Saving recoloured images (for retrospective feedback): time-consuming        
        if save_recoloured == 1:            
            destFile1 = os.path.join(destDir, experiment + "_" + element + ".tif")#.ome.tif
            
            image_recoloured = image_med.maplut(lut)
            image_recoloured.write_to_file(destFile1)                            
        
        #Building recoloured stack        
        stack_layers.append(image_med)
        stack_layers_log.append(image_log_med)

#Info table
out3 = np.array(out2) 
stats_table = pd.DataFrame(out3)
stats_table.columns =['min', 'max', 'sum', 
                      'sumOfSquares', 'mean', 'stddev', 
                      'min_x_coord', 'min_y_coord', 
                      'max_x_coord', 'max_y_coord']
stats_table.insert(0, "experiment", experiment_list, True)
stats_table.insert(1, "element", element_list, True)
stats_table.insert(2, "path", path_list, True)
#save
file_name1 = os.path.join(destDir, "descriptiveStats.csv")
stats_table.to_csv(file_name1, sep=',', 
                   encoding='utf-8', index=False, header=True)


#Saving recoloured stack (uint8)
#input for PCA analysis 'wsi_dimPCA_v1.m' and autoencoder 'DSA_images_v4.py' scripts

if save_stack_transformed == 1:  
    
    #Linear data    
    image_stack_recoloured = stack_layers[0].bandjoin(stack_layers[1:])

    #Medicine (fixing pyvips vertical flip)
    image_flipped = image_stack_recoloured.flipver()

    destDir2 = os.path.join(destDir, 'linear_pyramid')
    destDir2_files = destDir2 + '_files'
    try:
        remove(destDir2_files)
    except:
        print("Producing linear pyramid for the first time")

    image_flipped.dzsave(destDir2, suffix='.tif', 
                    skip_blanks=-1, background=0, 
                    depth='one', overlap=0, tile_size= 1024, 
                    layout='dz') #Tile overlap in pixels*2
    
    ###################

    #Log-transformed data
    image_stack_log_recoloured = stack_layers_log[0].bandjoin(stack_layers_log[1:])

    #Medicine (fixing pyvips vertical flip)
    image_log_flipped = image_stack_log_recoloured.flipver()

    destDir3 = os.path.join(destDir, 'log_pyramid')
    destDir3_files = destDir3 + '_files'
    try:
        remove(destDir3_files)
    except:
        print("Producing natural log pyramid for the first time")

    image_log_flipped.dzsave(destDir3, suffix='.tif', 
                    skip_blanks=-1, background=0, 
                    depth='one', overlap=0, tile_size= 1024, 
                    layout='dz') #Tile overlap in pixels*2    
    

#Saving original stack (used for 'ROIimageAnalysis_v7_wsi.m' script)

if save_stack == 1:
    pages = [pyvips.Image.new_from_file(path) for path in path_list]
    image_stack = pages[0].bandjoin(pages[1:])

    #Medicine (fixing pyvips vertical flip)
    image_flipped = image_stack.flipver()

    destDir4 = os.path.join(destDir, 'original_pyramid')
    destDir4_files = destDir4 + '_files'
    try:
        remove(destDir4_files)
    except:
        print("Producing original pyramid for the first time")

    image_flipped.dzsave(destDir4, suffix='.tif', 
                    skip_blanks=-1, background=0, 
                    depth='one', overlap=0, tile_size= 1024, 
                    layout='dz') #Tile overlap in pixels*2
