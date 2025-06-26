# AcevedoEtAl._2024b_autoencoder

This repository contains original image analysis scripts presented in [Acevedo Zamora et al. (2024)](https://doi.org/https://doi.org/10.1016/j.chemgeo.2024.121997) to study optical microscopy and micro-analytical maps of rock thin-sections (see Zenodo [repository](https://zenodo.org/records/10669251)). For example, on this oceanic gabbro (Figure 9):

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024b_autoencoder/assets/61703106/a116201f-b2eb-44d9-a826-236ead59b04e" width=100% height=100%>


The scripts are commented and a broad explanation of the workflow steps (see below) can be found in Supplementary Material 1:

**Acevedo Zamora, M. A., Kamber, B. S., Jones, M. W. M., Schrank, C. E., Ryan, C. G., Howard, D. L., Paterson, D. J., Ubide, T., & Murphy, D. T. (2024). Tracking element-mineral associations with unsupervised learning and dimensionality reduction in chemical and optical image stacks of thin sections. ***Chemical Geology***, 121997. https://doi.org/https://doi.org/10.1016/j.chemgeo.2024.121997** 

For the original dataset (for trialling the code), please, cite:

**Acevedo Zamora, M. A. (2024). Tracking element-mineral associations with unsupervised learning and dimensionality reduction in chemical and optical image stacks of thin sections: original datasets. In Chemical Geology (version 1). ***Zenodo***. https://doi.org/10.5281/zenodo.10669251**

The optical scans stacked with chemical element maps are in the [virtual microscope](https://qutrocks.qut.edu.au/) and can be accessed using (user ; password): QUTguest_paper3 ; vs200_paper3

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024b_autoencoder/assets/61703106/c106927a-3323-4899-8abd-0c5655166976" width=80% height=80%>

## Workflow

The main workflow combines the chemical images with optical microscopy images for enabling pixel classification (semantic segmentation) using QuPath software ([Bankhead et al., 2017](https://pubmed.ncbi.nlm.nih.gov/29203879/)). 
  1. For registering the image montages we use the routine explained in this [video playlist](https://youtu.be/YpxTobsB-RM) (following [Bogovic et al., 2016](https://ieeexplore.ieee.org/document/7493463)).
  2. The pixel-based segmentation is converted into a MatLab array using our previous [work](https://github.com/marcoaaz/Acevedo-Kamber/tree/main/QuPath_generatingMaps).

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024b_autoencoder/assets/61703106/595b94b7-9cf7-45ef-b63c-715586068d63" width=60% height=60%>

To be able to run the scripts, you require installing:

  +	MATLAB Version: 9.14.0.2239454 (R2023a) Update 1
  +	MATLAB License Number: 729365
  +	Operating System: Microsoft Windows 10 Enterprise Version 10.0 (Build 19045)
  +	Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode

, with toolboxes:

  +	Computer Vision Toolbox, Version 10.4
  +	Curve Fitting Toolbox, Version 3.9
  +	Deep Learning Toolbox, Version 14.6
  +	Fixed-Point Designer, Version 7.6
  +	Image Processing Toolbox, Version 11.7
  +	MATLAB Compiler, Version 8.6
  +	Mapping Toolbox, Version 5.5
  +	Parallel Computing Toolbox, Version 7.8
  +	Signal Processing Toolbox, Version 9.2
  +	Statistics and Machine Learning Toolbox, Version 12.5
  + Symbolic Math Toolbox, Version 9.3
  + Wavelet Toolbox, Version 6.3

Meanwhile, the Python script uses Python version 3.7.12. The autoencoder script used numpy 1.21.6, scipy 1.7.3, h5py 3.6.0, matplotlib 3.2.2, pytorch 1.12.1, pillow 9.0.1, and tqdm 4.64.1. An install of pyvips is required to enable whole-slide imaging.

## Video explaining the MatLab ROI Tool

The region of interest (ROI) tool shown in the figure above (panel E) is explained in this brief [video](https://youtu.be/poPmVhwMwbA).

## Unit conversion

Interested readers can convert the Synchrotron XFM maps provided in the data repository (available at [Zenodo](https://zenodo.org/records/10669251)) to ppm using the expression 10,000*(image pixels)/(mean density (g/cm<sup>3</sup>) *thickness (µm)).

## Update (21-Jun-2024)

In the paper we mentioned that the image analysis approach for tracking and segmentation were 'scalable'. Therefore, we have updated the original scripts to support working with whole-slide imaging using image pyramids (multi-gigapixel chemical image stacks). This was possible using the powerful [pyvips library](https://github.com/libvips/pyvips) and translating the code from MatLab to Python. The results should be stored in a folder output that looks like this:

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024b_autoencoder/assets/61703106/9ce0321a-7e41-44a9-92c3-55dd8fd97097" width=90% height=90%>

The process follows:
  1. To begin, the user can find the 'DSA autoencoder>tilingAndStacking_v4.py' script for producing the 1-level pyramid (at full resolution).
     + The script produces the folders 'linear_pyramid_files', 'log_pyramid_files', and 'original_pyramid_files' in a new directory.
     + It also saves the recoloured images (of each element, of each experiment) in the same new directory (for visual expert observation).
  3. Then, they can use 'wsi_dimPCA_v2.m' to produce the PCA image (of the whole-slide) and shorten the desired element list to elements of interest (with visible patterns). The implementation uses an [external package](https://au.mathworks.com/matlabcentral/fileexchange/88872-incrementalpca) for MatLab.
     + The script produces a folder called 'pca_results' that contains the image and processed metadata (montage_pcaInfo.mat) for each trial run.
  4. Next, they need to use 'DSA_images_training_v3.py' to train the autoencoder model (*.tar) that will produce the DSA image. The directories need to be changed manually to interpret the element list pre-selected in step (2).
     + The script produces a folder called 'dsa_results' that contains the pytorch model and processed metadata (sessionVariables.pckl) for each trial run.
  5. After, they must use 'DSA_images_predicting_v4.py' to build the DSA image (of the whole-slide) loading the autoencoder model trained in (3).
     + The script saves the image in 'dsa_results' for each trial run.
  6. Finally, they can (optional) use the new ROI Tool ('ROIimageAnalysis_v7_wsi.m') that allows navigating chemical image stacks (of the whole-slide) with interactive plots that load the quantitative data (float numbers) from the 'original_pyramid_files' folder.

With the update, there is no limit to how big the input images (geochemical maps) are for producing representation (DSA or PCA images) or deploying the ROI Tool (in MatLab):

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024b_autoencoder/assets/61703106/e1cf0558-6c83-4bc8-9fca-7544a8cdd075" width=100% height=100%>

If there is any missing function dependency that I have not uploaded, please, let me know. For smooth dependencies installation, I recommend using Visual Studio Code (and 'pip') rather than Anaconda-Spyder (and 'conda install') as an IDE.

## Future Updates

File management is done manually and the scripts do not have a Graphical User Interface. We expect interested users to fork our code into other software that is already adapted to perform petrological analysis (e.g., [XMapTools](https://github.com/xmaptools)). Also, we expect to increase the support for more file formats from different micro-analytical instruments.

Particularly after the last update, I expect to continue making improvements in performance and parallelising the approach. 

Please, feel free to contact us if there still are issues after running the code with your own data and getting familiar with the workflow. 

Thanks. Cordially,
Marco AAZ.
