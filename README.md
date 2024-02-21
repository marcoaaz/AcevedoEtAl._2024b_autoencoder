# AcevedoEtAl._2024b_autoencoder
This repository contains the original scripts to process the micro-analytical data of rock thin sections (available at [Zenodo](https://zenodo.org/records/10669251)). The full documentation of the code is in the submitted manuscript (Supplementary Material 1):

**Acevedo Zamora, M. A., Kamber, B. S., Jones, M. W. M., Schrank, C. E., Ryan, C. G., Howard, D. L., Paterson, D. J., Ubide, T., & Murphy, D. (in revision). Tracking element-mineral associations with unsupervised learning and dimensionality reduction in chemical and optical image stacks of thin sections. Chemical Geology.** 

The main workflow combines the chemical images with optical microscopy images for enabling semantic segmentation using QuPath software. 
  1. For registering the image montages we use the routine explained in this [video playlist](https://youtu.be/YpxTobsB-RM) (following [Bogovic et al., 2016](https://ieeexplore.ieee.org/document/7493463)).
  2. The pixel-based segmentation is converted into a MatLab array using our previous [work](https://github.com/marcoaaz/Acevedo-Kamber/tree/main/QuPath_generatingMaps).

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024b_autoencoder/assets/61703106/595b94b7-9cf7-45ef-b63c-715586068d63" width=80% height=80%>

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

The script Python version was 3.7.12. The autoencoder script used numpy 1.21.6, scipy 1.7.3, h5py 3.6.0, matplotlib 3.2.2, pytorch 1.12.1, pillow 9.0.1, and tqdm 4.64.1.

## Video explaining the MatLab ROI Tool

The region of interest (ROI) tool shown in the figure above (panel E) is explained in this brief [video](https://youtu.be/poPmVhwMwbA).

## Unit conversion

Interested readers can convert the Synchrotron XFM maps provided in the data repository (available at [Zenodo](https://zenodo.org/records/10669251)) to ppm using the expression 10,000*(image pixels)/(mean density (g/cm3) *thickness (µm)).

## Future Updates

We expect to increase the support for more file formats from different micro-analytical instruments. Also, file management is done manually and the scripts do not have a Graphical User Interface. 

Feel free to contact Marco if there still are issues after running the code with your own data. Thanks.
