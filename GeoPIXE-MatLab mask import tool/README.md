
# GeoPIXE-MatLab mask import tool

## Description

The code allows sub-setting a QuPath semantic segmentation map and saving the mask coordinates of the target minerals as 'region' files for GeoPIXE. This allows interrogating the Synchrotron X-ray fluorescence microscopy spectra directly within GeoPIXE (hosted in a Linux supercomputer in Australia). The spectra can then be re-fitted and re-quantified to newly export mineral-mask calibrated maps (fitting hidden peaks that are not apparent in the overall map area spectra). 

For example:

<img src="https://github.com/marcoaaz/AcevedoEtAl._2024b_autoencoder/assets/61703106/5b6cf7e7-859e-43bc-b98b-a382a296e3af" width=80% height=80%>

Description of scripts:

  + 'regionToGeoPIXE_v2.m' = script to convert a QuPath software mineral phase map 'region of interest' (ROI) into a GeoPIXE Q-vector file (readable as a region in GeoPIXE). 
  +	backgroundColourInverter.m = whitens the optical phase map background considering the foreground mask
  + load_csvTableExport.m = normalises spectra and compares the different regions and spectra exported from GeoPIXE
  + figureGeoPIXE_v3.m = script for generating spectral figures for Synchrotron XFM (for comparison). It supports GeoPIXE spectral output conventions. The script allows reproducing the figure above (see paper).

Authorship:
  + The MatLab implementation side (this folder) was done by Marco Acevedo (QUT).
  + The GeoPIXE implementation side to read 'region' files was done directly on IDL language by Dr. Chris Ryan (CSIRO).
  + When using the code, please cite: "Acevedo Zamora, M. A., Kamber, B. S., Jones, M. W. M., Schrank, C. E., Ryan, C. G., Howard, D. L., Paterson, D. J., Ubide Garralda, T., & Murphy, D. (under review). Tracking element-mineral associations with unsupervised learning and dimensionality reduction in chemical and optical image stacks of thin sections. Chemical Geology."
