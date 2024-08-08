# ECLiPSE
ECLiPSE analysis: Enhanced Classification of Localized Pointclouds by Shape Extraction

The current version of the ECLiPSE software (v2.0; April 10, 2024) to automatically describe and classify super-resolution fluorescence microscopy data (i.e., 2D & 3D localizations). 
This version only contains scripts/functions to do the analysis, but a GUI version of ECLiPSE is being prepared to provide a user-friendly experience. 

## System requirements
The files were tested in MATLAB R2022b (The MathWorks, USA), on a i7-12700H 2.30 GHz (32 GB RAM) laptop running Windows 64-bit. No errors were encountered (only some warnings sometimes related to the polyshape function of MATLAB, but these can safely be ignored). 

The code was not tested in another version of MATLAB, but should be working in any version of MATLAB where the polyshape.m function exists (i.e., MATLAB R2017b or newer).

We are aware that in very rare cases (<0.001% of the tested data), a small bug shows up where the software gets stuck in a loop while calculating the descriptors for the clusters. You should notice this by the waitbar that does not advance (Note: for clusters with sizes > 1μm, it may take a few seconds to calculate the descriptors. These are not to be confused with the small bug in the calculations).

## Installation and run guide
To install:
  - Download the full folder to your computer.
  - Add the full folder to your MATLAB path 
  ```
  Option 1: Navigate to the folder through the 'Current Folder' menu and right click -> Add To Path -> Selected Folders and Subfolders
  Option 2: Home tab in MATLAB -> Environment group: Set Path -> Add with Subfolders -> Select the folder in the input dialog -> Save -> Close
  ```
  - Run the 'ExampleScript.m' file
  - Run the 'ExampleScript3D.m' file
  ```
  In the Command Window, type: ExampleScript.m
  In the Command Window, type: ExampleScript3D.m
  ```
A typical "installation" should not take you longer than a minute.

2D example script: The run time on a i7-12700H 2.30 GHz (32 GB RAM) laptop with Windows 64-bit is ~50 minutes. This is to classify the example data that is included in this software (500 clusters of 5 different biological structures).
3D example script: The run time on a i7-12700H 2.30 GHz (32 GB RAM) laptop with Windows 64-bit is ~70 minutes. This is to classify the example data that is included in this software (500 clusters of 2 different biological structures).

**To take full advantage of the available algorithms implemented in this software, please install the PLS Toolbox (www.eigenvector.com). A free 45-days trial is available using an academic email address during signup. We refer to their website for help on the installation of that toolbox.**

If the PLS toolbox is not installed, the software can still be used (with less capabilities). See the User Guide for more information.

## More information
We refer to the User Guide for more information on the different steps of the software.

Please refer to: Hugelier et al. BioRxiv 2023 for more information about this work: DOI: https://doi.org/10.1101/2023.05.10.540077.
