# Multi-temporal function-based dasymetric model for mobile phone data

#### This repository provides supporting information to the following article:

 - Järv Olle, Tenkanen Henrikki & Toivonen Tuuli (2017). *Enhancing spatial accuracy of mobile
  phone data using multi-temporal dasymetric interpolation.* Published in **International Journal of 
  Geographical Information Science**.

Copyright (C) 2017.  Digital Geography Group / Accessibility Research Group, University of Helsinki (Järv, Tenkanen & Toivonen).
Programmed by: Henrikki Tenkanen, University of Helsinki, Finland.

#### Check the latest version!

This repository is the released version (1.0) that original one belonging to the article mentioned [above](README.md#this-repository-provides-supporting-information-to-the-following-article). However, as things might change and softwares tend to develop further, we highly recommend checking the [latest version from here](https://github.com/AccessibilityRG/MFD_Mobile).

## Purpose:
  
This repository and script(s) demonstrate in practice how multi-temporal dasymetric interpolation model can be used to interpolate population distributed in spatially varying coverage areas of mobile phone base stations into a desired layout of predefined statistical units using ancillary data sources. 

## Requirements:

Python 3 with following packages and their dependencies: pandas, geopandas.

### Installations

**We highly recommend using** [Anaconda](https://www.continuum.io/anaconda-overview) which is an open source distribution of the Python programming language for large-scale data processing, predictive analytics, and scientific computing, that aims to simplify package management and deployment. Pandas comes with the basic installation package of Anaconda:

 - Download [Anaconda installer (64 bit)](https://www.continuum.io/downloads) for Windows / Linux / Mac.
 
 - Install Geopandas and it's dependencies with [conda](http://conda.pydata.org/docs/using/using.html) in your terminal or command prompt:
 
    ```
    conda install -y -c ioos geopandas=0.2.1
    conda install -y -c ioos gdal=2.1.2
    ```

## Data:

This script is planned to work with data that is described in the main article, see Järv et al. (2017). In short, the script requires three datasets:
  
  1) a disaggregated physical surface layer, 
  2) time-dependent human activity data,
  3) mobile phone data (CDR), or similar. 
  
See more information about the datasets from [here]()

## Contact:

You can contact Olle Järv or Henrikki Tenkanen if you have any questions related to the model or the Python implementation of the model:
  
 - olle.jarv (a) helsinki.fi
 - henrikki.tenkanen (a) helsinki.fi
 - http://www.helsinki.fi/science/accessibility

## License:

mfd_interpolation.py by Accessibility Research Group (University of Helsinki) is licensed under a Creative Commons Attribution 4.0 International License. More information about license: http://creativecommons.org/licenses/by/4.0/


