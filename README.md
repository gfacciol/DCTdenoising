% Multiscale DCT image denoising.

# ABOUT

* Author : Nicola Pierazzo   <nicolapierazzo@gmail.com>
* Author : Gabriele Facciolo <gfacciol@gmail.com>
* Copyright : (C) 2017 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt
* Based on the 2010 implementation of DCT denoising by:
  Guoshen Yu <yu@cmap.polytechnique.fr> and Guillermo Sapiro <guille@umn.edu>
* Latest version available at: https://github.com/gfacciol/DCTdenoising

# OVERVIEW

This source code provides an implementation of the "Multiscale DCT denoising"
algorithm described in the IPOL article: http://www.ipol.im/pub/art/2017/201

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS. 

- Compilation. 
Automated compilation requires the Cmake and make.

- Dependencies.
This code requires the libpng, libtiff, libjpeg, and libfftw.

- Image formats. 
Only the PNG, JPEG, and TIFF (float) formats are supported. 
 
-------------------------------------------------------------------------
Usage:
1. Download the code package and extract it. Go to that directory. 

2. Compile the source code (on Unix/Linux/Mac OS). 

    mkdir build; cd build;
    cmake ..; make;

3. Runing DCT image denoising: parameters
 
    ./dctdenoising sigma [input [output]]   # noise std, noisy image, output  
       [-w patch_size (default 8)]   # DCT denoising patch size  
       [-1 | -2 guide]               # -1: only hard thresh., -2: provide guide  
       [-no_adaptive_aggregation]    # disable adaptive aggregation weights  
       [-n scales(4)]                # multiscale: number of scales  
       [-c factor(.5)]               # multiscale: recomposition factor  
       [-single output_singlescale]  # multiscale: save also one-scale result  


The flag -1 permits to run DCT denoising only with the hard thresholding step,
while -2 allows to specify the guide for the wiener filtering step.
When not set both steps are executed using the output of the first one as guide
for the second one.
Setting no_adaptive_aggregation disables the aggregation weights.


Example, run

    ./dctdenoising 15 ../noisy.tiff denoised.png


To visualize tiff (float) images use PVFLIP (https://github.com/gfacciol/pvflip) 
or ImageJ (https://imagej.nih.gov/ij/index.html)


# ABOUT THIS FILE
Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
