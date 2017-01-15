% Multiscale DCT image denoising.

# ABOUT

* Author : Nicola Pierazzo   <nicolapierazzo@gmail.com>
* Author : Gabriele Facciolo <gfacciol@gmail.com>
* Copyright : (C) 2017 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt
* Based on the 2010 implementation of DCT denoising by:
  Guoshen Yu <yu@cmap.polytechnique.fr> and Guillermo Sapiro <guille@umn.edu>

# OVERVIEW

This source code provides an implementation of the Multiscale DCT image denoising.

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

3. Run DCT image denoising
 
./dctdenoising sigma [input [output]] [-1 | -2 guide] [-w patch_size (default 16)] 
               [-c factor] [-n scales] [-single file] [-no_adaptive_aggregation]

Example, run
./dctdenoising 15 ../noisy.tiff denoised.png


To visualize tiff (float) images use PVFLIP (https://github.com/gfacciol/pvflip) 
or ImageJ (https://imagej.nih.gov/ij/index.html)


# ABOUT THIS FILE
Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
