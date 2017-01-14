% Multiscale DCT image denoising.

# ABOUT

* Modified version (2016) : Nicola Pierazzo   <nicolapierazzo@gmail.com>
* Modified version (2016) : Gabriele Facciolo <gfacciol@gmail.com>

* Original Author    : Guoshen Yu <yu@cmap.polytechnique.fr>
* Original Author    : Guillermo Sapiro <guille@umn.edu>
* Copyright : (C) 2017 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides an implementation of the DCT image denoising.

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

3. Run DCT image denoising.
./dctdenoising
 
Example, run
./dctdenoising sigma [input [output]] [-1 | -2 guide] [-w patch_size (default 16)] [-c factor] [-n scales] [-single file]

./dctdenoising ../noisy.tiff 15 denoised.png


Use PVFLIP (https://github.com/gfacciol/pvflip) to visualize tiff (float) images.


# ABOUT THIS FILE
Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
