% DCT image denoising.

# ABOUT

* Author    : Guoshen Yu <yu@cmap.polytechnique.fr>
* Author    : Guillermo Sapiro <guille@umn.edu>
* Copyright : (C) 2010 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

# OVERVIEW

This source code provides an implementation of the DCT image denoising.

# UNIX/LINUX/MAC USER GUIDE

The code is compilable on Unix/Linux and Mac OS. 

- Compilation. 
Automated compilation requires the make program.

- Library. 
This code requires the libpng library. You can automatically download, 
compile and include this library to the compiled program by adding the 
LOCAL_LIBS=1 option to the make commands.

- Image format. 
Only the PNG format is supported. 
 
-------------------------------------------------------------------------
Usage:
1. Download the code package and extract it. Go to that directory. 

2. Compile the source code (on Unix/Linux/Mac OS). 
There are two ways to compile the code. 
(1) RECOMMENDED, with Open Multi-Processing multithread parallelization 
(http://openmp.org/). Roughly speaking, it accelerates the program using the 
multiple processors in the computer. Run
make OMP=1

OR
(2) If the complier does not support OpenMp, run 
make

ATTENTION:
If libpng (the official PNG reference library) is not installed in your computer, 
an option LOCAL_LIBS=1 should be added after make. Example
make OMP=1 LOCAL_LIBS=1
The compilation will automatically download and compile libpng and zlib and include 
the library to the program.

3. Run DCT image denoising.
./demo_DCTdenoising 
 
Example, run
./demo_DCTdenoising cinput.png 10 ImNoisy.png ImDenoised.png ImDiff.png


# ABOUT THIS FILE

Copyright 2010 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
