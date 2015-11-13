/*
 * Copyright (c) 2015, Gabriele Facciolo <gfacciol@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string.h>

#include "DCTdenoising.h"
extern "C"{ 
#include "iio.h"
}

using namespace std;

/**
 * @brief write image
 *
 * @param name : path+name+extension of the image
 * @param img : vector which contains the image
 * @param width, height, chnls : size of the image
 * @param truncate8bit: clip the range to [0-255]
 *
 * @return EXIT_SUCCESS if the image has been saved, EXIT_FAILURE otherwise
 **/
int save_image(
    char* name
,   std::vector<float> &img
,   const unsigned width
,   const unsigned height
,   const unsigned chnls
,   bool truncate8bit
){
    //! Allocate Memory
    float* tmp = new float[width * height * chnls];

    //! Check for boundary problems
    for (unsigned k = 0; k < width * height * chnls; k++)
       if(truncate8bit)
        tmp[k] = (img[k] > 255.0f ? 255.0f : (img[k] < 0.0f ? 0.0f : img[k]));
       else
        tmp[k] = img[k];

    iio_save_image_float_split(name, tmp, width, height, chnls);


    //! Free Memory
    delete[] tmp;

    return EXIT_SUCCESS;
}

 /**
  * @brief Load image, check the number of channels
  *
  * @param name : name of the image to read
  * @param img : vector which will contain the image : R, G and B concatenated
  * @param width, height, chnls : size of the image
  *
  * @return EXIT_SUCCESS if the image has been loaded, EXIT_FAILURE otherwise
  **/
int load_image(
    char* name
,   vector<float> &img
,   unsigned * width
,   unsigned * height
,   unsigned * chnls
){
    //! read input image
	cout << endl << "Read input image...";
	size_t h, w, c;
	float *tmp = NULL;
   int ih, iw, ic;

   tmp = iio_read_image_float_split(name, &iw, &ih, &ic);
   w=iw; h=ih; c=ic;
	if (!tmp)
	{
		cout << "error :: " << name << " not found or not a correct image" << endl;
		return EXIT_FAILURE;
	}
	cout << "done." << endl;

	//! test if image is really a color image and exclude the alpha channel
	if (c > 2)
	{
	    unsigned k = 0;
	    while (k < w * h && tmp[k] == tmp[w * h + k] && tmp[k] == tmp[2 * w * h + k])
            k++;
        c = (k == w * h ? 1 : 3);
	}

	//! Some image informations
	cout << "image size :" << endl;
	cout << " - width          = " << w << endl;
	cout << " - height         = " << h << endl;
	cout << " - nb of channels = " << c << endl;

	//! Initializations
	*width  = w;
	*height = h;
	*chnls  = c;
	img.resize(w * h * c);
	for (unsigned k = 0; k < w * h * c; k++)
        img[k] = tmp[k];

    return EXIT_SUCCESS;
}

int nextPowOf2(int n) {
   int k = 1;
   while (k < n)  k *= 2;
   return k;
}


// c: pointer to original argc
// v: pointer to original argv
// o: option name after hyphen
// d: default value (if NULL, the option takes no argument)
const char *pick_option(int *c, char **v, const char *o, const char *d) {
  int id = d ? 1 : 0;
  for (int i = 0; i < *c - id; i++) {
    if (v[i][0] == '-' && 0 == strcmp(v[i] + 1, o)) {
      char *r = v[i + id] + 1 - id;
      for (int j = i; j < *c - id; j++)
        v[j] = v[j + id + 1];
      *c -= id + 1;
      return r;
    }
  }
  return d;
}

/**
 * @file   main.cpp
 * @brief  Main executable file. 
 *
 * @author Gabriele Facciolo
 */
int main(int argc, char **argv)
{
   //! Variables initialization
   unsigned int dctsz = atoi(pick_option(&argc, argv, "w", "16"));
//   //next power of two? 
//   dctsz = nextPowOf2(dctsz);
//   cout << "patch size: " << dctsz << endl;

   //! Check if there is the right call for the algorithm
   if (argc < 4) {
      cerr << "usage: " << argv[0] << " input sigma output\n\
         [-w {8,16} DCT window size 8x8, 16x16 (default)]" << endl;
      return EXIT_FAILURE;
   }

   //! Declarations
   vector<float> img_noisy, img_denoised;
   unsigned width, height, chnls;

   //! Load image
   if(load_image(argv[1], img_noisy, &width, &height, &chnls) != EXIT_SUCCESS)
      return EXIT_FAILURE;

   float fSigma = atof(argv[2]);

   //! Denoising
   img_denoised.resize(img_noisy.size());
   DCTdenoising(img_noisy, img_denoised, width, height, chnls, fSigma, dctsz);

   //! save noisy, denoised and differences images
   cout << endl << "Save images...";

   if (save_image(argv[3], img_denoised, width, height, chnls, false) != EXIT_SUCCESS)
      return EXIT_FAILURE;

   cout << "done." << endl;

   return EXIT_SUCCESS;
}
