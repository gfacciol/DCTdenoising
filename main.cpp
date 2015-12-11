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
	cout << "Read input image "<< name << "...";
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

	//! test if image is really a color image and exclude the alpha channel
	if (c > 2)
	{
	    unsigned k = 0;
	    while (k < w * h && tmp[k] == tmp[w * h + k] && tmp[k] == tmp[2 * w * h + k])
            k++;
        c = (k == w * h ? 1 : 3);
	}

	//! Some image informations
	cout << "(" << w <<  "x" << h << "x" << c << ")" << endl;

	//! Initializations
	*width  = w;
	*height = h;
	*chnls  = c;
	img.resize(w * h * c);
	for (unsigned k = 0; k < w * h * c; k++)
        img[k] = tmp[k];

    return EXIT_SUCCESS;
}

int sym(int x, int nx) {
   return  (x)<0 ?  \
          ( -(x)-1 >= nx  ? nx-1: -(x)-1 )  :  \
          ( (x) >= nx ? (-(x)+2*nx-1 < 0 ? 0 : -(x)+2*nx-1)   : (x)  ) ;
}

vector<float> padimage (vector<float> in, unsigned int *w, unsigned int *h, unsigned int *c, unsigned int pad) {
   int outw=*w+2*pad;
   int outh=*h+2*pad;
   int outc=*c;

   vector<float> out (outc*outw*outh);

   for(int k=0; k<outc; k++)
   for(int j=0; j<outh; j++)
   for(int i=0; i<outw; i++) {
      out[i + j*outw + k*outw*outh] = 
         in[ sym(i-pad,*w) + sym(j-pad,*h)**w + k**w**h];
   }
   *w=outw;
   *h=outh;
   return out;
}


vector<float> unpadimage (vector<float> in, unsigned int *w, unsigned int *h, unsigned int *c, unsigned int pad) {
   int outw=*w-2*pad;
   int outh=*h-2*pad;
   int outc=*c;

   if(outw<=0 || outh<=0) return in;
   vector<float> out (outc*outw*outh);

   for(int k=0; k<outc; k++)
   for(int j=0; j<outh; j++)
   for(int i=0; i<outw; i++) {
      out[i + j*outw + k*outw*outh] = 
         in[ (i+pad) + (j+pad)**w + k**w**h];
   }
   *w=outw;
   *h=outh;
   return out;
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
   bool only1step = pick_option(&argc, argv, "1", NULL) ? true : false;
   const char *guide2ndstep = pick_option(&argc, argv, "2", "");
   bool no1ststep = guide2ndstep[0] != '\0';
   if (no1ststep) only1step = false;

//   //next power of two? 
//   dctsz = nextPowOf2(dctsz);
//   cout << "patch size: " << dctsz << endl;

   //! Check if there is the right call for the algorithm
   if (argc < 4) {
      cerr << "usage: " << argv[0] << " input sigma output\n"
         "[-1         if set disable 2nd step (soft thresholding)]\n"
         "[-2 guide   use the guide image and disable the 1st step]\n"
         "[-w {8,16}  DCT window size 8x8, 16x16 (default)]\n";
      return EXIT_FAILURE;
   }

   //! Declarations
   vector<float> img_noisy, img_denoised, img_guide;
   unsigned width, height, chnls;

   //! Load image
   if(load_image(argv[1], img_noisy, &width, &height, &chnls) != EXIT_SUCCESS)
      return EXIT_FAILURE;
//   img_noisy = padimage(img_noisy, &width, &height, &chnls, dctsz);

   float fSigma = atof(argv[2]);

   //! Denoising
   img_guide.resize(img_noisy.size());
   img_denoised.resize(img_noisy.size());

   if(no1ststep) {
      if(load_image((char *) guide2ndstep, img_guide, &width, &height, &chnls) != EXIT_SUCCESS)
      return EXIT_FAILURE;
      cerr << "only step2 " << endl;
 //     img_guide = padimage(img_guide, &width, &height, &chnls, dctsz);
   } else {
      DCTdenoising(img_noisy, img_guide, width, height, chnls, fSigma, dctsz);
   }

   if(only1step)  {
      img_denoised = img_guide;
      cerr << "only step1 " << endl;
   }
   else {
      DCTdenoisingGuided(img_noisy, img_guide, img_denoised, width, height, chnls, fSigma, dctsz);
   }

   //! save noisy, denoised and differences images
   //cout << "Save images."<< endl;

  // img_denoised = unpadimage(img_denoised, &width, &height, &chnls, dctsz);
   if (save_image(argv[3], img_denoised, width, height, chnls, false) != EXIT_SUCCESS)
      return EXIT_FAILURE;

   return EXIT_SUCCESS;
}
