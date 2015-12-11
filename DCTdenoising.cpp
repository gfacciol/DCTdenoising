/*
 * Original code Copyright (c) 2010, Guoshen Yu <yu@cmap.polytechnique.fr>,
 *                                   Guillermo Sapiro <guille@umn.edu>
 * Modified code Copyright (c) 2015, Gabriele Facciolo <gfacciol@gmail.com>
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


/*--------------------------- DCTdenoising  -------------------------*/
// This code implements "DCT image denoising: a simple and effective image 
// denoising algorithm".
// http://www.ipol.im/pub/algo/ys_dct_denoising
// Copyright, Guoshen Yu, Guillermo Sapiro, 2010.
// Please report bugs and/or send comments to Guoshen Yu 
// yu@cmap.polytechnique.fr
/*---------------------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>
#include "DCTdenoising.h"

#include <fftw3.h>
#include <stdlib.h>

extern "C" {
#include "smapa.h"
}

SMART_PARAMETER_FLOAT(HARD_THRESHOLD,2.7) // hard threshold value from BM3D

#ifdef _OPENMP
#include <omp.h>
#endif


// Thread pool and intermeiary variables
fftwf_plan plan_forward[100];
fftwf_plan plan_backward[100];
float* dataspace[100];
float* datafreq[100];

int _DCT2D_psz = -1;
int _DCT2D_nch = -1;


void _DCT2D_init(int sz, int nch) {
   _DCT2D_psz = sz;
   _DCT2D_nch = nch;
   unsigned int w = sz, h = sz;
   unsigned int N = w*h*nch;

   int nthreads = 1; 

#ifdef _OPENMP 
   nthreads = omp_get_max_threads();
   if ( nthreads > 100 ) exit(1);
#endif

   for (int i=0; i<nthreads; i++) {
      dataspace[i] = (float*) fftwf_malloc(sizeof(float) * N);
      datafreq[i]  = (float*) fftwf_malloc(sizeof(float) * N);

      int n[] = {(int)w, (int)h};
      fftwf_r2r_kind dct2[] = {FFTW_REDFT10, FFTW_REDFT10};
      plan_forward[i] = fftwf_plan_many_r2r(2, n, nch, dataspace[i], NULL, 1, w * h,
            datafreq[i], NULL, 1, w * h,
            dct2, FFTW_ESTIMATE);

      fftwf_r2r_kind idct2[] = {FFTW_REDFT01, FFTW_REDFT01};
      plan_backward[i] = fftwf_plan_many_r2r(2, n, nch, datafreq[i], NULL, 1, w * h,
            dataspace[i], NULL, 1, w * h,
            idct2, FFTW_ESTIMATE);
   }
}


void _DCT2D_end() {
   int nthreads = 1; 

#ifdef _OPENMP 
   nthreads = omp_get_max_threads();
   if ( nthreads > 100 ) exit(1);
#endif

   for (int i=0; i<nthreads; i++) {
      free(dataspace[i]);
      free(datafreq[i]);
      fftwf_destroy_plan(plan_forward[i]);
      fftwf_destroy_plan(plan_backward[i]);
   }
}


void _DCT2D(vector< float > & patch, int invflag) {
   int tid = 0;
#ifdef _OPENMP 
   tid = omp_get_thread_num();
#endif

   int psz = _DCT2D_psz;
   int nch = _DCT2D_nch;
   int N = psz*psz*nch;

   fftwf_plan p;
   float *idata, *odata;
   if (invflag>0) { 
      p = plan_forward[tid];
      idata = dataspace[tid];
      odata = datafreq[tid];

      float norm = 1.0/(2*psz);
      float isqrt2 = 1.0/sqrt(2.0);

      for (int i=0; i<N; i++) 
         idata[i] = patch[i];

      fftwf_execute(p);

      for (int i=0; i<N; i++) 
         patch[i] = odata[i] * norm;

      for (int t=0; t<nch; t++) 
         for (int i=0; i<psz; i++) {
            patch[t*psz*psz + i    ] *= isqrt2;
            patch[t*psz*psz + i*psz] *= isqrt2;
         }
   } else {
      p = plan_backward[tid];
      idata = datafreq[tid];
      odata = dataspace[tid];

      float norm = 1.0/(2*psz);
      float sqrt2 = sqrt(2);

      for (int t=0; t<nch; t++) 
         for (int i=0; i<psz; i++) {
            patch[t*psz*psz + i    ] *= sqrt2;
            patch[t*psz*psz + i*psz] *= sqrt2;
         }

      for (int i=0; i<N; i++) 
         idata[i] = patch[i] * norm;

      fftwf_execute(p);

      for (int i=0; i<N; i++) 
         patch[i] = odata[i];
   }
}




#define ABS(x)    (((x) > 0) ? (x) : (-(x)))

// Define a 3x3 1D DCT basis (each ROW is a vector of the forward
// transform basis).
const float DCTbasis3x3[3][3] = {
    {   0.5773502588272094726562500000000000000000,
        0.5773502588272094726562500000000000000000,
        0.5773502588272094726562500000000000000000,     },

    {   0.7071067690849304199218750000000000000000,
      0.0000000000000000000000000000000000000000,
      -0.7071067690849304199218750000000000000000, },

    {
        0.4082483053207397460937500000000000000000,
        -0.8164966106414794921875000000000000000000,
        0.4082483053207397460937500000000000000000      }
};

void ColorTransform(vector<float>&, vector<float>&, int, int, int);


void extract_patch(vector<float> &im, int width, int height, int channel, 
      int i, int j, vector<float> &patch, int width_p, int height_p) {
   // loop over the pixels in the patch
   int size1 = width * height;
   for (int kp = 0; kp < channel; kp++)
      for (int jp = 0; jp < height_p; jp ++)
         for (int ip = 0; ip < width_p; ip ++) {
            patch[kp*width_p*height_p + jp*width_p + ip] = 
               im[kp*size1 + (j+jp)*width + i + ip];
         }
}

// Denoise an image with sliding DCT thresholding.
// ipixels, opixels: noisy and denoised images.
// width, height, channel: image width, height and number of channels.
// sigma: standard deviation of Gaussian white noise in ipixels.
void DCTdenoising(vector<float>& ipixels, vector<float>& opixels, int width,
      int height, int channel, float sigma, int dct_sz)
{
   // Threshold
   float Th = HARD_THRESHOLD() * sigma; 

   // DCT window size
   const int width_p=dct_sz, height_p=dct_sz;

   _DCT2D_init(dct_sz, channel);

   std::vector<float> tpixels = ipixels;
   if (channel == 3) {
      tpixels.resize(width*height*channel);
      // 3-point DCT transform in the color dimension
      ColorTransform(ipixels, tpixels, width, height, 1);
   } 

   const int size1 = width * height;
   const int size = size1 * channel;

   // clean the image
   std::vector<float> im(size);
   std::vector<float> im_weight(size);
   for (int i = 0; i < size; i ++)
      im[i] = im_weight[i] = 0;

   // Loop over the patch positions
#pragma omp parallel for 
   for (int j = 0; j < height - height_p + 1; j ++)
      for (int i = 0; i < width - width_p + 1; i ++) {

         vector< float >  patch(channel*height_p*width_p);

         // extract one patch
         extract_patch(tpixels, width, height, channel, i, j, patch, width_p, height_p);

         // 2D DCT forward
         _DCT2D(patch, 1);

         // Thresholding
         double patch_weight = 0;
         for (int kp = 0; kp < channel; kp ++)
            for (int jp = 0; jp < height_p; jp ++)
               for (int ip = 0; ip < width_p; ip ++) {
                  int idx = kp*width_p*height_p + jp*width_p + ip;
                  if ( ABS(patch[idx]) < Th )
                     patch[idx] = 0;
                  else
                     patch_weight ++;
               }
         // patch weights
         //patch_weight = 1.0/fmax(1, patch_weight);
         patch_weight = 1.0;

         // 2D DCT inverse
         _DCT2D(patch, -1);

         for (int kp = 0; kp < channel; kp ++)
            for (int jp = 0; jp < height_p; jp ++)
               for (int ip = 0; ip < width_p; ip ++) {
                  int idx = kp*width_p*height_p + jp*width_p + ip;
                  int idxim = kp*size1 + (j+jp)*width + i + ip;
                  im[idxim] += patch_weight * patch[idx];
                  //im_weight[idxim] ++;
                  im_weight[idxim] += patch_weight;
               }

      }


   // Normalize by the weight
   for (int i = 0; i < size; i ++)
      im[i] = im[i] / im_weight[i];


   // If the image is colored (3 channels), reverse the 3-point DCT transform
   // in the color dimension.
   if (channel == 3) {
      // inverse 3-point DCT transform in the color dimension
      ColorTransform(im, opixels, width, height, -1);
   } else {
      opixels = im;
   }

   _DCT2D_end();
}



// Denoise an image with sliding DCT thresholding.
// ipixels, gpixels, opixels: noisy, guide and denoised images.
// width, height, channel: image width, height and number of channels.
// sigma: standard deviation of Gaussian white noise in ipixels.
void DCTdenoisingGuided(vector<float>& ipixels, vector<float>& gpixels,  vector<float>& opixels, int width,
      int height, int channel, float sigma, int dct_sz)
{
   // Threshold
   float sigma2 = sigma * sigma;

   // DCT window size
   const int width_p=dct_sz, height_p=dct_sz;

   _DCT2D_init(dct_sz, channel);

   std::vector<float> tpixels = ipixels;
   std::vector<float> tgpixels = gpixels;
   if (channel == 3) {
      tpixels.resize(width*height*channel);
      tgpixels.resize(width*height*channel);
      // 3-point DCT transform in the color dimension
      ColorTransform(ipixels, tpixels, width, height, 1);
      ColorTransform(gpixels, tgpixels, width, height, 1);
   } 

   const int size1 = width * height;
   const int size = size1 * channel;

   // clean the image
   std::vector<float> im(size);
   std::vector<float> im_weight(size);
   for (int i = 0; i < size; i ++)
      im[i] = im_weight[i] = 0;

   // Loop over the patch positions
#pragma omp parallel for 
   for (int j = 0; j < height - height_p + 1; j ++)
      for (int i = 0; i < width - width_p + 1; i ++) {

         vector< float >  patch(channel*height_p*width_p);
         vector< float >  gpatch(channel*height_p*width_p);

         // extract one patch
         extract_patch(tpixels, width, height, channel, i, j, patch, width_p, height_p);
         extract_patch(tgpixels, width, height, channel, i, j, gpatch, width_p, height_p);

         // 2D DCT forward
         _DCT2D(patch, 1);
         _DCT2D(gpatch, 1);

         // Wiener
         double patch_weight = 0;
         for (int kp = 0; kp < channel; kp ++)
            for (int jp = 0; jp < height_p; jp ++)
               for (int ip = 0; ip < width_p; ip ++) {
                  int idx = kp*width_p*height_p + jp*width_p + ip;
                  float G2 = gpatch[idx] * gpatch[idx];
                  float w = G2 / ( G2 + sigma2 );
                  patch[idx] = patch[idx] * w;
                  patch_weight += w*w;
               }
         // patch weights
         //patch_weight = 1.0/patch_weight;
         patch_weight = 1.0;

         // 2D DCT inverse
         _DCT2D(patch, -1);

         for (int kp = 0; kp < channel; kp ++)
            for (int jp = 0; jp < height_p; jp ++)
               for (int ip = 0; ip < width_p; ip ++) {
                  int idx = kp*width_p*height_p + jp*width_p + ip;
                  int idxim = kp*size1 + (j+jp)*width + i + ip;
                  im[idxim] += patch_weight * patch[idx];
                  //im_weight[idxim] ++;
                  im_weight[idxim] += patch_weight;
               }

      }


   // Normalize by the weight
   for (int i = 0; i < size; i ++)
      im[i] = im[i] / im_weight[i];


   // If the image is colored (3 channels), reverse the 3-point DCT transform
   // in the color dimension.
   if (channel == 3) {
      // inverse 3-point DCT transform in the color dimension
      ColorTransform(im, opixels, width, height, -1);
   } else {
      opixels = im;
   }

   _DCT2D_end();
}




// Do a 3-point DCT transform in the image color dimension.
// flag=1/-1 --> forward/inverse transform
void ColorTransform(vector<float>& in, vector<float>& out, int width, 
                    int height, int flag)
{
    int size1 = width * height;

    int j = 0;

    // forward transform
    if ( flag == 1 ) {
#pragma omp parallel for private(j)
        for (j = 0; j < height; j ++)
            for (int i = 0; i < width; i ++) {
                int idx_pixel0 = j*width + i;
                int idx_pixel1 = 1*size1 + j*width + i;
                int idx_pixel2 = 2*size1 + j*width + i;
                out[idx_pixel0] =
                    (  in[idx_pixel0] * DCTbasis3x3[0][0]
                       + in[idx_pixel1] * DCTbasis3x3[0][1]
                       + in[idx_pixel2] * DCTbasis3x3[0][2] );

                out[idx_pixel1] =
                    (  in[idx_pixel0] * DCTbasis3x3[1][0]
                       + in[idx_pixel1] * DCTbasis3x3[1][1]
                       + in[idx_pixel2] * DCTbasis3x3[1][2] );

                out[idx_pixel2] =
                    (  in[idx_pixel0] * DCTbasis3x3[2][0]
                       + in[idx_pixel1] * DCTbasis3x3[2][1]
                       + in[idx_pixel2] * DCTbasis3x3[2][2] );
            }

    }
    // reverse transform
    else if (flag == -1) {
#pragma omp parallel for private(j)
        for (int j = 0; j < height; j ++)
            for (int i = 0; i < width; i ++) {
                int idx_pixel0 = j*width + i;
                int idx_pixel1 = 1*size1 + j*width + i;
                int idx_pixel2 = 2*size1 + j*width + i;
                out[idx_pixel0] =
                    (  in[idx_pixel0] * DCTbasis3x3[0][0]
                       + in[idx_pixel1] * DCTbasis3x3[1][0]
                       + in[idx_pixel2] * DCTbasis3x3[2][0] );

                out[idx_pixel1] =
                    (  in[idx_pixel0] * DCTbasis3x3[0][1]
                       + in[idx_pixel1] * DCTbasis3x3[1][1]
                       + in[idx_pixel2] * DCTbasis3x3[2][1] );

                out[idx_pixel2] =
                    (  in[idx_pixel0] * DCTbasis3x3[0][2]
                       + in[idx_pixel1] * DCTbasis3x3[1][2]
                       + in[idx_pixel2] * DCTbasis3x3[2][2] );
            }
    } else {
        printf ("Error: ColorTransform flag should be 1 (forward) or -1 \
 (inverse). \n");
        //exit (1);
    }
}
