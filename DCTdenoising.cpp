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


#ifdef _OPENMP
#include <omp.h>
#endif


// Thread pool and intermeiary variables
fftwf_plan plan_forward[100];
fftwf_plan plan_backward[100];
float* dataspace[100];
float* datafreq[100];
int psz = -1;


void _DCT2D_init(int sz, int nch) {
   unsigned int w=sz ,h=sz;
   unsigned int N = w*h*nch;
   psz = sz;

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


void _DCT2D(vector< vector< float > >& patch, int invflag) {
   int tid = 0;
#ifdef _OPENMP 
   tid = omp_get_thread_num();
#endif

   fftwf_plan p;
   float *idata, *odata;
   if (invflag>0) { 
      p = plan_forward[tid];
      idata = dataspace[tid];
      odata = datafreq[tid];

      float norm = 1.0/(2*psz);
      float isqrt2 = 1.0/sqrt(2.0);

      for (int j=0, t=0; j<psz; j++)
         for (int i=0; i<psz; i++, t++) 
            idata[t] = patch[j][i];
   
      fftwf_execute(p);
   
      for (int j=0, t=0; j<psz; j++) 
         for (int i=0; i<psz; i++, t++)
            patch[j][i] = odata[t] * norm;
      
      for (int i=0; i<psz; i++) {
            patch[0][i] *= isqrt2;
            patch[i][0] *= isqrt2;
      }
   } else {
      p = plan_backward[tid];
      idata = datafreq[tid];
      odata = dataspace[tid];

      float norm = 1.0/(2*psz);
      float isqrt2 = sqrt(2);

      for (int i=0; i<psz; i++) {
            patch[0][i] *= isqrt2;
            patch[i][0] *= isqrt2;
      }

      for (int j=0, t=0; j<psz; j++)
         for (int i=0; i<psz; i++, t++) 
            idata[t] = patch[j][i] * norm;
   
      fftwf_execute(p);
   
      for (int j=0, t=0; j<psz; j++) 
         for (int i=0; i<psz; i++, t++)
            patch[j][i] = odata[t];
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

void Image2Patches(vector<float>&, 
     vector< vector< vector< vector< float > > > >&, int, int, int, int, int);
void Patches2Image(vector<float>&, 
     vector< vector< vector< vector< float > > > >&, int, int, int, int, int);
void ColorTransform(vector<float>&, vector<float>&, int, int, int);

// Denoise an image with sliding DCT thresholding.
// ipixels, opixels: noisy and denoised images.
// width, height, channel: image width, height and number of channels.
// sigma: standard deviation of Gaussian white noise in ipixels.
void DCTdenoising(vector<float>& ipixels, vector<float>& opixels, int width,
                  int height, int channel, float sigma, int dct_sz)
{
    // Threshold
    float Th = 3 * sigma;

    // DCT window size
    int width_p=dct_sz, height_p=dct_sz;

    _DCT2D_init(dct_sz, 1);

    int num_patches = (width - width_p + 1) * (height - height_p + 1);

    std::vector< vector< vector< vector< float > > > > patches;
    patches.resize(num_patches);
    for (int p = 0; p < num_patches; p ++) {
        patches[p].resize(channel);
        for (int k = 0; k < channel; k ++) {
            patches[p][k].resize(height_p);
            for (int j = 0; j < height_p; j ++)
                patches[p][k][j].resize(width_p);
        }
    }

    // If the image is colored (3 channels), first decorrelate
    // the color channels with a 3-point DCT transform
    // in the color dimension.
    if (channel == 3) {
        std::vector<float> tpixels;
        tpixels.resize(width*height*channel);

        // 3-point DCT transform in the color dimension
        ColorTransform(ipixels, tpixels, width, height, 1);

        // Decompose the image into patches
        Image2Patches(tpixels, patches, width, height, channel, width_p,
                      height_p);
    } else {
        // Decompose the image into patches
        Image2Patches(ipixels, patches, width, height, channel, width_p,
                      height_p);
    }

    // 2D DCT forward
#pragma omp parallel for private(p)
    for (int p = 0; p < num_patches; p ++) {
        for (int k = 0; k < channel; k ++) {
            _DCT2D(patches[p][k], 1);
        }
    }

    // Thresholding
#pragma omp parallel for private(p)
    for (int p = 0; p < num_patches; p ++)
        for (int k = 0; k < channel; k ++)
            for (int j = 0; j < height_p; j ++)
                for (int i = 0; i < width_p; i ++) {
                    if ( ABS(patches[p][k][j][i]) < Th )
                        patches[p][k][j][i] = 0;
                }

    // 2D DCT inverse
#pragma omp parallel for private(p)
    for (int p = 0; p < num_patches; p ++) {
        for (int k = 0; k < channel; k ++) {
            _DCT2D(patches[p][k], -1);
        }

    }

    // If the image is colored (3 channels), reverse the 3-point DCT transform
    // in the color dimension.
    if (channel == 3) {
        std::vector<float> tpixels;
        tpixels.resize(width*height*channel);

        // Decompose the image into patches
        Patches2Image(tpixels, patches, width, height, channel, width_p,
                      height_p);

        // inverse 3-point DCT transform in the color dimension
        ColorTransform(tpixels, opixels, width, height, -1);
    } else {
        Patches2Image(opixels, patches, width, height, channel, width_p,
                      height_p);
    }
    _DCT2D_end();
}

// Transfer an image im of size width x height x channel to sliding patches of 
// size width_p x height_p xchannel.
// The patches are stored in patches, where each ROW is a patch after being 
// reshaped to a vector.
void Image2Patches(vector<float>& im, 
                   vector< vector< vector< vector< float > > > >& patches, 
                   int width, int height, int channel, int width_p, 
                   int height_p)
{
    int size1 = width * height;

    int counter_patch = 0;

    // Loop over the patch positions
    for (int j = 0; j < height - height_p + 1; j ++)
        for (int i = 0; i < width - width_p + 1; i ++) {
            int counter_pixel = 0;
            // loop over the pixels in the patch
            for (int kp = 0; kp < channel; kp++)
                for (int jp = 0; jp < height_p; jp ++)
                    for (int ip = 0; ip < width_p; ip ++) {
                        patches[counter_patch][kp][jp][ip] = 
                                         im[kp*size1 + (j+jp)*width + i + ip];
                        counter_pixel ++;
                    }
            counter_patch ++;
        }
}

// Transfer sliding patches of size width_p x height_p xchannel to an image im 
// of size width x height x channel.
// The patches are stored in patches, where each ROW is a patch after being 
// reshaped to a vector.
void Patches2Image(vector<float>& im, 
                   vector< vector< vector< vector< float > > > >& patches, 
                   int width, int height, int channel, int width_p, 
                   int height_p)
{
    int size1 = width * height;
    int size = size1 * channel;

    // clean the image
    for (int i = 0; i < size; i ++)
        im[i] = 0;

    // Store the weight
    std::vector<float> im_weight;
    im_weight.resize(size);
    for (int i = 0; i < size; i ++)
        im_weight[i] = 0;

    int counter_patch = 0;

    // Loop over the patch positions
    for (int j = 0; j < height - height_p + 1; j ++)
        for (int i = 0; i < width - width_p + 1; i ++) {
            int counter_pixel = 0;
            // loop over the pixels in the patch
            for (int kp = 0; kp < channel; kp++)
                for (int jp = 0; jp < height_p; jp ++)
                    for (int ip = 0; ip < width_p; ip ++) {
                        im[kp*size1 + (j+jp)*width + i + ip] += 
                                       patches[counter_patch][kp][jp][ip];
                        im_weight[kp*size1 + (j+jp)*width + i + ip] ++;
                        counter_pixel ++;
                    }
            counter_patch ++;
        }

    // Normalize by the weight
    for (int i = 0; i < size; i ++)
        im[i] = im[i] / im_weight[i];
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
