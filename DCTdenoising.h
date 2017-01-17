/*
 * Code Copyright (c) 2017, Nicola Pierazzo <nicolapierazzo@gmail.com>,
 *                          Gabriele Facciolo <gfacciol@gmail.com>
 * Based on the 2010 article by Guoshen Yu <yu@cmap.polytechnique.fr>,
 *                              Guillermo Sapiro <guille@umn.edu>
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


/*-----------------------  Multiscale DCTdenoising  -------------------------*/
// This code implements "Multiscale DCT denoising".
// http://www.ipol.im/pub/art/2017/201 
// Copyright, Nicola Pierazzo, Gabriele Facciolo, 2017.
// Please report bugs and/or send comments to G. Facciolo gfacciol@gmail.com
/*---------------------------------------------------------------------------*/

#ifndef DCTDENOISING_DCTDENOISING_HPP
#define DCTDENOISING_DCTDENOISING_HPP

#include "Image.hpp"

imgutils::Image DCTdenoising(const imgutils::Image &noisy, float sigma,
                             int dct_size, bool adaptive_aggregation = true, 
                             int nthreads = 0);
imgutils::Image DCTdenoisingGuided(const imgutils::Image &noisy,
                                   const imgutils::Image &guide,
                                   float sigma, int dct_size, 
                                   bool adaptive_aggregation = true, 
                                   int nthreads = 0);

#endif  // DCTDENOISING_DCTDENOISING_HPP
