/*
 * Copyright (c) 2016, Gabriele Facciolo <gfacciol@gmail.com>
 *                     Nicola Pierazzo <nicolapierazzo@gmail.com>
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
#include <cmath>
#include <cstring>

#include "DCTdenoising.h"
#include "utils.hpp"
#include "demoutils.hpp"

using imgutils::pick_option;
using imgutils::read_image;
using imgutils::save_image;
using imgutils::Image;
using std::cerr;
using std::endl;
using std::move;
using std::vector;

/**
 * @file   main.cpp
 * @brief  Main executable file. 
 *
 * @author Gabriele Facciolo
 * @author Nicola Pierazzo
 */
int main(int argc, char **argv) {
  // read DCT denoising options
  const bool  usage = static_cast<bool>(pick_option(&argc, argv, "h", nullptr));
  const int   dct_sz = atoi(pick_option(&argc, argv, "w", "8"));
  const char  *second_step_guide = pick_option(&argc, argv, "2", "");
  const bool  no_second_step = static_cast<bool>(pick_option(&argc, argv, "1", NULL));
  const bool  no_first_step = second_step_guide[0] != '\0';
  const bool  adaptive_aggregation 
    = ! static_cast<bool>(pick_option(&argc, argv, "no_adaptive_aggregation", NULL));
  // read multiscaler options
  const char  *out_single = pick_option(&argc, argv, "single", "");
  const int   scales = atoi(pick_option(&argc, argv, "n", "4"));
  const float recompose_factor
    = static_cast<float>(atof(pick_option(&argc, argv, "c", ".5")));

  //! Check if there is the right call for the algorithm
  if (usage || argc < 2) {
    cerr << "usage: " << argv[0] << " sigma [input [output]] [-1 | -2 guide] "
         << "[-w patch_size (default 8)] [-c factor(.5)] [-n scales(4)] "
         << "[-single output_singlescale] [-no_adaptive_aggregation]" << endl;
    return usage ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  if (no_second_step && no_first_step) {
    cerr << "You can't use -1 and -2 together." << endl;
    return EXIT_FAILURE;
  }

#ifndef _OPENMP
  cerr << "Warning: OpenMP not available. The algorithm will run in a single" <<
       " thread." << endl;
#endif

  // read input
  Image noisy = read_image(argc > 2 ? argv[2] : "-");
  const float sigma = static_cast<float>(atof(argv[1]));
  // generate the DCT pyramid 
  vector<Image> noisy_p = decompose(noisy, scales);
  vector<Image> guide_p, denoised_p;
  if (no_first_step) {
    Image guide = read_image(second_step_guide);
    guide_p = decompose(guide, scales);
  }

  // apply DCT denoising at each scale of the pyramid
  for (int layer = 0; layer < scales; ++layer) {
    // noise at the current scale is proportional to the number of pixels
    float s = sigma * sqrt(static_cast<float>(noisy_p[layer].pixels()) / noisy.pixels());
    if (!no_first_step) {
      Image guide = DCTdenoising(noisy_p[layer], s, dct_sz, adaptive_aggregation);
      guide_p.push_back(move(guide));
    }
    if (!no_second_step) {
      Image result =
          DCTdenoisingGuided(noisy_p[layer], guide_p[layer], s, dct_sz, adaptive_aggregation);
      denoised_p.push_back(move(result));
    } else {
      denoised_p.push_back(move(guide_p[layer]));
    }
  }

  // recompose pyramid
  if (strlen(out_single)) save_image(denoised_p[0], out_single);
  Image result = recompose(denoised_p, recompose_factor);

  save_image(result, argc > 3 ? argv[3] : "TIFF:-");

  return EXIT_SUCCESS;
}
