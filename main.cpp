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
#include <cmath>
#include <cstring>

#include "DCTdenoising.h"
#include "utils.hpp"

using imgutils::pick_option;
using imgutils::read_image;
using imgutils::save_image;
using imgutils::Image;
using std::cerr;
using std::endl;
using std::move;

/**
 * @file   main.cpp
 * @brief  Main executable file. 
 *
 * @author Gabriele Facciolo
 * @author Nicola Pierazzo
 */
int main(int argc, char **argv) {
  const bool usage = static_cast<bool>(pick_option(&argc, argv, "h", nullptr));
  const int dct_sz = atoi(pick_option(&argc, argv, "w", "16"));
  const bool no_second_step = static_cast<bool>(pick_option(&argc, argv, "1", NULL));
  const char *second_step_guide = pick_option(&argc, argv, "2", "");
  const bool no_first_step = second_step_guide[0] != '\0';

  //! Check if there is the right call for the algorithm
  if (usage || argc < 4) {
    cerr << "usage: " << argv[0] << " input sigma output [-1 | -2 guide] " <<
        "[-w patch_size (default 16)]" << endl;
    return usage ? EXIT_SUCCESS : EXIT_FAILURE;
  }

  if (no_second_step && no_first_step) {
    cerr << "You can't use -1 and -2 together." << endl;
    return EXIT_FAILURE;
  }

  Image noisy = read_image(argv[1]);
  Image guide, result;
  const float sigma = static_cast<float>(atof(argv[2]));

  if (!no_first_step) {
    guide = DCTdenoising(noisy, sigma, dct_sz);
  } else {
    guide = read_image(second_step_guide);
  }

  if (!no_second_step) {
    result = DCTdenoisingGuided(noisy, guide, sigma, dct_sz);
  } else {
    result = move(guide);
  }

  save_image(result, argv[3]);

  return EXIT_SUCCESS;
}
