/*
 * Copyright (c) 2010, Guoshen Yu <yu@cmap.polytechnique.fr>,
 *                     Guillermo Sapiro <guille@umn.edu>
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


/*----------------------------------- DCT2D ---------------------------------*/
// Detect corresponding points in two images with the ASIFT method.
// Copyright, Guoshen Yu, Guillermo Sapiro, 2010.
// Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
/*---------------------------------------------------------------------------*/

#include <vector>
#include <iostream>
using namespace std;

void DCT2D16x16(vector< vector< float > >&, int);
