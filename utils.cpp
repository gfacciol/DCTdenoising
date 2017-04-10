//
// Created by Nicola Pierazzo on 21/10/15.
//

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <utility>
#include <fftw3.h>
#include "utils.hpp"

using std::pair;
using std::sqrt;
using std::vector;
using std::move;
using std::max;
using std::min;

namespace imgutils {

inline int SymmetricCoordinate(int pos, int size) {
  if (pos < 0) pos = -pos - 1;
  if (pos >= 2 * size) pos %= 2 * size;
  if (pos >= size) pos = 2 * size - 1 - pos;
  return pos;
}

pair<int, int> ComputeTiling(int rows, int columns, int tiles) {
  float best_r = sqrt(static_cast<float>(tiles * rows) / columns);
  int r_low = static_cast<int>(best_r);
  int r_up = r_low + 1;
  if (r_low < 1) return {1, tiles};
  if (r_up > tiles) return {tiles, 1};
  while (tiles % r_low != 0) --r_low;
  while (tiles % r_up != 0) ++r_up;
  if (r_up * r_low * columns > tiles * rows) {
    return {r_low, tiles / r_low};
  } else {
    return {r_up, tiles / r_up};
  }
}

vector<Image> SplitTiles(const Image &src, int pad_before, int pad_after,
                         pair<int, int> tiling) {
  vector<Image> result;
  for (int tr = 0; tr < tiling.first; ++tr) {
    int rstart = src.rows() * tr / tiling.first - pad_before;
    int rend = src.rows() * (tr + 1) / tiling.first + pad_after;
    for (int tc = 0; tc < tiling.second; ++tc) {
      int cstart = src.columns() * tc / tiling.second - pad_before;
      int cend = src.columns() * (tc + 1) / tiling.second + pad_after;
      Image tile(rend - rstart, cend - cstart, src.channels());
      for (int ch = 0; ch < src.channels(); ++ch) {
        for (int row = rstart; row < rend; ++row) {
          for (int col = cstart; col < cend; ++col) {
            tile.val(col - cstart, row - rstart, ch) = src.val(
                SymmetricCoordinate(col, src.columns()),
                SymmetricCoordinate(row, src.rows()),
                ch);
          }
        }
      }
      result.push_back(move(tile));
    }
  }
  return result;
}

Image MergeTiles(const vector<pair<Image, Image>> &src, pair<int, int> shape,
                 int pad_before, int pad_after, pair<int, int> tiling) {
  int channels = src[0].first.channels();
  Image result(shape.first, shape.second, channels);
  Image weights(shape.first, shape.second);
  auto tile = src.begin();
  for (int tr = 0; tr < tiling.first; ++tr) {
    int rstart = shape.first * tr / tiling.first - pad_before;
    int rend = shape.first * (tr + 1) / tiling.first + pad_after;
    for (int tc = 0; tc < tiling.second; ++tc) {
      int cstart = shape.second * tc / tiling.second - pad_before;
      int cend = shape.second * (tc + 1) / tiling.second + pad_after;
      for (int ch = 0; ch < channels; ++ch) {
        for (int row = max(0, rstart); row < min(shape.first, rend); ++row) {
          for (int col = max(0, cstart); col < min(shape.second, cend); ++col) {
            result.val(col, row, ch) +=
                tile->first.val(col - cstart, row - rstart, ch);
          }
        }
      }
      for (int row = max(0, rstart); row < min(shape.first, rend); ++row) {
        for (int col = max(0, cstart); col < min(shape.second, cend); ++col) {
          weights.val(col, row) += tile->second.val(col - cstart, row - rstart);
        }
      }
      ++tile;
    }
  }
  for (int ch = 0; ch < channels; ++ch) {
    for (int row = 0; row < shape.first; ++row) {
      for (int col = 0; col < shape.second; ++col) {
        result.val(col, row, ch) /= weights.val(col, row);
      }
    }
  }
  return result;
}

void dct_inplace(Image &img) {
  int n[] = {img.rows(), img.columns()};
  fftwf_r2r_kind dct2[] = {FFTW_REDFT10, FFTW_REDFT10};
  fftwf_plan plan = fftwf_plan_many_r2r(2, n, img.channels(), img.data(), NULL,
                                        1, img.pixels(), img.data(), NULL, 1,
                                        img.pixels(), dct2, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

  // Normalization
  for (int i = 0; i < img.samples(); ++i) {
    img.val(i) /= 4 * img.pixels();
  }
}

void idct_inplace(Image &img) {
  int n[] = {img.rows(), img.columns()};
  fftwf_r2r_kind idct2[] = {FFTW_REDFT01, FFTW_REDFT01};
  fftwf_plan plan = fftwf_plan_many_r2r(2, n, img.channels(), img.data(), NULL,
                                        1, img.pixels(), img.data(), NULL, 1,
                                        img.pixels(), idct2, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);
}

/*! \brief Dyadic DCT pyramid decomposition.
 *  
 *  Returns a vector containing #levels images
 */
vector<Image> decompose(const Image &img, int levels) {
  vector<Image> pyramid;
  Image freq = img.copy();
  dct_inplace(freq);
  int h{freq.rows()}, w{freq.columns()};
  for (int i = 0; i < levels; ++i) {
    // Copy data
    Image layer(h, w, freq.channels());
    for (int ch = 0; ch < freq.channels(); ++ch) {
      for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
          layer.val(c, r, ch) = freq.val(c, r, ch);
        }
      }
    }
    // Inverse DCT
    idct_inplace(layer);
    w /= 2;
    h /= 2;
    pyramid.push_back(move(layer));
  }
  return pyramid;
}

/*! \brief Dyadic DCT pyramid conservative recomposition.
 *  
 *  Takes a vector containing the layers of the pyramid. Computes
 *  the transform of each layer and recomposes the frequencies in a 
 *  fine-to-coarse fashion, copying into the final result only the 
 *  lowest recompose_factor frequencies of each coarse layer.
 *  After inverting the DCT, returns the resulting image.
 */
Image recompose(const vector<Image> &pyramid, float recompose_factor) {
  // Use the bigger image to determine width, height and number of channels
  Image output = pyramid[0].copy();
  // Perform the DCT
  dct_inplace(output);

  for (int i = 1; i < static_cast<int>(pyramid.size()); ++i) {
    // Read level i of the pyramid
    Image layer = pyramid[i].copy();
    // Perform the DCT
    dct_inplace(layer);
    // Copy data (selected by recompose_factor)
    for (int ch = 0; ch < layer.channels(); ++ch) {
      for (int r = 0; r < layer.rows() * recompose_factor; ++r) {
        for (int c = 0; c < layer.columns() * recompose_factor; ++c) {
          output.val(c, r, ch) = layer.val(c, r, ch);
        }
      }
    }
  }
  // IDCT of the output image
  idct_inplace(output);
  return output;
}

}  // namespace imgutils
