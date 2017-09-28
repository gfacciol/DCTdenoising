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

// use isometric DCT for multiscale decompose/recompose
#define ISOMETRIC_DCT

inline int SymmetricCoordinate(int pos, int size) {
  if (pos < 0) pos = -pos - 1;
  if (pos >= 2 * size) pos %= 2 * size;
  if (pos >= size) pos = 2 * size - 1 - pos;
  return pos;
}

/*! \brief Compute best tiling of image in at most ntiles
 *
 *  Returns the pair: rows, columns
 */
pair<int, int> ComputeTiling(int rows, int columns, int ntiles) {
  // The objective is extract ntiles square-ish tiles
  // For a square image the optimal number of rows is sqrt(ntiles)
  // The ratio rows/columns permits to handle rectangular images
  float best_r = sqrt(static_cast<float>(ntiles * rows) / columns);
  int r_low = static_cast<int>(best_r);
  int r_up = r_low + 1;
  if (r_low < 1) return {1, ntiles};      // single row
  if (r_up > ntiles) return {ntiles, 1};  // single column
  // look for the nearest integer divisors of ntiles
  while (ntiles % r_low != 0) --r_low;
  while (ntiles % r_up != 0) ++r_up;
  // At this point there are two possible tilings:
  //   {r_low, ntiles / r_low} and {r_up, ntiles / r_up}. 
  // We need to select the best. 
  // To do that, we consider the shape of the tiles.
  // In the first case, the tiles are roughly 
  //   {rows / r_low, columns * r_low / ntiles} pixels.
  // In the second case, the tiles are 
  //   {rows / r_up, columns * r_up / ntiles} pixels.
  // Since r_low <= best_r <= r_up the first tile will have i
  // more rows than columns and vice-versa.
  //
  // To select the best case we consider the ratio between the 
  // lengths of the longer and the shorter edge of a tile. 
  // The closer this ratio is to 1, the "squarer" the tile will be. 
  // In other words, we select the first tiling if
  //   (rows / r_low) / (columns * r_low / ntiles) < 
  //        (columns * r_up / ntiles) / (rows / r_up)
  // That is equivalent to (all values are > 0): 
  //   rows * ntiles < r_up * r_low * columns
  if (r_up * r_low * columns > ntiles * rows) {
    return {r_low, ntiles / r_low};
  } else {
    return {r_up, ntiles / r_up};
  }
}

/*! \brief Split image in tiles
 *
 *  Returns a vector containing tiling.first x tiling.sencond images
 *  each padded by pad_* pixels. Tiles are stored in lexicographic order.
 *  Padding outside the image is done by symmetrization
 */
vector<Image> SplitTiles(const Image &src, int pad_before, int pad_after,
                         pair<int, int> tiling) {
  vector<Image> result;
  for (int tr = 0; tr < tiling.first; ++tr) {
    int rstart = src.rows() * tr / tiling.first - pad_before;
    int rend = src.rows() * (tr + 1) / tiling.first + pad_after;
    for (int tc = 0; tc < tiling.second; ++tc) {
      int cstart = src.columns() * tc / tiling.second - pad_before;
      int cend = src.columns() * (tc + 1) / tiling.second + pad_after;
      // copy image to tile using the above computed limits
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

/*! \brief Recompose tiles produced by SplitTiles
 *
 *  Returns an image resulting of recomposing the tiling
 *  padded margins are averaged in the result
 */
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
      // copy tile to image using the above computed limits
      for (int ch = 0; ch < channels; ++ch) {
        for (int row = max(0, rstart); row < min(shape.first, rend); ++row) {
          for (int col = max(0, cstart); col < min(shape.second, cend); ++col) {
            result.val(col, row, ch) +=
                tile->first.val(col - cstart, row - rstart, ch);
          }
        }
      }
      // image weight due to the padding
      for (int row = max(0, rstart); row < min(shape.first, rend); ++row) {
        for (int col = max(0, cstart); col < min(shape.second, cend); ++col) {
          weights.val(col, row) += tile->second.val(col - cstart, row - rstart);
        }
      }
      ++tile;
    }
  }
  // normalize by the weight
  for (int ch = 0; ch < channels; ++ch) {
    for (int row = 0; row < shape.first; ++row) {
      for (int col = 0; col < shape.second; ++col) {
        result.val(col, row, ch) /= weights.val(col, row);
      }
    }
  }
  return result;
}

/*! \brief 2D DCT transform of image (channel-wise)
 *
 * Operates over the img in-place
 * Computes the normalized but not orthogonal 2D DCT
 */
void dct_inplace(Image &img) {
  int n[] = {img.rows(), img.columns()};
  fftwf_r2r_kind dct2[] = {FFTW_REDFT10, FFTW_REDFT10};
  fftwf_plan plan = fftwf_plan_many_r2r(2, n, img.channels(), img.data(), NULL,
                                        1, img.pixels(), img.data(), NULL, 1,
                                        img.pixels(), dct2, FFTW_ESTIMATE);
  fftwf_execute(plan);
  fftwf_destroy_plan(plan);

#ifdef ISOMETRIC_DCT
  ////> isometric normalization
  // this normalization (and scaling) affects several other functions: 
  //     idct_inplace, decompose, recompose
  // but they can all be removed only by applying the Normalization below 
  double norm_factor = sqrt(.25f / (img.rows() * img.columns()));
  for (int ch = 0; ch < img.channels(); ++ch) {
    for (int row = 0; row < img.rows(); ++row) {
      img.val(0, row, ch) /= sqrt(2.f);
      for (int col = 0; col < img.columns(); ++col) {
        img.val(col, row, ch) *= norm_factor;
      }
    }
    for (int col = 0; col < img.columns(); ++col) {
      img.val(col, 0, ch) /= sqrt(2.f);
    }
  }
#else
  ////> Normalization
  for (int i = 0; i < img.samples(); ++i) {
    img.val(i) /= 4 * img.pixels();
  }
#endif
}

/*! \brief 2D inverse DCT transform of image (channel-wise)
 *
 * Operates over the img in-place
 */
void idct_inplace(Image &img) {
#ifdef ISOMETRIC_DCT
  ////> isometric normalization
  long double norm_factor = sqrt(.25f / (img.rows() * img.columns()));
  for (int ch = 0; ch < img.channels(); ++ch) {
    for (int row = 0; row < img.rows(); ++row) {
      img.val(0, row, ch) *= sqrt(2.f);
      for (int col = 0; col < img.columns(); ++col) {
        img.val(col, row, ch) *= norm_factor;
      }
    }
    for (int col = 0; col < img.columns(); ++col) {
      img.val(col, 0, ch) *= sqrt(2.f);
    }
  }
#endif

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
#ifdef ISOMETRIC_DCT
    ////> isometric normalization scaling
    const double scaling = std::sqrt((double)(w*h)/((double)(img.rows()*img.columns())));
#else
    const double scaling = 1.0;
#endif
    for (int ch = 0; ch < freq.channels(); ++ch) {
      for (int r = 0; r < h; ++r) {
        for (int c = 0; c < w; ++c) {
          layer.val(c, r, ch) = freq.val(c, r, ch) * scaling;
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
#ifdef ISOMETRIC_DCT
    ////> isometric normalization scaling
    const double scaling = std::sqrt((double)(output.rows()*output.columns())/((double)(layer.rows()*layer.columns())));
#else
    const double scaling = 1.0;
#endif
    // Copy data (selected by recompose_factor)
    for (int ch = 0; ch < layer.channels(); ++ch) {
      for (int r = 0; r < layer.rows() * recompose_factor; ++r) {
        for (int c = 0; c < layer.columns() * recompose_factor; ++c) {
          output.val(c, r, ch) = layer.val(c, r, ch) * scaling;
        }
      }
    }
  }
  // IDCT of the output image
  idct_inplace(output);
  return output;
}

}  // namespace imgutils
