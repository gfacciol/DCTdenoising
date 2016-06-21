//
// Created by Nicola Pierazzo on 21/10/15.
//

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <utility>
#include "utils.hpp"

extern "C" {
#include "iio.h"
}

using std::string;
using std::free;
using std::strcmp;
using std::pair;
using std::sqrt;
using std::vector;
using std::move;
using std::max;
using std::min;

namespace imgutils {

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

Image read_image(const string &filename) {
  int w, h, c;
  float *data = iio_read_image_float_split(filename.c_str(), &w, &h, &c);
  Image im(data, h, w, c);
  free(data);
  return im;
}

void save_image(const Image &image, const string &filename) {
  iio_save_image_float_split(const_cast<char *>(filename.c_str()),
                             const_cast<float *>(image.data()),
                             image.columns(), image.rows(), image.channels());
}

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

}  // namespace imgutils
