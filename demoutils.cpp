//
// Created by Gabriele Facciolo on 10/04/17.
//

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <utility>
#include "demoutils.hpp"

extern "C" {
#include "iio.h"
}

using std::string;
using std::free;
using std::strcmp;

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


}  // namespace imgutils
