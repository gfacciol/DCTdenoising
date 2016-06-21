//
// Created by Nicola Pierazzo on 21/10/15.
//

#ifndef IMGUTILS_UTILS_HPP
#define IMGUTILS_UTILS_HPP

#include "Image.hpp"
#include <string>

namespace imgutils {

const char *pick_option(int *c, char **v, const char *o, const char *d);
Image read_image(const std::string& filename);
void save_image(const Image& image, const std::string& filename);

std::pair<int, int> ComputeTiling(int rows, int columns, int tiles);
std::vector<Image> SplitTiles(const Image &src, int pad_before, int pad_after,
                              std::pair<int, int> tiling);
Image MergeTiles(const std::vector<std::pair<Image, Image>> &src,
                       std::pair<int, int> shape, int pad_before, int pad_after,
                       std::pair<int, int> tiling);

}  // namespace imgutils

#endif //IMGUTILS_UTILS_HPP
