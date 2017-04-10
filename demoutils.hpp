//
// Created by Nicola Pierazzo on 21/10/15.
//

#ifndef IMGUTILS_DEMOUTILS_HPP
#define IMGUTILS_DEMOUTILS_HPP

#include "Image.hpp"
#include <string>

namespace imgutils {

const char *pick_option(int *c, char **v, const char *o, const char *d);
Image read_image(const std::string& filename);
void save_image(const Image& image, const std::string& filename);

}  // namespace imgutils

#endif //IMGUTILS_DEMOUTILS_HPP
