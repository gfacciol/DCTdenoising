/*
 * Image.hpp
 *
 *  Created on: 14/gen/2015
 *      Author: nicola
 */

#ifndef IMAGE_HPP_
#define IMAGE_HPP_

#include <cassert>
#include <vector>
#include <utility>

namespace imgutils {

class Image {
 public:
  Image() = default;
  Image(int rows, int columns, int channels = 1, float val = 0.f);
  // construct from C array
  Image(const float *data, int rows, int columns, int channels = 1);

  // disable copy constructor
  Image(const Image&) = delete;
  Image& operator=(const Image&) = delete;
  // instead of the copy constructor, we want explicit copy
  Image copy() const;

  // default move constructor
  Image(Image&&) = default;
  Image& operator=(Image&&) = default;

  ~Image() = default;

  void Clear(float val = 0.f) { std::fill(data_.begin(), data_.end(), val); }

  const float& val(int col, int row, int chan = 0) const;
  float& val(int col, int row, int chan = 0);
  const float& val(int pos) const;
  float& val(int pos);

  int channels() const { return channels_; }
  int columns() const { return columns_; }
  int rows() const { return rows_; }
  int pixels() const { return columns_ * rows_; }
  int samples() const { return channels_ * columns_ * rows_; }
  float* data() { return data_.data(); }
  const float* data() const { return data_.data(); }
  std::pair<int, int> shape() const { return {rows_, columns_}; }
  std::vector<float>::iterator begin() { return data_.begin(); }
  std::vector<float>::const_iterator begin() const { return data_.begin(); }
  std::vector<float>::iterator end() { return data_.end(); }
  std::vector<float>::const_iterator end() const { return data_.end(); }

 protected:
  int rows_{0};
  int columns_{0};
  int channels_{0};
  std::vector<float> data_{};
};

inline Image::Image(int rows, int columns, int channels, float val)
    : rows_(rows), columns_(columns), channels_(channels),
      data_(rows * columns * channels, val) {}

inline Image::Image(const float *data, int rows, int columns, int channels)
    : rows_(rows), columns_(columns), channels_(channels),
      data_(data, data + rows * columns * channels) {}

inline Image Image::copy() const {
  return Image(data(), rows(), columns(), channels());
}

inline const float& Image::val(int col, int row, int chan) const {
  assert(0 <= col && col < columns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return data_[chan * rows_ * columns_ + row * columns_ + col];
}

inline float& Image::val(int col, int row, int chan) {
  assert(0 <= col && col < columns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return data_[chan * rows_ * columns_ + row * columns_ + col];
}

inline const float& Image::val(int pos) const {
  return data_[pos];
}

inline float& Image::val(int pos) {
  return data_[pos];
}

} /* namespace imgutils */

#endif  // IMAGE_HPP_
