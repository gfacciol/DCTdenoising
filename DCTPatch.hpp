/*
 * DCTPatch.hpp
 *
 *  Created on: 11/feb/2015
 *      Author: nicola
 */

#ifndef DCTDENOISING_DFTPATCH_HPP_
#define DCTDENOISING_DFTPATCH_HPP_

#include <fftw3.h>
#include <cassert>
#include <complex>

namespace imgutils {

/*! \brief DCTPatch object
 *
 *  DCTPatch contains a patch and its DCT transform,
 *  and allows to convert back and forth between the two
 */
class DCTPatch {
 public:
  DCTPatch(int rows, int columns, int channels = 1);

  // disable copy constructor
  DCTPatch(const DCTPatch&) = delete;
  DCTPatch& operator=(const DCTPatch&) = delete;

  // default move constructor
  DCTPatch(DCTPatch&&) = default;
  DCTPatch& operator=(DCTPatch&&) = default;

  ~DCTPatch();

  void ToFreq();
  void ToSpace();
  int rows() const { return rows_; }
  int columns() const { return columns_; }
  int channels() const { return channels_; }
  float& space(int col, int row, int chan = 0);
  float& freq(int col, int row, int chan = 0);

 private:
  float *space_;
  float *freq_;
  fftwf_plan plan_forward_;
  fftwf_plan plan_backward_;
  const int rows_, columns_, channels_;
  const float norm_factor_;
};

 
/*! \brief Access the spatial representation of the patch
 *
 *  Returns pixel value
 */
inline float& DCTPatch::space(int col, int row, int chan) {
  assert(0 <= col && col < columns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return space_[chan * rows_ * columns_ + row * columns_ + col];
}

/*! \brief Access the DCT representation of the patch
 *
 *  Returns freq coefficient
 */
inline float& DCTPatch::freq(int col, int row, int chan) {
  assert(0 <= col && col < columns_);
  assert(0 <= row && row < rows_);
  assert(0 <= chan && chan < channels_);
  return freq_[chan * rows_ * columns_ + row * columns_ + col];
}

inline DCTPatch::DCTPatch(int rows, int columns, int channels)
    : rows_(rows), columns_(columns), channels_(channels),
      norm_factor_(std::sqrt(.25f / (rows * columns))) {
  int N = rows * columns * channels;
  space_ = reinterpret_cast<float *>(fftwf_malloc(sizeof(float) * N));
  freq_ = reinterpret_cast<float *>(fftwf_malloc(sizeof(float) * N));
  int n[] = {rows, columns};
  fftwf_r2r_kind dct2[] = {FFTW_REDFT10, FFTW_REDFT10};
  fftwf_r2r_kind idct2[] = {FFTW_REDFT01, FFTW_REDFT01};
#pragma omp critical
  {
    plan_forward_ = fftwf_plan_many_r2r(2, n, channels, space_, NULL, 1,
                                        rows * columns, freq_, NULL, 1,
                                        rows * columns, dct2,
                                        FFTW_MEASURE | FFTW_DESTROY_INPUT);
    plan_backward_ = fftwf_plan_many_r2r(2, n, channels, freq_, NULL, 1,
                                         rows * columns, space_, NULL, 1,
                                         rows * columns, idct2,
                                         FFTW_MEASURE | FFTW_DESTROY_INPUT);
  }
}

inline DCTPatch::~DCTPatch() {
  fftwf_free(space_);
  fftwf_free(freq_);
  fftwf_destroy_plan(plan_forward_);
  fftwf_destroy_plan(plan_backward_);
}


/*! \brief Computes the isometric DCT(Type2) of space, stores result in freq
 */
inline void DCTPatch::ToFreq() {
  fftwf_execute(plan_forward_);
  // normalize coefficients
  for (int ch = 0; ch < channels_; ++ch) {
    for (int row = 0; row < rows_; ++row) {
      freq(0, row, ch) /= sqrt(2.f);
      for (int col = 0; col < columns_; ++col) {
        freq(col, row, ch) *= norm_factor_;
      }
    }
    for (int col = 0; col < columns_; ++col) {
      freq(col, 0, ch) /= sqrt(2.f);
    }
  }
}

/*! \brief Computes the isometric iDCT(Type3) of freq, stores result in space
 */
inline void DCTPatch::ToSpace() {
  // normalize coefficients 
  for (int ch = 0; ch < channels_; ++ch) {
    for (int row = 0; row < rows_; ++row) {
      freq(0, row, ch) *= sqrt(2.f);
      for (int col = 0; col < columns_; ++col) {
        freq(col, row, ch) *= norm_factor_;
      }
    }
    for (int col = 0; col < columns_; ++col) {
      freq(col, 0, ch) *= sqrt(2.f);
    }
  }
  fftwf_execute(plan_backward_);
}

} /* namespace imgutils */

#endif  // DCTDENOISING_DFTPATCH_HPP_
