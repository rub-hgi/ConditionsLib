/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_VECTORIAL_BOOLEAN_FUNCTION_HPP_
#define CONDITIONS_LIB_VECTORIAL_BOOLEAN_FUNCTION_HPP_
#include <cstdint>
#include <ostream>
#include <vector>

/// Class for storing a vectorial Boolean function
class VectorialBooleanFunction {
private:
  size_t input_size = 0, output_size = 0; /** Input and output dimension of the function. */
  std::vector<uint64_t> values; /**< Lookup table. */
  bool ddt_calculated = false; /**< True if ddt does not have to be re-calculated. */

protected:
  std::vector<uint64_t> &GetValuesMutable(); /**< Get a mutable reference to the lookup-table */
  void TruncateOutputSize(size_t new_output_size); /**< Change the number of output bits */

public:
  std::vector<std::vector<uint64_t>> ddt; /**< Field for saving the ddt */
  /**
   * Check if the look-up table, input dimension and output dimension for
   * this function and \p other are identical.
   */
  bool operator==(const VectorialBooleanFunction &other) const;
  friend std::ostream &operator<<(std::ostream &os,
                                  const VectorialBooleanFunction &vec);
  VectorialBooleanFunction(const std::initializer_list<uint64_t> &x,
                           size_t input_size, size_t output_size);

  template <class V>
  VectorialBooleanFunction(const V &x, size_t input_size, size_t output_size);
  /// Get read-only reference to look-up table.
  const std::vector<uint64_t> &GetValues() const;
  /// Dimension of input space.
  [[nodiscard]] size_t InputSize() const;
  /// Dimension of output space.
  [[nodiscard]] size_t OutputSize() const;
  /// Evaluate function at position \p x .
  uint64_t operator()(uint64_t x) const;

  void calculate_ddt() {
    if (ddt_calculated)
      return;
    ddt.resize(1u << input_size);
    for (auto &x : ddt)
      x.resize(1u << output_size);
    size_t beta;
    for (size_t i = 0; i < (1u << input_size); ++i) {
      for (size_t x = 0; x < (1u << input_size); ++x) {
        ++ddt[i][values[x] ^ values[x ^ i]];
      }
    }
    ddt_calculated = true;
  }

  size_t differential_uniformity() {
    calculate_ddt();
    size_t uni = 0;
    for (size_t i = 1; i < (1u << input_size); ++i) {
      for (size_t j = 0; j < (1u << output_size); ++j) {
        if (ddt[i][j] > uni)
          uni = ddt[i][j];
      }
    }
    return uni;
  }

  bool IsAPN() {
    calculate_ddt();
    bool apn = true;
    for (size_t i = 1; i < (1u << input_size) && apn; ++i) {
      for (size_t j = 0; j < (1u << output_size) && apn; ++j) {
        apn &= ddt[i][j] == 0 || ddt[i][j] == 2;
      }
    }
    return apn;
  }
  size_t &OutputSizeMutable();
};

template <class V>
VectorialBooleanFunction::VectorialBooleanFunction(const V &x,
                                                   size_t input_size,
                                                   size_t output_size)
    : values(std::move(x)), input_size(input_size), output_size(output_size) {}

/// Generates random boolean function with one output bit and fixed weight.
template <class RNG>
VectorialBooleanFunction RandomVBFWithWeight(uint64_t dim, uint64_t weight,
                                             RNG &rng) {
  std::vector<uint64_t> base;
  base.resize(1u << dim);
  for (size_t i = 0; i < weight; ++i) {
    base[i] = 1;
  }
  for (size_t i = weight; i < 1u << dim; ++i) {
    base[i] = 0;
  }
  std::shuffle(base.begin(), base.end(), rng);
  VectorialBooleanFunction fun(base, dim, 1);
  return fun;
}

/// Generates random boolean function with one output bit and weight 2^(dim -
/// 1).
template <class RNG>
VectorialBooleanFunction RandomBalancedVBF(uint64_t dim, RNG &rng) {
  return RandomVBFWithWeight(dim, 1u << (dim - 1u), rng);
}

/// Generates random boolean function with one output bit.
template <class RNG>
VectorialBooleanFunction RandomVBF(uint64_t dim, RNG &rng) {
  std::vector<uint64_t> base;
  base.resize(1u << dim);
  for (size_t i = 0; i < 1u << dim; ++i) {
    base[i] = rng() % 2;
  }
  VectorialBooleanFunction fun(base, dim, 1);
  return fun;
}
#endif
