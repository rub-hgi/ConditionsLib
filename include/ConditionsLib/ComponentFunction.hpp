/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_COMPONENT_FUNCTION_HPP_
#define CONDITIONS_LIB_COMPONENT_FUNCTION_HPP_
#include "VectorialBooleanFunction.hpp"
#include <bit>

class ComponentFunction : public VectorialBooleanFunction {
 public:
  /// Construct a function representing \p component * \p underlying_function.
  /**
   * Construct a function representing \p component * \p underlying_function.
   * Here * denotes the the inner product
   * std::popcount(component & underlying_function).
   */
  ComponentFunction(const VectorialBooleanFunction &underlying_function,
                    uint64_t component) :
      VectorialBooleanFunction(underlying_function) {
    for (auto &x: GetValuesMutable()) {
      x = std::popcount(x & component) % 2;
    }
    TruncateOutputSize(1);
  }
};

#endif
