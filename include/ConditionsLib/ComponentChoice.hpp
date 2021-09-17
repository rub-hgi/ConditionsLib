/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_COMPONENTCHOICE_HPP_
#define CONDITIONS_LIB_COMPONENTCHOICE_HPP_

#include <bit>
#include <vector>
#include "VectorialBooleanFunction.hpp"

class ComponentChoice: public VectorialBooleanFunction {
 public:
  /// Construct function from multiple components.
  /**
   * Construct function from multiple components.
   * \param components List of components of \p v which will make up the coordinates
   * of the new function.
   * \param v Original function.
   */
  ComponentChoice(
      const std::vector<uint64_t> &components,
      const VectorialBooleanFunction &v) : VectorialBooleanFunction(v.GetValues(), v.InputSize(), v.OutputSize()) {
    std::vector<uint64_t> function(1u << v.InputSize());
    size_t i = 0;
    for (auto &x: components) {
      for (uint64_t inp = 0; inp < (1u << v.InputSize()); ++inp) {
        function[inp] ^= (std::popcount(x & v(inp)) & 1) << i;
      }
      i++;
    }
    this->GetValuesMutable() = function;
    this->OutputSizeMutable() = i;
  }
};


#endif //CONDITIONS_SRC_COMPONENTCHOICE_HPP_
