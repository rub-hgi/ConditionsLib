/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <iostream>
#include <vector>

#include "ConditionsLib/enumeration_algorithm/MinLeaves.hpp"

std::vector<uint64_t> rectangle_sbox{
  0x6, 0x5, 0xC, 0xA, 0x1, 0xE, 0x7, 0x9, 0xB, 0x0, 0x3, 0xD, 0x8, 0xF, 0x4, 0x2
};

VectorialBooleanFunction prepare_g(const VectorialBooleanFunction & f,
                                   uint64_t delta, uint64_t Delta) {
  const auto n = f.InputSize();
  std::vector<uint64_t> underlying_function(1 << n);
  for (uint64_t x = 0; x < 1 << n; ++x) {
      if ((f(x) ^ f(x^delta)) != Delta)
        underlying_function[x] = 1;
  }
  return { underlying_function, n, 1 };
}

int main() {
  uint64_t delta = 6;
  uint64_t Delta = 2;
  VectorialBooleanFunction vbf(rectangle_sbox, 4, 4);
  auto g = prepare_g(vbf, delta, Delta);
  try {
    auto tree_ptr = MinLeaves::StartSearch(g, false);
    std::cout << *tree_ptr << std::endl;
  } catch(std::logic_error &e) {

  }
  return 0;
}