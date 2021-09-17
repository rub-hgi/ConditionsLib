/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <iostream>
#include <vector>
#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/ComponentFunction.hpp"
#include "ConditionsLib/SboxTools.hpp"
#include "ConditionsLib/enumeration_algorithm/MinLeaves.hpp"

std::vector<uint64_t> present_sbox{0xc, 0x5, 0x6, 0xb,
                                   0x9, 0x0, 0xa, 0xd,
                                   0x3, 0xe, 0xf, 0x8,
                                   0x4, 0x7, 0x1, 0x2};

std::vector<uint64_t> present_inverse(16);

void init() {
  for (uint64_t i = 0; i < 16; ++i)
    present_inverse[present_sbox[i]] = i;
}

std::vector<uint64_t> prepare_f_i(uint64_t delta) {
  std::vector<uint64_t> result(16);
  for (uint64_t i = 0; i < 16; ++i) {
    result[i] = present_sbox[present_inverse[i] ^ delta] ^ i;
  }
  return result;
}

int main() {
  init();
  VectorialBooleanFunction fun(present_sbox, 4, 4);
  for (uint64_t component = 1; component < 16; ++component) {
    auto c = component;
    ComponentFunction cfun(fun, c);
    try {
      std::cout << "-------------------\n";
      std::cout << "Component: " << std::hex << c << std::dec << "\n";
      auto gs = MinLeaves::StartSearch(cfun, true);
      auto ls = SboxTools::LinearStructures(cfun, 4);
      std::cout << "Result: \n" << *gs << std::endl;
      std::cout << "0-Linear Structures" << std::endl;
      for (auto x: ls.first) {
        std::cout << std::dec << x << std::endl;
      }
      std::cout << "1-Linear Structures" << std::endl;
      for (auto x: ls.second) {
        std::cout << std::dec << x << std::endl;
      }
    } catch (std::logic_error &e) {
      std::cout << "For component " << component << " nothing found.\n";
    }
  }

  return 0;
}