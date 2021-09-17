/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <fstream>
#include <vector>

#include "ConditionsLib/ComponentChoice.hpp"
#include "ConditionsLib/SboxTools.hpp"
#include "ConditionsLib/VectorialBooleanFunction.hpp"

std::vector<uint64_t> rectangle_sbox{0x9, 0x4, 0xF, 0xA, 0xE, 0x1, 0x0, 0x6,
                                     0xC, 0x7, 0x3, 0x8, 0x2, 0xB, 0x5, 0xD};

std::vector<uint64_t> rectangle_inverse(16);

std::vector<uint64_t> rectangle_new_sbox{0x6, 0x5, 0xC, 0xA, 0x1, 0xE,
                                         0x7, 0x9, 0xB, 0x0, 0x3, 0xD,
                                         0x8, 0xF, 0x4, 0x2};

std::vector<uint64_t> rectangle_new_inverse(16);

void init() {
  for (uint64_t i = 0; i < 16; ++i) {
    rectangle_inverse[rectangle_sbox[i]] = i;
    rectangle_new_inverse[rectangle_new_sbox[i]] = i;
  }
}

void analyse_pairs(std::ostream &out, VectorialBooleanFunction &v,
                   std::vector<uint64_t> outputs) {
  auto cfun = ComponentChoice(outputs, v);
  for (size_t i = 1; i < (1u << v.InputSize()); ++i) {
    for (size_t j = 1 + i; j < (1u << v.InputSize()); ++j) {
      out << "Components (dec): " << std::dec;
      for (auto x : outputs)
        out << x << ",";
      out << "\n";
      out << "Fixed bits: " << i << "," << j << "\n";
      auto stump =
          BruteForceOptimization::fix_bits(std::vector<uint64_t>{i, j});
      double bound = v.InputSize();
      auto gs = BruteForceOptimization::StartSearchWithFixedStump(
          cfun, stump, bound, 0, false);
      out << *gs << "\n------------------" << std::endl;
    }
  }
}

int main() {
  init();
  VectorialBooleanFunction fun(rectangle_new_sbox, 4, 4);
  VectorialBooleanFunction fun_inverse(rectangle_new_inverse, 4, 4);
  std::vector<std::vector<uint64_t>> list_of_outputs{{2}, {4}, {2, 4}, {8}};
  for (auto &outputs : list_of_outputs) {
    std::stringstream ss;
    for (auto &x : outputs) {
      ss << std::dec << x << ".";
    }
    std::ofstream out("rect_fixed_straight." + ss.str() + "txt");
    analyse_pairs(out, fun, outputs);
  }
  for (auto &outputs : list_of_outputs) {
    std::stringstream ss;
    for (auto &x : outputs) {
      ss << std::dec << x << ".";
    }
    std::ofstream out("rect_fixed_inverse." + ss.str() + "txt");
    analyse_pairs(out, fun_inverse, outputs);
  }
  return 0;
}
