/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/#include <fstream>
#include <vector>

#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/ComponentChoice.hpp"
#include "ConditionsLib/SboxTools.hpp"

std::vector<uint64_t> noekeon_sbox
{0x7, 0xa, 0x2, 0xc, 0x4, 0x8, 0xf, 0x0,
                                0x5, 0x9, 0x1, 0xe, 0x3, 0xd, 0xb, 0x6};

void analyse_pairs(std::ostream& out, VectorialBooleanFunction & v, std::vector<uint64_t> outputs) {
  auto cfun = ComponentChoice(outputs, v);
  for (size_t i = 1; i < (1u << v.InputSize()); ++i) {
    for (size_t j = 1+i; j < (1u << v.InputSize()); ++j) {
      out << "Components (dec): " << std::dec;
      for (auto x: outputs) out << x << ",";
      out << "\n";
      out << "Fixed bits: "  << i << "," << j << "\n";
      auto stump = BruteForceOptimization::fix_bits(std::vector<uint64_t> {i, j});
      double bound = 1u << v.InputSize();
      auto gs = MinLeaves::StartSearchWithFixedStump(cfun, stump, bound, false);
      out << *gs << "\n------------------" << std::endl;
    }
  }
}

int main() {
  VectorialBooleanFunction fun(noekeon_sbox, 4, 4);
  std::vector<std::vector<uint64_t>> list_of_outputs {};
  for (uint64_t x1 = 1; x1 < 16; ++x1) {
    for (uint64_t x2 = x1 + 1; x2 < 16; ++x2) {
      list_of_outputs.push_back(std::vector<uint64_t> {x1, x2});
    }
  }
  for (auto &outputs: list_of_outputs) {
    std::stringstream ss;
    for (auto &x: outputs) {
      ss << std::dec << x << ".";
    }
    std::ofstream out("noekeon_fixed." + ss.str() + "txt");
    analyse_pairs(out, fun, outputs);
  }
  return 0;
}

