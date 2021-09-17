/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <vector>
#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/SboxTools.hpp"
#include "ConditionsLib/ComponentChoice.hpp"

void Analyse(std::ostream& os, const VectorialBooleanFunction& v, int dim) {
  os << "Function: " << std::endl;
  for (auto x: v.GetValues()) {
    os << x << " ";
  }
  os << "\n";
  try {
    os << "-------------------\n";
    auto ls = SboxTools::LinearStructures(v, v.InputSize());
    os << "0-Linear Structures" << std::endl;
    for (auto x: ls.first) {
      os << std::dec << x << std::endl;
    }
    auto g = BruteForce::StartSearch(v, false);
      os << *g << "\n";
      os << "Costs (Average guessing): " << g->associated_cost_ << std::endl;
  } catch(std::logic_error &e) {
    os << "Nothing found for whole input.\n";
  }
  os << "###################\n";
}

void choice(VectorialBooleanFunction & v) {
  std::stringstream os2;
  std::stringstream os3;
  for (uint64_t x = 1; x < 1u<<v.OutputSize(); x*=2) {
    for (uint64_t y = x * 2; y < 1u << v.OutputSize(); y *= 2) {
      ComponentChoice cfun(std::vector<uint64_t>{x, y}, v);
      os2 << "{" << x << "," << y << "}\n";
      Analyse(os2, cfun, cfun.InputSize());
      for (uint64_t z = y * 2; z < 1u << v.OutputSize(); z *= 2) {
        ComponentChoice cfun(std::vector<uint64_t>{x, y, z}, v);
        os3 << "{" << x << "," << y << "," << z << "}\n";
        Analyse(os3, cfun, cfun.InputSize());
      }
    }
  }
  std::cout << os2.str() << std::endl;
  std::cout << os3.str() << std::endl;
};
std::vector<uint64_t> rectangle_sbox{
  0x9,0x4,0xF,0xA,0xE,0x1,0x0,0x6,0xC,0x7,0x3,0x8,0x2,0xB,0x5,0xD
};

std::vector<uint64_t> rectangle_inverse(16);

std::vector<uint64_t> new_rectangle_sbox{
  0x6, 0x5, 0xC, 0xA, 0x1, 0xE, 0x7, 0x9, 0xB, 0x0, 0x3, 0xD, 0x8, 0xF, 0x4, 0x2
};

std::vector<uint64_t> new_rectangle_inverse(16);

void init() {
  for (uint64_t i = 0; i < 16; ++i) {
    rectangle_inverse[rectangle_sbox[i]] = i;
    new_rectangle_inverse[new_rectangle_sbox[i]] = i;
  }
}

int main() {
  init();
  VectorialBooleanFunction sbox(rectangle_inverse, 4, 4);
  choice(sbox);
  return 0;
}