/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/#include <iostream>
#include <vector>

#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/ComponentFunction.hpp"
#include "ConditionsLib/SboxTools.hpp"
#include "ConditionsLib/enumeration_algorithm/MinLeaves.hpp"

std::vector<uint64_t> noekeon_sbox
{0x7, 0xa, 0x2, 0xc, 0x4, 0x8, 0xf, 0x0,
                                0x5, 0x9, 0x1, 0xe, 0x3, 0xd, 0xb, 0x6};

int main() {
  const auto dim = 4u;
  std::vector<uint64_t> &sbox = noekeon_sbox;
  for (auto x: sbox) {
    std::cout << x << " ";
  }
  std::cout << "\n";

  VectorialBooleanFunction fun(sbox, dim, dim);
  for (uint64_t component = 1; component < 16; ++component) {
    ComponentFunction cfun(fun, component);
    try {
      std::cout << "-------------------\n";
      std::cout << "Component: " << std::hex << component << std::dec << "\n";
      auto gs = MinLeaves::StartSearch(cfun);
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
    } catch(std::logic_error &e) {
      std::cout << "For component " << component << " nothing found.\n";
    }
  }

    return 0;
}