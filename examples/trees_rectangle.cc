/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <fstream>
#include <vector>

#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/ComponentFunction.hpp"
#include "ConditionsLib/SboxTools.hpp"

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
void analyse(std::ostream& os, std::vector<uint64_t> & sbox, int dim,
             uint64_t root) {
  for (auto x: sbox) {
    os << x << " ";
  }
  os << "\n";
  VectorialBooleanFunction fun(sbox, dim, dim);
  for (uint64_t component = 1; component < 16; ++component) {
    ComponentFunction cfun(fun, component);
    try {
      os << "-------------------\n";
      os << "Component: " << std::hex << component << std::dec << "\n";
      double bound = cfun.InputSize();
      auto gs = BruteForce::StartSearchWithFixedRoot(cfun, bound, root);
      os << "Result: \n" << *gs << std::endl;
      os << "Costs (Average guessing): " << gs->associated_cost_ << std::endl;
    } catch(std::logic_error &e) {
      os << "For component " << component << " nothing found.\n";
    }
  }
}
void analyse_overall(std::ostream& os, std::vector<uint64_t> & sbox, int dim,
             uint64_t root) {
  for (auto x: sbox) {
    os << x << " ";
  }
  os << "\n";
  VectorialBooleanFunction fun(sbox, dim, dim);
  for (uint64_t component = 1; component < 16; ++component) {
    ComponentFunction cfun(fun, component);
    try {
      os << "-------------------\n";
      os << "Component: " << std::hex << component << std::dec << "\n";
      double bound = cfun.InputSize();
      auto gs = BruteForce::StartSearch(cfun, bound, false);
      os << "Result: \n" << *gs << std::endl;
      os << "Costs (Average guessing): " << gs->associated_cost_ << std::endl;
    } catch(std::logic_error &e) {
      os << "For component " << component << " nothing found.\n";
    }
  }

}

void Rectangle() {

  const auto dim = 4;
  std::vector<uint64_t> possible_roots {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  for (auto root: possible_roots) {
    {
      std::ofstream ffile("Root" + std::to_string(root) + "_Rectangle.txt");
      analyse(ffile, rectangle_sbox, dim, root);
    }
    {
      std::ofstream ffile("Root" + std::to_string(root) + "_RectangleInverse.txt");
      analyse(ffile, rectangle_inverse, dim, root);
    }
  }
    {
      std::ofstream ffile("Rectangle.txt");
      analyse_overall(ffile, rectangle_sbox, dim, 0);
    }
    {
      std::ofstream ffile("RectangleInverse.txt");
      analyse_overall(ffile, rectangle_inverse, dim, 0);
    }
}

void NewRectangle() {

  const auto dim = 4;
  std::vector<uint64_t> possible_roots {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  for (auto root: possible_roots) {
    {
      std::ofstream ffile("Root" + std::to_string(root) + "_NewRectangle.txt");
      analyse(ffile, new_rectangle_sbox, dim, root);
    }
    {
      std::ofstream ffile("Root" + std::to_string(root) + "_NewRectangleInverse.txt");
      analyse(ffile, new_rectangle_inverse, dim, root);
    }
  }
  {
    std::ofstream ffile("NewRectangle.txt");
    analyse_overall(ffile, new_rectangle_sbox, dim, 0);
  }
  {
    std::ofstream ffile("NewRectangleInverse.txt");
    analyse_overall(ffile, new_rectangle_inverse, dim, 0);
  }
}

int main() {
  init();
  std::vector<uint64_t> possible_roots {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  Rectangle();
  NewRectangle();
  return 0;
}