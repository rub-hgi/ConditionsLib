/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <fstream>
#include <iostream>
#include <vector>
#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/SboxTools.hpp"
#include "ConditionsLib/Uint64Subspace.hpp"
#include "ConditionsLib/ComponentFunction.hpp"

struct csv_line {
 public:
  std::string function;
  uint64_t sbox_number;
  uint64_t component;
  uint64_t uniformity;
  uint64_t linearity;
  uint64_t card_ls0;
  uint64_t card_ls1;
  uint64_t degree;
  double costs_A;
  double costs_D;
  friend std::ostream& operator<<(std::ostream&, const csv_line&);
  static std::ostream& printHeader(std::ostream& os) {
    os << "function" << ",";
    os << "sbox_number" << ",";
    os << "component" << ",";
    os << "uniformity" << ",";
    os << "linearity" << ",";
    os << "card_ls0" << ",";
    os << "card_ls1" << ",";
    os << "degree" << ",";
    os << "leaves" << ",";
    os << "domsize" << ",";
    return os;
  }
};

std::ostream &operator<<(std::ostream &os, const csv_line& line) {
  os << line.function << ",";
  os << line.sbox_number<< ",";
  os << line.component << ",";
  os << line.uniformity << ",";
  os << line.linearity << ",";
  os << line.card_ls0 << ",";
  os << line.card_ls1 << ",";
  os << line.degree << ",";
  os << line.costs_A << ",";
  os << line.costs_D << ",";
  return os;
}

auto Analyse(const VectorialBooleanFunction& fun) {
  csv_line result;
  result.function += "\"";
  for (auto x: fun.GetValues()) {
    result.function += std::to_string(x);
  }
  result.function += "\"";
  result.uniformity = SboxTools::Uniformity(fun);
  result.linearity = SboxTools::Linearity(fun);
  try {
    auto gs = MinLeaves::StartSearch(fun, false);
    result.costs_A = gs->associated_cost_;
    std::vector<Uint64Subspace> alls;
    std::vector<uint64_t> current_basis;
    SboxTools::AllSubspaces(gs, alls, current_basis, fun.InputSize());
    auto in = SboxTools::Intersection(alls);
    result.costs_D = // ((double)fun.InputSize() - in.Dimension()) -
        SboxTools::Dom(gs, fun.InputSize()).Dimension();
  } catch(std::logic_error &e) {
    std::cerr << "Nothing found.\n";
  }
  result.degree = SboxTools::NaiveDegree(fun);
  auto ls = SboxTools::LinearStructures(fun, fun.InputSize());
  result.card_ls0 = ls.first.size();
  result.card_ls1 = ls.second.size();
  return result;
}

int main() {
  const std::string filename = "../FunctionData/sbox_classes.dat";
  const uint64_t dim = 4;
  std::ifstream file(filename);
  auto functions_raw = SboxTools::ReadFile(file, dim, dim);
  csv_line::printHeader(std::cout) << std::endl;
  auto sbox_number = 0;
  for (auto &x_raw: functions_raw) {
    VectorialBooleanFunction x(x_raw, dim, dim);
    for (uint64_t comp = 1; comp < (1u << dim); ++comp) {
      auto y = ComponentFunction(x, comp);
      auto line = Analyse(y);
      line.sbox_number = sbox_number;
      line.component = comp;
      std::cout << line << std::endl;
    }
    sbox_number++;
  }
}
