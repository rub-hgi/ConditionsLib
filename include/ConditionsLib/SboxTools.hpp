/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_SBOX_TOOLS_HPP_
#define CONDITIONS_LIB_SBOX_TOOLS_HPP_

#include <sstream>

#include "SboxTools.hpp"
#include "VectorialBooleanFunction.hpp"
#include "enumeration_algorithm/BruteForce.hpp"
#include "enumeration_algorithm/BruteForceOptimizations.hpp"
#include "enumeration_algorithm/MinLeaves.hpp"

/// Convenience functions for reading Sboxes from files and analyzing them.
namespace SboxTools {

std::vector<uint64_t> ReadFunction(std::istream &stream, size_t input_size,
                                   size_t output_size);

std::vector<std::vector<uint64_t>>
ReadFile(std::ifstream &stream, size_t input_size, size_t output_size);

std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
LinearStructures(const VectorialBooleanFunction &cfun, size_t input_size);

void AllSubspaces(BruteForce::BinaryDecisionTreePointer &x,
                  std::vector<Uint64Subspace> &list,
                  std::vector<uint64_t> &current_basis, size_t dimension);

Uint64Subspace Intersection(std::vector<Uint64Subspace> &sub);
int64_t SignedLinearity(const VectorialBooleanFunction &fun);

int64_t Linearity(const VectorialBooleanFunction &fun);

uint64_t Uniformity(const VectorialBooleanFunction &fun);

uint64_t NaiveDegree(const VectorialBooleanFunction &fun);

Uint64Subspace Dom(const BruteForce::BinaryDecisionTreePointer &x,
                   uint64_t input_size);

struct csv_line {
public:
  std::string function;
  uint64_t uniformity;
  uint64_t linearity;
  uint64_t card_ls0;
  uint64_t card_ls1;
  uint64_t degree;
  double costs_A;
  double costs_D;
  friend std::ostream &operator<<(std::ostream &, const csv_line &);
  static std::ostream &PrintHeader(std::ostream &os);
};
std::ostream &operator<<(std::ostream &os, const csv_line &line);

template <class CSVLine>
CSVLine Analyse(const VectorialBooleanFunction &fun, bool print = false) {
  CSVLine result;
  result.function += "\"";
  for (auto x : fun.GetValues()) {
    result.function += std::to_string(x);
  }
  result.function += "\"";
  result.uniformity = SboxTools::Uniformity(fun);
  result.linearity = SboxTools::Linearity(fun);
  try {
    auto gs = MinLeaves::StartSearch(fun, print);
    /*auto test_gs = BruteForceOptimization::StartSearch(fun, print);
    if(test_gs->leaves() != gs->leaves()) {
      std::cout << *gs << std::endl;
      std::cout << "avgcost: " << std::dec<< gs->AveragePathLength() <<
    std::endl; std::cout << "leaves: " <<std::dec << gs->leaves() << std::endl;
      std::cout << *test_gs << std::endl;
      std::cout << "Test avgcost: " << std::dec << test_gs->AveragePathLength()
    << std::endl; std::cout << "Test leaves: " << std::dec << test_gs->leaves()
    << std::endl;
    }*/
    result.costs_A = gs->associated_cost_;
    std::vector<Uint64Subspace> alls;
    std::vector<uint64_t> current_basis;
    AllSubspaces(gs, alls, current_basis, fun.InputSize());
    auto in = Intersection(alls);
    result.costs_D = // ((double)fun.InputSize() - in.Dimension()) -
        Dom(gs, fun.InputSize()).Dimension();
  } catch (std::logic_error &e) {
    std::cerr << "Nothing found.\n";
  }
  result.degree = NaiveDegree(fun);
  auto ls = LinearStructures(fun, fun.InputSize());
  result.card_ls0 = ls.first.size();
  result.card_ls1 = ls.second.size();
  return result;
}

std::vector<uint64_t> ReadFunctionNodelim(std::istream &stream,
                                          size_t input_size,
                                          size_t output_size);
std::vector<VectorialBooleanFunction>
ReadFileNodelim(std::ifstream &stream, size_t input_size, size_t output_size,
                ssize_t frst_line = 0, ssize_t last_line = -1);
} // namespace SboxTools

#endif
