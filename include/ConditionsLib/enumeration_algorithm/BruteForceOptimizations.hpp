/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_ENUMERATION_ALGORITHM_BRUTEFORCE_OPTIMIZATION_HPP_
#define CONDITIONS_LIB_ENUMERATION_ALGORITHM_BRUTEFORCE_OPTIMIZATION_HPP_

#include <iostream>
#include <memory>
#include <optional>
#include <set>

#include "BruteForce.hpp"
#include "../Uint64Subspace.hpp"
#include "../VectorialBooleanFunction.hpp"

/// Same as BruteForce, only reducing the search by considering linear
/// structures.
namespace BruteForceOptimization {
typedef BruteForce::BinaryDecisionTreePointer BinaryDecisionTreePointer;
typedef BruteForce::BinaryDecisionTree BinaryDecisionTree;

BinaryDecisionTreePointer TreeSearch(const VectorialBooleanFunction &fun,
                                     const NTL::mat_GF2 &vectors_on_path,
                                     const NTL::vec_GF2 &choices,
                                     int last_choice,
                                     const NTL::vec_GF2 &coset_so_far,
                                     double &bound,
                                     int level, bool print= true);

BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun, bool print=true);
BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun,
                                      double &bound, bool print=true);
BinaryDecisionTreePointer StartSearchWithFixedRoot(const VectorialBooleanFunction &fun,
                                                   double &bound,
                                                   uint64_t root, bool print=true);
BinaryDecisionTreePointer StartSearchWithFixedStump(const VectorialBooleanFunction &fun,
                                                   const BinaryDecisionTreePointer &stump,
                                                   double &bound,
                                                   uint64_t root, bool print=true);
void RecursiveStump(const VectorialBooleanFunction &fun,
                                         const BinaryDecisionTreePointer &stump,
                                         const NTL::mat_GF2 &vectors_on_path,
                                         const NTL::vec_GF2 &choices,
                                         int last_choice,
                                         const NTL::vec_GF2 &coset_so_far,
                                         double &bound,
                                         int level,
                                         bool print);

  BinaryDecisionTreePointer generate_stump(std::vector<uint64_t>::iterator fixed_bits_begin,
                                           std::vector<uint64_t>::iterator fixed_bits_end);

  BinaryDecisionTreePointer fix_bits(std::vector<uint64_t> fixed_bits);
std::set<uint64_t> zero_linear_structures(const VectorialBooleanFunction &fun,
                                          const Uint64Subspace &space,
                                          uint64_t offset);
}

#endif
