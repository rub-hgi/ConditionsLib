/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_ENUMERATION_ALGORITHM_EXHAUSTIVEBF_HPP_
#define CONDITIONS_LIB_ENUMERATION_ALGORITHM_EXHAUSTIVEBF_HPP_

#include "BruteForce.hpp"

/// Find all trees for a function (which have no redundancies on any path).
namespace ExhaustiveBF {
using BinaryDecisionTree = BruteForce::BinaryDecisionTree;
using BinaryDecisionTreePointer = BruteForce::BinaryDecisionTreePointer;
typedef std::vector<BinaryDecisionTreePointer> BinaryDecisionTreeVector;

BinaryDecisionTreeVector TreeSearch(const VectorialBooleanFunction &fun,
                                    const NTL::mat_GF2 &vectors_on_path,
                                    const NTL::vec_GF2 &choices,
                                    int last_choice,
                                    const NTL::vec_GF2 &coset_so_far,
                                    double &bound,
                                    int level, bool print = true);

BinaryDecisionTreeVector StartSearch(const VectorialBooleanFunction &fun, bool print=true);
BinaryDecisionTreeVector StartSearch(const VectorialBooleanFunction &fun,
                                     double &bound, bool print=true);

};

#endif
