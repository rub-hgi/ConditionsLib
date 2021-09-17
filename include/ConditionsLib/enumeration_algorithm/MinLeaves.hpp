/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_ENUMERATION_ALGORITHM_MIN_LEAVES_HPP_
#define CONDITIONS_LIB_ENUMERATION_ALGORITHM_MIN_LEAVES_HPP_

#include <iostream>
#include <memory>
#include <optional>

#include "../Uint64Subspace.hpp"
#include "../VectorialBooleanFunction.hpp"
#include "BruteForce.hpp"

/// Implements the optimized algorithm for finding trees with minimal size.
namespace MinLeaves {
typedef BruteForce::BinaryDecisionTreePointer BinaryDecisionTreePointer;
typedef BruteForce::BinaryDecisionTree BinaryDecisionTree;

BinaryDecisionTreePointer
TreeSearch(const VectorialBooleanFunction &fun,
           const NTL::mat_GF2 &vectors_on_path, const NTL::vec_GF2 &choices,
           int last_choice, const NTL::vec_GF2 &coset_so_far, double &bound,
           int level, bool print = true);

/// Starts search for size-minimal tree.
/** Finds optimal (lowest average path length) tree representing a given
 * function.
 * \param fun Target function.
 * \param print If set to true, outputs the trees improving
 * the internal size bound.
 */
BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun,
                                      bool print = true);

/// Starts search for best tree whose size is below a certain
/// bound.
/** Finds tree better than a certain bound and saves its size.
 * \param fun Target function.
 * \param bound Pointer to upper bound on
 * the number of leaves. Is updated to contain the number of leaves of the best
 * tree found.
 * \param print If set to true, outputs the trees improving the
 * internal size bound.
 */
BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun,
                                      double &bound, bool print = true);

/// Starts search for best tree whose size is below a certain
/// bound with a fixed root.
/**
 * Finds tree better than a certain bound and saves the number of its leaves.
 * The root label will be fixed though. This is useful for parallelizing the
 * search.
 * \param fun Target function.
 * \param bound Pointer to upper bound on
 * the number of leaves. Is updated to contain the number of leaves of the best
 * tree found.
 * \param root Label of the trees root.
 * \param print If set to true,
 * outputs the trees improving the internal size bound.
 */
BinaryDecisionTreePointer
StartSearchWithFixedRoot(const VectorialBooleanFunction &fun, double &bound,
                         uint64_t root, bool print = true);

/// Start search with search with a fixed stump.
BinaryDecisionTreePointer
StartSearchWithFixedStump(const VectorialBooleanFunction &fun,
                          const BinaryDecisionTreePointer &stump, double &bound,
                          bool print = true);
void RecursiveStump(const VectorialBooleanFunction &fun,
                    const BinaryDecisionTreePointer &stump,
                    const NTL::mat_GF2 &vectors_on_path,
                    const NTL::vec_GF2 &choices, int last_choice,
                    const NTL::vec_GF2 &coset_so_far, double &bound, int level,
                    bool print);

} // namespace MinLeaves

#endif
