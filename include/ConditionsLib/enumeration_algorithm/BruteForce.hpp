/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_ENUMERATION_ALGORITHM_BRUTEFORCE_HPP_
#define CONDITIONS_LIB_ENUMERATION_ALGORITHM_BRUTEFORCE_HPP_

#include <iostream>
#include <memory>
#include <optional>

#include "../Uint64Subspace.hpp"
#include "../VectorialBooleanFunction.hpp"

/// Implements naive algorithm for finding trees.
namespace BruteForce {

/// Checks if function is constant on affine subspace.
/** Checks for a vectorial Boolean function f whether
 * f|A is constant, where A is an affine subspace
 * of GF(2)^n, n = f.InputSize().
 * \param fun   Function f
 * \param space Underlying vector space of A.
 * \param coset Displacement of A.
 */
bool IsConstantOnSubspace(const VectorialBooleanFunction &fun,
                          const Uint64Subspace &space, uint64_t coset);

/// Internal representation of affine decision trees.
/**
 * Affine subtrees are represented recursively as a value and
 * the two subtrees.
 */
class BinaryDecisionTree {
  typedef std::unique_ptr<BinaryDecisionTree> Subtree;
  static void print_spaces(std::ostream &os, int i);

  void recursive_print(std::ostream &os, int i);

public:
  Subtree left_ = nullptr;  ///< Pointer to the left subtree
  Subtree right_ = nullptr; ///< Pointer to the right subtree

  /// Prints tree in a (somewhat) human-readable form.
  friend std::ostream &operator<<(std::ostream &os, BinaryDecisionTree &x) {
    os << "Binary tree: \n";
    x.recursive_print(os, 0);
    return os;
  }

  uint64_t value_ =
      0; ///< Either label, if inner node, or function value if leaf.

  double associated_cost_ =
      0; ///< Field for saving associated data during searches.

  explicit BinaryDecisionTree() = default;

  /// Constructor for building a leaf.
  explicit BinaryDecisionTree(uint64_t value, double associated_cost = 0) {
    left_ = {};
    right_ = {};
    value_ = value;
    associated_cost_ = associated_cost;
  }

  /// Constructor for making a rooted tree.
  /**
   * \param value Node label.
   * \param left std::unique_ptr to left subtree. Will be moved.
   * \param right std::unique_ptr to right subtree. Will be moved.
   */
  BinaryDecisionTree(uint64_t value, std::unique_ptr<BinaryDecisionTree> left,
                     std::unique_ptr<BinaryDecisionTree> right)
      : value_(value), left_(std::move(left)), right_(std::move(right)) {}

  // Might drop some information, but not value_, the subtrees, and associated
  // cost.
  [[nodiscard]] BinaryDecisionTree DeepImperfectCopy() const {
    if (IsLeaf()) {
      auto result = BinaryDecisionTree(value_);
      result.associated_cost_ = 0;
      return result;
    } else {
      auto result = BinaryDecisionTree(
          value_,
          left_
              ? std::make_unique<BinaryDecisionTree>(left_->DeepImperfectCopy())
              : nullptr,
          right_ ? std::make_unique<BinaryDecisionTree>(
                       right_->DeepImperfectCopy())
                 : nullptr);
      result.associated_cost_ = associated_cost_;
      return result;
    }
  }

  [[nodiscard]] bool IsLeaf() const;

  [[nodiscard]] size_t Depth() const;

  [[nodiscard]] double AveragePathLength() const;

  /// Evaluates the function calculated by the tree.
  /** For a vector x, in an inner node with label a,
   * goes to the left subtree if a*x is zero and to the right subtree otherwise
   * until a terminal node is reached. Its value is the value of the underlying
   * function at position x.
   */
  [[nodiscard]] uint64_t EvaluateAt(uint64_t x) const;

  /// Calculates a lower bound on the input bit size of the underlying function.
  [[nodiscard]] size_t BitSize() const;

  /// Calculates a lower bound on the output bit size of the underlying
  /// function.
  [[nodiscard]] size_t OutputBitSize() const;

  /// Reconstructs the underlying function.
  /**
   * Reconstructs the underlying function. Might not be able to infer the
   * correct input and output sizes and by default (size_hint = 0,
   * output_size_hint = 0) uses the lowest possible values for that.
   */
  [[nodiscard]] VectorialBooleanFunction
  UnderlyingFunction(size_t size_hint = 0, size_t output_size_hint = 0) const;

  /// Checks if trees have equal labels (for inner and terminal nodes).
  bool HasSameLabelsAs(const BinaryDecisionTree &other) {
    if (IsLeaf() != other.IsLeaf()) {
      return false;
    } else if (value_ == other.value_ && !IsLeaf()) {
      return left_->HasSameLabelsAs(*other.left_) &&
             right_->HasSameLabelsAs(*other.right_);
    } else if (value_ == other.value_) {
      return true;
    }
    return false;
  }

  /// Returns number of leaves.
  int leaves() {
    if (IsLeaf())
      return 1;
    int total = 0;
    if (right_) {
      total += right_->leaves();
    }
    if (left_) {
      total += left_->leaves();
    }
    return total;
  }
};

typedef std::unique_ptr<BinaryDecisionTree> BinaryDecisionTreePointer;

NTL::vec_GF2 ComplementIn(const NTL::mat_GF2 &parity_check_1,
                          const NTL::vec_GF2 &new_vector);

BinaryDecisionTreePointer
TreeSearch(const VectorialBooleanFunction &fun,
           const NTL::mat_GF2 &vectors_on_path, const NTL::vec_GF2 &choices,
           int last_choice, const NTL::vec_GF2 &coset_so_far, double &bound,
           int level, bool print = true);

/// Starts search for optimal tree.
/** Finds optimal (lowest average path length) tree representing a given
 * function.
 * \param fun Target function.
 * \param print If set to true, outputs the trees improving
 * the internal size bound.
 */
BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun,
                                      bool print = true);
/// Starts search for best tree whose average path length is below a certain
/// bound.
/** Finds tree better than a certain bound and saves the lowest average path
 * length.
 * \param fun Target function.
 * \param bound Pointer to upper bound on
 * the number of leaves. Is updated to contain the number of leaves of the best
 * tree found.
 * \param print If set to true, outputs the trees improving the
 * internal size bound.
 */
BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun,
                                      double &bound, bool print = true);

/// Starts search for best tree whose average path length is below a certain
/// bound with a fixed root.
/**
 * Finds tree better than a certain bound and saves the average path length.
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

/// Standard report on a fixed component of a vectorial Boolean function.
void analyse_component(const VectorialBooleanFunction &function,
                       uint64_t component);
} // namespace BruteForce

#endif
