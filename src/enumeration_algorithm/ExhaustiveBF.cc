#include "ConditionsLib/enumeration_algorithm/ExhaustiveBF.hpp"

ExhaustiveBF::BinaryDecisionTreeVector
ExhaustiveBF::TreeSearch(const VectorialBooleanFunction &fun,
                         const NTL::mat_GF2 &vectors_on_path,
                         const NTL::vec_GF2 &choices, int last_choice,
                         const NTL::vec_GF2 &coset_so_far,
                         double &bound, // unclear what to do with it
                         int level, bool print) {
  /*if (last_choice > int(fun.InputSize()) - 1) {
    // function has degree at most n.
    return BinaryDecisionTreeVector();
    // Return an empty list of binary decision trees
  }*/
  uint64_t coset = basic_linear_algebra::ConvertToUint(coset_so_far);
  // in space we calculate the dual of vectors_on_path
  NTL::mat_GF2 space;
  NTL::transpose(space, vectors_on_path);
  NTL::kernel(space, space);
  BinaryDecisionTreeVector result;
  if (BruteForce::IsConstantOnSubspace(
          fun, Uint64Subspace(space, space.NumCols()), coset)) {
    result.push_back(
        std::move(std::make_unique<BinaryDecisionTree>(fun(coset))));
    return result;
  }
  auto possible_masks(vectors_on_path);
  possible_masks = basic_linear_algebra::ComplementSpace(possible_masks);
  const uint64_t num = 1u << possible_masks.NumRows();
  // It's crucial that possible_masks has linearly independent rows
  // Otherwise 0 will be enumerated
  auto new_vectors_on_path(vectors_on_path);
  new_vectors_on_path.SetDims(vectors_on_path.NumRows() + 1,
                              vectors_on_path.NumCols());
  auto new_choices(choices);
  for (uint64_t k = 1; k < num; ++k) {
    auto splitting_alpha = basic_linear_algebra::ConvertToUint(
        basic_linear_algebra::Embed(k, possible_masks));
    new_choices[last_choice + 1] = 0;
    auto last = basic_linear_algebra::ConvertToNtl(
        splitting_alpha, new_vectors_on_path.NumCols());
    new_vectors_on_path[new_vectors_on_path.NumRows() - 1] = last;
    auto q = BruteForce::ComplementIn(space, last);
    // The following calculates which coset is U0 and which is U1:
    NTL::vec_GF2 new_coset_0 = coset_so_far + q;
    NTL::vec_GF2 new_coset_1 = coset_so_far;
    if (NTL::IsZero(last * coset_so_far)) {
      new_coset_0.swap(new_coset_1);
    }
    // Build the subtrees:
    auto new_bound = std::min(
        2 * bound - 2, double(fun.InputSize() - new_vectors_on_path.NumRows()));
    auto left_trees =
        TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                   new_coset_0, new_bound, 1 + level, print);
    // If bound cannot be fulfilled, next alpha.
    if (left_trees.empty()) {
      continue;
    }
    for (const auto &left_tree : left_trees) {
      new_choices[last_choice + 1] = 1;
      // the vector of right trees is only dependent
      // on left_tree.associated_cost_ and bound
      // If we are really only looking for
      // depth, I just have to do this ones and not for every tree.
      new_bound =
          std::min(2 * bound - 2 - left_tree->associated_cost_,
                   double(fun.InputSize() - new_vectors_on_path.NumRows()));
      auto right_trees =
          TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                     new_coset_1, new_bound, 1 + level, print);
      if (right_trees.empty()) {
        continue;
      }
      for (const auto &right_tree : right_trees) {
        const auto total_costs = 1 + 0.5 * (left_tree->associated_cost_ +
                                            right_tree->associated_cost_);
        if (total_costs > bound)
          continue;
        // This copying is annoying, but I think unfortunately necessary
        // Until I have rewritten the trees with shared pointers
        // Not necessary for right_tree, perhaps.
        BruteForce::BinaryDecisionTreePointer left_tree_copy =
            std::make_unique<BruteForce::BinaryDecisionTree>(
                left_tree->DeepImperfectCopy());
        BruteForce::BinaryDecisionTreePointer right_tree_copy =
            std::make_unique<BruteForce::BinaryDecisionTree>(
                right_tree->DeepImperfectCopy());
        auto candidate = std::make_unique<BinaryDecisionTree>(
            splitting_alpha, std::move(left_tree_copy),
            std::move(right_tree_copy));
        candidate->associated_cost_ = total_costs;
        result.push_back(std::move(candidate));
        // bound = total_costs;
        if (level == 0 && print) {
          std::cout << "Tree \n"
                    << *result.back()
                    << "\nwith costs: " << result.back()->associated_cost_
                    << std::endl;
        }
      }
    }
  }
  return result;
}
ExhaustiveBF::BinaryDecisionTreeVector
ExhaustiveBF::StartSearch(const VectorialBooleanFunction &fun, bool print) {
  double bound = fun.InputSize() - 0.000001;
  return StartSearch(fun, bound, print);
}
ExhaustiveBF::BinaryDecisionTreeVector
ExhaustiveBF::StartSearch(const VectorialBooleanFunction &fun, double &bound,
                          bool print) {
  NTL::mat_GF2 vectors_on_path;
  vectors_on_path.SetDims(0, fun.InputSize());
  assert(vectors_on_path.NumCols() == fun.InputSize());
  NTL::vec_GF2 choices;
  choices.SetLength(fun.InputSize()); // I hope zero initializes
  return TreeSearch(fun, vectors_on_path, choices, -1, choices, bound, 0,
                    print);
}
