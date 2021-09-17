#include "ConditionsLib/enumeration_algorithm/MinLeaves.hpp"
#include "ConditionsLib/enumeration_algorithm/BruteForceOptimizations.hpp"
#include <random>
#include <set>

MinLeaves::BinaryDecisionTreePointer MinLeaves::TreeSearch(
    const VectorialBooleanFunction &fun, const NTL::mat_GF2 &vectors_on_path,
    const NTL::vec_GF2 &choices, int last_choice,
    const NTL::vec_GF2 &coset_so_far, double &bound, int level, bool print) {
  if (last_choice > int(fun.InputSize()) - 1) {
    // function has degree at most n.
    return nullptr;
  }
  uint64_t coset = basic_linear_algebra::ConvertToUint(coset_so_far);
  // in space we calculate the dual of vectors_on_path
  NTL::mat_GF2 space;
  NTL::transpose(space, vectors_on_path);
  NTL::kernel(space, space);
  auto linear_domain = Uint64Subspace(space, space.NumCols());
  if (BruteForce::IsConstantOnSubspace(fun, linear_domain, coset)) {
    return std::make_unique<MinLeaves::BinaryDecisionTree>(fun(coset), 1.0);
    // associated cost has to be set to 1 here, the tree with only
    // one node has one leaf
  }
  auto ls0 =
      BruteForceOptimization::zero_linear_structures(fun, linear_domain, coset);
  auto possible_masks(vectors_on_path);
  possible_masks = basic_linear_algebra::ComplementSpace(possible_masks);
  BinaryDecisionTreePointer result = nullptr;
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
    // (A representation of) the difference between cosets of the next level is
    // q. If two sister nodes only differ by a linear structure, there is a
    // better tree (with fewer levels). q might be a bad choice of
    // representation though
    auto q_num = basic_linear_algebra::ConvertToUint(q);
    if (ls0.contains(q_num)) {
      continue;
    }
    bool differs_by_linear_structure = false;
    for (auto &r : ls0) {
      auto r2 = basic_linear_algebra::ConvertToNtl(r, fun.InputSize());
      if (NTL::IsZero(new_vectors_on_path * (r2 + q))) {
        differs_by_linear_structure = true;
        break;
      }
    }
    if (differs_by_linear_structure)
      continue;
    // The following calculates which coset is U0 and which is U1:
    NTL::vec_GF2 new_coset_0 = coset_so_far + q;
    NTL::vec_GF2 new_coset_1 = coset_so_far;
    if (NTL::IsZero(last * coset_so_far)) {
      new_coset_0.swap(new_coset_1);
    }

    // Build the subtrees:
    auto new_bound = std::min(
        bound, std::pow(2.0, fun.InputSize() - new_vectors_on_path.NumRows()));
    // A subtree with more than 2^{n - d} leaves is complete
    // TODO: Encode this optimization somewhere else without the need to
    // TODO: Use std::pow
    auto left_tree =
        TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                   new_coset_0, new_bound, 1 + level, print);
    // If bound cannot be fulfilled, next alpha.
    if (!left_tree) {
      continue;
    }
    new_choices[last_choice + 1] = 1;
    new_bound = std::min(
        bound - left_tree->associated_cost_,
        std::pow(2.0, fun.InputSize() - new_vectors_on_path.NumRows()));
    auto right_tree =
        TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                   new_coset_1, new_bound, 1 + level, print);
    if (!right_tree) {
      continue;
    }
    const auto total_costs =
        left_tree->associated_cost_ + right_tree->associated_cost_;
    if (total_costs > bound)
      continue;
    if (!result) {
      result = std::make_unique<MinLeaves::BinaryDecisionTree>(
          splitting_alpha, std::move(left_tree), std::move(right_tree));
      result->associated_cost_ = total_costs;
      bound = total_costs;
      if (level == 0 && print) {
        std::cout << "Best tree so far\n"
                  << *result << "\nwith costs: " << result->associated_cost_
                  << std::endl;
      }
    } else {
      auto candidate = std::make_unique<BinaryDecisionTree>(
          splitting_alpha, std::move(left_tree), std::move(right_tree));
      if (candidate && (total_costs < result->associated_cost_)) {
        result = std::move(candidate);
        result->associated_cost_ = total_costs;
        bound = total_costs;
        if (level == 0 && print) {
          std::cout << "Best tree so far\n"
                    << *result << "\nwith costs: " << result->associated_cost_
                    << std::endl;
        }
      }
    }
  }
  return result;
}

MinLeaves::BinaryDecisionTreePointer
MinLeaves::StartSearchWithFixedRoot(const VectorialBooleanFunction &fun,
                                    double &bound, uint64_t root, bool print) {
  NTL::mat_GF2 vectors_on_path;
  vectors_on_path.SetDims(1, fun.InputSize());
  assert(vectors_on_path.NumCols() == fun.InputSize());
  NTL::vec_GF2 choices;
  choices.SetLength(fun.InputSize()); // I hope zero initializes
  NTL::vec_GF2 coset;
  coset.SetLength(fun.InputSize());
  if (root == 0) {
    throw std::logic_error("Root is 0");
  }
  // We need to calculate a coset
  uint16_t i = 0;
  uint64_t candidate = 0;
  std::random_device rd;
  std::mt19937_64 mt(rd());
  std::uniform_int_distribution<uint64_t> dist(0, 1u << fun.InputSize());
  for (i = 1; i; ++i) {
    candidate = dist(mt);
    if (std::popcount(candidate & root) & 1)
      break;
  }
  if (i == 0) {
    throw std::logic_error(
        "Bad luck, broken random device or insanely high dimension.");
  }
  vectors_on_path[0] =
      basic_linear_algebra::ConvertToNtl(root, vectors_on_path.NumCols());
  double left_bound = std::min(bound, std::pow(2.0, fun.InputSize() - 1));
  auto left_tree =
      TreeSearch(fun, vectors_on_path, choices, 0, coset, left_bound, 1, print);
  if (!left_tree) {
    throw std::logic_error("Nothing found.");
  }
  double right_bound = std::min(bound - left_tree->associated_cost_,
                                std::pow(2.0, fun.InputSize() - 1));
  choices[0] = 1;
  coset = basic_linear_algebra::ConvertToNtl(candidate, fun.InputSize());
  auto right_tree = TreeSearch(fun, vectors_on_path, choices, 0, coset,
                               right_bound, 1, print);
  if (!right_tree) {
    throw std::logic_error("Nothing found.");
  }
  auto total_cost = left_tree->associated_cost_ + right_tree->associated_cost_;
  auto result = std::make_unique<BinaryDecisionTree>(root, std::move(left_tree),
                                                     std::move(right_tree));
  result->associated_cost_ = total_cost;
  return result;
}

MinLeaves::BinaryDecisionTreePointer
MinLeaves::StartSearch(const VectorialBooleanFunction &fun, double &bound,
                       bool print) {
  NTL::mat_GF2 vectors_on_path;
  vectors_on_path.SetDims(0, fun.InputSize());
  assert(vectors_on_path.NumCols() == fun.InputSize());
  NTL::vec_GF2 choices;
  choices.SetLength(fun.InputSize());
  return TreeSearch(fun, vectors_on_path, choices, -1, choices, bound, 0,
                    print);
}

MinLeaves::BinaryDecisionTreePointer
MinLeaves::StartSearch(const VectorialBooleanFunction &fun, bool print) {
  double bound = std::pow(2.0, fun.InputSize());
  return StartSearch(fun, bound, print);
}

void MinLeaves::RecursiveStump(
    const VectorialBooleanFunction &fun,
    const MinLeaves::BinaryDecisionTreePointer &stump,
    const NTL::mat_GF2 &vectors_on_path, const NTL::vec_GF2 &choices,
    int last_choice, const NTL::vec_GF2 &coset_so_far, double &bound, int level,
    bool print) {
  if (last_choice > int(fun.InputSize()) - 1) {
    // function has degree at most n.
    throw std::logic_error("No no no no");
  }
  uint64_t coset = basic_linear_algebra::ConvertToUint(coset_so_far);
  // in space we calculate the dual of vectors_on_path
  NTL::mat_GF2 space;
  NTL::transpose(space, vectors_on_path);
  NTL::kernel(space, space);
  auto linear_domain = Uint64Subspace(space, space.NumCols());
  if (BruteForce::IsConstantOnSubspace(fun, linear_domain, coset)) {
    stump->left_ = nullptr;
    stump->right_ = nullptr;
    stump->value_ = fun(coset);
  } // I also cut
  auto ls0 =
      BruteForceOptimization::zero_linear_structures(fun, linear_domain, coset);
  auto possible_masks(vectors_on_path);
  possible_masks = basic_linear_algebra::ComplementSpace(possible_masks);
  BinaryDecisionTreePointer result = nullptr;
  const uint64_t num = 1u << possible_masks.NumRows();
  // It's crucial that possible_masks has linearly independent rows
  // Otherwise 0 will be enumerated
  auto new_vectors_on_path(vectors_on_path);
  new_vectors_on_path.SetDims(vectors_on_path.NumRows() + 1,
                              vectors_on_path.NumCols());
  auto new_choices(choices);
  auto splitting_alpha = stump->value_;
  new_choices[last_choice + 1] = 0;
  auto last = basic_linear_algebra::ConvertToNtl(splitting_alpha,
                                                 new_vectors_on_path.NumCols());
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
      bound, std::pow(2.0, fun.InputSize() - new_vectors_on_path.NumRows()));
  BinaryDecisionTreePointer left_tree;
  if (stump->left_) {
    RecursiveStump(fun, stump->left_, new_vectors_on_path, new_choices,
                   last_choice + 1, new_coset_0, new_bound, 1 + level, print);
  } else {
    stump->left_ =
        TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                   new_coset_0, new_bound, 1 + level, print);
  }
  // If bound cannot be fulfilled, next alpha.
  if (!stump->left_) {
    throw std::logic_error("This is not supposed to happen.");
  }
  new_choices[last_choice + 1] = 1;
  new_bound =
      std::min(bound - stump->left_->associated_cost_,
               std::pow(2.0, fun.InputSize() - new_vectors_on_path.NumRows()));
  if (stump->right_) {
    RecursiveStump(fun, stump->right_, new_vectors_on_path, new_choices,
                   last_choice + 1, new_coset_1, new_bound, 1 + level, print);
  } else {
    stump->right_ =
        TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                   new_coset_1, new_bound, 1 + level, print);
  }
  if (!stump->right_) {
    throw std::logic_error("This is not supposed to happen! bad");
  }
  auto total_cost =
      stump->left_->associated_cost_ + stump->right_->associated_cost_;
  stump->associated_cost_ = total_cost;
}

MinLeaves::BinaryDecisionTreePointer MinLeaves::StartSearchWithFixedStump(
    const VectorialBooleanFunction &fun,
    const MinLeaves::BinaryDecisionTreePointer &stump, double &bound,
    bool print) {
  MinLeaves::BinaryDecisionTreePointer temp =
      std::make_unique<MinLeaves::BinaryDecisionTree>(
          stump->DeepImperfectCopy());
  // Not efficient, but I don't care
  // Just go through the leaves in-order and start the optimal search with that
  // bound everywhere.
  NTL::mat_GF2 vectors_on_path;
  vectors_on_path.SetDims(0, fun.InputSize());
  assert(vectors_on_path.NumCols() == fun.InputSize());
  NTL::vec_GF2 choices;
  choices.SetLength(fun.InputSize()); // I hope zero initializes
  RecursiveStump(fun, temp, vectors_on_path, choices, -1, choices, bound, 0,
                 print);
  return temp;
}
