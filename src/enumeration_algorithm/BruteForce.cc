#include <random>

#include "ConditionsLib/ComponentFunction.hpp"
#include "ConditionsLib/enumeration_algorithm/BruteForce.hpp"

bool BruteForce::IsConstantOnSubspace(const VectorialBooleanFunction &fun,
                                      const Uint64Subspace &space,
                                      const uint64_t coset) {
  const uint64_t f0 = fun(coset);
  const uint64_t num = 1u << space.Dimension();
  for (uint64_t i = 1; i < num; ++i) {
    if (fun(space.ElementK(i) ^ coset) != f0) {
      return false;
    }
  }
  return true;
}

NTL::vec_GF2 BruteForce::ComplementIn(const NTL::mat_GF2 &parity_check_1,
                                      const NTL::vec_GF2 &new_vector) {
  ssize_t certainty_level = 128;
  NTL::vec_GF2 result;
  result.SetLength(new_vector.length());
  NTL::vec_GF2 random_vector;
  random_vector.SetLength(parity_check_1.NumRows());
  do {
    certainty_level--;
    NTL::random(random_vector, parity_check_1.NumRows());
    result = random_vector * parity_check_1;
  } while (certainty_level >= 0 && NTL::IsZero(new_vector * result));
  if (certainty_level < 0) {
    std::cout << parity_check_1 << std::endl;
    std::cout << new_vector << std::endl;
    throw std::logic_error("Subroutine used for determininng the cosets failed."
                           "Is the function constant?");
  }
  return result;
}

BruteForce::BinaryDecisionTreePointer BruteForce::TreeSearch(
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
  if (IsConstantOnSubspace(fun, Uint64Subspace(space, space.NumCols()),
                           coset)) {
    return std::make_unique<BinaryDecisionTree>(fun(coset));
  }
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
    auto q = ComplementIn(space, last);
    // The following calculates which coset is U0 and which is U1:
    NTL::vec_GF2 new_coset_0 = coset_so_far + q;
    NTL::vec_GF2 new_coset_1 = coset_so_far;
    if (NTL::IsZero(last * coset_so_far)) {
      new_coset_0.swap(new_coset_1);
    }
    // Build the subtrees:
    auto new_bound = std::min(
        2 * bound - 2, double(fun.InputSize() - new_vectors_on_path.NumRows()));
    auto left_tree =
        TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                   new_coset_0, new_bound, 1 + level, print);
    // If bound cannot be fulfilled, next alpha.
    if (!left_tree) {
      continue;
    }
    new_choices[last_choice + 1] = 1;
    new_bound =
        std::min(2 * bound - 2 - left_tree->associated_cost_,
                 double(fun.InputSize() - new_vectors_on_path.NumRows()));
    auto right_tree =
        TreeSearch(fun, new_vectors_on_path, new_choices, last_choice + 1,
                   new_coset_1, new_bound, 1 + level, print);
    if (!right_tree) {
      continue;
    }
    const auto total_costs =
        1 + 0.5 * (left_tree->associated_cost_ + right_tree->associated_cost_);
    if (total_costs > bound)
      continue;
    if (!result) {
      result = std::make_unique<BinaryDecisionTree>(
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

BruteForce::BinaryDecisionTreePointer
BruteForce::StartSearchWithFixedRoot(const VectorialBooleanFunction &fun,
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
  double left_bound = std::min(2 * bound - 2, double(fun.InputSize() - 1));
  auto left_tree =
      TreeSearch(fun, vectors_on_path, choices, 0, coset, left_bound, 1, print);
  if (!left_tree) {
    throw std::logic_error("Nothing found.");
  }
  double right_bound = std::min(2 * bound - left_tree->associated_cost_ - 2,
                                double(fun.InputSize() - 1));
  choices[0] = 1;
  coset = basic_linear_algebra::ConvertToNtl(candidate, fun.InputSize());
  auto right_tree = TreeSearch(fun, vectors_on_path, choices, 0, coset,
                               right_bound, 1, print);
  if (!right_tree) {
    throw std::logic_error("Nothing found.");
  }
  auto total_cost =
      1 + 0.5 * (left_tree->associated_cost_ + right_tree->associated_cost_);
  auto result = std::make_unique<BinaryDecisionTree>(root, std::move(left_tree),
                                                     std::move(right_tree));
  result->associated_cost_ = total_cost;
  return result;
}

BruteForce::BinaryDecisionTreePointer
BruteForce::StartSearch(const VectorialBooleanFunction &fun, double &bound,
                        bool print) {
  NTL::mat_GF2 vectors_on_path;
  vectors_on_path.SetDims(0, fun.InputSize());
  assert(vectors_on_path.NumCols() == fun.InputSize());
  NTL::vec_GF2 choices;
  choices.SetLength(fun.InputSize()); // I hope zero initializes
  return TreeSearch(fun, vectors_on_path, choices, -1, choices, bound, 0,
                    print);
}

BruteForce::BinaryDecisionTreePointer
BruteForce::StartSearch(const VectorialBooleanFunction &fun, bool print) {
  double bound = fun.InputSize();
  return StartSearch(fun, bound, print);
}

void BruteForce::analyse_component(const VectorialBooleanFunction &function,
                                   uint64_t component) {
  std::cout << "############################\n";
  ComponentFunction cfun(function, component);
  std::cout << function << '\n';
  std::cout << "Component: " << component << '\n';
  std::cout << cfun << '\n';
  try {
    auto gs = StartSearch(cfun);
    std::cout << *gs << '\n';
    std::cout << "Reconstructed function "
              << gs->UnderlyingFunction(function.InputSize()) << '\n';
    std::cout << "Average Path Length: " << gs->AveragePathLength() << '\n';
    if (cfun == gs->UnderlyingFunction(function.InputSize())) {
      std::cout << "Calculates the correct function!" << std::endl;
    } else {
      std::cerr << "The solution is wrong, unfortunately." << std::endl;
    }
  } catch (std::logic_error &e) {
    std::cout << e.what() << std::endl;
  }
}
void BruteForce::BinaryDecisionTree::print_spaces(std::ostream &os, int i) {
  for (int j = 0; j < i; ++j) {
    os << "|   ";
  }
}
void BruteForce::BinaryDecisionTree::recursive_print(std::ostream &os, int i) {
  print_spaces(os, i);
  if (IsLeaf())
    os << "Constant = ";
  os << std::hex << value_ << "\n";
  if (left_)
    left_->recursive_print(os, i + 1);
  if (right_)
    right_->recursive_print(os, i + 1);
}
bool BruteForce::BinaryDecisionTree::IsLeaf() const {
  return !left_ && !right_;
}
size_t BruteForce::BinaryDecisionTree::Depth() const {
  size_t result = 0;
  if (left_)
    result = left_->Depth() + 1;
  if (right_)
    result = std::max(result, right_->Depth() + 1);
  return result;
}
double BruteForce::BinaryDecisionTree::AveragePathLength() const {
  double result = 0;
  if (left_)
    result += left_->AveragePathLength() + 1;
  if (right_)
    result += right_->AveragePathLength() + 1;
  if (left_ && right_)
    result /= 2.0;
  return result;
}
uint64_t BruteForce::BinaryDecisionTree::EvaluateAt(uint64_t x) const {
  if (IsLeaf()) {
    return value_;
  }
  if ((std::popcount(x & value_) % 2) == 0) {
    return left_->EvaluateAt(x);
  }
  return right_->EvaluateAt(x);
}
size_t BruteForce::BinaryDecisionTree::BitSize() const {
  if (IsLeaf())
    return 0;
  return std::max<size_t>(
      {64u - std::countl_zero(value_), left_->BitSize(), right_->BitSize()});
}
size_t BruteForce::BinaryDecisionTree::OutputBitSize() const {
  if (IsLeaf())
    return 64u - std::countl_zero(value_);
  return std::max<size_t>(left_->OutputBitSize(), right_->OutputBitSize());
}
VectorialBooleanFunction BruteForce::BinaryDecisionTree::UnderlyingFunction(
    size_t size_hint, size_t output_size_hint) const {
  size_t bit_size = size_hint ? size_hint : BitSize();
  size_t output_bit_size =
      output_size_hint ? output_size_hint : OutputBitSize();
  std::vector<uint64_t> base(1u << bit_size);
  std::iota(base.begin(), base.end(), 0);
  std::transform(base.begin(), base.end(), base.begin(),
                 [&](auto x) { return EvaluateAt(x); });
  return VectorialBooleanFunction(base, bit_size, output_bit_size);
}
