/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <iostream>
#include <fstream>

#include "ConditionsLib/enumeration_algorithm/ExhaustiveBF.hpp"
#include "ConditionsLib/ComponentFunction.hpp"
#include "ConditionsLib/SboxTools.hpp"

std::vector<uint64_t> rectangle_sbox{
  0x6, 0x5, 0xC, 0xA, 0x1, 0xE, 0x7, 0x9, 0xB, 0x0, 0x3, 0xD, 0x8, 0xF, 0x4, 0x2
};
// old Rectangle:
/*
{
  0x9,0x4,0xF,0xA,0xE,0x1,0x0,0x6,0xC,0x7,0x3,0x8,0x2,0xB,0x5,0xD
};*/

std::vector<uint64_t> rectangle_inverse(16);

void init() {
  for (uint64_t i = 0; i < 16; ++i)
    rectangle_inverse[rectangle_sbox[i]] = i;

}

std::vector<uint64_t> prepare_f_i(uint64_t delta) {
  std::vector<uint64_t> result(16);
  for (uint64_t i = 0; i < 16;++i) {
    result[i] = rectangle_sbox[rectangle_inverse[i]^delta] ;
  }
  return result;
}
std::vector<uint64_t> prepare_f_i_plus_x(uint64_t delta) {
  std::vector<uint64_t> result(16);
  for (uint64_t i = 0; i < 16;++i) {
    result[i] = rectangle_sbox[rectangle_inverse[i]^delta] ^ i ;
  }
  return result;
}
std::vector<uint64_t> prepare_f_i_straight(uint64_t delta) {
  std::vector<uint64_t> result(16);
  for (uint64_t i = 0; i < 16;++i) {
    result[i] = rectangle_inverse[rectangle_sbox[i]^delta] ;
  }
  return result;
}
std::vector<uint64_t> prepare_f_i_plus_x_straight(uint64_t delta) {
  std::vector<uint64_t> result(16);
  for (uint64_t i = 0; i < 16;++i) {
    result[i] = rectangle_inverse[rectangle_sbox[i]^delta] ^ i ;
  }
  return result;
}

void analyse_old(std::ostream& os, std::vector<uint64_t> & sbox, int dim) {

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
      auto gs = BruteForce::StartSearch(cfun);
      auto ls = SboxTools::LinearStructures(cfun, 4);
      os << "Result: \n" << *gs << std::endl;
      os << "Costs (Average guessing): " << gs->associated_cost_ << std::endl;
      os << "0-Linear Structures" << std::endl;
      for (auto x: ls.first) {
        os << std::dec << x << std::endl;
      }
      os << "1-Linear Structures" << std::endl;
      for (auto x: ls.second) {
        os << std::dec << x << std::endl;
      }
    } catch(std::logic_error &e) {
      os << "For component " << component << " nothing found.\n";
    }
  }
}


std::vector<std::vector<uint64_t>> GetAllPaths(
    const BruteForce::BinaryDecisionTree &tree
    ) {
  if (tree.IsLeaf()) {
      return std::vector<std::vector<uint64_t>>{
          std::vector<uint64_t>(0)
      };
  }
  auto left_paths = GetAllPaths(*tree.left_);
  for (auto &x: left_paths) {
    x.push_back(tree.value_);
  }
  auto right_paths = GetAllPaths(*tree.right_);
  for (auto &x: right_paths) {
    x.push_back(tree.value_);
  }
  std::vector<std::vector<uint64_t>> result(0);
  result.reserve(left_paths.size() + right_paths.size());
  result.insert(result.end(), left_paths.begin(), left_paths.end());
  result.insert(result.end(), right_paths.begin(), right_paths.end());
  return result;
}

bool pathCondition(const std::vector<uint64_t> &path) {
  if (path.size() > 3) return false;
  // Find first non_zero value
  uint64_t uniq = 0;
  for (; uniq < path.size(); ++uniq) {
    if ((path[uniq]&3) != 0) {
      break;
    }
  }
  for (uint64_t i = uniq+1; i < path.size(); ++i) {
    if ((path[i]&3) == 0) continue;
    if ((path[i]&3) != (path[uniq]&3)) {
      return true; // Do not keep
    }
  }
  return false;
}

void filter(ExhaustiveBF::BinaryDecisionTreeVector &to_be_filtered) {
  for (auto ptr = to_be_filtered.begin(); ptr != to_be_filtered.end(); ) {
    bool cond1 = ptr->get()->left_->Depth() == 3;
    bool cond2 = ptr->get()->right_->Depth() == 3;
    bool cond3 = ptr->get()->left_->left_->Depth() == 2;
    bool cond4 = ptr->get()->left_->right_->Depth() == 2;
    bool cond5 = ptr->get()->right_->left_->Depth() == 2;
    bool cond6 = ptr->get()->right_->right_->Depth() == 2;
    auto non_trivial_paths = GetAllPaths(*ptr->get());
    auto cond7 =
        std::any_of(non_trivial_paths.begin(),
                    non_trivial_paths.end(),
                    pathCondition);
    if (cond1 && cond2 || cond7) {
    //if (cond3 && cond4 && cond5 && cond6 || cond7) {
      ptr = to_be_filtered.erase(ptr);
    } else {
      // What to keep
      ptr++;
    }
  }
}

void Analyse(std::ostream& os, std::vector<uint64_t> & sbox, int dim) {

  for (auto x: sbox) {
    os << x << " ";
  }
  os << "\n";
  VectorialBooleanFunction fun(sbox, dim, dim);
    try {
      os << "-------------------\n";
      auto gs = BruteForceOptimization::StartSearch(fun, false);
      //filter(gs);
      auto ls = SboxTools::LinearStructures(fun, 4);
      os << "0-Linear Structures" << std::endl;
      for (auto x: ls.first) {
        os << std::dec << x << std::endl;
      }
      os << "1-Linear Structures" << std::endl;
      for (auto x: ls.second) {
        os << std::dec << x << std::endl;
      }
      //os << "Results: \n" << gs.size() << std::endl;
      os << *gs << "\n";
      os << "Costs (Average guessing): " << gs->associated_cost_ << std::endl;
      os << "Costs (Domsize): " << std::log2((1u << fun.InputSize()) / ls.first.size()) << std::endl;
      /*for (const auto &g: gs) {
        os << *g << "\n";
        os << "Costs (Average guessing): " << g->associated_cost_ << std::endl;
      }*/
    } catch(std::logic_error &e) {
      os << "Nothing found for whole input.\n";
    }
}

int main() {
  init();
  const auto dim = 4u;
  for (uint64_t delta = 1; delta < 16; ++delta) {
    {
      std::vector<uint64_t> sbox = prepare_f_i_straight(delta);
      std::ofstream ffile("NewRectStraightF" + std::to_string(delta) + ".txt");
      analyse_old(ffile, sbox, dim);
    }
    {
      std::vector<uint64_t> sbox = prepare_f_i_plus_x_straight(delta);
      std::ofstream ffile("NewRectStraightF" + std::to_string(delta) + "PlusX_FOUR.txt");
      analyse_old(ffile, sbox, dim);
    }
  }
  for (uint64_t delta = 1; delta < 16; ++delta) {
    {
      std::vector<uint64_t> sbox = prepare_f_i(delta);
      std::ofstream ffile("NewRectF" + std::to_string(delta) + ".txt");
      analyse_old(ffile, sbox, dim);
    }
    {
      std::vector<uint64_t> sbox = prepare_f_i_plus_x(delta);
      std::ofstream ffile("NewRectF" + std::to_string(delta) + "PlusX_TWO.txt");
      analyse_old(ffile, sbox, dim);
    }
  }
  for (uint64_t delta = 1; delta < 16; ++delta) {
    {
      std::vector<uint64_t> sbox = prepare_f_i_straight(delta);
      std::ofstream ffile("NewRectAllStraightF" + std::to_string(delta) + ".txt");
      Analyse(ffile, sbox, dim);
    }
    {
      std::vector<uint64_t> sbox = prepare_f_i_plus_x_straight(delta);
      std::ofstream ffile("NewRectAllStraightF" + std::to_string(delta) + "PlusX.txt");
      Analyse(ffile, sbox, dim);
    }
  }
  for (uint64_t delta = 1; delta < 16; ++delta) {
    {
      std::vector<uint64_t> sbox = prepare_f_i(delta);
      std::ofstream ffile("NewRectAllF" + std::to_string(delta) + ".txt");
      Analyse(ffile, sbox, dim);
    }
    {
      std::vector<uint64_t> sbox = prepare_f_i_plus_x(delta);
      std::ofstream ffile("NewRectAllF" + std::to_string(delta) + "PlusX.txt");
      Analyse(ffile, sbox, dim);
    }
  }
  return 0;
}

