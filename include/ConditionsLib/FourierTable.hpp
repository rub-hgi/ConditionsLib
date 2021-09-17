/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_FOURIER_TABLE_HPP_
#define CONDITIONS_LIB_FOURIER_TABLE_HPP_

#include <bit>
#include <cstdint>
#include <vector>
#include <iostream>
#include <random>
#include <utility>
#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>
#include <numeric>
#include "VectorialBooleanFunction.hpp"
#include "FourierTable.hpp"
#include "basic_linear_algebra.hpp"

namespace FourierTableAux {
typedef uint64_t(*popcount_like)(uint64_t);
uint64_t MyPopcount(uint64_t x);
}

template<FourierTableAux::popcount_like po = FourierTableAux::MyPopcount>
class FourierTable {
 private:
  std::vector<std::vector<int64_t>> table;
 public:
  template<class Function>
  explicit FourierTable(const Function &f);

  int64_t operator()(uint64_t alpha, uint64_t beta) const;
};

template<FourierTableAux::popcount_like po>
template<class Function>
FourierTable<po>::FourierTable(const Function &f) {
  table.resize(1u << f.InputSize());
  for (size_t alpha = 0; alpha < 1u << f.InputSize(); ++alpha) {
    table[alpha].resize(1u << f.OutputSize());
    std::fill(table[alpha].begin(),
              table[alpha].end(), int64_t(0));
    for (size_t beta = 0; beta < 1u << f.OutputSize(); ++beta) {
      for (size_t x = 0; x < 1u << f.InputSize(); ++x) {
        table[alpha][beta] += (po(x & alpha)
            ^ po(f(x) & beta)) & 1u ? -1 : 1;
      }
    }
  }
}
template<FourierTableAux::popcount_like po>
int64_t FourierTable<po>::operator()(uint64_t alpha, uint64_t beta) const {
  return table[alpha][beta];
}

#endif //CONDITIONS__FOURIERTABLE_HPP_
