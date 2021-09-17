/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <chrono>
#include <iostream>
#include <random>

#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/SboxTools.hpp"
const auto dim = 6u;
const auto iterations = 100u;
const bool balanced = false;

int main() {
  std::cout << "Dimension: " << dim << std::endl;
  std::cout << "Iterations: " << iterations << std::endl;
  std::random_device re;
  std::vector<VectorialBooleanFunction> random_functions;
  for (size_t i = 0; i < iterations; ++i) {
    random_functions.emplace_back(balanced ? RandomBalancedVBF(dim, re)
                                           : RandomVBF(dim, re));
  }
  std::cout << "Start timing ..." << std::endl;
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (const auto &random_function: random_functions) {
    MinLeaves::StartSearch(random_function, false);
  }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  const auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
  const auto et_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
  std::cout << "Elapsed time: " << elapsed_time << " microseconds" << std::endl;
  std::cout << "Time/iteration: " << (double(elapsed_time)/iterations) << " microseconds" << std::endl;
  std::cout << "Elapsed time: " << et_seconds << " seconds" << std::endl;
  std::cout << "Time/iteration: " << (double(et_seconds)/iterations) << " seconds" << std::endl;

  return 0;
}

