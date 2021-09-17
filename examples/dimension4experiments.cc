/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <fstream>
#include <iostream>
#include <omp.h>

#include "ConditionsLib/SboxTools.hpp"

int main() {
  const size_t outdim = 1;
  const size_t dim = 4;
  std::string filename = "../FunctionData/F4b.txt";
  std::ifstream file(filename);
  auto functions=
      SboxTools::ReadFile(file, dim, outdim);
  SboxTools::csv_line::PrintHeader(std::cout);
  for (const auto &x_raw: functions) {
    VectorialBooleanFunction x(x_raw, dim, outdim);
    auto line = SboxTools::Analyse<SboxTools::csv_line>(x);
    {
      std::cout << line;
    }
  }
}
