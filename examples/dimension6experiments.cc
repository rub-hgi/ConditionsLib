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

int main(int argc, char** argv) {
  const size_t outdim = 1;
  const size_t dim = 6;
  if (argc < 2) {
    // First command line argument has to be a file.
    std::cerr << "First argument has to be a filename." << std::endl;
    return 1;
  }
  std::string filename;
  filename = argv[1];
  std::ifstream file(filename);
  ssize_t frst_line = 0; // start numbering from 0
  ssize_t stop_line = -1; // exclusive
  if (argc >= 3) {
    frst_line = atoll(argv[2]);
  }
  if (argc >= 4) {
    stop_line = atoll(argv[3]);
  }
  auto functions =
      SboxTools::ReadFileNodelim(file, dim, outdim, frst_line, stop_line);
  for (const auto &x: functions) {
    auto line = SboxTools::Analyse<SboxTools::csv_line>(x);
    {
      std::cout << line;
    }
  }
}
