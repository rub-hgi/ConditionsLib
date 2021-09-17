/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#include <sstream>
#include <omp.h>
#include <vector>
#include <fstream>
#include "ConditionsLib/VectorialBooleanFunction.hpp"
#include "ConditionsLib/enumeration_algorithm/ExhaustiveBF.hpp"
#include "ConditionsLib/SboxTools.hpp"
#include "ConditionsLib/ComponentChoice.hpp"
#include "ConditionsLib/ComponentFunction.hpp"

std::vector<std::vector<uint64_t>> serpent_sb {
    {0x3,0x8,0xf,0x1,0xa,0x6,0x5,0xb,0xe,0xd,0x4,0x2,0x7,0x0,0x9,0xc},
    {0xf,0xc,0x2,0x7,0x9,0x0,0x5,0xa,0x1,0xb,0xe,0x8,0x6,0xd,0x3,0x4},
    {0x8,0x6,0x7,0x9,0x3,0xc,0xa,0xf,0xd,0x1,0xe,0x4,0x0,0xb,0x5,0x2},
    {0x0,0xf,0xb,0x8,0xc,0x9,0x6,0x3,0xd,0x1,0x2,0x4,0xa,0x7,0x5,0xe},
    {0x1,0xf,0x8,0x3,0xc,0x0,0xb,0x6,0x2,0x5,0x4,0xa,0x9,0xe,0x7,0xd},
    {0xf,0x5,0x2,0xb,0x4,0xa,0x9,0xc,0x0,0x3,0xe,0x8,0xd,0x6,0x7,0x1},
    {0x7,0x2,0xc,0x5,0x8,0x4,0x6,0xb,0xe,0x9,0x1,0xf,0xd,0x3,0xa,0x0},
    {0x1,0xd,0xf,0x0,0xe,0x8,0x2,0xb,0x7,0x4,0xc,0xa,0x9,0x3,0x5,0x6},
};

std::vector<uint64_t> Inverse(const std::vector<uint64_t> & s) {
  std::vector<uint64_t > result(s.size());
  for (size_t i = 0; i < s.size(); ++i) {
    result[s[i]] = i;
  }
  return result;
}

std::vector<VectorialBooleanFunction> vbf_serpent_sb;
std::vector<VectorialBooleanFunction> vbf_serpent_inv;
void init() {
  for (const auto &sb: serpent_sb) {
    vbf_serpent_sb.emplace_back(sb, 4, 4);
    vbf_serpent_inv.emplace_back(std::move(Inverse(sb)), 4, 4);
  }
}

void filter(std::vector<ExhaustiveBF::BinaryDecisionTreePointer> & vec) {
  for (auto ptr = vec.begin(); ptr != vec.end();) {
    bool filter = false;
    filter |= SboxTools::Dom(*ptr, 4).Dimension() >= 4;
    if (filter) {
      ptr = vec.erase(ptr);
    } else {
      ptr++;
    }
  }
}

std::vector<uint64_t> prepare_f_i(std::vector<uint64_t> sbox,
                                  std::vector<uint64_t> sbox_inv,
                                  uint64_t delta) {
  std::vector<uint64_t> result(16);
  for (uint64_t i = 0; i < 16;++i) {
    result[i] = sbox_inv[sbox[i]^delta] ;
  }
  return result;
}

std::vector<uint64_t> prepare_f_i_plus_i(std::vector<uint64_t> sbox,
                                         std::vector<uint64_t> sbox_inv,
                                         uint64_t delta) {
  std::vector<uint64_t> result = prepare_f_i(sbox, sbox_inv, delta);
  for (uint64_t i = 0; i < 16;++i) {
    result[i] = result[i]^i ;
  }
  return result;
}



void analyse_sbox_components(const VectorialBooleanFunction &f,
                                    std::vector<uint64_t> components, std::string filename,
                                    bool domtriv = false) {
  ComponentChoice cfun(components, f);
  std::stringstream ss;
  bool first = true;
  for (auto c: components) {
    if (! first) ss << ".";
    first = false;
    ss << std::dec << std::to_string(c);
  }
  double bound = f.InputSize();
  auto opt_tree = BruteForceOptimization::StartSearch(cfun, bound, false);
  auto gs = ExhaustiveBF::StartSearch(cfun, bound, false);
  if (gs.empty()) {
    std::cout << "Only trivial trees for " << ss.str() << std::endl;
    auto ls = SboxTools::LinearStructures(cfun, 4).first;
    std::cout << "# linear structures: " << ls.size() << std::endl;
    for (auto &x: ls) {
      std::cout << x << "," ;
    }
    std::cout << std::endl;
    return;
  }
  const auto dom = SboxTools::Dom(gs.back(), 4);
  if (!domtriv)
    filter(gs);
  if (gs.empty()) {
    std::cout << "Only trivial trees for " << ss.str() << std::endl;
    auto ls = SboxTools::LinearStructures(cfun, 4).first;
    std::cout << "# linear structures: " << ls.size() << std::endl;
    for (auto &x: ls) {
      std::cout << x << "," ;
    }
    std::cout << std::endl;
    return;
  }
  std::ofstream outfile(filename +  ".component" + ss.str());
  for (auto &t: gs) {
    outfile << *t << "\n";
    outfile << "Domain: \n";
    outfile << dom << "\n";
  }
}
void analyse_sbox_component(const VectorialBooleanFunction &f,
                            uint64_t component, std::string filename,
                            bool domtriv = false) {
  analyse_sbox_components(f, std::vector<uint64_t> {component}, filename, domtriv);
}

void analyse_Si_inv(int i, bool domtriv=false) {
  std::string filename ="serpent_inv" + std::to_string(i);
  std::cout << "sbox_inv " << i << std::endl;
  for (size_t c = 1; c < 16; ++c) {
    analyse_sbox_component(vbf_serpent_inv[i], c, filename, domtriv);
  }
}
void analyse_fi(int i, int delta, bool domtriv=false) {
  std::string filename ="serpent_fi" + std::to_string(i) + ".delta" + std::to_string(delta);
  std::cout << "sbox_fi " << i << std::endl;
  auto vals = prepare_f_i(vbf_serpent_sb[i].GetValues(),
                         vbf_serpent_inv[i].GetValues(),
                         delta);
  VectorialBooleanFunction fun(vals, 4, 4);
  for (size_t c = 1; c < 16; ++c) {
    analyse_sbox_component(fun, c, filename, domtriv);
  }
}
void analyse_fi_plus_x(int i, int delta, bool domtriv = false) {
  std::string filename ="serpent_fi_plus_x" + std::to_string(i) + ".delta" + std::to_string(delta);
  std::cout << "sbox_fi " << i << std::endl;
  auto vals = prepare_f_i_plus_i(vbf_serpent_sb[i].GetValues(),
                         vbf_serpent_inv[i].GetValues(),
                         delta);
  VectorialBooleanFunction fun(vals, 4, 4);
  for (size_t c = 1; c < 16; ++c) {
    analyse_sbox_component(fun, c, filename, domtriv);
  }
}
void analyse_fi_full(int i, int delta, bool domtriv=false) {
  std::string filename ="serpent_fi" + std::to_string(i) + ".delta" + std::to_string(delta);
  std::cout << "sbox_fi " << i << std::endl;
  auto vals = prepare_f_i(vbf_serpent_sb[i].GetValues(),
                          vbf_serpent_inv[i].GetValues(),
                          delta);
  VectorialBooleanFunction fun(vals, 4, 4);
  analyse_sbox_components(fun, std::vector<uint64_t>{1, 2, 4, 8}, filename, domtriv);
}
void analyse_fi_plus_x_full(int i, int delta, bool domtriv = false) {
  std::string filename ="serpent_fi_plus_x" + std::to_string(i) + ".delta" + std::to_string(delta);
  std::cout << "sbox_fi " << i << std::endl;
  auto vals = prepare_f_i_plus_i(vbf_serpent_sb[i].GetValues(),
                                 vbf_serpent_inv[i].GetValues(),
                                 delta);
  VectorialBooleanFunction fun(vals, 4, 4);
  analyse_sbox_components(fun, std::vector<uint64_t> {1, 2, 4, 8}, filename, domtriv);
}

void analyse_Si_straight(int i, bool domtriv = false) {
  std::string filename ="serpent_straight" + std::to_string(i);
  std::cout << "sbox_straight " << i << std::endl;
  for (size_t c = 1; c < 16; ++c) {
    analyse_sbox_component(vbf_serpent_sb[i], c, filename, domtriv);
  }
}
void analyse_Si_straight_comps(int i, uint64_t components, bool domtriv = false) {
  std::string filename ="serpent_straight" + std::to_string(i);
  std::cout << "sbox_straight" << i << std::endl;
  std::vector<uint64_t> cs;
  for (size_t j = 0; j < 4; ++j) {
    if ((components >> j) & 1)
      cs.push_back(1 << j);
  }
  analyse_sbox_components(vbf_serpent_sb[i], cs, filename, domtriv);
}


int main() {
  init();
  std::cout << "For g" << std::endl;
  ComponentFunction cfun(vbf_serpent_sb[0], 8);
  auto g = BruteForceOptimization::StartSearch(cfun, false);
  std::cout << *g << std::endl;

#pragma omp parallel for default(none) collapse(2)
  for (size_t i = 0; i < 2; ++i) {
    for (uint64_t delta = 1; delta < 16; ++delta) {
      analyse_fi_full(i, delta, true);
      analyse_fi_plus_x_full(i, delta, true);
    }
  }
  for (size_t i = 0; i < vbf_serpent_sb.size(); ++i) {
    std::ofstream file("ddt" + std::to_string(i));
    vbf_serpent_sb[i].calculate_ddt();
    auto&ddt = vbf_serpent_sb[i].ddt;
    for (auto delta = 0; delta < ddt.size(); ++delta) {
      for (auto nabla = 0; nabla < ddt.size(); ++nabla) {
        file << ddt[delta][nabla] << ",";
      }
      file << std::endl;
    }
  }

#pragma omp parallel for default(none)
  for (size_t i = 2; i < 4; ++i) {
    analyse_Si_inv(i);
  }

#pragma omp parallel for default(none) shared(std::cout)
  for (size_t i = 0; i < 2; ++i) {
    analyse_Si_straight(i, true);
  }

  #pragma omp parallel for default(none) collapse(2)
  for (size_t i = 0; i < 2; ++i) {
    for (size_t delta = 1; delta < 16; ++delta) {
      analyse_fi(i, delta, true);
    }
  }
#pragma omp parallel for default(none) collapse(2)
  for (size_t i = 0; i < 2; ++i) {
    for (size_t delta = 1; delta < 16; ++delta) {
      analyse_fi_plus_x(i, delta, true);
    }
  }

std::vector<uint64_t> output_comps {5, 6, 7, 9, 0xA, 0xC, 0xB};
#pragma omp parallel for default(none) shared(output_comps)
  for (size_t i = 0; i < output_comps.size(); ++i) {
    analyse_Si_straight_comps(0, output_comps[i], true);
  }
  return 0;
}
