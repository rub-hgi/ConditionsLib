#include "ConditionsLib/SboxTools.hpp"
#include "ConditionsLib/FourierTable.hpp"
#include <fstream>
#include <vector>

std::vector<uint64_t> SboxTools::ReadFunction(std::istream &stream,
                                              size_t input_size,
                                              size_t output_size) {
  std::string current_value;
  uint64_t current_int;
  std::vector<uint64_t> values(1u << input_size);
  uint64_t i = 0;
  while (std::getline(stream, current_value, ',') && i < (1u << input_size)) {
    current_int = std::stoull(current_value, nullptr, 16);
    values[i++] = current_int % (1u << output_size);
  }
  return values;
}

std::vector<uint64_t> SboxTools::ReadFunctionNodelim(std::istream &stream,
                                                     size_t input_size,
                                                     size_t output_size) {
  char current_value;
  uint64_t current_int;
  std::vector<uint64_t> values(1u << input_size);
  uint64_t i = 0;
  while (i < 1u << input_size && stream.get(current_value) &&
         i < (1u << input_size)) {
    if (!isdigit(current_value))
      break; // ignore non-decimal digits
    // sufficient here
    current_int = std::stoull(std::to_string(current_value));
    values[i++] = current_int % (1u << output_size);
  }
  return values;
}

std::vector<std::vector<uint64_t>> SboxTools::ReadFile(std::ifstream &stream,
                                                       size_t input_size,
                                                       size_t output_size) {
  std::string line;
  std::vector<std::vector<uint64_t>> result;
  while (std::getline(stream, line)) {
    std::istringstream iss(line);
    result.push_back(ReadFunction(iss, input_size, output_size));
  }
  return result;
}
std::vector<VectorialBooleanFunction>
SboxTools::ReadFileNodelim(std::ifstream &stream, size_t input_size,
                           size_t output_size, ssize_t frst_line,
                           ssize_t last_line) {
  std::string line;
  ssize_t line_number = -1;
  std::vector<VectorialBooleanFunction> result;
  while (std::getline(stream, line)) {
    line_number++;
    if (line_number < frst_line)
      continue;
    if (last_line >= 0 && line_number >= last_line)
      break;
    std::istringstream iss(line);
    result.emplace_back(ReadFunctionNodelim(iss, input_size, output_size),
                        input_size, output_size);
  }
  return result;
}
std::ostream &SboxTools::csv_line::PrintHeader(std::ostream &os) {
  os << "function"
     << ",";
  os << "uniformity"
     << ",";
  os << "linearity"
     << ",";
  os << "card_ls0"
     << ",";
  os << "card_ls1"
     << ",";
  os << "degree"
     << ",";
  os << "leaves"
     << ",";
  os << "domsize"
     << ",";
  os << std::endl;
  return os;
}
std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
SboxTools::LinearStructures(const VectorialBooleanFunction &cfun,
                            size_t input_size) {
  std::vector<uint64_t> result_0;
  std::vector<uint64_t> result_1;
  for (uint64_t u = 0; u < (1u << input_size); ++u) {
    bool is_ls = true;
    uint64_t x = 1;
    for (x = 1; x < (1u << input_size); ++x) {
      if (cfun(x) ^ cfun(x ^ u) ^ cfun(0) ^ cfun(u)) {
        is_ls = false;
        break;
      }
    }
    if (!is_ls)
      continue;
    if (cfun(0) ^ cfun(u)) {
      result_1.push_back(u);
    } else {
      result_0.push_back(u);
    }
  }
  return std::make_pair(result_0, result_1);
}
void SboxTools::AllSubspaces(BruteForce::BinaryDecisionTreePointer &x,
                             std::vector<Uint64Subspace> &list,
                             std::vector<uint64_t> &current_basis,
                             size_t dimension) {
  if (x->IsLeaf()) {
    list.emplace_back(current_basis, dimension);
    list.back() = list.back().OrthogonalComplement();
    return;
  }
  current_basis.push_back(x->value_);
  AllSubspaces(x->left_, list, current_basis, dimension);
  AllSubspaces(x->right_, list, current_basis, dimension);
  current_basis.pop_back();
}
Uint64Subspace SboxTools::Dom(const BruteForce::BinaryDecisionTreePointer &x,
                              uint64_t input_size) {
  if (x->IsLeaf()) {
    return Uint64Subspace(std::vector<uint64_t>{}, input_size);
  }
  std::vector<uint64_t> current_label{x->value_};
  return Dom(x->left_, input_size) + Dom(x->right_, input_size) +
         Uint64Subspace(current_label, input_size);
}
uint64_t SboxTools::NaiveDegree(const VectorialBooleanFunction &fun) {
  uint64_t current_coefficient;
  uint64_t result = 0;
  for (uint64_t u = 0; u < (1u << fun.InputSize()); ++u) {
    current_coefficient = 0;
    for (uint64_t x = 0; x < (1u << fun.InputSize()); ++x) {
      if (!((x & u) == x && (x | u) == u))
        continue;
      current_coefficient ^= fun(x);
    }
    if (current_coefficient != 0 && std::popcount(u) > result) {
      result = std::popcount(u);
    }
  }
  return result;
}
uint64_t SboxTools::Uniformity(const VectorialBooleanFunction &fun) {
  std::vector<uint64_t> ddt_line((1u << fun.InputSize()));
  uint64_t result = 0;
  for (auto inp = 1; inp < (1u << fun.InputSize()); ++inp) {
    std::fill(ddt_line.begin(), ddt_line.end(), 0);
    for (auto x = 0; x < (1u << fun.InputSize()); x++) {
      const auto out = fun(x ^ inp) ^ fun(x);
      ddt_line[out]++;
    }
    for (auto u : ddt_line) {
      if (u > result)
        result = u;
    }
  }
  return result;
}
Uint64Subspace SboxTools::Intersection(std::vector<Uint64Subspace> &sub) {
  auto result = sub[0];
  for (size_t i = 0; i < sub.size(); ++i) {
    result = sub[i] & result;
  }
  return result;
}

int64_t SboxTools::SignedLinearity(const VectorialBooleanFunction &fun) {
  FourierTable ftbl(fun);
  int64_t result = 0;
  for (uint64_t inp = 0; inp < (1u << fun.InputSize()); ++inp) {
    for (uint64_t out = 1; out < (1u << fun.OutputSize()); ++out) {
      if (std::abs(ftbl(inp, out)) > std::abs(result))
        result = ftbl(inp, out);
    }
  }
  return result;
}
std::ostream &SboxTools::operator<<(std::ostream &os,
                                    const SboxTools::csv_line &line) {
  os << line.function << ",";
  os << std::dec << line.uniformity << ",";
  os << line.linearity << ",";
  os << line.card_ls0 << ",";
  os << line.card_ls1 << ",";
  os << line.degree << ",";
  os << line.costs_A << ",";
  os << line.costs_D << ",";
  os << std::endl;
  return os;
}
int64_t SboxTools::Linearity(const VectorialBooleanFunction &fun) {
  return std::abs(SignedLinearity(fun));
}
