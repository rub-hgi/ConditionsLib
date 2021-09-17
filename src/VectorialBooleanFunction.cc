#include <ostream>
#include "ConditionsLib/VectorialBooleanFunction.hpp"
VectorialBooleanFunction::VectorialBooleanFunction(const std::initializer_list<
    uint64_t> &x,
                                                   size_t input_size,
                                                   size_t output_size) :
    values(x), input_size(input_size),
    output_size(output_size) {
}
size_t VectorialBooleanFunction::InputSize() const {
  return input_size;
}
size_t VectorialBooleanFunction::OutputSize() const {
  return output_size;
}
size_t& VectorialBooleanFunction::OutputSizeMutable() {
  return output_size;
}
uint64_t VectorialBooleanFunction::operator()(uint64_t x) const {
  return values[x] & ((uint64_t(1) << OutputSize()) - 1);
}
const std::vector<uint64_t> &VectorialBooleanFunction::GetValues() const {
  return values;
}
std::vector<uint64_t> &VectorialBooleanFunction::GetValuesMutable() {
  return values;
}
void VectorialBooleanFunction::TruncateOutputSize(size_t new_output_size) {
  output_size = new_output_size;
}
std::ostream &operator<<(std::ostream &os,
                         const VectorialBooleanFunction &vec) {
  os << "Vectorial Boolean function from " << vec.InputSize()
     << " to " << vec.OutputSize() << " with values: \n";
  os << vec.GetValues()[0];
  for (size_t i = 1; i < vec.GetValues().size(); ++i) {
    os << " " << vec.GetValues()[i];
  }
  return os;
}
bool VectorialBooleanFunction::operator==(const VectorialBooleanFunction &other) const {
  if (other.InputSize() != InputSize()) return false;
  if (other.OutputSize() != OutputSize()) return false;
  for (uint64_t i = 0; i < 1u << InputSize(); ++i) {
    if (other(i) != operator()(i))
      return false;
  }
  return true;
}
