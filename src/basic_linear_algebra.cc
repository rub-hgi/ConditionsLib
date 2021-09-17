#include "ConditionsLib/basic_linear_algebra.hpp"

NTL::vec_GF2 basic_linear_algebra::ConvertToNtl(uint64_t x, size_t length) {
  NTL::vec_GF2 bit_vector;
  bit_vector.SetLength(length);
  for (size_t i = 0; i < length; ++i) {
    bit_vector[i] = (x >> i) & 1;
  }
  return bit_vector;
}

uint64_t basic_linear_algebra::ConvertToUint(const NTL::vec_GF2 &x) {
  uint64_t result = 0;
  /*uint64_t test = 0;
  for (uint64_t i = 0; i < x.length(); ++i) {
    test ^= (NTL::conv<uint64_t>(x[i]) << i);
  }*/
  for (long i = 0; i < (x.length() + NTL_BITS_PER_LONG - 1) / NTL_BITS_PER_LONG;
       ++i) {
    result ^= x.rep.elts()[i] << (NTL_BITS_PER_LONG * i);
  }
  // assert(test == result);
  return result;
}

NTL::mat_GF2
basic_linear_algebra::OrthogonalComplement(const NTL::mat_GF2 &matrix) {
  NTL::mat_GF2 complement;
  NTL::kernel(complement, NTL::transpose(matrix));
  return complement;
}

NTL::vec_GF2 basic_linear_algebra::Embed(uint64_t i,
                                         const NTL::mat_GF2 &matrix) {
  NTL::vec_GF2 sum;
  sum.SetLength(matrix.NumCols());
  for (uint64_t j = 0; j < matrix.NumRows(); ++j) {
    if ((i >> j) & 1u) {
      sum += matrix[j];
    }
  }
  return sum;
}

long basic_linear_algebra::Rre(NTL::mat_GF2 &mat) {
  long rank = NTL::gauss(mat);
  long j = 0;
  for (long i = 0; i < rank; ++i) {
    // Find pivot:
    while (NTL::IsZero(mat[i][j]))
      ++j;
    // remove all other ones from this column
    // by adding the whole row to other rows appropriately.
    for (long r = 0; r < i; ++r) {
      if (!NTL::IsZero(mat[r][j]))
        mat[r] += mat[i];
    }
  }
  return rank;
}

NTL::mat_GF2 basic_linear_algebra::ComplementSpace(NTL::mat_GF2 &factor_space,
                                                   bool already_rre) {
  long j = 0, new_line = 0;
  long rank = factor_space.NumRows();
  if (!already_rre)
    rank = basic_linear_algebra::Rre(factor_space);
  NTL::mat_GF2 result;
  result.SetDims(factor_space.NumCols() - rank, factor_space.NumCols());
  for (long i = 0; i < rank; ++i) {
    // Find pivot:
    while (NTL::IsZero(factor_space[i][j])) {
      result[new_line++][j++] = 1;
    };
    j++;
  }
  // j++; // we stayed on a leading one in the last iteration.
  while (new_line < factor_space.NumCols() - rank) {
    result[new_line++][j++] = 1;
  }
  return result;
}

// Test for equality: image(A) = image(B)
void ListOfVectorSpacesRec(std::vector<std::vector<NTL::mat_GF2>> &spaces_bins,
                           NTL::mat_GF2 &accumulator, long dim,
                           long max_subspace_dim, long row, uint64_t min) {
  // Adding to appropriate bin if necessary
  // and leaving the routine if the current space is redundant.
  auto copy(accumulator);
  long rank = basic_linear_algebra::Rre(copy);
  if (rank == row) {
    // We assume that only Rre matrices have been written to spaces_bins
    // Comparison with == is correct then.
    auto it =
        std::find_if(spaces_bins[row].begin(), spaces_bins[row].end(),
                     [copy](NTL::mat_GF2 &other) { return other == copy; });
    if (it == spaces_bins[row].end()) {
      spaces_bins[row].push_back(copy);
    } else
      return;
  } else
    return;
  if (row == max_subspace_dim)
    return;
  // enumerate all possible continuations and apply recursively
  for (uint64_t vec = min; vec < (1u << dim); ++vec) {
    accumulator[row] = basic_linear_algebra::ConvertToNtl(vec, dim);
    ListOfVectorSpacesRec(spaces_bins, accumulator, dim, max_subspace_dim,
                          row + 1, vec + 1);
  }
  NTL::clear(accumulator[row]);
}

std::vector<std::vector<NTL::mat_GF2>>
basic_linear_algebra::ListOfVectorSpaces(long dim, long max_subspace_dim) {
  std::vector<std::vector<NTL::mat_GF2>> spaces_bins(max_subspace_dim + 1);
  NTL::mat_GF2 accumulator;
  accumulator.SetDims(max_subspace_dim, dim);
  long row = 0;
  uint64_t min = 1;
  ListOfVectorSpacesRec(spaces_bins, accumulator, dim, max_subspace_dim, row,
                        min);
  return spaces_bins;
}
bool basic_linear_algebra::RreLt(const NTL::mat_GF2 &mat,
                                 const NTL::mat_GF2 &other) {
  auto a = mat.NumRows();
  if (mat.NumRows() > other.NumRows())
    a = other.NumRows();
  for (long i = 0; i < a; ++i) {
    if (VecGf2Lt(mat[i], other[i]))
      return true;
    if (VecGf2Lt(other[i], mat[i]))
      return false;
  }
  if (mat.NumRows() != other.NumRows()) {
    return mat.NumRows() < other.NumRows();
  }
  return false;
}
bool basic_linear_algebra::VecGf2Lt(const NTL::vec_GF2 &x,
                                    const NTL::vec_GF2 &y) {
  return ConvertToUint(x) < ConvertToUint(y);
}
NTL::mat_GF2 &basic_linear_algebra::AppendMatrix(NTL::mat_GF2 &target,
                                                 const NTL::mat_GF2 &other) {
  auto olde = target.NumRows();
  assert(target.NumCols() == other.NumCols());
  target.SetDims(other.NumRows() + target.NumRows(), target.NumCols());
  for (long i = olde; i < target.NumRows(); ++i) {
    target[i] = other[i - olde];
  }
  return target;
}
