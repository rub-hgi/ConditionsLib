/*
Copyright 2021 Marek Broll
This file is part of ConditionsLib published under the GPL version 3 or any
later version.
Please refer to the official repository for the full licence text:
https://github.com/rub-hgi/ConditionsLib
*/
#ifndef CONDITIONS_LIB_BASIC_LINEAR_ALGEBRA_HPP_
#define CONDITIONS_LIB_BASIC_LINEAR_ALGEBRA_HPP_

#include <bit>
#include <cassert>
#include <vector>

#include <NTL/mat_GF2.h>

/// Linear algebra routines for NTL::mat_GF2
namespace basic_linear_algebra {

// Passing from uint to vec_GF2 and back
/// Calculate the binary represenation of \p x and
/// Store it in a vector of length \p length.
NTL::vec_GF2 ConvertToNtl(uint64_t x, size_t length);

/// Natural number having \p x as its binary representation
uint64_t ConvertToUint(const NTL::vec_GF2 &x);

/// Use binary representation of \p x to select a linear combination of the rows
/// of \p matrix.
NTL::vec_GF2 Embed(uint64_t x, const NTL::mat_GF2 &matrix);

/// Add rows of \p other to \p target.
NTL::mat_GF2 &AppendMatrix(NTL::mat_GF2 &target, const NTL::mat_GF2 &other);

// Calculation of related spaces
/// Direct complement
NTL::mat_GF2 ComplementSpace(NTL::mat_GF2 &factor_space,
                             bool already_rre = false);
/// Orthogonal complement (or equivalently, parity check matrix).
NTL::mat_GF2 OrthogonalComplement(const NTL::mat_GF2 &matrix);

// Reduced row echelon form and related stuff
/// Bring matrix into reduced row echelon form.
/**
 * Takes a reference to a matrix and uses Gauss transformations
 * to bring it into reduced row echelon form.
 * \param mat Matrix to be braught in reduced row echelon form.
 * \return Rank of \p mat .
 */
long Rre(NTL::mat_GF2 &mat);

/// Total ordering for binary vectors.
/**
 * Checks if the natural number with binary representation
 * \p x is less than  the natural number with binary represenation \p y .
 * \return true if x <= y, where x and y are regarded as natural numbers.
 */
bool VecGf2Lt(const NTL::vec_GF2 &x, const NTL::vec_GF2 &y);

/// Row-wise lexicographic comparison of matrices.
/**
 * Checks if the \p mat is lexicographically less than \p other
 * where all lines are compared using \ref basic_linear_algebra::VecGg2Lt.
 * \return true if x <= y, where x and y are regarded as natural numbers.
 */
bool RreLt(const NTL::mat_GF2 &mat, const NTL::mat_GF2 &other);

/// List of subspaces of GF(2)^n
/** List subspaces of GF(2)^(\p dim) of dimension at most \p max_subspace_dim
 *
 * \param dim Dimension of ambient space.
 * \param max_subspace_dim  Maximal dimension of subspaces generated.
 */
std::vector<std::vector<NTL::mat_GF2>>
ListOfVectorSpaces(long dim, long max_subspace_dim);

template <class VecType>
void walsh_hadamard_inplace(VecType &truth_table, size_t bitlength) {
  size_t window_size = ((size_t)(1u) << bitlength);
  size_t number_of_windows = 1;
  for (size_t i = 0; i < bitlength; ++i) {
    window_size >>= 1u;
    for (size_t j = 0; j < number_of_windows; ++j) {
      for (size_t k = 0; k < window_size; ++k) {
        truth_table[(2 * j) * window_size + k] +=
            truth_table[(2 * j + 1) * window_size + k];
        truth_table[(2 * j + 1) * window_size + k] *= -2;
        truth_table[(2 * j + 1) * window_size + k] +=
            truth_table[(2 * j) * window_size + k];
      }
    }
    number_of_windows <<= 1u;
  }
}

}; // namespace basic_linear_algebra

#endif
