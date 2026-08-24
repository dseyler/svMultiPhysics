/* Copyright (c) Stanford University, The Regents of the University of
 * California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "../test_common.h"
#include "mat_fun.h"

#include <array>
#include <cmath>
#include <random>

using namespace mat_fun;

namespace {

// Deterministic fill so a failure is reproducible.
template <int nsd>
Tensor<nsd> random_tensor(unsigned seed) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  Tensor<nsd> T;
  for (int i = 0; i < T.size(); ++i) {
    T.data()[i] = dist(gen);
  }
  return T;
}

// C_ijmn = sum_kl A_ij(dimsA) * B_mn(dimsB), written out directly. This is the
// reference the fast path must reproduce; it depends on no Eigen machinery.
template <int nsd>
Tensor<nsd> contract_by_hand_trailing(const Tensor<nsd>& A, const Tensor<nsd>& B) {
  Tensor<nsd> C;
  C.setZero();
  for (int i = 0; i < nsd; ++i) {
    for (int j = 0; j < nsd; ++j) {
      for (int m = 0; m < nsd; ++m) {
        for (int n = 0; n < nsd; ++n) {
          double sum = 0.0;
          for (int k = 0; k < nsd; ++k) {
            for (int l = 0; l < nsd; ++l) {
              sum += A(i, j, k, l) * B(m, n, k, l);
            }
          }
          C(i, j, m, n) = sum;
        }
      }
    }
  }
  return C;
}

// Eigen's contract(), i.e. the path double_dot_product uses for dimension
// pairs outside the fast path.
template <int nsd>
Tensor<nsd> contract_with_eigen(const Tensor<nsd>& A, const std::array<int, 2>& dimsA,
                                const Tensor<nsd>& B, const std::array<int, 2>& dimsB) {
  Eigen::array<Eigen::IndexPair<int>, 2> dims = {
      Eigen::IndexPair<int>(dimsA[0], dimsB[0]),
      Eigen::IndexPair<int>(dimsA[1], dimsB[1])
  };
  return A.contract(B, dims);
}

template <int nsd>
void expect_tensors_near(const Tensor<nsd>& actual, const Tensor<nsd>& expected,
                         double tol) {
  ASSERT_EQ(actual.size(), expected.size());
  for (int i = 0; i < expected.size(); ++i) {
    const double scale = std::fmax(1.0, std::fabs(expected.data()[i]));
    EXPECT_NEAR(actual.data()[i], expected.data()[i], tol * scale)
        << "  at flat index " << i;
  }
}

}  // namespace

// The reshape the fast path relies on is only valid for column-major storage.
TEST(TensorOps, TensorStorageIsColumnMajor) {
  EXPECT_EQ(static_cast<int>(Tensor<3>::Options), 0) << "Tensor<nsd> is no longer "
      "column major; the double_dot_product reshape is invalid.";
  EXPECT_EQ(static_cast<int>(Tensor<2>::Options), 0);
}

// The fast path must reproduce the contraction itself, not merely agree with
// whatever Eigen happens to do.
TEST(TensorOps, DoubleDotTrailingDimsMatchesExplicitSum3D) {
  const auto A = random_tensor<3>(12345);
  const auto B = random_tensor<3>(67890);
  expect_tensors_near<3>(double_dot_product<3>(A, {2, 3}, B, {2, 3}),
                         contract_by_hand_trailing<3>(A, B), 1e-13);
}

TEST(TensorOps, DoubleDotTrailingDimsMatchesExplicitSum2D) {
  const auto A = random_tensor<2>(2468);
  const auto B = random_tensor<2>(1357);
  expect_tensors_near<2>(double_dot_product<2>(A, {2, 3}, B, {2, 3}),
                         contract_by_hand_trailing<2>(A, B), 1e-13);
}

// The fast path replaces Eigen's contract() for this dimension pair, so the
// two must still agree.
TEST(TensorOps, DoubleDotTrailingDimsMatchesEigenContract) {
  const auto A = random_tensor<3>(11111);
  const auto B = random_tensor<3>(22222);
  expect_tensors_near<3>(double_dot_product<3>(A, {2, 3}, B, {2, 3}),
                         contract_with_eigen<3>(A, {2, 3}, B, {2, 3}), 1e-13);
}

// Dimension pairs outside the fast path must be left on the contract() route
// untouched. These are the active-strain calls in mat_models.cpp.
TEST(TensorOps, DoubleDotOtherDimPairsUnchanged) {
  const auto A = random_tensor<3>(33333);
  const auto B = random_tensor<3>(44444);

  const std::array<std::pair<std::array<int, 2>, std::array<int, 2>>, 3> cases = {{
      {{2, 3}, {1, 3}},   // mat_models.cpp: CC contracted with CC_bar
      {{1, 3}, {0, 1}},   // mat_models.cpp: CC_bar contracted with CC
      {{0, 1}, {0, 1}},
  }};

  for (const auto& [dimsA, dimsB] : cases) {
    expect_tensors_near<3>(double_dot_product<3>(A, dimsA, B, dimsB),
                           contract_with_eigen<3>(A, dimsA, B, dimsB), 1e-13);
  }
}

// The two hot call sites in bar_to_iso pass the same tensor as both operands
// on some paths, and mat_models.cpp:675-676 self-assign the result.
TEST(TensorOps, DoubleDotHandlesAliasedOperands) {
  const auto A = random_tensor<3>(55555);
  expect_tensors_near<3>(double_dot_product<3>(A, {2, 3}, A, {2, 3}),
                         contract_by_hand_trailing<3>(A, A), 1e-13);
}
