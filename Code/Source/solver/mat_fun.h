// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAT_FUN_H 
#define MAT_FUN_H 
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"
#include <stdexcept>

#include "Array.h"
#include "Tensor4.h"
#include "Vector.h"
#include "FE/Common/FEException.h"

/// @brief The classes defined here duplicate the data structures in the 
/// Fortran MATFUN module defined in MATFUN.f. 
///
/// This module defines data structures for generally performed matrix and tensor operations.
///
/// \todo [TODO:DaveP] this should just be a namespace?
//
namespace mat_fun {
    // Define templated type aliases for Eigen matrices and tensors for convenience
    template<size_t nsd>
    using Matrix = Eigen::Matrix<double, nsd, nsd>;

    template<size_t nsd>
    using Tensor = Eigen::TensorFixedSize<double, Eigen::Sizes<nsd, nsd, nsd, nsd>>;

    // Function to convert Array<double> to Eigen::Matrix
    template <typename MatrixType>
    MatrixType convert_to_eigen_matrix(const Array<double>& src) {
        MatrixType mat;
        for (int i = 0; i < mat.rows(); ++i)
            for (int j = 0; j < mat.cols(); ++j)
                mat(i, j) = src(i, j);
        return mat;
    }

    // Function to convert Eigen::Matrix to Array<double>
    template <typename MatrixType>
    void convert_to_array(const MatrixType& mat, Array<double>& dest) {
        for (int i = 0; i < mat.rows(); ++i)
            for (int j = 0; j < mat.cols(); ++j)
                dest(i, j) = mat(i, j);
    }

    // Function to convert a higher-dimensional array like Dm
    template <typename MatrixType>
    void copy_Dm(const MatrixType& mat, Array<double>& dest) {
        if ((mat.rows() != dest.nrows()) || (mat.cols() != dest.ncols())) {
          const std::string mat_dims = "(" + std::to_string(mat.rows()) + "x" + std::to_string(mat.cols()) + ")";
          const std::string dest_dims = "(" + std::to_string(dest.nrows()) + "x" + std::to_string(dest.ncols()) + ")";
          svmp::raise<svmp::FE::InvalidArgumentException>(
              "The 'mat" + mat_dims + "' and 'dest" + dest_dims +
              "' arrays have incompatible sizes.");
        }

        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                dest(i, j) = mat(i, j);
            }
        }
    }

    template <int nsd>
    Eigen::Matrix<double, nsd, 1> cross_product(const Eigen::Matrix<double, nsd, 1>& u, const Eigen::Matrix<double, nsd, 1>& v) {
        if constexpr (nsd == 2) {
            return Eigen::Matrix<double, 2, 1>(v(1), - v(0));
        }
        else if constexpr (nsd == 3) {
            return u.cross(v);
        }
        else {
            throw std::runtime_error("[cross_product] Invalid number of spatial dimensions '" + std::to_string(nsd) + "'. Valid dimensions are 2 or 3.");
        }
    }

    double mat_ddot(const Array<double>& A, const Array<double>& B, const int nd);
    
    template <int nsd>
    double double_dot_product(const Matrix<nsd>& A, const Matrix<nsd>& B) {
        return A.cwiseProduct(B).sum();
    }
    
    double mat_det(const Array<double>& A, const int nd);
    Array<double> mat_dev(const Array<double>& A, const int nd);

    Array<double> mat_dyad_prod(const Vector<double>& u, const Vector<double>& v, const int nd);

    Array<double> mat_id(const int nsd);
    Array<double> mat_inv(const Array<double>& A, const int nd, bool debug = false);
    Array<double> mat_inv_ge(const Array<double>& A, const int nd, bool debug = false);
    Array<double> mat_inv_lp(const Array<double>& A, const int nd);

    /// @brief Multiply a matrix by a vector, returning A*v.
    ///
    /// @param[in] A matrix with as many columns as v has entries.
    /// @param[in] v vector.
    /// @return the product, of size rows(A).
    ///
    /// Throws InvalidArgumentException if the sizes are incompatible.
    Vector<double> mat_mul(const Array<double>& A, const Vector<double>& v);

    /// @brief Multiply two matrices, returning A*B.
    ///
    /// @param[in] A left operand.
    /// @param[in] B right operand, with as many rows as A has columns.
    /// @return the product, of size rows(A) by cols(B).
    ///
    /// Throws InvalidArgumentException if the sizes are incompatible. The
    /// result is freshly allocated, so the operands may alias it, as in
    /// A = mat_mul(A, B).
    Array<double> mat_mul(const Array<double>& A, const Array<double>& B);

    /// @brief Multiply two matrices, writing A*B into an existing result.
    ///
    /// @param[in]  A left operand.
    /// @param[in]  B right operand, with as many rows as A has columns.
    /// @param[out] result the product. The caller sizes it rows(A) by cols(B).
    ///
    /// Throws InvalidArgumentException if the sizes are incompatible. Preferred
    /// in loops, where it reuses the caller's storage instead of allocating a
    /// result on every call. The result must not alias A or B.
    void mat_mul(const Array<double>& A, const Array<double>& B, Array<double>& result);

    /// @brief Matrix product with the operand shape supplied at compile time.
    ///
    /// Overloads of mat_mul rather than differently named helpers, so a call
    /// site states the shape and otherwise reads exactly as before:
    ///
    /// @code
    ///   mat_mul(Dm, Bm.rslice(b), DBm);        // runtime shape check
    ///   mat_mul<6, 6, 3>(Dm, Bm.rslice(b), DBm);  // no check, same arguments
    /// @endcode
    ///
    /// The generic mat_mul overload above will dispatch to this overload
    /// when the shapes are known at compile time.
    ///
    /// @tparam M rows of A and of the result
    /// @tparam K columns of A and rows of B, the contracted dimension
    /// @tparam N columns of B and of the result
    template <int M, int K, int N>
    void mat_mul(const Array<double>& A, const Array<double>& B,
                 Array<double>& C)
    {
      Eigen::Map<const Eigen::Matrix<double, M, K>> a(A.data());
      Eigen::Map<const Eigen::Matrix<double, K, N>> b(B.data());
      Eigen::Map<Eigen::Matrix<double, M, N>>       c(C.data());

      c.noalias() = a * b;
    }

    /// @brief As above, but with the column count known only at run time.
    ///
    /// For operands with one column per element node, where the width depends on
    /// the element type. The row counts are still compile-time, which is where
    /// most of the benefit comes from.
    template <int M, int K>
    void mat_mul(const Array<double>& A, const Array<double>& B,
                 Array<double>& C)
    {
      Eigen::Map<const Eigen::Matrix<double, M, K>> a(A.data());
      Eigen::Map<const Eigen::Matrix<double, K, Eigen::Dynamic>> b(B.data(), K, B.ncols());
      Eigen::Map<Eigen::Matrix<double, M, Eigen::Dynamic>>       c(C.data(), M, C.ncols());

      c.noalias() = a * b;
    }

    Array<double> mat_symm(const Array<double>& A, const int nd);

    /**
     * @brief Symmetric part of a 2nd order tensor, 0.5 * (A + A^T).
     *
     * Fixed-size overload for the Eigen matrices used by the element kernels.
     *
     * @tparam nsd, the number of spatial dimensions
     */
    template <int nsd>
    Matrix<nsd> mat_symm(const Matrix<nsd>& A) {
        return 0.5 * (A + A.transpose());
    }

    /**
     * @brief Deviatoric part of a 2nd order tensor, A - tr(A)/nsd * I.
     *
     * Fixed-size overload for the Eigen matrices used by the element kernels.
     *
     * @tparam nsd, the number of spatial dimensions
     */
    template <int nsd>
    Matrix<nsd> mat_dev(const Matrix<nsd>& A) {
        return A - (A.trace() / static_cast<double>(nsd)) * Matrix<nsd>::Identity();
    }

    Array<double> mat_symm_prod(const Vector<double>& u, const Vector<double>& v, const int nd);

    double mat_trace(const Array<double>& A, const int nd);

    Tensor4<double> ten_asym_prod12(const Array<double>& A, const Array<double>& B, const int nd);
    Tensor4<double> ten_ddot(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);
    Tensor4<double> ten_ddot_2412(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);
    Tensor4<double> ten_ddot_3424(const Tensor4<double>& A, const Tensor4<double>& B, const int nd);

    /**
     * @brief Contracts two 4th order tensors A and B over two dimensions, 
     * 
     */
    template <int nsd>
    Tensor<nsd>
    double_dot_product(const Tensor<nsd>& A, const std::array<int, 2>& dimsA,
                        const Tensor<nsd>& B, const std::array<int, 2>& dimsB) {

        // Fast path for contraction over the two trailing dimensions,
        // C_ijmn = A_ijkl * B_mnkl, which is every call in the element loops.
        //
        // Tensor<nsd> is column major, so the linear index of (i,j,k,l) is
        // i + n*j + n^2*k + n^3*l. Viewed as an NxN column-major matrix with
        // N = nsd^2, row i + n*j and column k + n*l give that same index, so
        // the free pair indexes rows and the contracted pair indexes columns
        // with no copy or permutation. The contraction is then A * B^T.
        //
        // contract() reaches the same result, and both routes end in Eigen's
        // blocked GEMM. The difference is how the operands are read while
        // being packed: contract() goes through TensorContractionSubMapper,
        // which indirects through the tensor evaluator on every element, while
        // the matrix path uses a plain strided pointer mapper. Measured at
        // 1.8x for nsd = 3, with bitwise identical results.
        if (dimsA[0] == 2 && dimsA[1] == 3 && dimsB[0] == 2 && dimsB[1] == 3) {
            constexpr int N = nsd * nsd;
            Tensor<nsd> C;
            Eigen::Map<const Eigen::Matrix<double, N, N>> a(A.data());
            Eigen::Map<const Eigen::Matrix<double, N, N>> b(B.data());
            Eigen::Map<Eigen::Matrix<double, N, N>> c(C.data());
            c.noalias() = a * b.transpose();
            return C;
        }

        // Define the contraction dimensions
        Eigen::array<Eigen::IndexPair<int>, 2> contractionDims = {
            Eigen::IndexPair<int>(dimsA[0], dimsB[0]), // Contract A's dimsA[0] with B's dimsB[0]
            Eigen::IndexPair<int>(dimsA[1], dimsB[1])  // Contract A's dimsA[1] with B's dimsB[1]
        };

        // Return the double dot product
        return A.contract(B, contractionDims);

        // For some reason, in this case the Eigen::Tensor contract function is
        // faster than a for loop implementation.
    }

    /**
     * @brief Applies the isochoric projection PP : CC_bar : PP^T, where
     * PP = S - (1/nsd) * dyad(Ci, C) and S is the symmetric fourth order identity.
     *
     * PP differs from S by a single rank one term, and a hyperelastic tangent is
     * minor symmetric, so S : CC_bar : S is CC_bar and the projection collapses
     * to three rank one updates driven by one matrix vector product. PP is never
     * formed and the two nsd^2 x nsd^2 products disappear.
     *
     * Requires CC_bar to carry the minor and major symmetries of a hyperelastic
     * tangent. Instrumenting bar_to_iso across the unit test suite showed all
     * three holding exactly (max deviation 0.0 over 315802 calls, covering every
     * isochoric model). A model that produced an unsymmetric tangent would
     * silently get a wrong answer here.
     *
     * @tparam nsd, the number of spatial dimensions
     * @param[in] CC_bar, the fictitious elasticity tensor
     * @param[in] Ci, the inverse of the right Cauchy-Green deformation tensor
     * @param[in] C, the right Cauchy-Green deformation tensor
     * @return the isochoric elasticity tensor
     */
    template <int nsd>
    Tensor<nsd> isochoric_projection(const Tensor<nsd>& CC_bar,
                                     const Matrix<nsd>& Ci, const Matrix<nsd>& C) {
        constexpr int N = nsd * nsd;
        using VecN = Eigen::Matrix<double, N, 1>;

        // Column-major storage lets the 4th order tensor be read as an NxN
        // matrix and the 2nd order tensors as N-vectors, with no copy. See
        // double_dot_product for the index correspondence.
        Eigen::Map<const Eigen::Matrix<double, N, N>> cc(CC_bar.data());
        Eigen::Map<const VecN> ci(Ci.data());
        Eigen::Map<const VecN> c(C.data());

        const VecN w = cc * c;
        const double cw = c.dot(w);

        Tensor<nsd> CC_iso;
        Eigen::Map<Eigen::Matrix<double, N, N>> out(CC_iso.data());
        out.noalias() = cc
                      - (1.0 / nsd) * w * ci.transpose()
                      - (1.0 / nsd) * ci * w.transpose()
                      + (cw / (nsd * nsd)) * ci * ci.transpose();
        return CC_iso;
    }

    Tensor4<double> ten_dyad_prod(const Array<double>& A, const Array<double>& B, const int nd);
    
    /**
     * @brief Compute the dyadic product of two 2nd order tensors A and B, C_ijkl = A_ij * B_kl
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the first 2nd order tensor
     * @param B, the second 2nd order tensor
     * @return Tensor<nsd>
     */
    template <int nsd>
    Tensor<nsd>
    dyadic_product(const Matrix<nsd>& A, const Matrix<nsd>& B) {
        constexpr int N = nsd * nsd;

        // C_ijkl = A_ij * B_kl is an outer product in the NxN view of the
        // tensor: with column-major storage the index pair (i,j) is the row
        // and (k,l) the column, so the second order tensors read as N-vectors.
        //
        // Written as an index loop the innermost index is l, whose stride is
        // nsd^3, so every innermost write lands on a different cache line and
        // nothing vectorises. This form writes down columns instead.
        Tensor<nsd> C;
        Eigen::Map<const Eigen::Matrix<double, N, 1>> a(A.data());
        Eigen::Map<const Eigen::Matrix<double, N, 1>> b(B.data());
        Eigen::Map<Eigen::Matrix<double, N, N>> c(C.data());
        c.noalias() = a * b.transpose();
        return C;
    }

    Tensor4<double> ten_ids(const int nd);

    /**
     * @brief Create a 4th order identity tensor:
     * I_ijkl = 0.5 * (δ_ik * δ_jl + δ_il * δ_jk)
     * 
     * @tparam nsd, the number of spatial dimensions
     * @return Tensor<nsd> 
     */
    template <int nsd>
    Tensor<nsd>
    fourth_order_identity() {
        // Initialize as zero
        Tensor<nsd> I;
        I.setZero();

        // Set only non-zero entries
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                I(i,j,i,j) += 0.5;
                I(i,j,j,i) += 0.5;
            }
        }

        return I;
    }

    Array<double> ten_mddot(const Tensor4<double>& A, const Array<double>& B, const int nd);

    Tensor4<double> ten_symm_prod(const Array<double>& A, const Array<double>& B, const int nd);
    
    /// @brief Create a 4th order tensor from symmetric outer product of two matrices: C_ijkl = 0.5 * (A_ik * B_jl + A_il * B_jk)
    ///
    /// Reproduces 'FUNCTION TEN_SYMMPROD(A, B, nd) RESULT(C)'.
    //
    template <int nsd>
    Tensor<nsd>
    symmetric_dyadic_product(const Matrix<nsd>& A, const Matrix<nsd>& B) {
        
        // Initialize the result tensor
        Tensor<nsd> C;

        // C_ijkl = 0.5 * (A_ik * B_jl + A_il * B_jk).
        //
        // Same ordering point as dyadic_product: with column-major storage the
        // flat index is i + nsd*j + nsd^2*k + nsd^3*l, so iterating l innermost
        // strides by nsd^3. Running i innermost instead makes the writes
        // contiguous and the inner statement a vectorisable axpy, since only
        // A(i,k) and A(i,l) vary with i.
        double* c = C.data();
        for (int l = 0; l < nsd; ++l) {
            for (int k = 0; k < nsd; ++k) {
                const double* a_k = A.data() + nsd * k;   // A(:,k)
                const double* a_l = A.data() + nsd * l;   // A(:,l)
                for (int j = 0; j < nsd; ++j) {
                    const double b_jl = B(j, l);
                    const double b_jk = B(j, k);
                    double* out = c + nsd * j + nsd * nsd * k + nsd * nsd * nsd * l;
                    for (int i = 0; i < nsd; ++i) {
                        out[i] = 0.5 * (a_k[i] * b_jl + a_l[i] * b_jk);
                    }
                }
            }
        }

        // Return the symmetric product
        return C;
    }

    Tensor4<double> ten_transpose(const Tensor4<double>& A, const int nd);

    /**
     * @brief Performs a tensor transpose operation on a 4th order tensor A, B_ijkl = A_klij
     * 
     * @tparam nsd, the number of spatial dimensions
     * @param A, the input 4th order tensor
     * @return Tensor<nsd>
     */
    template <int nsd>
    Tensor<nsd>
    transpose(const Tensor<nsd>& A) {

        // Initialize the result tensor
        Tensor<nsd> B;

        // Permute the tensor indices to perform the transpose operation
        for (int i = 0; i < nsd; ++i) {
            for (int j = 0; j < nsd; ++j) {
                for (int k = 0; k < nsd; ++k) {
                    for (int l = 0; l < nsd; ++l) {
                        B(i,j,k,l) = A(k,l,i,j);
                    }
                }
            }
        }

        return B;
    }

    Array<double> transpose(const Array<double>& A);

    void ten_init(const int nd);

};

#endif
