// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PRESTRAIN_H
#define PRESTRAIN_H

#include "Array.h"
#include "CmMod.h"
#include "ComMod.h"
#include "mat_fun.h"

#include <cmath>
#include <string>

/// @brief Prestrain of a body imaged in a loaded configuration.
///
/// The mesh is the imaged, loaded geometry. The prestrain is a deformation
/// gradient F0 at every Gauss point that maps an unknown stress-free
/// configuration to the mesh. It need not be the gradient of any displacement
/// field: the stress-free configuration may be incompatible, as with residual
/// stresses. A displacement u from the mesh then gives the deformation
/// gradient from the stress-free configuration
///
///     F = (I + Grad u) F0,
///
/// with Grad taken on the mesh. The solid kernels integrate on the stress-free
/// configuration: where the assembly loops compute the shape function
/// gradients, they pull them back in place, F0^T Grad N, together with the
/// volume, dV / det F0, and the fibers, which are given on the mesh and pulled
/// back as material lines with F0^{-1}. The kernels then form F = F0 + Grad_0 u
/// and keep their total Lagrangian form. The formulation is objective:
/// rotating the deformed body rigidly rotates F and the internal forces with it.
///
/// Like the shape function gradients, F0 is constant within an element whose
/// shape functions are linear, so those elements are pulled back once. A
/// prestrain run produces such fields, and read() checks fields from file.
///
/// F0 is found by pseudo-transient continuation (see
/// PseudoTransientContinuation.h). Every pseudo-step starts from rest on the
/// imaged geometry under the full load; the displacement u it reaches is
/// accumulated as F0 <- (I + Grad u) F0 and discarded. At the fixed point the
/// mesh carries the load in Cauchy stress equilibrium without moving. This is
/// the incremental scheme of Gee, Forster and Wall, Int. J. Numer. Meth.
/// Biomed. Engng. 26 (2010) 52-72, section 3, with the increments composed on
/// the imaged mesh, so that equilibrium holds on the imaged geometry rather
/// than on a virtually deformed one.
///
/// F0 is stored per mesh as Array<double>(nsd*nsd*nG, nEl): one column per
/// element holding the nG Gauss point tensors back to back, each in the
/// column-major order of Matrix<nsd>. Written to and read from VTU cell arrays
/// named Prestrain_F_g<g>, one per Gauss point.
///
namespace prestrain {

/// @brief Deformation gradient from the stress-free configuration to the mesh
/// at Gauss point g of element e, or the identity when the mesh carries none.
template <int nsd>
mat_fun::Matrix<nsd> deformation_gradient(const mshType& lM, const int e, const int g)
{
  if (lM.F0.size() == 0) {
    return mat_fun::Matrix<nsd>::Identity();
  }
  return Eigen::Map<const mat_fun::Matrix<nsd>>(lM.F0.data() + e*lM.F0.nrows() + g*nsd*nsd);
}

/// @brief Element average of the deformation gradient from the stress-free
/// configuration, or the identity when the mesh carries none. For output whose
/// quadrature differs from the mesh's.
template <int nsd>
mat_fun::Matrix<nsd> mean_deformation_gradient(const mshType& lM, const int e)
{
  mat_fun::Matrix<nsd> F0 = mat_fun::Matrix<nsd>::Zero();
  for (int g = 0; g < lM.nG; g++) {
    F0 += deformation_gradient<nsd>(lM, e, g);
  }
  return F0 / lM.nG;
}

/// @brief The same, for callers whose number of spatial dimensions is not
/// known at compile time.
Array<double> deformation_gradient(const mshType& lM, int e, int g, int nsd);

/// @brief Pull the shape function gradients Nx of Gauss point g of element e
/// back to the stress-free configuration, in place: Nx <- F0^T Nx. Returns
/// det F0, the ratio of mesh to stress-free volume, by which the integration
/// Jacobian is divided.
double pull_back(const mshType& lM, int e, int g, Array<double>& Nx);

/// @brief Element e's fibers, read from the mesh and pulled back to the
/// stress-free configuration that F0 maps to the mesh, written to fN. Each
/// direction is a material line, pulled back with F0^{-1} and rescaled to its
/// given length; a second direction given perpendicular to the first stays
/// perpendicular to it, within the pulled-back plane of the two. Reading from
/// the mesh each time, repeated calls do not compound. Does nothing when the
/// mesh has no fibers.
template <int nsd>
void pull_back_fibers(const mat_fun::Matrix<nsd>& F0, const mshType& lM, const int e, Array<double>& fN)
{
  if (lM.fN.size() == 0) {
    return;
  }
  using Direction = Eigen::Matrix<double, nsd, 1>;
  const mat_fun::Matrix<nsd> F0_inv = F0.inverse();
  const Eigen::Map<const Eigen::Matrix<double, nsd, Eigen::Dynamic>> f(lM.fN.data() + e*lM.fN.nrows(), nsd, lM.fN.nrows() / nsd);
  auto f0 = mat_fun::eigen_view_mutable(fN);

  for (int k = 0; k < f.cols(); k++) {
    const Direction d = F0_inv * f.col(k);
    const double length = d.norm();
    f0.col(k) = (length > 0.0) ? Direction(d * (f.col(k).norm() / length)) : Direction::Zero();
  }

  if (f.cols() >= 2) {
    const double lf = f.col(0).norm();
    const double ls = f.col(1).norm();
    constexpr double perpendicular_tol = 1.0e-6;
    if (lf > 0.0 && ls > 0.0 && std::abs(f.col(0).dot(f.col(1))) <= perpendicular_tol * lf * ls) {
      const Direction a = f0.col(0) / lf;
      const Direction s = f0.col(1) - a.dot(f0.col(1)) * a;
      f0.col(1) = s * (ls / s.norm());
    }
  }
}

/// @brief The same at Gauss point g of element e, for callers whose number of
/// spatial dimensions is not known at compile time.
void pull_back_fibers(const mshType& lM, int e, int g, Array<double>& fN);

/// @brief Make a solid kernel's tangent that of a prestrain step.
///
/// A prestrain step composes the displacement u it reaches into F0 and starts
/// again from rest, so the residual the next step sees is that of the Cauchy
/// stress on the mesh, which does not move. The kernel's tangent also carries
/// the change of the current volume and of the current shape function gradients
/// with u; this subtracts that part,
///
///     w (tau Grad_x N_a (x) Grad_x N_b - tau Grad_x N_b (x) Grad_x N_a),
///
/// with tau = P F^T the Kirchhoff stress, Grad_x N = F^-T Grad_0 N, and w the
/// stress-free weight. Without it each step is an inexact Newton step whose
/// error grows with the stress, and the prestrain can stop converging once the
/// pseudo time step is large. The fibers' dependence on F0 is not included.
///
/// @param Nx0 shape function gradients on the stress-free configuration.
/// @param w_afu the weight times the stiffness scaling of lK.
/// @param lK the element tangent, lK(i*dof + j, a, b).
template <int nsd, class Gradients>
void correct_step_tangent(const mat_fun::Matrix<nsd>& F, const mat_fun::Matrix<nsd>& P, const Gradients& Nx0,
                          const double w_afu, const int dof, Array3<double>& lK)
{
  const mat_fun::NodalMatrix<nsd> Nxs = F.inverse().transpose() * Nx0;
  const mat_fun::NodalMatrix<nsd> tau_Nxs = (P * F.transpose()) * Nxs;
  for (int b = 0; b < Nxs.cols(); b++) {
    for (int a = 0; a < Nxs.cols(); a++) {
      for (int i = 0; i < nsd; i++) {
        for (int j = 0; j < nsd; j++) {
          lK(i*dof + j, a, b) -= w_afu * (tau_Nxs(i,a) * Nxs(j,b) - tau_Nxs(i,b) * Nxs(j,a));
        }
      }
    }
  }
}

/// @brief Read the prestrain of a mesh from the Prestrain_F_g<g> cell arrays
/// of a VTU file, in the mesh's original element order. Called on the master
/// before the mesh is partitioned.
void read(const std::string& file_name, mshType& mesh, int nsd);

/// @brief Validate the settings and the prestrains read from file, and seed
/// the identity where a prestrain run has none. Called once the meshes are
/// partitioned.
void init(ComMod& com_mod);

/// @brief Accumulate a pseudo-transient step into the prestrain,
/// F0 <- (I + Grad u) F0 at every Gauss point. Registered with the
/// Integrator's pseudo-transient continuation, which passes the step's
/// displacement u as the full tDof x tnNo array; only the solid equation's
/// displacement rows are read.
void accumulate(ComMod& com_mod, const Array<double>& update);

/// @brief Gather a mesh's prestrain onto the master, in the mesh's original
/// element order, for writing. All ranks must call it.
void gather(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, Array<double>& gF0);

};

#endif
