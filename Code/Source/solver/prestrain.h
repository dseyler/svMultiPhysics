// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PRESTRAIN_H
#define PRESTRAIN_H

#include "Array.h"
#include "CmMod.h"
#include "ComMod.h"
#include "mat_fun.h"

#include <string>

/// @brief Prestrain by an imprinted deformation gradient.
///
/// Implements the modified updated Lagrangian formulation (MULF) of Gee,
/// Forster and Wall, Int. J. Numer. Meth. Biomed. Engng. 26 (2010) 52-72,
/// section 3. The mesh is taken to be an imaged, already loaded configuration.
/// Its geometry never changes; instead each converged load step imprints the
/// deformation gradient it reached, F0, into every Gauss point, and the next
/// step starts from it:
///
///     F = F0 + Grad(u)
///
/// with Grad taken on the mesh as usual. This is the paper's F_{t+1} F~ with the
/// increment's gradient re-expressed on the mesh (Eqs. 39, 40 and 44), so the
/// element kernels change only in what they add Grad(u) to. During the
/// prestress phase the displacement is reset to zero every step, so the
/// solution converges to a state that carries the load without deforming.
/// A forward simulation then loads F0 from file and accumulates displacement
/// from the imaged geometry as usual.
///
/// F0 is stored per mesh as Array<double>(nsd*nsd*nG, nEl): one column per
/// element holding the nG Gauss point tensors back to back, each in the
/// column-major order of Matrix<nsd>. Written to and read from VTU cell
/// arrays named Prestrain_F_g<g>, one per Gauss point.
///
namespace prestrain {

/// @brief Imprinted deformation gradient at Gauss point g of element e, or the
/// identity when the mesh carries none.
template <int nsd>
mat_fun::Matrix<nsd> deformation_gradient(const mshType& lM, const int e, const int g)
{
  if (lM.F0.size() == 0) {
    return mat_fun::Matrix<nsd>::Identity();
  }
  return Eigen::Map<const mat_fun::Matrix<nsd>>(lM.F0.data() + e*lM.F0.nrows() + g*nsd*nsd);
}

/// @brief The same, for callers whose number of spatial dimensions is not
/// known at compile time.
Array<double> deformation_gradient(const mshType& lM, const int e, const int g, const int nsd);

/// @brief Read the imprinted deformation gradient of a mesh from the
/// Prestrain_F_g<g> cell arrays of a VTU file, in the mesh's original element
/// order. Called on the master before the mesh is partitioned.
void read(const std::string& file_name, mshType& mesh, const int nsd);

/// @brief Validate the imprinted deformation gradients read from file and seed
/// the identity where a prestrain run has none. Called once the meshes are
/// partitioned.
void init(ComMod& com_mod);

/// @brief Start a prestrain step: accumulate the deformation gradient the
/// previous step reached, F0 <- F0 + Grad(u), at every Gauss point.
///
/// u is the displacement the kernels saw, the generalized-alpha level
/// displacement Dg, which is the state the residual was driven to zero at.
/// Accumulating the end-of-step displacement instead would overshoot by 1/alpha_f.
void accumulate(ComMod& com_mod, const Array<double>& Dg);

/// @brief Set the pseudo time step of the next prestrain step from the finished
/// one (switched evolution relaxation): dt <- dt * R_prev / R_last, capped by
/// dt_max, where R is the first Newton residual of a step, the out-of-balance
/// force of the accumulated state. A step that did not converge halves dt
/// instead, and its displacement is not accumulated. Time is fictitious in a
/// prestrain run, so this requires steady loads. Called before the time step
/// counter advances; does nothing unless adaptive stepping was requested.
void adapt_time_step(ComMod& com_mod);

/// @brief Gather a mesh's imprinted deformation gradients onto the master, in
/// the mesh's original element order, for writing. All ranks must call it.
void gather(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, Array<double>& gF0);

/// @brief Print the largest nodal displacement of the prestrain step. It tends
/// to zero as the imprinted state converges.
void print_displacement_norm(const ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Dn);

};

#endif
