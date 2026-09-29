// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PRESTRAIN_H
#define PRESTRAIN_H

#include "Array.h"
#include "CmMod.h"
#include "ComMod.h"
#include "consts.h"

#include <string>

/// @brief Prestrain by an imprinted deformation gradient.
///
/// Implements the modified updated Lagrangian formulation (MULF) of Gee,
/// Forster and Wall, Int. J. Numer. Meth. Biomed. Engng. 26 (2010) 52-72,
/// section 3. The mesh is taken to be an imaged, already loaded configuration.
/// Its geometry never changes; instead each load step accumulates the
/// displacement it reached into a nodal field U that the geometry would have
/// to move by to reach the paper's virtually deformed configuration, and the
/// next step starts from it. The deformation gradient is
///
///     F = I + Grad(U + u)
///
/// with Grad taken on the mesh as usual. This is the paper's F_{t+1} F~ with
/// the increment's gradient re-expressed on the mesh (Eqs. 39, 40 and 44).
/// Every increment is the gradient of a nodal displacement, so the imprinted
/// deformation gradient is always compatible, and U is all that needs to be
/// stored: the element kernels see U added to the element displacement they
/// take the gradient of, and nothing else. During the prestrain phase the
/// displacement is reset to zero every step, so the solution converges to a
/// state that carries the load without deforming. A forward simulation then
/// loads U from file and accumulates displacement from the imaged geometry
/// as usual.
///
/// U is held in com_mod.prestrainU and written to and read from the VTU point
/// array Prestrain_displacement, in the mesh file's units like Displacement.
///
namespace prestrain {

/// @brief Add the prestrain displacement of element e's nodes to rows
/// s..s+nsd-1 of the element displacement dl, the rows a solid kernel takes
/// the deformation gradient from. Does nothing for other physics or when the
/// run carries no prestrain.
inline void add_to_element(const ComMod& com_mod, const mshType& lM, const int e, const int s,
                           const consts::EquationType phys, Array<double>& dl)
{
  const auto& U = com_mod.prestrainU;
  if (U.size() == 0) {
    return;
  }
  if (phys != consts::EquationType::phys_struct && phys != consts::EquationType::phys_ustruct) {
    return;
  }
  for (int a = 0; a < lM.eNoN; a++) {
    const int Ac = lM.IEN(a,e);
    for (int i = 0; i < com_mod.nsd; i++) {
      dl(s+i,a) += U(i,Ac);
    }
  }
}

/// @brief Read the prestrain displacement of a mesh from the
/// Prestrain_displacement point array of a VTU file into com_mod.prestrainU.
/// Called on the master before the mesh is partitioned.
void read(const std::string& file_name, ComMod& com_mod, mshType& mesh);

/// @brief Validate the settings and seed U = 0 for a prestrain run that
/// starts from nothing. Called once the meshes are partitioned.
void init(ComMod& com_mod);

/// @brief Start a prestrain step: add the displacement the previous step
/// reached to U.
///
/// Dn is the displacement the previous step ended with. The kernels saw the
/// generalized-alpha level displacement alpha_f*Dn (the step started from
/// rest), which is the state the Newton update moved to, so that is what is
/// added. Adding Dn itself would overshoot by 1/alpha_f.
void accumulate(ComMod& com_mod, const Array<double>& Dn);

/// @brief Set the pseudo time step of the next prestrain step from the finished
/// one (switched evolution relaxation): dt <- dt * R_prev / R_last, capped by
/// dt_max, where R is the first Newton residual of a step, the out-of-balance
/// force of the accumulated state. A step that did not converge halves dt
/// instead, and its displacement is not accumulated. Time is fictitious in a
/// prestrain run, so this requires steady loads. Called before the time step
/// counter advances; does nothing unless adaptive stepping was requested.
void adapt_time_step(ComMod& com_mod);

/// @brief Gather a mesh's prestrain displacement onto the master in the
/// mesh's original node order and the mesh file's units, for writing.
/// All ranks must call it.
void global_displacement(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, Array<double>& gU);

/// @brief Print the largest nodal displacement of the prestrain step. It tends
/// to zero as the prestrained state converges.
void print_displacement_norm(const ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Dn);

};

#endif
