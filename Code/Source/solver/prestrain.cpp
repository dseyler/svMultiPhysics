// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "prestrain.h"

#include "VtkData.h"
#include "all_fun.h"
#include "consts.h"
#include "vtk_xml.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "mpi.h"

namespace prestrain {

/// @brief Name of the VTU point array holding the prestrain displacement.
static const std::string displacement_array = "Prestrain_displacement";

void read(const std::string& file_name, ComMod& com_mod, mshType& mesh)
{
  const int nsd = com_mod.nsd;

  {
    std::unique_ptr<VtkData> vtk_data(VtkData::create_reader(file_name));
    if (vtk_data == nullptr) {
      throw std::runtime_error("Failed to read the prestrain file '" + file_name + "'.");
    }
    if (!vtk_data->has_point_data(displacement_array)) {
      throw std::runtime_error("No point array named '" + displacement_array + "' in the prestrain file '" +
          file_name + "' for the mesh named '" + mesh.name + "'.");
    }
  }

  // Into mesh.x, in the mesh's original node order and the file's units.
  mesh.x.resize(nsd, mesh.gnNo);
  vtk_xml::read_vtu_pdata(file_name, displacement_array, nsd, nsd, 0, mesh);

  auto& U = com_mod.prestrainU;
  if (U.size() == 0) {
    U.resize(nsd, com_mod.gtnNo);
    U = 0.0;
  }
  for (int a = 0; a < mesh.gnNo; a++) {
    const int Ac = mesh.gN[a];
    for (int i = 0; i < nsd; i++) {
      U(i,Ac) = mesh.x(i,a) * mesh.scF;
    }
  }
  mesh.x.clear();
}

void init(ComMod& com_mod)
{
  using namespace consts;
  const int nsd = com_mod.nsd;

  if (com_mod.prestrainEq) {
    if (com_mod.pstEq) {
      throw std::runtime_error("Prestress and Prestrain can not both be set.");
    }
    if (com_mod.stFileFlag) {
      throw std::runtime_error("A Prestrain run can not be restarted from a .bin file. "
          "Point Prestrain_file_path at its last VTU file instead.");
    }

    bool solid = false;
    for (auto& eq : com_mod.eq) {
      if (eq.phys == EquationType::phys_struct || eq.phys == EquationType::phys_ustruct) {
        solid = true;
      }
    }
    if (!solid) {
      throw std::runtime_error("Prestrain requires a struct or ustruct equation.");
    }
  }

  // A prestrain run that starts from nothing starts from the identity.
  auto& U = com_mod.prestrainU;
  if (com_mod.prestrainEq && U.size() == 0) {
    U.resize(nsd, com_mod.tnNo);
    U = 0.0;
  }
}

/// @brief The solid equation of a prestrain run, or nullptr.
static const eqType* solid_equation(const ComMod& com_mod)
{
  using namespace consts;
  for (auto& eq : com_mod.eq) {
    if (eq.phys == EquationType::phys_struct || eq.phys == EquationType::phys_ustruct) {
      return &eq;
    }
  }
  return nullptr;
}

void adapt_time_step(ComMod& com_mod)
{
  auto& ser = com_mod.prestrainDt;
  if (!com_mod.prestrainEq || !ser.adaptive) {
    return;
  }
  const eqType* eq = solid_equation(com_mod);
  if (eq == nullptr || eq->itr == 0) {
    return;                       // no step has run yet
  }

  double& dt = com_mod.dt;
  if (ser.dt_max <= 0.0) {
    ser.dt_max = 100.0 * dt;
  }

  // The finished step: its first Newton residual (pNorm is relative to iNorm,
  // the first residual of the run) and whether it converged rather than
  // running out of iterations.
  //
  // With Max_iterations 1 each step is a single linear solve from rest,
  // which is pseudo-transient continuation (Kelley & Keyes 1998): the
  // step is always accepted, and a residual that grows shrinks the next
  // time step through the same ratio that grows it when the residual falls.
  const double R_last = eq->pNorm * eq->iNorm;
  const double r_final = eq->FSILS.RI.iNorm / eq->iNorm;
  const bool single_solve = (eq->maxItr == 1);
  ser.converged = single_solve || (eq->itr < eq->maxItr) || (r_final <= eq->tol) || (r_final <= eq->tol * eq->pNorm);

  if (!ser.converged) {
    dt = 0.5 * dt;
  } else if (!(R_last > 0.0) || !std::isfinite(R_last)) {
    // The residual has reached round-off: the state is converged and the
    // ratio is meaningless, so leave dt where it is.
  } else if (ser.residual > 0.0) {
    dt = std::min(ser.dt_max, dt * ser.residual / R_last);
    ser.residual = R_last;
  } else {
    ser.residual = R_last;        // first completed step: nothing to compare with yet
  }
}

void accumulate(ComMod& com_mod, const Array<double>& Dn)
{
  using namespace consts;
  const int nsd = com_mod.nsd;

  // A step that did not converge is repeated at a smaller time step.
  if (!com_mod.prestrainDt.converged) {
    return;
  }

  // The displacement dofs of the solid equation and its alpha_f.
  int s = -1;
  double af = 1.0;
  for (auto& eq : com_mod.eq) {
    if (eq.phys == EquationType::phys_struct || eq.phys == EquationType::phys_ustruct) {
      s = eq.s;
      af = eq.af;
      break;
    }
  }
  if (s < 0) {
    return;
  }

  auto& U = com_mod.prestrainU;
  for (int a = 0; a < com_mod.tnNo; a++) {
    for (int i = 0; i < nsd; i++) {
      U(i,a) += af * Dn(s+i,a);
    }
  }
}

void global_displacement(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, Array<double>& gU)
{
  const int nsd = com_mod.nsd;
  const auto& U = com_mod.prestrainU;

  // On the mesh's own nodes, in the mesh file's units like Displacement.
  Array<double> uM(nsd, lM.nNo);
  for (int a = 0; a < lM.nNo; a++) {
    const int Ac = lM.gN(a);
    for (int i = 0; i < nsd; i++) {
      uM(i,a) = U(i,Ac) / lM.scF;
    }
  }
  gU = all_fun::global(com_mod, cm_mod, lM, uM);
}

void print_displacement_norm(const ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Dn)
{
  using namespace consts;
  const auto& cm = com_mod.cm;
  const int nsd = com_mod.nsd;

  // The displacement dofs of the solid equation.
  int s = -1;
  for (auto& eq : com_mod.eq) {
    if (eq.phys == EquationType::phys_struct || eq.phys == EquationType::phys_ustruct) {
      s = eq.s;
      break;
    }
  }
  if (s < 0) {
    return;
  }

  double dmax = 0.0;
  for (int a = 0; a < com_mod.tnNo; a++) {
    double d2 = 0.0;
    for (int i = 0; i < nsd; i++) {
      d2 += Dn(s+i,a) * Dn(s+i,a);
    }
    // Let a NaN through rather than have max() hide it.
    dmax = (d2 > dmax || std::isnan(d2)) ? d2 : dmax;
  }
  dmax = std::sqrt(cm.reduce(cm_mod, dmax, MPI_MAX));

  if (cm.mas(cm_mod)) {
    std::cout << " Prestrain: max nodal displacement = " << dmax;
    if (const eqType* eq = solid_equation(com_mod)) {
      std::cout << "  (dt = " << com_mod.dt << ", first Newton residual = " << eq->pNorm * eq->iNorm << ")";
    }
    std::cout << std::endl;
  }
}

};
