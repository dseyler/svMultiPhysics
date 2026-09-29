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
    if (!com_mod.pseudoTransient.enabled) {
      throw std::runtime_error("Prestrain is solved by pseudo-transient continuation, which is not enabled.");
    }
  }

  // A prestrain run that starts from nothing starts from the identity.
  auto& U = com_mod.prestrainU;
  if (com_mod.prestrainEq && U.size() == 0) {
    U.resize(nsd, com_mod.tnNo);
    U = 0.0;
  }
}

void accumulate(ComMod& com_mod, const Array<double>& update)
{
  using namespace consts;
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

  auto& U = com_mod.prestrainU;
  for (int a = 0; a < com_mod.tnNo; a++) {
    for (int i = 0; i < nsd; i++) {
      U(i,a) += update(s+i,a);
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

};
