// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "prestrain.h"

#include "VtkData.h"
#include "consts.h"
#include "nn.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "mpi.h"

namespace prestrain {

/// @brief Name of the VTU cell array holding Gauss point g.
static std::string array_name(const int g)
{
  return "Prestrain_F_g" + std::to_string(g);
}

Array<double> deformation_gradient(const mshType& lM, const int e, const int g, const int nsd)
{
  auto F0 = mat_fun::mat_id(nsd);
  if (lM.F0.size() == 0) {
    return F0;
  }

  // Column-major nsd x nsd block at this Gauss point, as Matrix<nsd> lays it out.
  const double* block = lM.F0.data() + e*lM.F0.nrows() + g*nsd*nsd;
  for (int j = 0; j < nsd; j++) {
    for (int i = 0; i < nsd; i++) {
      F0(i,j) = block[i + j*nsd];
    }
  }
  return F0;
}

void read(const std::string& file_name, mshType& mesh, const int nsd)
{
  std::unique_ptr<VtkData> vtk_data(VtkData::create_reader(file_name));
  if (vtk_data == nullptr) {
    throw std::runtime_error("Failed to read the prestrain file '" + file_name + "'.");
  }

  if (vtk_data->num_elems() != mesh.gnEl) {
    throw std::runtime_error("The number of elements (" + std::to_string(vtk_data->num_elems()) +
        ") in the prestrain file '" + file_name +
        "' is not equal to the number of elements (" + std::to_string(mesh.gnEl) +
        ") for the mesh named '" + mesh.name + "'.");
  }

  const int ncomp = nsd*nsd;
  mesh.F0 = Array<double>(ncomp*mesh.nG, mesh.gnEl);

  Array<double> Fg(ncomp, mesh.gnEl);
  for (int g = 0; g < mesh.nG; g++) {
    const auto name = array_name(g);
    if (!vtk_data->has_cell_data(name)) {
      throw std::runtime_error("No cell array named '" + name + "' in the prestrain file '" +
          file_name + "': the mesh named '" + mesh.name + "' integrates with " + std::to_string(mesh.nG) +
          " Gauss points per element.");
    }
    vtk_data->copy_cell_data(name, Fg);

    for (int e = 0; e < mesh.gnEl; e++) {
      for (int c = 0; c < ncomp; c++) {
        mesh.F0(g*ncomp + c, e) = Fg(c,e);
      }
    }
  }
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

  for (auto& msh : com_mod.msh) {
    const int rows = nsd*nsd*msh.nG;

    if (msh.F0.size() != 0 && msh.F0.nrows() != rows) {
      throw std::runtime_error("The prestrain of the mesh named '" + msh.name +
          "' holds " + std::to_string(msh.F0.nrows() / (nsd*nsd)) + " Gauss points per element but the mesh integrates with " +
          std::to_string(msh.nG) + ".");
    }

    if (!com_mod.prestrainEq) {
      continue;
    }

    // A prestrain run that starts from nothing starts from the identity.
    if (msh.F0.size() == 0) {
      msh.F0 = Array<double>(rows, msh.nEl);
      msh.F0 = 0.0;
      for (int e = 0; e < msh.nEl; e++) {
        for (int g = 0; g < msh.nG; g++) {
          for (int k = 0; k < nsd; k++) {
            msh.F0(g*nsd*nsd + k*(nsd+1), e) = 1.0;
          }
        }
      }
    }

  }
}

/// @brief F0 <- F0 + Grad(u) at Gauss point g of element e, in place.
template <int nsd>
static void add_gradient(Array<double>& F0, const int e, const int g,
                         const Array<double>& dl, const Array<double>& Nx)
{
  Eigen::Map<mat_fun::Matrix<nsd>> F(F0.data() + e*F0.nrows() + g*nsd*nsd);
  F += mat_fun::eigen_view_rows<nsd>(dl, 0) * mat_fun::eigen_view<nsd>(Nx).transpose();
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

  for (auto& msh : com_mod.msh) {
    if (msh.F0.size() == 0) {
      continue;
    }
    const int eNoN = msh.eNoN;
    Array<double> xl(nsd,eNoN), dl(nsd,eNoN), Nx(nsd,eNoN), ksix(nsd,nsd);
    double Jac = 0.0;

    for (int e = 0; e < msh.nEl; e++) {
      for (int a = 0; a < eNoN; a++) {
        const int Ac = msh.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          xl(i,a) = com_mod.x(i,Ac);
          dl(i,a) = af * Dn(s+i,Ac);
        }
      }

      for (int g = 0; g < msh.nG; g++) {
        // Shape function gradients on the mesh, constant within linear elements
        if (g == 0 || !msh.lShpF) {
          auto Nx_g = msh.Nx.slice(g);
          nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
        }
        if (nsd == 3) {
          add_gradient<3>(msh.F0, e, g, dl, Nx);
        } else {
          add_gradient<2>(msh.F0, e, g, dl, Nx);
        }
      }
    }
  }
}

void gather(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, Array<double>& gF0)
{
  const auto& cm = com_mod.cm;

  if (cm.seq()) {
    gF0 = lM.F0;
    return;
  }

  // Every rank holds the same row count; the columns are its own elements.
  const int rows = lM.F0.nrows();
  const int np = cm.np();

  Vector<int> sCount(np), disps(np);
  for (int i = 0; i < np; i++) {
    disps(i) = lM.eDist(i) * rows;
    sCount(i) = lM.eDist(i+1) * rows - disps(i);
  }

  Array<double> tmp;
  if (cm.mas(cm_mod)) {
    tmp.resize(rows, lM.gnEl);
  }

  MPI_Gatherv(lM.F0.data(), rows*lM.nEl, cm_mod::mpreal, tmp.data(), sCount.data(), disps.data(),
      cm_mod::mpreal, cm_mod.master, cm.com());

  // Back to the original element order, the inverse of the partitioning.
  if (cm.mas(cm_mod)) {
    gF0.resize(rows, lM.gnEl);
    for (int e = 0; e < lM.gnEl; e++) {
      const int Ec = lM.otnIEN(e);
      gF0.set_col(e, tmp.col(Ec));
    }
  }
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
