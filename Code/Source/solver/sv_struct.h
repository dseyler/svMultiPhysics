// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef STRUCT_H 
#define STRUCT_H 

#include "ComMod.h"
#include "SolutionStates.h"

namespace struct_ns {

/// @brief The element nodal data the solid kernels read, as one argument.
///
/// Holds references to arrays the caller already owns; it copies nothing and
/// must not outlive them. Callers assemble these arrays for several kernel
/// families, so ownership stays with the caller and this only groups the
/// subset the solid kernels need.
struct SolidElementInput
{
  /// @brief Number of element nodes.
  int eNoN;
  /// @brief Number of fiber directions.
  int nFn;

  /// @brief Nodal acceleration, velocity and displacement, tDof x eNoN.
  const Array<double> &al, &yl, &dl;
  /// @brief Nodal body force, nsd x eNoN.
  const Array<double> &bfl;
  /// @brief Fiber directions, nsd x nFn.
  const Array<double> &fN;
  /// @brief Nodal prestress in Voigt form, nsymd x eNoN.
  const Array<double> &pS0l;
  /// @brief Nodal active tension along fiber, sheet and normal, eNoN each.
  const Vector<double> &ya_l_f, &ya_l_s, &ya_l_n;
};

/// @brief Working storage for the solid kernels, reused for every Gauss point.
///
/// Owned by the caller and sized once per mesh, so the kernels do not allocate
/// on each call.
struct SolidScratch
{
  /// @brief Deformation gradient, prestress and velocity gradient, nsd x nsd.
  Array<double> F, S0, vx;
  /// @brief 2nd Piola-Kirchhoff stress, and material stiffness in Voigt form.
  Array<double> S, Dm;
  /// @brief Viscous stress and its tangent contributions.
  Array<double> Svis;
  Array3<double> Kvis_u, Kvis_v;
  /// @brief 1st Piola-Kirchhoff stress, strain-displacement matrix, and Dm*Bm.
  Array<double> P, DBm;
  Array3<double> Bm;
  /// @brief Body force acceleration term and body force vector, nsd.
  Vector<double> ud, fb;

  /// @brief Size every member for the given mesh dimensions.
  void resize(const int nsd, const int num_element_nodes);
};

void b_struct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void b_struct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void construct_dsolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

/// @brief Assemble the 2D solid element residual and stiffness at one Gauss point.
///
/// @param[in] element Nodal data for the element being assembled.
/// @param[in] w Gauss weight times the Jacobian.
/// @param[in] N,Nx Shape functions and their spatial derivatives at this point.
/// @param[in,out] scratch Caller-owned working storage, sized for this mesh.
/// @param[out] pSl Second Piola-Kirchhoff stress in Voigt form at this point.
/// @param[in,out] lR,lK Element residual and stiffness, accumulated into.
void struct_2d(ComMod &com_mod, CepMod &cep_mod,
               const SolidElementInput &element, const double w,
               const Vector<double> &N, const Array<double> &Nx,
               SolidScratch &scratch, Vector<double> &pSl,
               Array<double> &lR, Array3<double> &lK);

/// @brief Assemble the 3D solid element residual and stiffness at one Gauss point.
///
/// @param[in] element Nodal data for the element being assembled.
/// @param[in] w Gauss weight times the Jacobian.
/// @param[in] N,Nx Shape functions and their spatial derivatives at this point.
/// @param[in,out] scratch Caller-owned working storage, sized for this mesh.
/// @param[out] pSl Second Piola-Kirchhoff stress in Voigt form at this point.
/// @param[in,out] lR,lK Element residual and stiffness, accumulated into.
void struct_3d(ComMod &com_mod, CepMod &cep_mod,
               const SolidElementInput &element, const double w,
               const Vector<double> &N, const Array<double> &Nx,
               SolidScratch &scratch, Vector<double> &pSl,
               Array<double> &lR, Array3<double> &lK);
};

#endif

