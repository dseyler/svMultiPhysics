// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef STRUCT_H 
#define STRUCT_H 

#include "ComMod.h"
#include "SolutionStates.h"

namespace struct_ns {


void b_struct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void b_struct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void construct_dsolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

void struct_2d(ComMod &com_mod, CepMod &cep_mod, const int eNoN, const int nFn,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Array<double> &fN, const Array<double> &pS0l,
               Vector<double> &pSl, const Vector<double> &ya_l_f,
               const Vector<double> &ya_l_s, const Vector<double> &ya_l_n,
               Array<double> &lR, Array3<double> &lK);

/// @brief Assemble a 3D solid element in one pass, for linear shape functions.
///
/// Valid only where the shape function derivatives are constant over the element
/// (lM.lShpF). The quantities that depend on Nx alone are then the same at every
/// quadrature point, and the remaining terms are linear in the stress, so the
/// quadrature sum is formed analytically rather than by looping.
///
/// The constitutive model is evaluated at the weighted-average activation, which
/// is exact only where the stress is affine in it: true for the active-stress
/// models, not for active strain.
///
/// @param[in] eNoN Number of element nodes.
/// @param[in] nFn Number of fiber directions.
/// @param[in] nG Number of quadrature points.
/// @param[in] wg,Ng Quadrature weights and shape functions at those points.
/// @param[in] Jac Element Jacobian, constant for a linear element.
/// @param[in] Nx Shape function spatial derivatives.
/// @param[in] al,yl,dl,bfl Nodal acceleration, velocity, displacement, body force.
/// @param[in] fN Fiber directions.
/// @param[in] pS0l Nodal prestress.
/// @param[in] ya_l_f,ya_l_s,ya_l_n Nodal active tension along fiber, sheet, normal.
/// @param[in,out] lR,lK Element residual and stiffness, accumulated into.
void struct_3d_linear(ComMod &com_mod, CepMod &cep_mod, const int eNoN,
                      const int nFn, const int nG, const Vector<double> &wg,
                      const Array<double> &Ng, const double Jac,
                      const Array<double> &Nx, const Array<double> &al,
                      const Array<double> &yl, const Array<double> &dl,
                      const Array<double> &bfl, const Array<double> &fN,
                      const Array<double> &pS0l, const Vector<double> &ya_l_f,
                      const Vector<double> &ya_l_s, const Vector<double> &ya_l_n,
                      Array<double> &lR, Array3<double> &lK);

/// @brief Assemble the 3D solid element residual and stiffness at one Gauss point.
///
/// @param[in] eNoN Number of element nodes.
/// @param[in] nFn Number of fiber directions.
/// @param[in] w Gauss weight times the Jacobian.
/// @param[in] N,Nx Shape functions and their spatial derivatives.
/// @param[in] al,yl,dl,bfl Nodal acceleration, velocity, displacement, body force.
/// @param[in] fN Fiber directions.
/// @param[in] pS0l Nodal prestress.
/// @param[in] ya_l_f,ya_l_s,ya_l_n Nodal active tension along fiber, sheet, normal.
/// @param[out] pSl Second Piola-Kirchhoff stress in Voigt form at this point.
/// @param[in,out] lR,lK Element residual and stiffness, accumulated into.
void struct_3d(ComMod &com_mod, CepMod &cep_mod, const int eNoN, const int nFn,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Array<double> &fN, const Array<double> &pS0l,
               Vector<double> &pSl, const Vector<double> &ya_l_f,
               const Vector<double> &ya_l_s, const Vector<double> &ya_l_n,
               Array<double> &lR, Array3<double> &lK);
};

#endif

