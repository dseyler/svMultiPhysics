// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef STRUCT_H 
#define STRUCT_H 

#include "ComMod.h"
#include "SolutionStates.h"

namespace struct_ns {

/// @brief Per-Gauss-point scratch buffers for struct_3d.
///
/// struct_3d is called once per element per Gauss point -- 28 M times for ten
/// steps of the 108k-node cardiac case -- and used to construct all of these
/// on entry and destroy them on exit. Array::allocate() is a new[] followed by
/// a memset, so that was 13 allocations and 474 zeroed doubles per call, of
/// which only vx actually depended on the zeroing.
///
/// Sizes depend only on nsd and eNoN, both fixed for a mesh, so the caller
/// builds one of these outside the element loop and hands it in.
struct Struct3dScratch {
  Array<double> F, S0, vx, S, Dm, Svis, P, DBm;
  Array3<double> Kvis_u, Kvis_v, Bm;
  Vector<double> ud, fb;

  Struct3dScratch(const int nsd, const int eNoN) :
    F(nsd,nsd), S0(nsd,nsd), vx(nsd,nsd), S(nsd,nsd),
    Dm(3*(nsd-1), 3*(nsd-1)), Svis(nsd,nsd), P(nsd,nsd), DBm(3*(nsd-1), nsd),
    Kvis_u(nsd*nsd, eNoN, eNoN), Kvis_v(nsd*nsd, eNoN, eNoN),
    Bm(3*(nsd-1), nsd, eNoN),
    ud(nsd), fb(nsd) {}
};

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

void struct_3d(ComMod &com_mod, CepMod &cep_mod, const int eNoN, const int nFn,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const Array<double> &bfl,
               const Array<double> &fN, const Array<double> &pS0l,
               Vector<double> &pSl, const Vector<double> &ya_l_f,
               const Vector<double> &ya_l_s, const Vector<double> &ya_l_n,
               Array<double> &lR, Array3<double> &lK,
               Struct3dScratch &scr, const bool update_elem_invariants);
};

#endif

