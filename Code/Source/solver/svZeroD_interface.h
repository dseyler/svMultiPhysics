// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SVZEROD_H
#define SVZEROD_H 

#include "Simulation.h"
#include "consts.h"
#include "svZeroD_interface/LPNSolverInterface.h"
#include <vector>

#include <string>

namespace svZeroD {

void get_coupled_QP(ComMod& com_mod, double QCoupled[], double QnCoupled[], double PCoupled[], double PnCoupled[]);

void print_svZeroD(int* nSrfs, const std::vector<int>& surfID, double Q[], double P[]);

void init_svZeroD(ComMod& com_mod, const CmMod& cm_mod);

void calc_svZeroD(ComMod& com_mod, const CmMod& cm_mod, char BCFlag);

/// Accessor for the active LPNSolverInterface; used by the
/// cplBC_I_analytical branch in set_bc.cpp::calc_der_cpl_bc to call
/// get_coupling_jacobian directly. Returns nullptr if no interface has
/// been initialized.
LPNSolverInterface* get_interface();

};

#endif
