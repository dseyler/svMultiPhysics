// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "PseudoTransientContinuation.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "mpi.h"

PseudoTransientContinuation::PseudoTransientContinuation(const PseudoTransientSettings& settings)
  : settings_(settings)
{
}

void PseudoTransientContinuation::finish_step(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions)
{
  const auto& eq = com_mod.eq[settings_.equation];
  double& dt = com_mod.dt;

  // pNorm is the step's first residual relative to iNorm, the run's first.
  const double residual = eq.pNorm * eq.iNorm;
  const double umax = max_update(com_mod, cm_mod, solutions.current.get_displacement() * eq.af);

  if (com_mod.cm.mas(cm_mod)) {
    std::cout << " Pseudo-transient: max update = " << umax << "  (dt = " << dt
              << ", first residual = " << residual << ", R/R0 = " << eq.pNorm << ")" << std::endl;
  }

  if (settings_.adaptive_dt) {
    if (settings_.max_dt <= 0.0) {
      settings_.max_dt = 100.0 * dt;
    }
    dt = next_time_step(dt, residual_prev_, residual);
  }
  residual_prev_ = residual;
}

// Switched evolution relaxation (Mulder & van Leer 1985), the controller analysed by Kelley & Keyes (1998).
double PseudoTransientContinuation::next_time_step(const double dt, const double residual_prev, const double residual) const
{
  // Nothing to compare with on the first step; leave dt where it is once the
  // residual has reached round-off.
  if (residual_prev <= 0.0 || !(residual > 0.0) || !std::isfinite(residual)) {
    return dt;
  }
  return std::min(settings_.max_dt, dt * residual_prev / residual);
}

double PseudoTransientContinuation::max_update(const ComMod& com_mod, const CmMod& cm_mod, const Array<double>& update) const
{
  const auto& eq = com_mod.eq[settings_.equation];
  const int s = eq.s;
  const int n = std::min(eq.dof, com_mod.nsd);

  double umax = 0.0;
  for (int a = 0; a < com_mod.tnNo; a++) {
    double u2 = 0.0;
    for (int i = 0; i < n; i++) {
      u2 += update(s+i,a) * update(s+i,a);
    }
    // Let a NaN through rather than have max() hide it.
    umax = (u2 > umax || std::isnan(u2)) ? u2 : umax;
  }
  return std::sqrt(com_mod.cm.reduce(cm_mod, umax, MPI_MAX));
}
