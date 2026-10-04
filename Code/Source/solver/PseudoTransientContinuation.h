// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PSEUDO_TRANSIENT_CONTINUATION_H
#define PSEUDO_TRANSIENT_CONTINUATION_H

#include "Array.h"
#include "CmMod.h"
#include "ComMod.h"
#include "SolutionStates.h"

/// @brief Pseudo-transient continuation: solve a steady nonlinear problem by
/// time stepping with a pseudo time step that grows as the residual falls, so
/// that the iteration turns into Newton's method (Kelley & Keyes, SIAM J.
/// Numer. Anal. 35 (1998) 508-523).
///
/// Every pseudo-step takes one Newton iteration. The pseudo time step of the
/// next step is set from the ratio of successive first residuals, capped by
/// max_dt; a time step so large that the inertia term vanishes makes the step
/// a plain Newton iteration. Time is fictitious, so the loads must be steady.
///
/// Owned and driven by the Integrator. Enabled by Pseudo_transient on an
/// equation, whose residual is the one monitored. Prestrain is the only
/// feature that uses it so far: each of its steps starts from rest and the
/// displacement it reaches is composed into the prestrain (see prestrain.h).
class PseudoTransientContinuation {

  public:
    PseudoTransientContinuation() = default;
    explicit PseudoTransientContinuation(const PseudoTransientSettings& settings);

    bool enabled() const { return settings_.enabled; }

    /// @brief Finish a pseudo-step: report the update it reached and set the
    /// next pseudo time step from its first residual. Called once the Newton
    /// loop has exited.
    void finish_step(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);

  private:
    /// @brief The pseudo time step for the next step, from the first residuals
    /// of the last two.
    double next_time_step(const double dt, const double residual_prev, const double residual) const;

    /// @brief Largest nodal update of the monitored equation's displacement
    /// dofs, over all ranks.
    double max_update(const ComMod& com_mod, const CmMod& cm_mod, const Array<double>& update) const;

    PseudoTransientSettings settings_;

    /// @brief First Newton residual of the previous pseudo-step
    double residual_prev_ = 0.0;
};

#endif
