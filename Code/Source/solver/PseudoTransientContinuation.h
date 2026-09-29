// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PSEUDO_TRANSIENT_CONTINUATION_H
#define PSEUDO_TRANSIENT_CONTINUATION_H

#include "Array.h"
#include "CmMod.h"
#include "ComMod.h"
#include "SolutionStates.h"

#include <functional>

/// @brief Pseudo-transient continuation: solve a steady nonlinear problem by
/// time stepping from rest with a pseudo time step that grows as the residual
/// falls, so that the iteration turns into Newton's method (Kelley & Keyes,
/// SIAM J. Numer. Anal. 35 (1998) 508-523).
///
/// Every pseudo-step starts from rest under the full load and takes one Newton
/// iteration. The state it reaches is handed to a client, which folds it into
/// whatever the run is really solving for, and the solution is then reset to
/// zero, so the geometry never moves. The pseudo time step of the next step is
/// set from the ratio of successive first residuals, capped by max_dt; a time
/// step so large that the inertia term vanishes makes the step a plain Newton
/// iteration. Time is fictitious, so the loads must be steady.
///
/// Owned and driven by the Integrator. Enabled by Pseudo_transient on an
/// equation, whose residual is the one monitored. A feature registers what
/// to do with each step's update through set_state_update; prestrain is the
/// first such client.
class PseudoTransientContinuation {

  public:
    /// @brief Receives the displacement update a pseudo-step reached, as the
    /// full tDof x tnNo array; the client reads the rows it owns.
    using StateUpdate = std::function<void(const Array<double>& update)>;

    PseudoTransientContinuation() = default;
    explicit PseudoTransientContinuation(const PseudoTransientSettings& settings);

    bool enabled() const { return settings_.enabled; }

    /// @brief Register the client that accumulates each pseudo-step's update.
    void set_state_update(StateUpdate update);

    /// @brief Start a pseudo-step: hand the state the previous step reached to
    /// the client and reset the solution to rest. Called from the predictor.
    ///
    /// The kernels saw the generalized-alpha level displacement alpha_f*Dn,
    /// which is the state the Newton update moved to, so that is the update;
    /// Dn itself would overshoot by 1/alpha_f. At the predictor Dn has been
    /// copied into Do.
    void start_step(ComMod& com_mod, SolutionStates& solutions);

    /// @brief Finish a pseudo-step: report it and set the next pseudo time
    /// step from its first residual. Called once the Newton loop has exited.
    void finish_step(ComMod& com_mod, const CmMod& cm_mod, const SolutionStates& solutions);

  private:
    /// @brief The pseudo time step for the next step, from the first residuals
    /// of the last two.
    double next_time_step(const double dt, const double residual_prev, const double residual) const;

    /// @brief Largest nodal update of the monitored equation's displacement
    /// dofs, over all ranks.
    double max_update(const ComMod& com_mod, const CmMod& cm_mod, const Array<double>& update) const;

    PseudoTransientSettings settings_;
    StateUpdate state_update_;

    /// @brief First Newton residual of the previous pseudo-step
    double residual_prev_ = 0.0;
};

#endif
