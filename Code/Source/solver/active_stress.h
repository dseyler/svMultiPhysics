#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"

namespace active_stress {

// Row indices into the (N_ACTIVE_STRESS_PARAMS x nElem) per-element array carried
// on mshType::active_stress_params (and, before partitioning, on
// fibStrsType::elemental_distribution).
// To add a spatially-varying parameter: add a field to ElementActiveStressParams,
// add a row index here, bump N_ACTIVE_STRESS_PARAMS, read/fill/validate it in the
// VTU loader + read_files, and consume it where it applies.
enum AsParamRow : int {
  AS_IDX_ETA_F = 0,
  AS_IDX_ETA_S = 1,
  AS_IDX_ETA_N = 2,
  AS_IDX_DELAY = 3,
};

static constexpr int N_ACTIVE_STRESS_PARAMS = 4;

struct ElementActiveStressParams
{
  double eta_f = 1.0;
  double eta_s = 0.0;
  double eta_n = 0.0;
  double delay = 0.0;

  // double scale = 1.0;   // future: spatially-varying active-stress magnitude
};

/// @brief Resolve the per-element active-stress parameters for a single element.
///
/// When a per-element distribution is available (params != nullptr and elem_id
/// is in range) every field is read from the local (N_ACTIVE_STRESS_PARAMS x nElem)
/// column elem_id. Otherwise `defaults` (the uniform domain values) are returned.
/// No heap allocation or virtual dispatch: safe to call in the Gauss-point loop.
inline ElementActiveStressParams resolve_active_stress_params(
    const Array<double>* params, const int elem_id, const ElementActiveStressParams& defaults)
{
  ElementActiveStressParams p = defaults;

  if (params != nullptr && params->size() != 0 && elem_id >= 0 && elem_id < params->ncols()) {
    p.eta_f = (*params)(AS_IDX_ETA_F, elem_id);
    p.eta_s = (*params)(AS_IDX_ETA_S, elem_id);
    p.eta_n = (*params)(AS_IDX_ETA_N, elem_id);
    p.delay = (*params)(AS_IDX_DELAY, elem_id);
  }

  return p;
}

} // namespace active_stress

#endif
