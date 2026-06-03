#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"

namespace active_stress {

// Number of per-element active-stress parameters currently stored per element
// in the (N_ACTIVE_STRESS_PARAMS x nElem) array carried on mshType.
// Row layout: row 0 -> eta_f, row 1 -> eta_s, row 2 -> eta_n.
// Bump to 5 and fill rows 3/4 to add spatially-varying scale/delay.
static constexpr int N_ACTIVE_STRESS_PARAMS = 3;

struct ElementActiveStressParams
{
  double eta_f = 1.0;
  double eta_s = 0.0;
  double eta_n = 0.0;

  // Reserved for future spatially-varying magnitude / time delay.
  double scale = 1.0;
  double delay = 0.0;
};

/// @brief Resolve the active-stress directional parameters for a single element.
///
/// When a per-element distribution is available (params != nullptr and elem_id
/// is in range) the eta fractions are read from the local (N_ACTIVE_STRESS_PARAMS
/// x nElem) column elem_id. Otherwise the uniform domain defaults are returned.
/// No heap allocation or virtual dispatch: safe to call in the Gauss-point loop.
inline ElementActiveStressParams resolve_active_stress_params(
    const Array<double>* params, const int elem_id,
    const double def_eta_f, const double def_eta_s, const double def_eta_n)
{
  ElementActiveStressParams p;
  p.eta_f = def_eta_f;
  p.eta_s = def_eta_s;
  p.eta_n = def_eta_n;

  if (params != nullptr && params->size() != 0 && elem_id >= 0 && elem_id < params->ncols()) {
    p.eta_f = (*params)(0, elem_id);
    p.eta_s = (*params)(1, elem_id);
    p.eta_n = (*params)(2, elem_id);
    // future: p.scale = (*params)(3, elem_id); p.delay = (*params)(4, elem_id);
  }

  return p;
}

} // namespace active_stress

#endif
