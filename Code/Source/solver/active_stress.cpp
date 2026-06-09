#include "active_stress.h"

#include "consts.h"
#include "vtk_xml_parser.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace active_stress {

//------------------
// read_from_vtu
//------------------
// Load per-element active-stress parameters (directional fractions and/or
// activation delay/magnitude) from a VTU file. eta, delay and scale are
// independently optional; unsupplied eta rows are filled with the uniform values.
//
void ActiveStressField::read_from_vtu(const std::string& file_path)
{
  vtk_xml_parser::load_active_stress_directional_distribution_vtu(
      file_path, data_, vtu_has_eta_, vtu_has_delay_, vtu_has_scale_);

  spatially_variable_ = vtu_has_eta_ || vtu_has_delay_ || vtu_has_scale_;
  if (!spatially_variable_) {
    data_.clear();
    return;
  }

  const int ncols = data_.ncols();

  // Fill the eta rows: VTU values when supplied, else the uniform eta in every
  // column. (The delay and scale rows are always written by the loader.)
  if (!vtu_has_eta_) {
    for (int e = 0; e < ncols; e++) {
      data_(AS_IDX_ETA_F, e) = uniform_.eta_f;
      data_(AS_IDX_ETA_S, e) = uniform_.eta_s;
      data_(AS_IDX_ETA_N, e) = uniform_.eta_n;
    }
  }

  // Validate eta only when supplied by the VTU (uniform defaults are validated
  // at parse time): non-negative and summing to 1 per element.
  if (vtu_has_eta_) {
    const double tol = 1.0e-10;
    for (int e = 0; e < ncols; e++) {
      double eta_f = data_(AS_IDX_ETA_F, e);
      double eta_s = data_(AS_IDX_ETA_S, e);
      double eta_n = data_(AS_IDX_ETA_N, e);

      if (eta_f < 0.0 || eta_s < 0.0 || eta_n < 0.0) {
        throw std::runtime_error("Elemental directional distribution contains negative eta values at element index " +
            std::to_string(e+1) + ".");
      }

      double eta_sum = eta_f + eta_s + eta_n;
      if (std::abs(eta_sum - 1.0) > tol) {
        throw std::runtime_error("Elemental directional distribution fractions must sum to 1.0 at element index " +
            std::to_string(e+1) + ". Sum is " + std::to_string(eta_sum) + ".");
      }
    }
  }

  // Validate delay only when supplied: non-negative.
  if (vtu_has_delay_) {
    for (int e = 0; e < ncols; e++) {
      double delay = data_(AS_IDX_DELAY, e);
      if (delay < 0.0) {
        throw std::runtime_error("Active-stress delay must be >= 0; got " + std::to_string(delay) +
            " at element index " + std::to_string(e+1) + ".");
      }
    }
  }
}

//------------
// validate
//------------
// Setup-time checks that need the constitutive model and the steady/unsteady
// flag: a per-element delay only applies to an unsteady curve, and eta_s/eta_n
// are only supported by Guccione/HO/HO-ma.
//
void ActiveStressField::validate(consts::ConstitutiveModelType isoType, bool is_unsteady) const
{
  // A per-element delay requires a temporal (unsteady) stress curve.
  if (vtu_has_delay_ && !is_unsteady) {
    throw std::runtime_error("A per-element active-stress 'delay' field was provided, but the fiber "
        "reinforcement stress is not Unsteady. Delay requires a temporal stress curve "
        "(Fiber_reinforcement_stress type=\"Unsteady\").");
  }

  // Sheet / sheet-normal active stress is only supported by Guccione, HO, HO-ma.
  // Checked against the configured fractions (uniform and per-element) so a
  // misconfiguration fails immediately and deterministically at setup.
  bool supports_directional_distribution = (isoType == consts::ConstitutiveModelType::stIso_Gucci ||
                                            isoType == consts::ConstitutiveModelType::stIso_HO ||
                                            isoType == consts::ConstitutiveModelType::stIso_HO_ma);
  if (!supports_directional_distribution) {
    double max_eta_s = uniform_.eta_s;
    double max_eta_n = uniform_.eta_n;
    if (spatially_variable_) {
      for (int e = 0; e < data_.ncols(); e++) {
        max_eta_s = std::max(max_eta_s, data_(AS_IDX_ETA_S, e));
        max_eta_n = std::max(max_eta_n, data_(AS_IDX_ETA_N, e));
      }
    }
    if (max_eta_s > 0.0 || max_eta_n > 0.0) {
      throw std::runtime_error("Directional distribution of active stress (eta_s > 0 or eta_n > 0) "
        "is only supported for Guccione, Holzapfel-Ogden (HO), and Holzapfel-Ogden Modified Anisotropy (HO-ma) models. "
        "Current model does not support sheet or sheet-normal stress contributions. "
        "Set Fiber_direction=1.0, Sheet_direction=0.0, Sheet_normal_direction=0.0.");
    }
  }
}

} // namespace active_stress
