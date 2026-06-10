#include "active_stress.h"

#include "consts.h"
#include "VtkData.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace active_stress {

//------------------
// read_from_vtu
//------------------
// Load per-element active-stress parameters from a VTU file (CellData arrays keyed
// by Int32 GlobalElementID). eta_f/eta_s/eta_n (all-three-or-none), delay and scale
// are independently optional; unsupplied eta rows are filled with the uniform eta,
// delay defaults to 0 and scale to 1. Per the codebase idiom (VtkVtuData /
// copy_cell_data) a wrong-typed CellData array silently does not load. All value
// validation is deferred to validate().
//
void ActiveStressField::read_from_vtu(const std::string& file_path)
{
  VtkVtuData vtu(file_path);   // reads the file; throws on a missing/unreadable file
  const int n = vtu.num_elems();
  if (n <= 0) {
    throw std::runtime_error("Failed reading active stress VTU file '" + file_path + "'.");
  }

  bool has_eta_f = vtu.has_cell_data("eta_f");
  bool has_eta_s = vtu.has_cell_data("eta_s");
  bool has_eta_n = vtu.has_cell_data("eta_n");
  bool any_eta = has_eta_f || has_eta_s || has_eta_n;
  vtu_has_eta_   = has_eta_f && has_eta_s && has_eta_n;
  vtu_has_delay_ = vtu.has_cell_data("delay");
  vtu_has_scale_ = vtu.has_cell_data("scale");

  if (any_eta && !vtu_has_eta_) {
    throw std::runtime_error("Active stress directional distribution file '" + file_path +
        "' must either provide all three CellData arrays (eta_f, eta_s, eta_n) or none of them.");
  }

  spatially_variable_ = vtu_has_eta_ || vtu_has_delay_ || vtu_has_scale_;
  if (!spatially_variable_) {
    data_.clear();
    return;
  }

  if (!vtu.has_cell_data("GlobalElementID")) {
    throw std::runtime_error("Active stress directional distribution file '" + file_path +
        "' does not contain required CellData Int32 array 'GlobalElementID'.");
  }
  Vector<int> gid(n);
  vtu.copy_cell_data("GlobalElementID", gid);

  // Read the supplied arrays into cell-order temporaries (pre-sized; zero-init).
  Vector<double> eta_f, eta_s, eta_n, delay, scale;
  if (vtu_has_eta_) {
    eta_f.resize(n); eta_s.resize(n); eta_n.resize(n);
    vtu.copy_cell_data("eta_f", eta_f);
    vtu.copy_cell_data("eta_s", eta_s);
    vtu.copy_cell_data("eta_n", eta_n);
  }
  if (vtu_has_delay_) { delay.resize(n); vtu.copy_cell_data("delay", delay); }
  if (vtu_has_scale_) { scale.resize(n); vtu.copy_cell_data("scale", scale); }

  // Build the (N x nElem) array, indexed by GlobalElementID-1. Each row is filled
  // directly: VTU value where supplied, else the uniform eta / delay 0 / scale 1.
  data_.resize(N_ACTIVE_STRESS_PARAMS, n);
  std::vector<bool> seen(n, false);
  for (int cell = 0; cell < n; cell++) {
    int global_elem_id = gid(cell);
    if (global_elem_id <= 0 || global_elem_id > n) {
      throw std::runtime_error("GlobalElementID value '" + std::to_string(global_elem_id) +
          "' in active stress directional distribution file '" + file_path +
          "' is out of valid range [1, " + std::to_string(n) + "].");
    }

    int idx = global_elem_id - 1;
    if (seen[idx]) {
      throw std::runtime_error("Duplicate GlobalElementID value '" + std::to_string(global_elem_id) +
          "' in active stress directional distribution file '" + file_path + "'.");
    }
    seen[idx] = true;

    data_(AS_IDX_ETA_F, idx) = vtu_has_eta_   ? eta_f(cell) : uniform_.eta_f;
    data_(AS_IDX_ETA_S, idx) = vtu_has_eta_   ? eta_s(cell) : uniform_.eta_s;
    data_(AS_IDX_ETA_N, idx) = vtu_has_eta_   ? eta_n(cell) : uniform_.eta_n;
    data_(AS_IDX_DELAY, idx) = vtu_has_delay_ ? delay(cell) : 0.0;
    data_(AS_IDX_SCALE, idx) = vtu_has_scale_ ? scale(cell) : 1.0;
  }

  for (int i = 0; i < n; i++) {
    if (!seen[i]) {
      throw std::runtime_error("Missing GlobalElementID '" + std::to_string(i+1) +
          "' in active stress directional distribution file '" + file_path + "'.");
    }
  }
}

//------------
// validate
//------------
// All setup-time validation for the active-stress distribution: supplied
// per-element eta sum to 1 and are non-negative; supplied delay is non-negative
// and requires an unsteady curve; and eta_s/eta_n are only supported by
// Guccione/HO/HO-ma. (Uniform eta is validated at parse time.)
//
void ActiveStressField::validate(consts::ConstitutiveModelType isoType, bool is_unsteady) const
{
  const int ncols = data_.ncols();

  // Per-element eta supplied by the VTU: non-negative and summing to 1 (eta filled
  // from the already-validated uniform values needs no re-check).
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

  // Per-element delay supplied by the VTU: non-negative, and only meaningful for a
  // temporal (unsteady) stress curve.
  if (vtu_has_delay_) {
    for (int e = 0; e < ncols; e++) {
      double delay = data_(AS_IDX_DELAY, e);
      if (delay < 0.0) {
        throw std::runtime_error("Active-stress delay must be >= 0; got " + std::to_string(delay) +
            " at element index " + std::to_string(e+1) + ".");
      }
    }

    if (!is_unsteady) {
      throw std::runtime_error("A per-element active-stress 'delay' field was provided, but the fiber "
          "reinforcement stress is not Unsteady. Delay requires a temporal stress curve "
          "(Fiber_reinforcement_stress type=\"Unsteady\").");
    }
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
