#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"

#include <string>

// Forward declaration so the header doesn't pull in consts.h (which would create
// an include cycle via ComMod.h). The enumerators are only used in active_stress.cpp.
namespace consts { enum class ConstitutiveModelType; }

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
  AS_IDX_SCALE = 4,
};

static constexpr int N_ACTIVE_STRESS_PARAMS = 5;

struct ElementActiveStressParams
{
  double eta_f = 1.0;
  double eta_s = 0.0;
  double eta_n = 0.0;
  double delay = 0.0;
  double scale = 1.0;   // per-element active-stress magnitude multiplier
};

/// @brief Active stress distributed among the fiber / sheet / sheet-normal
/// directions (the output of compute_fib_stress). Member order matches the
/// structured-binding call site: auto [Tfa, Tsa, Tna] = compute_fib_stress(...).
struct DirectionalActiveStress
{
  double Tfa = 0.0;   // fiber direction
  double Tsa = 0.0;   // sheet direction
  double Tna = 0.0;   // sheet-normal direction
};

/// @brief Encapsulates the spatial active-stress directional distribution:
/// the uniform fallback (eta_f/s/n, delay, scale), the optional per-element
/// distribution read from a VTU, and the read/validation/resolution logic.
///
/// One instance lives on the domain material (stM.Tf). `data_` is the global
/// read/distribution buffer (master only); during partitioning the bridge copies
/// it onto the owning mesh (mshType::active_stress_params) and `clear_data()`s it.
/// At runtime `params()` resolves against the caller-supplied LOCAL mesh array
/// (the scattered slice), falling back to the uniform defaults.
class ActiveStressField
{
  public:
    ActiveStressField() = default;

    // --- setup (master, from read_mat_model) ---

    /// @brief Set the uniform directional fractions (from the XML
    /// Directional_distribution block). Uniform delay/scale stay 0/1.
    void set_uniform_eta(double eta_f, double eta_s, double eta_n)
    {
      uniform_.eta_f = eta_f;
      uniform_.eta_s = eta_s;
      uniform_.eta_n = eta_n;
    }

    /// @brief Load the per-element distribution from a VTU and fill unsupplied eta
    /// rows from the uniform eta. Records which groups were supplied and sets the
    /// spatially-variable flag. All validation is deferred to validate().
    void read_from_vtu(const std::string& file_path);

    /// @brief Setup-time validation (all of it): supplied per-element eta sum to 1
    /// and are non-negative; supplied delay is non-negative and requires an unsteady
    /// curve; and eta_s/eta_n are only supported by Guccione/HO/HO-ma. Uses the
    /// uniform and (if present) per-element fractions.
    void validate(consts::ConstitutiveModelType isoType, bool is_unsteady) const;

    // --- distribution support (used by the bridge in distribute.cpp) ---
    bool spatially_variable() const { return spatially_variable_; }
    Array<double>& data() { return data_; }
    void clear_data() { data_.clear(); spatially_variable_ = false; }
    ElementActiveStressParams& uniform() { return uniform_; }

    // --- hot path (per Gauss point) ---
    /// @brief Resolve the per-element params for one element against the caller-
    /// supplied LOCAL array (the mesh's scattered slice): read its column if present
    /// and in range, else fall back to the uniform defaults. Allocation-free,
    /// non-virtual.
    ElementActiveStressParams params(const Array<double>* local_array, int elem_id) const
    {
      ElementActiveStressParams p = uniform_;
      if (local_array != nullptr && local_array->size() != 0 &&
          elem_id >= 0 && elem_id < local_array->ncols()) {
        p.eta_f = (*local_array)(AS_IDX_ETA_F, elem_id);
        p.eta_s = (*local_array)(AS_IDX_ETA_S, elem_id);
        p.eta_n = (*local_array)(AS_IDX_ETA_N, elem_id);
        p.delay = (*local_array)(AS_IDX_DELAY, elem_id);
        p.scale = (*local_array)(AS_IDX_SCALE, elem_id);
      }
      return p;
    }

  private:
    bool spatially_variable_ = false;
    ElementActiveStressParams uniform_;   // uniform defaults (eta from XML; delay 0 / scale 1)
    Array<double> data_;                  // (N x nElem) read/distribution buffer, global order
    bool vtu_has_eta_ = false;
    bool vtu_has_delay_ = false;
    bool vtu_has_scale_ = false;
};

} // namespace active_stress

#endif
