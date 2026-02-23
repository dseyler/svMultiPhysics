#include "active_stress.h"

namespace active_stress {

UniformSpatialField::UniformSpatialField(const ElementActiveStressParams& params)
  : params_(params)
{
}

ElementActiveStressParams UniformSpatialField::params(const int elem_id) const
{
  return params_;
}

ElementalDirectionalSpatialField::ElementalDirectionalSpatialField(
    const Array<double>& elemental_distribution,
    const ElementActiveStressParams& fallback_params)
  : elemental_distribution_(elemental_distribution)
  , fallback_params_(fallback_params)
{
}

ElementActiveStressParams ElementalDirectionalSpatialField::params(const int elem_id) const
{
  if (elem_id < 0 || elemental_distribution_.size() == 0 || elem_id >= elemental_distribution_.ncols()) {
    return fallback_params_;
  }

  ElementActiveStressParams elem_params = fallback_params_;
  elem_params.eta_f = elemental_distribution_(0, elem_id);
  elem_params.eta_s = elemental_distribution_(1, elem_id);
  elem_params.eta_n = elemental_distribution_(2, elem_id);
  return elem_params;
}

PrescribedActiveStressModel::PrescribedActiveStressModel(std::unique_ptr<ActiveStressSpatialField> spatial_field)
  : spatial_field_(std::move(spatial_field))
{
}

void PrescribedActiveStressModel::directional_stress(const double total_active_stress, const int elem_id, double& Tfa, double& Tsa, double& Tna) const
{
  auto elem_params = spatial_field_->params(elem_id);
  Tfa = elem_params.eta_f * total_active_stress;
  Tsa = elem_params.eta_s * total_active_stress;
  Tna = elem_params.eta_n * total_active_stress;
}

} // namespace active_stress
