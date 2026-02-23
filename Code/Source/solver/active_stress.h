#ifndef ACTIVE_STRESS_H
#define ACTIVE_STRESS_H

#include "Array.h"

#include <memory>

namespace active_stress {

struct ElementActiveStressParams
{
  double eta_f = 1.0;
  double eta_s = 0.0;
  double eta_n = 0.0;

  // Stage-2 placeholders.
  double scale = 1.0;
  double delay = 0.0;
};

class ActiveStressSpatialField
{
  public:
    virtual ~ActiveStressSpatialField() = default;
    virtual ElementActiveStressParams params(const int elem_id) const = 0;
};

class UniformSpatialField : public ActiveStressSpatialField
{
  public:
    explicit UniformSpatialField(const ElementActiveStressParams& params);
    ElementActiveStressParams params(const int elem_id) const override;

  private:
    ElementActiveStressParams params_;
};

class ElementalDirectionalSpatialField : public ActiveStressSpatialField
{
  public:
    ElementalDirectionalSpatialField(const Array<double>& elemental_distribution, const ElementActiveStressParams& fallback_params);
    ElementActiveStressParams params(const int elem_id) const override;

  private:
    const Array<double>& elemental_distribution_;
    ElementActiveStressParams fallback_params_;
};

class PrescribedActiveStressModel
{
  public:
    explicit PrescribedActiveStressModel(std::unique_ptr<ActiveStressSpatialField> spatial_field);
    void directional_stress(const double total_active_stress, const int elem_id, double& Tfa, double& Tsa, double& Tna) const;

  private:
    std::unique_ptr<ActiveStressSpatialField> spatial_field_;
};

} // namespace active_stress

#endif
