#ifndef GAUSSIAN_PROPAGATION_LOSS_MODEL_H
#define GAUSSIAN_PROPAGATION_LOSS_MODEL_H

#include "three-gpp-propagation-loss-model.h"
#include "ns3/random-variable-stream.h"
#include "ns3/mobility-model.h"
#include "ns3/core-module.h"

namespace ns3 {

/**
 * @brief Gaussian propagation loss model.
 *
 * Inherits from ThreeGppPropagationLossModel to follow the same template.
 * Pathloss = FreeSpaceLoss(3D distance, freq) + GaussianOffset
 * NLOS adds an extra penalty (attribute).
 */
class GaussianPropagationLossModel : public ThreeGppPropagationLossModel
{
public:
  static TypeId GetTypeId();

  GaussianPropagationLossModel();
  ~GaussianPropagationLossModel() override;

  // delete copy/assign (follow base)
  GaussianPropagationLossModel(const GaussianPropagationLossModel&) = delete;
  GaussianPropagationLossModel& operator=(const GaussianPropagationLossModel&) = delete;

private:
  // Implement abstract virtuals from ThreeGppPropagationLossModel
  double GetLossLos(Ptr<MobilityModel> a, Ptr<MobilityModel> b) const override;
  double GetLossNlos(Ptr<MobilityModel> a, Ptr<MobilityModel> b) const override;
  double GetO2iDistance2dIn() const override;
  double GetShadowingStd(Ptr<MobilityModel> a,
                         Ptr<MobilityModel> b,
                         ChannelCondition::LosConditionValue cond) const override;
  double GetShadowingCorrelationDistance(ChannelCondition::LosConditionValue cond) const override;

  int64_t DoAssignStreams(int64_t stream) override;

  // helper
  static double FreeSpaceLoss(double frequencyHz, double distanceMeters);

private:
  // Parameters (attributes)
  double m_gaussMean;         //!< mean of Gaussian offset in dB
  double m_gaussStd;          //!< std dev of Gaussian offset in dB
  double m_nlosPenalty;       //!< extra penalty for NLOS in dB
  double m_shadowingStd;      //!< shadowing std dev (dB)
  double m_shadowingCorrDist; //!< shadowing correlation distance (m)
  double m_o2iMin;            //!< min for O2I 2d-in uniform draws
  double m_o2iMax;            //!< max for O2I 2d-in uniform draws

  Ptr<NormalRandomVariable> m_gaussVar;  //!< gaussian RV for per-link offset
  Ptr<UniformRandomVariable> m_o2iVar1;  //!< uniform rv for O2I (draw 1)
  Ptr<UniformRandomVariable> m_o2iVar2;  //!< uniform rv for O2I (draw 2)
};

} // namespace ns3

#endif // GAUSSIAN_PROPAGATION_LOSS_MODEL_H
