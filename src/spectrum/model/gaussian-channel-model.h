#ifndef GAUSSIAN_CHANNEL_MODEL_H
#define GAUSSIAN_CHANNEL_MODEL_H

#include "three-gpp-channel-model.h"



#include "ns3/double.h"
#include "ns3/integer.h"
#include "ns3/log.h"

#include <string>

namespace ns3
{

/**
 * @brief GaussianChannelModel
 *
 * Follows the ThreeGppChannelModel template and lifecycle but implements a very
 * simple fading model: per-cluster, per-(tx,rx) complex Gaussian samples
 * (i.i.d. CN(0, sigma^2)). The code uses the same maps / keys / update logic
 * as the ThreeGppChannelModel and stores params in the same ThreeGppChannelParams
 * structure (minimal fields filled).
 */
class GaussianChannelModel : public ThreeGppChannelModel
{
  public:
    static TypeId GetTypeId();

    GaussianChannelModel();
    ~GaussianChannelModel() override;

    // assign streams to internal RNGs (same name as ThreeGppChannelModel)
    int64_t AssignStreams(int64_t stream);

  protected:
    // We override the two extension points used by the base GetChannel flow:
    Ptr<const ParamsTable> GetThreeGppTable(const Ptr<const MobilityModel> aMob,
                                            const Ptr<const MobilityModel> bMob,
                                            Ptr<const ChannelCondition> channelCondition) const override;

    Ptr<ThreeGppChannelParams> GenerateChannelParameters(
        const Ptr<const ChannelCondition> channelCondition,
        const Ptr<const ParamsTable> table3gpp,
        const Ptr<const MobilityModel> aMob,
        const Ptr<const MobilityModel> bMob) const ;

    // This generates the actual matrix (replaces 3GPP cluster/ray math).
    Ptr<ChannelMatrix> GetNewChannel(Ptr<const ThreeGppChannelParams> channelParams,
                                    Ptr<const ParamsTable> table3gpp,
                                    const Ptr<const MobilityModel> sMob,
                                    const Ptr<const MobilityModel> uMob,
                                    Ptr<const PhasedArrayModel> sAntenna,
                                    Ptr<const PhasedArrayModel> uAntenna) const override;

    void DoDispose() override;

  private:
    // Model attributes
    double m_sigmaLinear;   //!< linear variance used for CN(0, sigma^2) entries
    uint8_t m_clusters;     //!< how many clusters (keeps interface: 1 by default)
    uint8_t m_raysPerCluster; //!< rays per cluster (keeps interface)

    // RNGs
    Ptr<NormalRandomVariable> m_normalRv; //!< for real/imag Gaussian draws
    Ptr<UniformRandomVariable> m_uniformRv; //!< used for simple draws (if needed)
};

} // namespace ns3

#endif // GAUSSIAN_CHANNEL_MODEL_H

