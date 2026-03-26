#include "gaussian-channel-model.h"

#include "ns3/log.h"
#include "ns3/mobility-model.h"
#include "ns3/ptr.h"
#include "ns3/string.h"
#include "ns3/uinteger.h"
#include "ns3/vector.h"
#include "ns3/simulator.h"

#include <complex>

namespace ns3
{

NS_LOG_COMPONENT_DEFINE("GaussianChannelModel");
NS_OBJECT_ENSURE_REGISTERED(GaussianChannelModel);

TypeId
GaussianChannelModel::GetTypeId()
{
  static TypeId tid = TypeId("ns3::GaussianChannelModel")
                          .SetParent<ThreeGppChannelModel>()
                          .SetGroupName("Spectrum")
                          .AddConstructor<GaussianChannelModel>()
                          .AddAttribute("SigmaLinear",
                                        "Linear sigma used to scale the CN(0,1) samples. Channel gain variance = sigma^2.",
                                        DoubleValue(1.0),
                                        MakeDoubleAccessor(&GaussianChannelModel::m_sigmaLinear),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("Clusters",
                                        "Number of clusters to synthesize (keeps compatibility with 3GPP tables).",
                                        UintegerValue(1),
                                        MakeUintegerAccessor(&GaussianChannelModel::m_clusters),
                                        MakeUintegerChecker<uint8_t>())
                          .AddAttribute("RaysPerCluster",
                                        "Rays per cluster (compatibility switch, typically 1).",
                                        UintegerValue(1),
                                        MakeUintegerAccessor(&GaussianChannelModel::m_raysPerCluster),
                                        MakeUintegerChecker<uint8_t>());

  return tid;
}

GaussianChannelModel::GaussianChannelModel()
    : m_sigmaLinear(1.0), m_clusters(1), m_raysPerCluster(1)
{
  NS_LOG_FUNCTION(this);

  // Create and configure RNGs here (consistent with ThreeGppChannelModel style)
  m_normalRv = CreateObject<NormalRandomVariable>();
  m_normalRv->SetAttribute("Mean", DoubleValue(0.0));
  // By convention ns-3 NormalRandomVariable wants Variance attribute for variance
  m_normalRv->SetAttribute("Variance", DoubleValue(1.0));

  m_uniformRv = CreateObject<UniformRandomVariable>();
  m_uniformRv->SetAttribute("Min", DoubleValue(0.0));
  m_uniformRv->SetAttribute("Max", DoubleValue(1.0));

  // Leave other ThreeGppChannelModel members intact (we inherit them)
}

GaussianChannelModel::~GaussianChannelModel() = default;

void
GaussianChannelModel::DoDispose()
{
  NS_LOG_FUNCTION(this);
  m_normalRv = nullptr;
  m_uniformRv = nullptr;
  ThreeGppChannelModel::DoDispose();
}

/**
 * AssignStreams:
 * assign stream indices to RNGs to allow deterministic parallel runs.
 */
int64_t
GaussianChannelModel::AssignStreams(int64_t stream)
{
  NS_LOG_FUNCTION(this << stream);
  int64_t current = stream;
  if (m_normalRv)
  {
    m_normalRv->SetStream(current);
    current++;
  }
  if (m_uniformRv)
  {
    m_uniformRv->SetStream(current);
    current++;
  }
  // ThreeGppChannelModel (base) has its own RNGs; do not attempt to reassign them here.
  return (current - stream);
}

/**
 * GetThreeGppTable:
 * We provide a minimal ParamsTable that matches the interface expected by the
 * ThreeGppChannelModel::GetChannel() control flow. The table values are not used
 * for 3GPP equations (we will ignore them), but some code paths expect a non-null
 * Ptr<ParamsTable>.
 */
Ptr<const ThreeGppChannelModel::ParamsTable>
GaussianChannelModel::GetThreeGppTable(const Ptr<const MobilityModel> aMob,
                                       const Ptr<const MobilityModel> bMob,
                                       Ptr<const ChannelCondition> channelCondition) const
{
  NS_LOG_FUNCTION(this << aMob << bMob);

  Ptr<ThreeGppChannelModel::ParamsTable> table = Create<ThreeGppChannelModel::ParamsTable>();
  // Minimal compat defaults:
  table->m_numOfCluster = m_clusters;
  table->m_raysPerCluster = m_raysPerCluster;
  // other fields left as defaults (they won't be used)
  return table;
}

/**
 * GenerateChannelParameters:
 * Minimal parameter generator: fills a ThreeGppChannelParams with
 * - distances (2D/3D)
 * - los condition (from channelCondition if provided)
 * - cluster power (single cluster with power = 1.0)
 *
 * The base GetChannel() logic will cache these params in m_channelParamsMap.
 */
Ptr<ThreeGppChannelModel::ThreeGppChannelParams>
GaussianChannelModel::GenerateChannelParameters(const Ptr<const ChannelCondition> channelCondition,
                                                const Ptr<const ThreeGppChannelModel::ParamsTable> table3gpp,
                                                const Ptr<const MobilityModel> aMob,
                                                const Ptr<const MobilityModel> bMob) const
{
  NS_LOG_FUNCTION(this << aMob << bMob);

  Ptr<ThreeGppChannelModel::ThreeGppChannelParams> params =
      Create<ThreeGppChannelModel::ThreeGppChannelParams>();

  // LOS/NLOS if available from channelCondition
  if (channelCondition)
  {
    params->m_losCondition = channelCondition->GetLosCondition();
    params->m_o2iCondition = channelCondition->GetO2iCondition();
  }
  else
  {
    // fallback: assume NLOS
    params->m_losCondition = ChannelCondition::LOS;
    //params->m_o2iCondition = ChannelCondition::NOT_O2I;
  }

  // positions
  Vector pa = aMob->GetPosition();
  Vector pb = bMob->GetPosition();
  double dx = pb.x - pa.x;
  double dy = pb.y - pa.y;
  double dz = pb.z - pa.z;
  params->m_dis2D = std::sqrt(dx * dx + dy * dy);
  params->m_dis3D = std::sqrt(dx * dx + dy * dy + dz * dz);

  // Minimal cluster setup: single cluster with unit power (keeps signatures compatible)
  params->m_clusterPower.clear();
  for (uint8_t c = 0; c < m_clusters; ++c)
  {
    params->m_clusterPower.push_back(1.0 / static_cast<double>(m_clusters)); // normalized
  }
  params->m_reducedClusterNumber = static_cast<uint8_t>(params->m_clusterPower.size());

  // minimal other defaults
  params->m_DS = 0.0;
  params->m_K_factor = 0.0;
  params->m_cluster1st = 0;
  params->m_cluster2nd = 0;

  return params;
}

/**
 * GetNewChannel:
 * Build a ChannelMatrix where each (rx,tx,cluster) coefficient is drawn as
 * h = (1/sqrt(2)) * sigma * (N(0,1) + j*N(0,1))  which is CN(0, sigma^2).
 *
 * We respect antenna element counts and cluster counts so the ChannelMatrix shape
 * is compatible with the rest of the code using ThreeGppChannelModel.
 */
Ptr<MatrixBasedChannelModel::ChannelMatrix>
GaussianChannelModel::GetNewChannel(Ptr<const ThreeGppChannelModel::ThreeGppChannelParams> channelParams,
                                   Ptr<const ThreeGppChannelModel::ParamsTable> table3gpp,
                                   const Ptr<const MobilityModel> sMob,
                                   const Ptr<const MobilityModel> uMob,
                                   Ptr<const PhasedArrayModel> sAntenna,
                                   Ptr<const PhasedArrayModel> uAntenna) const
{
  NS_LOG_FUNCTION(this << channelParams << sMob << uMob);

  // Create channel matrix object (same factory used by ThreeGppChannelModel)
  Ptr<MatrixBasedChannelModel::ChannelMatrix> channelMatrix = Create<MatrixBasedChannelModel::ChannelMatrix>();

  // Number of antenna elements
  size_t sSize = sAntenna->GetNumElems();
  size_t uSize = uAntenna->GetNumElems();

  // Determine number of clusters/pages. Keep compatibility with ThreeGpp logic:
  uint16_t numOverallCluster = static_cast<uint16_t>(channelParams->m_reducedClusterNumber);

  if (numOverallCluster == 0)
  {
    numOverallCluster = 1;
  }

  // If the 3GPP code expects some offset (e.g., +2/+4), we keep it minimal and safe:
  // (we produce exactly m_reducedClusterNumber pages)
  Complex3DVector hUsn(uSize, sSize, numOverallCluster);

  // Fill each (u,s,n) with complex Gaussian CN(0, sigma^2)
  double sigma = m_sigmaLinear;
  double scale = sigma / std::sqrt(2.0); // for real and imag parts

  for (size_t uIndex = 0; uIndex < uSize; ++uIndex)
  {
    for (size_t sIndex = 0; sIndex < sSize; ++sIndex)
    {
      for (uint16_t nIndex = 0; nIndex < numOverallCluster; ++nIndex)
      {
        double re = m_normalRv->GetValue();
        double im = m_normalRv->GetValue();
        std::complex<double> val(re * scale, im * scale);

        // Optionally scale by cluster power if provided
        double clusterPower = 1.0;
        if (nIndex < channelParams->m_clusterPower.size())
        {
          clusterPower = channelParams->m_clusterPower[nIndex];
        }
        // Apply cluster power as sqrt(power)
        val *= std::sqrt(clusterPower);

        hUsn(uIndex, sIndex, nIndex) = val;
      }
    }
  }

  // assign channel matrix fields compatible with ThreeGppChannelModel expectations
  channelMatrix->m_channel = hUsn;
  channelMatrix->m_antennaPair = std::make_pair(sAntenna->GetId(), uAntenna->GetId());

  return channelMatrix;
}

} // namespace ns3

