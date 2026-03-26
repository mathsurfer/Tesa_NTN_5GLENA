#include "gaussian-propagation-loss-model.h"
#include "ns3/log.h"
#include "ns3/double.h"
#include "ns3/boolean.h"
#include <cmath>
#include <fstream>
#include <streambuf>

// >>> ADD
#include <map>
// <<< ADD

static std::ofstream g_propLog;


namespace ns3 {

NS_LOG_COMPONENT_DEFINE("GaussianPropagationLossModel");
NS_OBJECT_ENSURE_REGISTERED(GaussianPropagationLossModel);

TypeId
GaussianPropagationLossModel::GetTypeId()
{
  static TypeId tid = TypeId("ns3::GaussianPropagationLossModel")
                          .SetParent<ThreeGppPropagationLossModel>()
                          .SetGroupName("Propagation")
                          .AddConstructor<GaussianPropagationLossModel>()
                          .AddAttribute("GaussianMeanDb",
                                        "Mean (dB) of the Gaussian offset added to the pathloss.",
                                        DoubleValue(0.0),
                                        MakeDoubleAccessor(&GaussianPropagationLossModel::m_gaussMean),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("GaussianStdDb",
                                        "Std deviation (dB) of the Gaussian offset added to the pathloss.",
                                        DoubleValue(2.0),
                                        MakeDoubleAccessor(&GaussianPropagationLossModel::m_gaussStd),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("NlosPenaltyDb",
                                        "Additional dB penalty to add for NLOS conditions.",
                                        DoubleValue(20.0),
                                        MakeDoubleAccessor(&GaussianPropagationLossModel::m_nlosPenalty),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("ShadowingStdDb",
                                        "Shadowing standard deviation in dB (used by GetShadowing()).",
                                        DoubleValue(4.0),
                                        MakeDoubleAccessor(&GaussianPropagationLossModel::m_shadowingStd),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("ShadowingCorrelationDistance",
                                        "Shadowing correlation distance in meters.",
                                        DoubleValue(50.0),
                                        MakeDoubleAccessor(&GaussianPropagationLossModel::m_shadowingCorrDist),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("O2iDistance2dInMin",
                                        "Minimum for the O2I 2d-in uniform draws (meters).",
                                        DoubleValue(0.0),
                                        MakeDoubleAccessor(&GaussianPropagationLossModel::m_o2iMin),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("O2iDistance2dInMax",
                                        "Maximum for the O2I 2d-in uniform draws (meters).",
                                        DoubleValue(25.0),
                                        MakeDoubleAccessor(&GaussianPropagationLossModel::m_o2iMax),
                                        MakeDoubleChecker<double>());

  return tid;
}

GaussianPropagationLossModel::GaussianPropagationLossModel()
    : m_gaussMean(0.0),
      m_gaussStd(2.0),
      m_nlosPenalty(20.0),
      m_shadowingStd(0.0),
      m_shadowingCorrDist(50.0),
      m_o2iMin(0.0),
      m_o2iMax(25.0)
{
  NS_LOG_FUNCTION(this);
  if (!g_propLog.is_open())
  {
    g_propLog.open("propagation_log.txt", std::ios::out);
    g_propLog << "time(s)\tdistance(m)\tfspl(dB)\tgauss(dB)\tpathloss(dB)\n";
  }
  m_gaussVar = CreateObject<NormalRandomVariable>();
  m_gaussVar->SetAttribute("Mean", DoubleValue(m_gaussMean));
  m_gaussVar->SetAttribute("Variance", DoubleValue(m_gaussStd * m_gaussStd));
  // Note: ns-3 NormalRandomVariable accepts Mean and Variance attributes.

  m_o2iVar1 = CreateObject<UniformRandomVariable>();
  m_o2iVar2 = CreateObject<UniformRandomVariable>();
  m_o2iVar1->SetAttribute("Min", DoubleValue(m_o2iMin));
  m_o2iVar1->SetAttribute("Max", DoubleValue(m_o2iMax));
  m_o2iVar2->SetAttribute("Min", DoubleValue(m_o2iMin));
  m_o2iVar2->SetAttribute("Max", DoubleValue(m_o2iMax));
}

GaussianPropagationLossModel::~GaussianPropagationLossModel() = default;

int64_t
GaussianPropagationLossModel::DoAssignStreams(int64_t stream)
{
  NS_LOG_FUNCTION(this << stream);
  // assign streams to our random variables
  if (m_gaussVar)
  {
    m_gaussVar->SetStream(stream);
    stream++;
  }
  if (m_o2iVar1)
  {
    m_o2iVar1->SetStream(stream);
    stream++;
  }
  if (m_o2iVar2)
  {
    m_o2iVar2->SetStream(stream);
    stream++;
  }
  return stream;
}

/** Helper: free-space loss in dB, using frequency (Hz) and distance (m).
 *  returns large number if distance == 0 (caller should handle).
 */
double
GaussianPropagationLossModel::FreeSpaceLoss(double frequencyHz, double distanceMeters)
{
  if (distanceMeters <= 0.0)
  {
    return 0.0;
  }
  const double c = 299792458.0; // m/s
  double lambda = c / frequencyHz;
  double fspl = 20.0 * std::log10(4.0 * M_PI * distanceMeters / lambda);
  return fspl;
}

double
GaussianPropagationLossModel::GetLossLos(Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{
  NS_LOG_FUNCTION(this << a << b);
  // compute 3D distance
  Vector pa = a->GetPosition();
Vector pb = b->GetPosition();

double dx = pb.x - pa.x;
double dy = pb.y - pa.y;
double dz = pb.z - pa.z;

double d3 = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (d3 < 1e-9)
  {
    return 0.0; // zero pathloss (same position) -- base will handle tx power -> rx power
  }

  double freq = GetFrequency(); // from base class
  double fspl = FreeSpaceLoss(freq, d3);

  double gaussDb = 0.0;
  if (m_gaussVar)
  {
    // NormalRandomVariable::GetValue() uses its internal mean/variance attributes
    gaussDb = m_gaussVar->GetValue();
  }

  //double pathloss = fspl + gaussDb;
  double pathloss = fspl ;
  // ================= DEBUG RX POWER =================
double txPowerDbm = 30.0; 
double rxPowerDbm = txPowerDbm - pathloss;

NS_LOG_UNCOND(
  "[PROPAGATION] LOS"
  << " d3=" << d3 << " m"
  << " FSPL=" << fspl << " dB"
  << " Gaussian=" << gaussDb << " dB"
  << " Pathloss=" << pathloss << " dB"
  << " RxPower=" << rxPowerDbm << " dBm"
);
// ==================================================


  return pathloss;
}

double
GaussianPropagationLossModel::GetLossNlos(Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{
  NS_LOG_FUNCTION(this << a << b);
  double losLoss = GetLossLos(a, b);
  double pathloss = losLoss + m_nlosPenalty;
  return pathloss;
}

double
GaussianPropagationLossModel::GetO2iDistance2dIn() const
{
  NS_LOG_FUNCTION(this);
  // follow the pattern in the header you provided: generate two uniforms and return min
  double v1 = m_o2iVar1->GetValue();
  double v2 = m_o2iVar2->GetValue();
  return std::min(v1, v2);
}

double
GaussianPropagationLossModel::GetShadowingStd(Ptr<MobilityModel> a,
                                              Ptr<MobilityModel> b,
                                              ChannelCondition::LosConditionValue cond) const
{
  NS_LOG_FUNCTION(this << a << b << cond);
  // simple model returning configured shadowing std
  //return m_shadowingStd;
  return 0.0;
}

double
GaussianPropagationLossModel::GetShadowingCorrelationDistance(ChannelCondition::LosConditionValue cond) const
{
  NS_LOG_FUNCTION(this << cond);
  return m_shadowingCorrDist;
}

} // namespace ns3
