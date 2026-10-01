// Unit tests for nuclear densities, fluctuations and UPC screening
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <array>
#include <catch.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Nuclear/MBreakup.h"
#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MFluctuation.h"
#include "Graniitti/Nuclear/MGlauber.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MPhoto.h"
#include "Graniitti/Nuclear/MPhoton.h"
#include "Graniitti/Nuclear/MRecord.h"
#include "Graniitti/Nuclear/MShadow.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Tech/MLHE.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "json.hpp"
#include "support/nuclear_test_support.hh"

using gra::aux::indices;

namespace {

constexpr int    pb208_pdg  = 1000822080;
constexpr double pb208_mass = 193.687;

// Construct the default lead nucleus with reduced unit-test quadratures
gra::nuclear::MNucleus MakeLead() {
  gra::nuclear::GeometryParam geometry;
  geometry.radius_scale = 1.2;
  geometry.skin         = 0.55;
  geometry.nodes        = 96;
  geometry.tail_skin    = 18.0;
  geometry.cdf_nodes    = 2048;
  geometry.form_q_max   = 3.0;
  geometry.form_nodes   = 4097;
  geometry.form_abs_tol = 1.0e-7;
  geometry.nucleus.push_back({208, 82, {6.624, 0.549}, {6.620, 0.546}});
  auto param         = gra::nuclear::DefaultNucleusParam(pb208_pdg, pb208_mass, geometry);
  param.charge.nodes = 96;
  param.matter.nodes = 96;
  return gra::nuclear::MNucleus(param);
}

// Construct a light isoscalar nucleus for configuration and screening tests
gra::nuclear::MNucleus MakeOxygen() {
  gra::nuclear::NucleusParam param;
  param.pdg    = gra::nuclear::EncodeNuclearPDG(16, 8);
  param.mass   = 14.899;
  param.charge = {2.608, 0.513, 11.842, 64, 2048, 3.0, 4097, 1.0e-7};
  param.matter = param.charge;
  return gra::nuclear::MNucleus(param);
}

// Construct a neutron-rich nucleus with distinct charge and matter radii
gra::nuclear::MNucleus MakeNeutronSkin() {
  gra::nuclear::NucleusParam param;
  param.pdg    = gra::nuclear::EncodeNuclearPDG(200, 40);
  param.mass   = 186.0;
  param.charge = {5.4, 0.45, 13.5, 96, 2048, 3.0, 4097, 1.0e-7};
  param.matter = {6.0, 0.55, 15.9, 96, 2048, 3.0, 4097, 1.0e-7};
  return gra::nuclear::MNucleus(param);
}

// Compute complete compact configuration controls for unit tests
gra::nuclear::ConfigParam ConfigControls(const std::size_t count = 8) {
  gra::nuclear::ConfigParam param;
  param.count             = count;
  param.max_trials        = 10000;
  param.sweeps            = 32;
  param.neutron_nodes     = 4096;
  param.density_rel_tol   = 1.0e-6;
  param.negative_norm_tol = 1.0e-7;
  param.norm_tol          = 2.0e-5;
  return param;
}

// Sample one reproducible configuration bank with an explicit test RNG
gra::nuclear::MConfigBank SampleBank(const gra::nuclear::MNucleus &nucleus, const gra::nuclear::ConfigParam &param,
                                     const std::uint32_t seed) {
  const gra::nuclear::MConfigSampler sampler(nucleus, param);
  gra::MRandom                       random;
  random.SetSeed(seed);
  return gra::nuclear::MConfigBank(sampler, random);
}

// Sample one shared reproducible configuration bank with an explicit test RNG
std::shared_ptr<const gra::nuclear::MConfigBank> SampleBankPtr(const gra::nuclear::MNucleus    &nucleus,
                                                               const gra::nuclear::ConfigParam &param,
                                                               const std::uint32_t              seed) {
  return std::make_shared<const gra::nuclear::MConfigBank>(SampleBank(nucleus, param, seed));
}

// Compute a nuclear charge-current bank or the unit source of a nonnuclear leg
std::vector<std::complex<double>> ChargeBank(const gra::nuclear::MConfigBank *bank, const gra::M3Vec &q) {
  if (bank == nullptr) { return {1.0}; }
  std::vector<std::complex<double>> current(bank->Size());
  for (const auto &i : indices(current)) { current[i] = bank->At(i).ChargeCurrent(q[0], q[1], q[2]); }
  return current;
}

// Integrate a smooth screening transform with an independent fixed polar rule
std::vector<gra::nuclear::LoopNode> FixedNodes(const gra::nuclear::MUPC &upc) {
  const auto                          radial  = gra::math::PolarRadialRule(upc.Param().loop);
  const auto                          azimuth = gra::math::PolarAzimuthRule(upc.Param().loop);
  std::vector<gra::nuclear::LoopNode> nodes;
  for (const auto &i : indices(radial.node)) {
    const double k      = radial.node[i];
    const double weight = -radial.measure[i] * upc.ProfileTransform(k) / (4.0 * gra::math::PIPI);
    for (const auto &j : indices(azimuth.node)) {
      const double phi = 2.0 * gra::math::PI * (j + 0.5) / azimuth.node.size();
      nodes.push_back({k, phi, k * std::cos(phi), k * std::sin(phi), weight * azimuth.measure[j], {}});
    }
  }
  return nodes;
}

// Compute complete compact photonuclear attenuation controls
gra::nuclear::PhotoParam PhotoControls() {
  gra::nuclear::PhotoParam param;
  gra::test::Photo(param);
  param.b_nodes = 64;
  param.z_nodes = 64;
  return param;
}

// Compute complete compact UPC controls for unit tests
gra::nuclear::UPCParam UPCControls() {
  gra::nuclear::UPCParam param;
  param.structure               = gra::nuclear::StructureType::Nucleon;
  param.loop.radial_integrator  = "GL";
  param.loop.azimuth_integrator = "Trap";
  param.loop.radial_map         = gra::math::RadialMap::Square;
  gra::test::Glauber(param.glauber, 0.0);
  param.glauber.b_max              = 20.0;
  param.glauber.q_max              = 3.0;
  param.glauber.b_nodes            = 64;
  param.glauber.q_nodes            = 64;
  param.glauber.profile_nodes      = 256;
  param.convolution.b_max          = 20.0;
  param.convolution.smooth_b_nodes = 64;
  param.convolution.sample_b_nodes = 64;
  gra::test::Convolution(param.convolution);
  param.config                    = ConfigControls();
  for (auto &photo : param.photo) {
    gra::test::Photo(photo);
    photo.b_nodes = 32;
    photo.z_nodes = 32;
  }
  param.loop.r_min              = 0.0;
  param.loop.r_max              = 0.2;
  param.loop.radial_intervals   = 8;
  param.loop.azimuth_nodes      = 4;
  param.convolution.b_phi_nodes = 5;
  return param;
}

// Compute the direct ordered current used before compact proton indexing
std::complex<double> DirectConfigCurrent(const gra::nuclear::MConfig &config, const gra::M3Vec &q,
                                         const bool charge_only) {
  constexpr double     hbarc   = gra::PDG::GeV2fm;
  std::complex<double> current = 0.0;
  for (const auto &nucleon : config.Nucleons()) {
    if (charge_only && nucleon.type != gra::nuclear::NucleonType::Proton) { continue; }
    current += std::polar(1.0, gra::BilinearProduct(q, nucleon.x) / hbarc);
  }
  return current;
}

// Construct synthetic electromagnetic breakup controls for response and quadrature checks
gra::nuclear::BreakupParam BreakupControls(const unsigned int a, const unsigned int z, const double z_emit,
                                           const double gamma, const unsigned int omega_nodes = 96) {
  gra::nuclear::BreakupParam param;
  param.a       = a;
  param.z       = z;
  param.z_emit  = z_emit;
  param.gamma   = gamma;
  param.emitter = gra::nuclear::EmitterType::Point;
  gra::test::BreakupNumerics(param);
  param.photo.energy_min      = 6.0;
  param.photo.gdr.systematics = {27.47, 22.063, 0.0277, 1.9, 60.0, 1.222};
  param.photo.gdr.isotope     = {{208, 82, 13.373, 3.938, 1.33716}};
  param.photo.qd.deuteron     = {2.224, 61.2};
  param.photo.qd.levinger     = 6.5;
  param.photo.qd.pauli        = {20.0, 140.0, 73.3, {8.3714e-2, -9.8343e-3, 4.1222e-4, -3.4762e-6, 9.3537e-9}, 24.2};
  param.photo.resonance.threshold = 100.0;
  param.photo.resonance.isotope   = {{208, 82, 350.0, 110.0, 23000.0}};
  param.photo.continuum           = {a, z, 400.0, 600.0, 80000.0, 15.2, 0.06, 0.09, 100.0, 580.0};
  param.transfer.E0               = 50.0;
  param.transfer.match            = 200.0;
  param.response                  = {omega_nodes, 1.0e-2, 1.0e-3};
  param.profile.nodes             = 2048;
  param.profile.b_min             = 1.0e-6;
  return param;
}

// Compute the sampled charge and matter form-factor root-mean-square error
double ConfigFormError(const gra::nuclear::MConfigBank &bank) {
  const auto                     &nucleus  = bank.Nucleus();
  constexpr std::array<double, 5> momentum = {0.025, 0.050, 0.075, 0.100, 0.125};
  double                          error    = 0.0;
  for (const double q : momentum) {
    const std::complex<double> charge =
        (bank.ChargeStat(q, 0.0, 0.0).mean + bank.ChargeStat(0.0, q, 0.0).mean + bank.ChargeStat(0.0, 0.0, q).mean) /
        (3.0 * static_cast<double>(nucleus.Z()));
    const std::complex<double> matter =
        (bank.MatterStat(q, 0.0, 0.0).mean + bank.MatterStat(0.0, q, 0.0).mean + bank.MatterStat(0.0, 0.0, q).mean) /
        (3.0 * static_cast<double>(nucleus.A()));
    error += std::norm(charge - nucleus.ChargeDensity().Form(q));
    error += std::norm(matter - nucleus.MatterDensity().Form(q));
  }
  return std::sqrt(error / (2.0 * static_cast<double>(momentum.size())));
}

// Integrate one scaled piecewise-linear nucleon probability in polar space
double NNCrossSection(const gra::nuclear::MGlauber &glauber, const double sigma) {
  const auto      &profile       = glauber.Param().profile;
  const double     scale         = std::sqrt(sigma / profile.sigma);
  constexpr double inverse_sqrt3 = 0.57735026918962576451;
  double           integral      = 0.0;
  for (std::size_t i = 1; i < profile.b_node.size(); ++i) {
    const double lower      = scale * profile.b_node[i - 1];
    const double upper      = scale * profile.b_node[i];
    const double center     = 0.5 * (lower + upper);
    const double half_width = 0.5 * (upper - lower);
    for (const double sign : {-1.0, 1.0}) {
      const double b = center + sign * half_width * inverse_sqrt3;
      integral += half_width * b * (1.0 - gra::math::pow2(glauber.NNAmp(b, sigma)));
    }
  }
  return 2.0 * gra::math::PI * integral;
}

// Integrate one normalized transverse thickness over the plane
double ThicknessNorm(const gra::nuclear::MDensity &density) {
  const auto rule     = gra::math::GaussLegendreRule(96, 0.0, density.Param().r_max);
  double     integral = 0.0;
  for (const auto &i : indices(rule.first)) {
    integral += rule.second[i] * rule.first[i] * density.Thick(rule.first[i]);
  }
  return 2.0 * gra::math::PI * integral;
}

// Compute compact numerical controls for Glauber unit tests
gra::nuclear::GlauberParam GlauberControls(const double omega) {
  gra::nuclear::GlauberParam param;
  gra::test::Glauber(param, omega);
  param.b_max         = 20.0;
  param.q_max         = 3.0;
  param.b_nodes       = 64;
  param.q_nodes       = 64;
  param.profile_nodes = 256;
  return param;
}

// Require one real probability to be physical within roundoff
void RequireProbability(const double probability) {
  REQUIRE(std::isfinite(probability));
  REQUIRE(probability >= -1.0e-14);
  REQUIRE(probability <= 1.0 + 1.0e-14);
}

// Require normalized coherent means and incoherent variances
void RequireCurrentMoments(const gra::nuclear::MUPC::SectorRatios &ratio, const std::size_t count) {
  REQUIRE(ratio[0].size() == count);
  REQUIRE(ratio[1].size() == count);
  std::complex<double> coherent_mean     = 0.0;
  std::complex<double> incoherent_mean   = 0.0;
  double               incoherent_second = 0.0;
  for (std::size_t sample = 0; sample < count; ++sample) {
    coherent_mean += ratio[0][sample];
    incoherent_mean += ratio[1][sample];
    incoherent_second += std::norm(ratio[1][sample]);
  }
  coherent_mean /= static_cast<double>(count);
  incoherent_mean /= static_cast<double>(count);
  incoherent_second /= static_cast<double>(count);
  REQUIRE(coherent_mean.real() == Approx(1.0).epsilon(2.0e-12));
  REQUIRE(coherent_mean.imag() == Approx(0.0).margin(2.0e-12));
  REQUIRE(std::abs(incoherent_mean) < 2.0e-12);
  REQUIRE(incoherent_second == Approx(1.0).epsilon(2.0e-12));
}

// Store survival and screened Good-Walker observables for convergence tests
struct UPCMoment {
  double                survival = 0.0;
  std::array<double, 4> sector{};
};

// Store the resolved Good-Walker survival sectors at fixed impact parameter
struct SurvivalMoment {
  gra::nuclear::NuclearGoodWalkerWeights sector{};
  double                                 probability = 0.0;
  double                                 mean_second = 0.0;
};

// Store one-body and pair observables of a sampled nuclear geometry
struct GeometryMoment {
  double radius2 = 0.0;
  double pair12  = 0.0;
  double pair20  = 0.0;
};

// Compute normalized geometry observables from one configuration bank
GeometryMoment ConfigGeometry(const gra::nuclear::MConfigBank &bank) {
  GeometryMoment moment;
  double         nucleon_count = 0.0;
  double         pair_count    = 0.0;
  for (std::size_t sample = 0; sample < bank.Size(); ++sample) {
    const auto &nucleon = bank.At(sample).Nucleons();
    for (const auto &i : indices(nucleon)) {
      moment.radius2 += gra::SquaredNorm(nucleon[i].x);
      nucleon_count += 1.0;
      for (std::size_t j = i + 1; j < nucleon.size(); ++j) {
        const double distance2 = gra::SquaredNorm(gra::Subtract(nucleon[i].x, nucleon[j].x));
        moment.pair12 += distance2 < 1.44 ? 1.0 : 0.0;
        moment.pair20 += distance2 < 4.00 ? 1.0 : 0.0;
        pair_count += 1.0;
      }
    }
  }
  moment.radius2 /= nucleon_count;
  moment.pair12 /= pair_count;
  moment.pair20 /= pair_count;
  return moment;
}

// Project one deterministic Cartesian survival amplitude ensemble
SurvivalMoment SurvivalSectors(const gra::nuclear::MGlauber &glauber, const gra::nuclear::MNucleus &nucleus,
                               const std::size_t count, const std::uint32_t seed, const double b) {
  auto upper = ConfigControls(count);
  // Match the production minimum nucleon separation
  upper.d_min                                  = 0.8;
  auto                              lower      = upper;
  const auto                        upper_bank = SampleBank(nucleus, upper, seed);
  const auto                        lower_bank = SampleBank(nucleus, lower, seed ^ 0x9e3779b9U);
  std::vector<std::complex<double>> amplitude(count * count, 0.0);
  SurvivalMoment                    moment;
  for (const auto &sample : indices(amplitude)) {
    amplitude[sample] = glauber.SampleAmp(b, &upper_bank, &lower_bank, sample);
    moment.mean_second += std::norm(amplitude[sample]);
  }
  moment.mean_second /= static_cast<double>(amplitude.size());
  moment.sector      = gra::nuclear::NuclearGoodWalkerProject(amplitude, {count, count});
  moment.probability = glauber.ConfigProb(b, &upper_bank, &lower_bank);
  return moment;
}

// Construct a configuration UPC model with independent loop and transform quadratures
gra::nuclear::MUPC ConfigUPC(const std::size_t count, const std::uint32_t seed, const unsigned int k_nodes,
                             const unsigned int phi_nodes, const unsigned int transform_nodes = 256,
                             const unsigned int impact_nodes = 48, const unsigned int impact_angles = 5) {
  auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto upper  = ConfigControls(count);
  // Match the production minimum nucleon separation
  upper.d_min           = 0.8;
  auto       lower      = upper;
  const auto upper_bank = SampleBankPtr(*oxygen, upper, seed);
  const auto lower_bank = SampleBankPtr(*oxygen, lower, seed ^ 0x9e3779b9U);

  auto param                        = UPCControls();
  param.survival                    = gra::nuclear::SurvivalType::MCGGCF;
  param.config                      = upper;
  param.emission                    = {gra::nuclear::CoherenceType::Inclusive, gra::nuclear::CoherenceType::Inclusive};
  param.glauber                     = GlauberControls(0.15);
  param.glauber.b_nodes             = 48;
  param.convolution.smooth_b_nodes  = 48;
  param.convolution.sample_b_nodes  = impact_nodes;
  param.convolution.b_phi_nodes = impact_angles;
  param.convolution.smooth_kt_nodes = transform_nodes;
  param.convolution.sample_kt_nodes = transform_nodes;
  param.loop.r_min                  = 0.0;
  param.loop.r_max                  = 2.0;
  param.loop.radial_intervals       = k_nodes;
  param.loop.azimuth_nodes          = phi_nodes;
  return gra::nuclear::MUPC({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus}, {oxygen, oxygen}, param,
                            {upper_bank, lower_bank});
}

// Compute the exact coherent mean and orthogonal fluctuation norm of one charge current
std::array<std::complex<double>, 2> ChargeSectorAmps(const gra::nuclear::MUPC &upc, const int leg,
                                                     const gra::M3Vec &transfer) {
  const auto *bank = upc.Bank(leg);
  if (bank == nullptr || bank->Size() == 0) { throw std::logic_error("ChargeSectorAmps: missing configuration bank"); }
  std::vector<std::complex<double>> current(bank->Size(), 0.0);
  for (const auto &sample : indices(current)) {
    current[sample] = bank->At(sample).ChargeCurrent(transfer[0], transfer[1], transfer[2]);
  }
  const auto stat = gra::statistics::ComplexMoments(current);
  if (!(stat.variance >= 0.0) || !std::isfinite(stat.variance)) {
    throw std::runtime_error("ChargeSectorAmps: invalid charge-current variance");
  }
  return {stat.mean, std::sqrt(stat.variance)};
}

// Compute a smooth hard amplitude with an exact Good-Walker current decomposition
std::vector<std::complex<double>> FusionAmplitude(const gra::nuclear::MUPC        &upc,
                                                  const std::array<gra::M3Vec, 2> &transfer) {
  constexpr double slope = 6.0;
  const double     norm2 = gra::SquaredNorm(transfer[0]) + gra::SquaredNorm(transfer[1]);
  const auto       hard  = std::polar(std::exp(-slope * norm2), 0.7 * transfer[0][0] - 0.4 * transfer[1][1]);
  const auto       upper = ChargeSectorAmps(upc, 1, transfer[0]);
  const auto       lower = ChargeSectorAmps(upc, 2, transfer[1]);
  std::vector<std::complex<double>> amplitude;
  amplitude.reserve(4);
  for (const auto source1 : upper) {
    for (const auto source2 : lower) { amplitude.push_back(hard * source1 * source2); }
  }
  return amplitude;
}

// Compute one circular transverse component v_lambda = v_x + i lambda v_y
std::complex<double> Circular(const gra::M3Vec &vector, const int helicity) {
  return {vector[0], static_cast<double>(helicity) * vector[1]};
}

// Compute the exact impact-space closure of two displaced Gaussian EPA currents
std::array<std::complex<double>, 2> GaussianCurrentClosure(const gra::nuclear::MUPC        &upc,
                                                           const std::array<gra::M3Vec, 2> &transfer,
                                                           const std::array<double, 2> &slope, const int helicity) {
  constexpr double hbarc    = gra::PDG::GeV2fm;
  const double     total    = slope[0] + slope[1];
  const double     cx       = (slope[1] * transfer[1][0] - slope[0] * transfer[0][0]) / total;
  const double     cy       = (slope[1] * transfer[1][1] - slope[0] * transfer[0][1]) / total;
  const double     c        = std::hypot(cx, cy);
  const double     q1t2     = gra::math::pow2(transfer[0][0]) + gra::math::pow2(transfer[0][1]);
  const double     q2t2     = gra::math::pow2(transfer[1][0]) + gra::math::pow2(transfer[1][1]);
  const double     exponent = -slope[0] * q1t2 - slope[1] * q2t2 + total * c * c;
  const double     norm     = std::exp(exponent);

  // G(k) = norm exp[-total |k-c|^2]
  // M0 = int d^2b S(b) g(b) / hbarc^2
  // M1_lambda is the transform of (kx + i lambda ky) G(k)
  double       scalar         = 0.0;
  double       vector         = 0.0;
  const double gaussian_limit = 2.0 * hbarc * std::sqrt(45.0 * total);
  const double b_max          = upc.Param().convolution.b_max;
  // Integrate the tabulated survival region and its exact unit-survival tail separately
  const auto integrate = [&](const double lower, const double upper, const bool screened, const unsigned int nodes) {
    if (!(upper > lower)) { return; }
    const auto [b_node, b_weight] = gra::math::GaussLegendreRule(nodes, lower, upper);
    for (const auto &i : indices(b_node)) {
      const double b        = b_node[i];
      const double survival = screened ? upc.SurvivalAmp(b) : 1.0;
      const double radial   = b_weight[i] * b * survival * std::exp(-b * b / (4.0 * total * hbarc * hbarc));
      const double phase    = c * b / hbarc;
      scalar += radial * std::cyl_bessel_j(0.0, phase);
      vector +=
          radial * (c * std::cyl_bessel_j(0.0, phase) - b * std::cyl_bessel_j(1.0, phase) / (2.0 * total * hbarc));
    }
  };
  integrate(0.0, b_max, true, 384);
  integrate(b_max, std::max(b_max, gaussian_limit), false, 128);

  const double               prefactor = norm / (2.0 * total * hbarc * hbarc);
  const std::complex<double> direction =
      c > 0.0 ? std::complex<double>(cx, static_cast<double>(helicity) * cy) / c : std::complex<double>(0.0, 0.0);
  const std::complex<double> moment0 = prefactor * scalar;
  const std::complex<double> moment1 = prefactor * vector * direction;
  // A1_scr = q1_lambda M0 + M1_lambda and A2_scr = q2_lambda M0 - M1_lambda
  return {Circular(transfer[0], helicity) * moment0 + moment1, Circular(transfer[1], helicity) * moment0 - moment1};
}

// Evaluate displaced Gaussian transverse currents through the public screening API
std::vector<std::complex<double>> GaussianCurrentScreen(const gra::nuclear::MUPC        &upc,
                                                        const std::array<gra::M3Vec, 2> &transfer,
                                                        const std::array<double, 2>     &slope) {
  const auto amplitude = [&](const std::array<gra::M3Vec, 2> &q) {
    const double form = std::exp(-slope[0] * (gra::math::pow2(q[0][0]) + gra::math::pow2(q[0][1])) -
                                 slope[1] * (gra::math::pow2(q[1][0]) + gra::math::pow2(q[1][1])));
    return std::vector<std::complex<double>>{form * Circular(q[0], 1), form * Circular(q[1], 1),
                                             form * Circular(q[0], -1), form * Circular(q[1], -1)};
  };

  gra::nuclear::ScreenLayout layout;
  layout.type        = gra::nuclear::ScreenType::Scalar;
  layout.fixed.valid = true;
  gra::nuclear::ScreenPoint point;
  point.transfer  = transfer;
  point.amplitude = amplitude(transfer);
  gra::nuclear::MUPCScreen screen(upc, layout, point);
  for (const auto &node : upc.Nodes(transfer)) {
    auto shifted = transfer;
    shifted[0][0] += node.kx;
    shifted[0][1] += node.ky;
    shifted[1][0] -= node.kx;
    shifted[1][1] -= node.ky;
    point.transfer  = shifted;
    point.amplitude = amplitude(shifted);
    screen.Add(node, point);
  }
  return screen.Result().amplitude;
}

// Evaluate both Gaussian EPA currents with their exact propagator denominators
std::vector<std::complex<double>> GaussianEPAScreen(const gra::nuclear::MUPC        &upc,
                                                    const std::array<gra::M3Vec, 2> &transfer,
                                                    const std::array<double, 2>     &slope) {
  // J1_lambda J2_lambda' = q1_lambda q2_lambda' exp(-a1 q1T^2-a2 q2T^2)/(Q1^2 Q2^2)
  const auto amplitude = [&](const std::array<gra::M3Vec, 2> &q) {
    const double qt1  = gra::math::pow2(q[0][0]) + gra::math::pow2(q[0][1]);
    const double qt2  = gra::math::pow2(q[1][0]) + gra::math::pow2(q[1][1]);
    const double den  = (qt1 + gra::math::pow2(q[0][2])) * (qt2 + gra::math::pow2(q[1][2]));
    const double form = std::exp(-slope[0] * qt1 - slope[1] * qt2) / den;
    return std::vector<std::complex<double>>{
        form * Circular(q[0], 1) * Circular(q[1], 1), form * Circular(q[0], 1) * Circular(q[1], -1),
        form * Circular(q[0], -1) * Circular(q[1], 1), form * Circular(q[0], -1) * Circular(q[1], -1)};
  };

  gra::nuclear::ScreenLayout layout;
  layout.type        = gra::nuclear::ScreenType::Scalar;
  layout.fixed.valid = true;
  gra::nuclear::ScreenPoint point;
  point.transfer  = transfer;
  point.amplitude = amplitude(transfer);
  gra::nuclear::MUPCScreen screen(upc, layout, point);
  for (const auto &node : upc.Nodes(transfer)) {
    auto shifted = transfer;
    shifted[0][0] += node.kx;
    shifted[0][1] += node.ky;
    shifted[1][0] -= node.kx;
    shifted[1][1] -= node.ky;
    point.transfer  = shifted;
    point.amplitude = amplitude(shifted);
    screen.Add(node, point);
  }
  return screen.Result().amplitude;
}

// Evaluate real survival and screened Good-Walker moments
UPCMoment ScreenMoment(const gra::nuclear::MUPC &upc) {
  gra::nuclear::ScreenLayout layout;
  layout.type             = gra::nuclear::ScreenType::Fusion;
  layout.fusion.rows      = 1;
  layout.fusion.sector[0] = gra::nuclear::CoherenceSectors(upc.Param().emission[0]);
  layout.fusion.sector[1] = gra::nuclear::CoherenceSectors(upc.Param().emission[1]);

  gra::nuclear::ScreenPoint born;
  born.transfer  = {{{0.08, -0.05, 0.02}, {-0.06, 0.04, -0.015}}};
  born.amplitude = FusionAmplitude(upc, born.transfer);
  gra::nuclear::MUPCScreen screen(upc, layout, born);
  for (const auto &node : upc.Nodes(born.transfer)) {
    auto point = born;
    point.transfer[0][0] += node.kx;
    point.transfer[0][1] += node.ky;
    point.transfer[1][0] -= node.kx;
    point.transfer[1][1] -= node.ky;
    point.amplitude = FusionAmplitude(upc, point.transfer);
    screen.Add(node, point);
  }
  const auto result = screen.Result();
  UPCMoment  moment;
  if (result.helicity_norm.size() != moment.sector.size()) {
    throw std::runtime_error("ScreenMoment: invalid fusion sector count");
  }
  moment.survival = gra::math::pow2(upc.SurvivalAmp(6.0));
  std::copy(result.helicity_norm.cbegin(), result.helicity_norm.cend(), moment.sector.begin());
  return moment;
}

// Compute the relative RMS spread of one seed ensemble
double SeedSpread(const std::vector<double> &value) {
  double mean = 0.0;
  for (const double item : value) { mean += item; }
  mean /= static_cast<double>(value.size());
  double variance = 0.0;
  for (const double item : value) { variance += gra::math::pow2(item - mean); }
  variance /= static_cast<double>(value.size());
  return std::sqrt(variance) / std::abs(mean);
}

// Compute the configuration ensemble RMS error against a nested reference
double EnsembleError(const std::vector<double> &value, const std::vector<double> &reference) {
  double mean_reference = 0.0;
  double error          = 0.0;
  for (const auto &i : indices(value)) {
    mean_reference += reference[i];
    error += gra::math::pow2(value[i] - reference[i]);
  }
  mean_reference /= static_cast<double>(value.size());
  error /= static_cast<double>(value.size());
  return std::sqrt(error) / std::abs(mean_reference);
}

// Compute the relative L1 distance between two resolved screening results
double ScreenDistance(const UPCMoment &value, const UPCMoment &reference) {
  double distance = 0.0;
  double norm     = 0.0;
  for (const auto &i : indices(value.sector)) {
    distance += std::abs(value.sector[i] - reference.sector[i]);
    norm += reference.sector[i];
  }
  return distance / norm;
}

// Require finite physical survival and resolved Good-Walker probabilities
void RequireScreenMoment(const UPCMoment &moment) {
  RequireProbability(moment.survival);
  for (const double sector : moment.sector) {
    REQUIRE(std::isfinite(sector));
    REQUIRE(sector >= 0.0);
  }
  REQUIRE(gra::Sum(moment.sector) > 0.0);
}

}  // namespace

TEST_CASE("Nuclear PDG identities round trip and reject invalid encodings", "[gra::nuclear][identity]") {
  const int encoded = gra::nuclear::EncodeNuclearPDG(208, 82);
  REQUIRE(encoded == pb208_pdg);
  REQUIRE(gra::nuclear::IsNuclearPDG(encoded));

  const auto lead = gra::nuclear::DecodeNuclearPDG(encoded);
  REQUIRE(lead.a == 208);
  REQUIRE(lead.z == 82);
  REQUIRE(lead.lambda == 0);
  REQUIRE(lead.isomer == 0);
  REQUIRE_FALSE(lead.anti);

  const int  anti_code = gra::nuclear::EncodeNuclearPDG(208, 82, 0, 0, true);
  const auto anti_lead = gra::nuclear::DecodeNuclearPDG(anti_code);
  REQUIRE(anti_lead.anti);
  REQUIRE(anti_lead.a == lead.a);
  REQUIRE(anti_lead.z == lead.z);

  REQUIRE_FALSE(gra::nuclear::IsNuclearPDG(2212));
  REQUIRE_FALSE(gra::nuclear::IsNuclearPDG(1000830820));
  REQUIRE_FALSE(gra::nuclear::IsNuclearPDG(1000820000));
  REQUIRE_THROWS_AS(gra::nuclear::DecodeNuclearPDG(1000830820), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::nuclear::EncodeNuclearPDG(82, 83), std::invalid_argument);

  const auto nucleus = MakeLead();
  REQUIRE(nucleus.A() == 208);
  REQUIRE(nucleus.Z() == 82);
  REQUIRE(nucleus.N() == 126);
  REQUIRE(nucleus.Charge() == 82);
  REQUIRE(nucleus.Mass() == Approx(pb208_mass).epsilon(1.0e-14));
}

TEST_CASE("Lead charge and matter densities have normalized low-q limits", "[gra::nuclear][density]") {
  const auto  lead   = MakeLead();
  const auto &charge = lead.ChargeDensity();
  const auto &matter = lead.MatterDensity();

  REQUIRE(charge.Form(0.0) == Approx(1.0).margin(2.0e-14));
  REQUIRE(matter.Form(0.0) == Approx(1.0).margin(2.0e-14));
  REQUIRE(ThicknessNorm(charge) == Approx(1.0).epsilon(2.0e-5));
  REQUIRE(ThicknessNorm(matter) == Approx(1.0).epsilon(2.0e-5));
  REQUIRE(charge.Rms() > 5.0);
  REQUIRE(charge.Rms() < 6.5);

  constexpr double hbarc = gra::PDG::GeV2fm;
  constexpr double q     = 1.0e-3;
  const double     low_q = 1.0 - q * q * charge.Rms() * charge.Rms() / (6.0 * hbarc * hbarc);
  REQUIRE(charge.Form(q) == Approx(low_q).margin(2.0e-7));
  REQUIRE(std::abs(charge.Form(0.2)) < 0.2);
  REQUIRE(charge.Rho(0.0) > charge.Rho(charge.Param().radius));

  // Resolve C(b) / (pi b^2) approaching the central transverse thickness
  for (const double b : {1.0e-6, 1.0e-8, 1.0e-10}) {
    REQUIRE(charge.Cylinder(b) / (gra::math::PI * b * b) == Approx(charge.Thick(0.0)).epsilon(2.0e-8));
  }
}

// Check the tuned table against a finer grid, including signed diffraction minima
TEST_CASE("Lead form factors resolve the oscillatory high-q tail", "[gra::nuclear][density][form-factor]") {
  const auto numerics =
      nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("modeldata/TUNE0/NUMERICS.json")))
          .at("NUMERICS_NUCLEAR")
          .at("STRUCTURE");
  auto param                = MakeLead().Param();
  param.charge.r_max        = param.charge.radius + 16.0 * param.charge.skin;
  param.charge.form_q_max   = numerics.at("form_q_max").get<double>();
  param.charge.form_nodes   = numerics.at("form_nodes").get<std::size_t>();
  param.charge.form_abs_tol = numerics.at("form_abs_tol").get<double>();
  param.matter              = param.charge;
  const gra::nuclear::MNucleus lead(param);
  param.charge.form_nodes = 4 * (param.charge.form_nodes - 1) + 1;
  param.charge.form_abs_tol *= 0.01;
  param.matter = param.charge;
  const gra::nuclear::MNucleus reference(param);

  const auto &density = lead.ChargeDensity();
  REQUIRE(density.Form(4.0) == Approx(1.4754e-11).margin(2.0e-13));
  REQUIRE(density.Form(8.0) == Approx(8.7524e-12).margin(2.0e-13));
  double max_error = 0.0;
  for (const double q : gra::math::linspace(0.0, 0.5, 2001)) {
    max_error = std::max(max_error, std::abs(density.Form(q) - reference.ChargeDensity().Form(q)));
  }
  REQUIRE(max_error < 1.0e-7);
  const gra::nuclear::MPhoton photon(lead), fine(reference);
  for (const double q : {0.0013, 0.05, 0.12, 0.131, 0.132, 0.133, 0.18, 0.3, 0.8, 2.1}) {
    const double source = photon.CoherentAmp(1.0e-4, -q * q, q);
    const double truth  = fine.CoherentAmp(1.0e-4, -q * q, q);
    for (const int helicity : {-1, 1}) {
      const auto amplitude = gra::nuclear::PhotonAmp(source, 1.0e-4, lead.Z(), 0.37, helicity);
      const auto expected  = gra::nuclear::PhotonAmp(truth, 1.0e-4, lead.Z(), 0.37, helicity);
      REQUIRE(std::abs(amplitude - expected) < 1.0e-6 + 1.0e-4 * std::abs(expected));
    }
  }
}

// Check the finite-density Fourier integral across the cache and large recoil momenta
TEST_CASE("Nuclear form factors cover large recoil momenta", "[gra::nuclear][density][form-factor]") {
  auto param         = MakeOxygen().Param().charge;
  param.form_nodes   = 16385;
  param.form_abs_tol = 1.0e-10;
  const auto rule    = gra::math::GaussLegendreRule(128, 0.0, 1.0);
  for (const double radius : {param.r_max, param.radius + param.skin}) {
    param.r_max = radius;
    const gra::nuclear::MDensity density(param);
    for (const double q : {2.9999, 3.0001, 9.0, 20.0, 70.0, 200.0}) {
      double       integral = 0.0;
      const double step     = radius / 128.0;
      for (std::size_t j = 0; j < 128; ++j) {
        for (const auto &i : indices(rule.first)) {
          const double r = (j + rule.first[i]) * step, phase = q * r / gra::PDG::GeV2fm;
          integral += step * rule.second[i] * r * r * density.Rho(r) * std::sin(phase) / phase;
        }
      }
      const double expected = 4.0 * gra::math::PI * integral;
      CAPTURE(radius, q, expected);
      CHECK(density.Form(q) == Approx(expected).margin(param.form_abs_tol));
      CHECK(density.Form(-q) == Approx(expected).margin(param.form_abs_tol));
    }
  }
}

// Check the complex impulse amplitude and incoherent limit under rotations and beam reversal
TEST_CASE("Impulse photoproduction covers large recoil momenta", "[gra::nuclear][photo][form-factor]") {
  using namespace gra::nuclear;
  const auto   nucleus = MakeLead();
  const MPhoto photo(nucleus, PhotoControls(), PhotoModel::Impulse);
  const auto   profile = gra::test::PhotoProfile(0.0, 0.0);
  for (const double q : {70.0, 200.0, 1.0e4}) {
    const double form = nucleus.MatterDensity().Form(q);
    for (const auto &transfer : std::array<gra::M3Vec, 3>{{{q, 0.0, 0.0}, {0.0, q, 0.0}, {0.0, 0.0, -q}}}) {
      for (const auto direction : {PhotonDirection::PositiveZ, PhotonDirection::NegativeZ}) {
        const auto value = photo.Factors(profile, transfer[0], transfer[1], transfer[2], direction);
        CHECK(std::abs(value.coherent - std::complex<double>(nucleus.A() * form, 0.0)) < 1.0e-12);
        CHECK(value.incoherent * value.incoherent == Approx(nucleus.A()).epsilon(1.0e-10));
      }
    }
  }
}

TEST_CASE("Configuration currents obey the Good-Walker moment identity", "[gra::nuclear][config][good-walker]") {
  const auto oxygen   = MakeOxygen();
  auto       param    = ConfigControls(32);
  param.d_min         = 0.35;
  param.max_trials    = 20000;
  const auto bank     = SampleBank(oxygen, param, 1701);
  const auto repeated = SampleBank(oxygen, param, 1701);

  REQUIRE(bank.Size() == param.count);
  REQUIRE(bank.Param().sweeps == param.sweeps);
  for (std::size_t sample = 0; sample < bank.Size(); ++sample) {
    const auto &nucleon = bank.At(sample).Nucleons();
    for (std::size_t i = 0; i < nucleon.size(); ++i) {
      for (std::size_t j = i + 1; j < nucleon.size(); ++j) {
        REQUIRE(std::sqrt(gra::SquaredNorm(gra::Subtract(nucleon[i].x, nucleon[j].x))) + 2.0e-14 >= param.d_min);
      }
    }
  }
  const auto origin = bank.ChargeStat(0.0, 0.0);
  REQUIRE(origin.mean.real() == Approx(8.0).margin(2.0e-14));
  REQUIRE(std::abs(origin.mean.imag()) < 2.0e-14);
  REQUIRE(origin.second == Approx(64.0).margin(2.0e-13));
  REQUIRE(origin.variance < 2.0e-13);

  const auto finite_q = bank.ChargeStat(0.18, -0.07);
  REQUIRE(finite_q.second == Approx(std::norm(finite_q.mean) + finite_q.variance).margin(2.0e-13));
  REQUIRE(finite_q.variance > 0.0);
  REQUIRE(finite_q.variance < finite_q.second);

  const auto deterministic = repeated.ChargeStat(0.18, -0.07);
  REQUIRE(deterministic.mean.real() == Approx(finite_q.mean.real()).margin(2.0e-14));
  REQUIRE(deterministic.mean.imag() == Approx(finite_q.mean.imag()).margin(2.0e-14));
  REQUIRE(deterministic.variance == Approx(finite_q.variance).margin(2.0e-14));
}

TEST_CASE("UPC nuclear configurations are sampled from the event RNG", "[gra::nuclear][config][event]") {
  const auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param  = UPCControls();
  param.survival    = gra::nuclear::SurvivalType::Optical;
  param.config      = ConfigControls(4);
  const gra::nuclear::MUPC model({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus}, {nullptr, oxygen},
                                 param);
  REQUIRE_FALSE(model.HasSamples());
  REQUIRE(model.Bank(2) == nullptr);

  gra::MRandom random;
  random.SetSeed(271828);
  const auto first  = model.Sample(random);
  const auto second = model.Sample(random);
  REQUIRE(first->HasSamples());
  REQUIRE(second->HasSamples());
  REQUIRE(first->Bank(2) != nullptr);
  REQUIRE(second->Bank(2) != nullptr);
  const auto &first_center  = first->Bank(2)->At(0).Nucleons()[0].x;
  const auto &second_center = second->Bank(2)->At(0).Nucleons()[0].x;
  REQUIRE(gra::SquaredNorm(gra::Subtract(first_center, second_center)) > 1.0e-12);

  gra::MRandom repeated_random;
  repeated_random.SetSeed(271828);
  const auto  repeated        = model.Sample(repeated_random);
  const auto &repeated_center = repeated->Bank(2)->At(0).Nucleons()[0].x;
  REQUIRE(gra::SquaredNorm(gra::Subtract(first_center, repeated_center)) < 1.0e-28);
}

TEST_CASE("Compact proton currents equal the direct ordered phase sum", "[gra::nuclear][config][current]") {
  using gra::nuclear::Nucleon;
  using gra::nuclear::NucleonType;
  const gra::nuclear::MConfig     config(7, 3,
                                         std::vector<Nucleon>{{{0.4, -0.7, 1.1}, NucleonType::Neutron},
                                                              {{-0.2, 0.9, -1.4}, NucleonType::Proton},
                                                              {{1.3, 0.1, 0.8}, NucleonType::Neutron},
                                                              {{-1.1, -0.5, 0.3}, NucleonType::Proton},
                                                              {{0.7, 1.5, -0.6}, NucleonType::Neutron},
                                                              {{0.0, -1.2, 1.7}, NucleonType::Proton},
                                                              {{-0.8, 0.6, -0.9}, NucleonType::Neutron}});
  const std::array<gra::M3Vec, 4> momentum = {gra::M3Vec{0.0, 0.0, 0.0}, gra::M3Vec{0.11, -0.04, 0.07},
                                              gra::M3Vec{-0.31, 0.28, -0.19}, gra::M3Vec{1.7, -0.9, 0.4}};
  const auto                     &order    = config.TransverseOrder();
  std::vector<bool>               seen(order.size(), false);
  REQUIRE(order.size() == config.Nucleons().size());
  for (const auto &i : indices(order)) {
    REQUIRE(order[i] < order.size());
    REQUIRE_FALSE(seen[order[i]]);
    seen[order[i]] = true;
    if (i > 0) { REQUIRE(config.Nucleons()[order[i - 1]].x[0] <= config.Nucleons()[order[i]].x[0]); }
  }
  for (const auto &q : momentum) {
    const auto charge = config.ChargeCurrent(q[0], q[1], q[2]);
    const auto matter = config.MatterCurrent(q[0], q[1], q[2]);
    REQUIRE(std::abs(charge - DirectConfigCurrent(config, q, true)) < 1.0e-15);
    REQUIRE(std::abs(matter - DirectConfigCurrent(config, q, false)) < 1.0e-15);
  }
}

TEST_CASE("Nuclear currents use the full rest-frame momentum transfer", "[gra::nuclear][config][kinematics]") {
  constexpr double mass     = 14.899;
  constexpr double energy   = 100.0;
  const double     momentum = std::sqrt(energy * energy - mass * mass);
  const gra::M4Vec upper_beam(0.0, 0.0, momentum, energy);
  const gra::M4Vec lower_beam(0.0, 0.0, -momentum, energy);
  const gra::M4Vec upper_q(0.08, -0.03, 0.70, 0.65);
  const gra::M4Vec lower_q(0.08, -0.03, -0.70, 0.65);

  const gra::M3Vec upper = gra::nuclear::RestTransfer(upper_beam, upper_q, mass);
  const gra::M3Vec lower = gra::nuclear::RestTransfer(lower_beam, lower_q, mass);
  REQUIRE(upper[0] == Approx(upper_q.Px()).margin(2.0e-14));
  REQUIRE(upper[1] == Approx(upper_q.Py()).margin(2.0e-14));
  REQUIRE(lower[0] == Approx(upper[0]).margin(2.0e-14));
  REQUIRE(lower[1] == Approx(upper[1]).margin(2.0e-14));
  REQUIRE(lower[2] == Approx(-upper[2]).epsilon(2.0e-13));

  const double rest_energy    = (upper_beam * upper_q) / mass;
  const double expected_norm2 = rest_energy * rest_energy - upper_q.M2();
  REQUIRE(gra::SquaredNorm(upper) == Approx(expected_norm2).epsilon(2.0e-13));
  REQUIRE(std::abs(upper[2] - upper_q.Pz()) > 1.0e-3);
  const auto unboosted = gra::nuclear::RestTransfer(gra::M4Vec(0.0, 0.0, 0.0, mass), upper_q, mass);
  REQUIRE(unboosted[0] == Approx(upper_q.Px()).margin(2.0e-14));
  REQUIRE(unboosted[1] == Approx(upper_q.Py()).margin(2.0e-14));
  REQUIRE(unboosted[2] == Approx(upper_q.Pz()).margin(2.0e-14));
  REQUIRE_THROWS_AS(gra::nuclear::RestTransfer(gra::M4Vec(0.0, 0.0, 1.0, 1.0), upper_q, mass), std::invalid_argument);

  const auto oxygen     = MakeOxygen();
  auto       param      = ConfigControls(8);
  const auto bank       = SampleBank(oxygen, param, 1701);
  const auto transverse = bank.ChargeStat(upper[0], upper[1], 0.0);
  const auto full       = bank.ChargeStat(upper[0], upper[1], upper[2]);
  REQUIRE(std::abs(full.mean - transverse.mean) > 1.0e-6);
}

// Check transfer closure after a low-energy recoil from a high-energy ion
TEST_CASE("Nuclear transfers remove beam-subtraction roundoff", "[gra::nuclear][kinematics][closure]") {
  constexpr double proton_energy = 6500.0;
  constexpr double ion_energy    = 532480.0;
  const gra::M4Vec beam1(0.0, 0.0, std::sqrt(proton_energy * proton_energy - gra::PDG::mp * gra::PDG::mp),
                         proton_energy);
  const gra::M4Vec beam2(0.0, 0.0, -std::sqrt(ion_energy * ion_energy - pb208_mass * pb208_mass), ion_energy);
  const gra::M4Vec q1(-0.0100395, -0.0106910, 17.7710, 17.7710);
  const gra::M4Vec q2(0.00399188, 0.00726120, -0.0267177, 0.0267177);
  const std::array<gra::M4Vec, 2> beam    = {beam1, beam2};
  const std::array<gra::M4Vec, 2> forward = {beam1 - q1, beam2 - q2};
  const gra::M4Vec                central = q1 + q2;
  const std::array<gra::M4Vec, 2> raw     = {beam[0] - forward[0], beam[1] - forward[1]};

  std::array<gra::M4Vec, 2> transfer;
  REQUIRE(gra::nuclear::CloseTransfers({gra::PDG::PDG_p, pb208_pdg}, beam, forward, central, transfer));
  CHECK(gra::math::IsExactEqual(transfer[0].E(), raw[0].E()));
  CHECK(gra::math::IsExactEqual(transfer[0].Pz(), raw[0].Pz()));
  REQUIRE(gra::math::CheckEMC(transfer[0] + transfer[1] - central, 1.0e-12));

  gra::M4Vec inconsistent = central;
  inconsistent.SetE(inconsistent.E() + 1.0e-3);
  REQUIRE_FALSE(gra::nuclear::CloseTransfers({gra::PDG::PDG_p, pb208_pdg}, beam, forward, inconsistent, transfer));
}

TEST_CASE("Sampled nuclear current variance removes the finite-bank bias",
          "[gra::nuclear][config][good-walker][statistics]") {
  // A large coherent current must not erase the centered incoherent fluctuations
  for (const double offset : {0.0, 1073741824.0}) {
    for (const std::complex<double> phase : {std::complex<double>{1.0, 0.0}, {0.0, 1.0}}) {
      std::vector<std::complex<double>> current = {{1.0, 0.0}, {-1.0, 0.0}, {0.0, 2.0}, {0.0, -2.0}};
      const std::complex<double>        mean    = phase * std::complex<double>{offset, -offset};
      for (auto &value : current) { value = phase * value + mean; }
      const auto stat = gra::nuclear::CurrentMoments(current);
      REQUIRE(std::abs(stat.mean - mean) < 1.0e-14);
      REQUIRE(stat.variance == Approx(10.0 / 3.0).epsilon(1.0e-12));
      REQUIRE(stat.second == Approx(std::norm(mean) + stat.variance));
    }
  }
  REQUIRE(gra::nuclear::CurrentMoments({{1.0, 2.0}}).variance == Approx(0.0).margin(1.0e-15));
  REQUIRE_THROWS_AS(gra::nuclear::CurrentMoments({}), std::invalid_argument);
}

TEST_CASE("Configuration sampling reproduces charge and matter form factors", "[gra::nuclear][config][density]") {
  const auto nucleus      = MakeNeutronSkin();
  auto       coarse_param = ConfigControls(16);
  coarse_param.d_min      = 0.0;
  auto fine_param         = coarse_param;
  fine_param.count        = 1024;

  const auto   coarse       = SampleBank(nucleus, coarse_param, 9157);
  const auto   fine         = SampleBank(nucleus, fine_param, 9157);
  const double coarse_error = ConfigFormError(coarse);
  const double fine_error   = ConfigFormError(fine);
  CAPTURE(coarse_error, fine_error);
  REQUIRE(fine_error < 0.5 * coarse_error);
  REQUIRE(fine_error < 0.02);

  const unsigned int         n = nucleus.N();
  constexpr double           q = 0.075;
  const std::complex<double> sampled_neutron =
      (fine.MatterStat(q, 0.0).mean - fine.ChargeStat(q, 0.0).mean) / static_cast<double>(n);
  const double expected_neutron = (static_cast<double>(nucleus.A()) * nucleus.MatterDensity().Form(q) -
                                   static_cast<double>(nucleus.Z()) * nucleus.ChargeDensity().Form(q)) /
                                  static_cast<double>(n);
  REQUIRE(sampled_neutron.real() == Approx(expected_neutron).margin(0.025));
  REQUIRE(std::abs(sampled_neutron.imag()) < 0.025);
}

TEST_CASE("Configuration sampling rejects a negative derived neutron density", "[gra::nuclear][config][density]") {
  gra::nuclear::NucleusParam nucleus_param;
  nucleus_param.pdg    = gra::nuclear::EncodeNuclearPDG(20, 19);
  nucleus_param.mass   = 18.6;
  nucleus_param.charge = {6.0, 0.8, 20.4, 64, 2048, 3.0, 4097, 1.0e-7};
  nucleus_param.matter = {2.8, 0.4, 10.0, 64, 2048, 3.0, 4097, 1.0e-7};
  const gra::nuclear::MNucleus nucleus(nucleus_param);
  auto                         config_param = ConfigControls(1);
  REQUIRE_THROWS_AS(gra::nuclear::MConfigSampler(nucleus, config_param), std::invalid_argument);
}

TEST_CASE("Nuclear recentering preserves separate proton and neutron radii", "[gra::nuclear][config][density]") {
  auto param   = MakeOxygen().Param();
  param.pdg    = gra::nuclear::EncodeNuclearPDG(8, 2);
  param.mass   = 7.48;
  param.charge = {2.0, 0.4, 12.0, 96, 4096, 3.0, 4097, 1.0e-7};
  param.matter = {2.8, 0.4, 12.0, 96, 4096, 3.0, 4097, 1.0e-7};
  const gra::nuclear::MNucleus nucleus(param);
  auto                         controls = ConfigControls(1);
  controls.d_min                        = 0.0;
  const gra::nuclear::MConfigSampler sampler(nucleus, controls);
  gra::MRandom                       random;
  random.SetSeed(94621);
  constexpr std::size_t samples = 32768;
  constexpr double      q       = 1.0e-3;
  std::array<double, 2> radius2{};
  std::array<double, 2> charge_form{};
  double                center2 = 0.0;
  for (std::size_t sample = 0; sample < samples; ++sample) {
    const auto            config = sampler.Draw(random);
    std::array<double, 3> center{};
    for (const auto &nucleon : config.Nucleons()) {
      gra::AddScaled(center, nucleon.x, 1.0);
      radius2[nucleon.type == gra::nuclear::NucleonType::Proton ? 0 : 1] += gra::SquaredNorm(nucleon.x);
    }
    center2 = std::max(center2, gra::SquaredNorm(center));
    charge_form[0] += config.ChargeCurrent(q, 0.0, 0.0).real();
    charge_form[1] += config.ChargeCurrent(0.0, 0.0, q).real();
  }
  radius2[0] /= static_cast<double>(samples * nucleus.Z());
  radius2[1] /= static_cast<double>(samples * nucleus.N());
  const double charge2  = gra::math::pow2(nucleus.ChargeDensity().Rms());
  const double matter2  = gra::math::pow2(nucleus.MatterDensity().Rms());
  const double neutron2 = (nucleus.A() * matter2 - nucleus.Z() * charge2) / nucleus.N();
  CAPTURE(radius2[0], radius2[1], charge2, neutron2);
  REQUIRE(center2 < 1.0e-24);
  REQUIRE(radius2[0] == Approx(charge2).epsilon(0.015));
  REQUIRE(radius2[1] == Approx(neutron2).epsilon(0.015));
  REQUIRE((nucleus.Z() * radius2[0] + nucleus.N() * radius2[1]) / nucleus.A() == Approx(matter2).epsilon(0.015));
  // The low-q charge form factor measures the same radius along either axis
  for (const double sum : charge_form) {
    const double form          = sum / static_cast<double>(samples * nucleus.Z());
    const double slope_radius2 = 6.0 * gra::math::pow2(gra::PDG::GeV2fm) * (1.0 - form) / (q * q);
    REQUIRE(slope_radius2 == Approx(charge2).epsilon(0.02));
  }

  // Equal-mass two-body recentering cannot produce distinct constituent radii
  param.pdg  = gra::nuclear::EncodeNuclearPDG(2, 1);
  param.mass = 1.8756;
  const gra::nuclear::MNucleus asymmetric_pair(param);
  REQUIRE_THROWS_AS(gra::nuclear::MConfigSampler(asymmetric_pair, controls), std::invalid_argument);
  param.matter = param.charge;
  const gra::nuclear::MNucleus       symmetric_pair(param);
  const gra::nuclear::MConfigSampler pair_sampler(symmetric_pair, controls);
  const auto                         pair = pair_sampler.Draw(random);
  REQUIRE(gra::SquaredNorm(gra::Add(pair.Nucleons()[0].x, pair.Nucleons()[1].x)) < 1.0e-24);
}

TEST_CASE("Coherent and incoherent nuclear photon densities are orthogonal", "[gra::nuclear][photon][good-walker]") {
  const auto                  lead = MakeLead();
  const gra::nuclear::MPhoton photon(lead);
  constexpr double            xi = 1.0e-4;
  constexpr double            pt = 0.05;
  const auto density = [&](const gra::nuclear::CoherenceType type, const double x, const double t, const double q) {
    return photon.Density(type, x, t, q);
  };

  const auto coherent   = density(gra::nuclear::CoherenceType::Coherent, xi, -pt * pt, pt);
  const auto incoherent = density(gra::nuclear::CoherenceType::Incoherent, xi, -pt * pt, pt);
  const auto inclusive  = density(gra::nuclear::CoherenceType::Inclusive, xi, -pt * pt, pt);
  REQUIRE(std::isfinite(coherent.Trace()));
  REQUIRE(std::isfinite(incoherent.Trace()));
  REQUIRE(coherent.parallel > 0.0);
  REQUIRE(std::abs(coherent.perpendicular) < 1.0e-14);
  REQUIRE(incoherent.parallel > 0.0);
  REQUIRE(incoherent.perpendicular > 0.0);
  REQUIRE(inclusive.parallel == Approx(coherent.parallel + incoherent.parallel).epsilon(1.0e-13));
  REQUIRE(inclusive.perpendicular == Approx(coherent.perpendicular + incoherent.perpendicular).epsilon(1.0e-13));

  const auto forward_incoherent = density(gra::nuclear::CoherenceType::Incoherent, xi, 0.0, pt);
  REQUIRE(std::abs(forward_incoherent.Trace()) < 1.0e-12);
  REQUIRE(density(gra::nuclear::CoherenceType::Coherent, -xi, -pt * pt, pt).Trace() < 1.0e-14);
  REQUIRE(density(gra::nuclear::CoherenceType::Incoherent, 1.0 / 208.0, -pt * pt, pt).Trace() < 1.0e-14);

  auto half_charge_param = lead.Param();
  half_charge_param.pdg  = gra::nuclear::EncodeNuclearPDG(208, 41);
  const gra::nuclear::MPhoton half_charge{gra::nuclear::MNucleus(half_charge_param)};
  REQUIRE(coherent.Trace() / half_charge.Density(gra::nuclear::CoherenceType::Coherent, xi, -pt * pt, pt).Trace() ==
          Approx(4.0).epsilon(1.0e-12));
  REQUIRE(incoherent.Trace() / half_charge.Density(gra::nuclear::CoherenceType::Incoherent, xi, -pt * pt, pt).Trace() ==
          Approx(2.0).epsilon(1.0e-12));

  constexpr double proton_mass      = gra::PDG::mp;
  const double     proton_xi        = static_cast<double>(lead.A()) * xi;
  const double     proton_q2        = (pt * pt + proton_xi * proton_xi * proton_mass * proton_mass) / (1.0 - proton_xi);
  const double     ge               = gra::form::G_E(proton_q2);
  const double     gm               = gra::form::G_M(proton_q2);
  const double     mass2            = proton_mass * proton_mass;
  const double     electric_part    = (4.0 * mass2 * ge * ge + proton_q2 * gm * gm) / (4.0 * mass2 + proton_q2);
  const double     delta            = pt * pt / (pt * pt + proton_xi * proton_xi * mass2);
  const double     common           = 16.0 * gra::math::PI * gra::qed::alpha_QED() / (proton_xi * pt * pt);
  const double     electric         = common * (1.0 - proton_xi) * delta * delta * electric_part;
  const double     magnetic         = common * proton_xi * proton_xi * delta * gm * gm / 4.0;
  const double     form             = lead.ChargeDensity().Form(pt);
  const double     current_variance = static_cast<double>(lead.Z()) * (1.0 - form * form);
  const double     scale            = static_cast<double>(lead.A()) * current_variance;
  REQUIRE(incoherent.parallel == Approx(scale * (electric + magnetic)).epsilon(1.0e-12));
  REQUIRE(incoherent.perpendicular == Approx(scale * magnetic).epsilon(1.0e-12));
}

// Check the finite proton magnetic response at zero transverse photon momentum
TEST_CASE("Incoherent nuclear photon density retains its forward magnetic limit", "[gra::nuclear][photon][forward]") {
  const auto                  lead = MakeLead();
  const gra::nuclear::MPhoton photon(lead);
  constexpr double            xi        = 1.0e-4;
  const double                proton_xi = lead.A() * xi;
  const double                q2        = gra::math::pow2(xi * lead.Mass()) / (1.0 - xi);
  const double                proton_q2 = gra::math::pow2(proton_xi * gra::PDG::mp) / (1.0 - proton_xi);
  const double                gm        = gra::form::G_M(proton_q2);
  const double                form      = lead.ChargeDensity().Form(std::sqrt(q2));
  const double                variance  = lead.Z() * (1.0 - form * form);
  const double                expected =
      4.0 * gra::math::PI * gra::qed::alpha_QED() * variance * gm * gm / (xi * gra::PDG::mp * gra::PDG::mp);
  REQUIRE(expected > 0.0);
  for (const double pt : {0.0, 1.0e-8, 1.0e-12, 1.0e-200}) {
    const auto density = photon.Density(gra::nuclear::CoherenceType::Incoherent, xi, -q2, pt);
    REQUIRE(density.parallel == Approx(expected).epsilon(1.0e-8));
    REQUIRE(density.perpendicular == Approx(expected).epsilon(1.0e-8));
  }
  REQUIRE(photon.Density(gra::nuclear::CoherenceType::Coherent, xi, -q2, 0.0).Trace() == Approx(0.0));
}

TEST_CASE("Coherent nuclear photon current retains diffraction signs", "[gra::nuclear][photon][amplitude]") {
  const gra::nuclear::MPhoton photon(MakeLead());
  constexpr double            xi       = 1.0e-4;
  constexpr double            q_before = 0.05;
  constexpr double            q_after  = 0.18;
  REQUIRE(photon.Nucleus().ChargeDensity().Form(q_before) > 0.0);
  REQUIRE(photon.Nucleus().ChargeDensity().Form(q_after) < 0.0);
  const double before = photon.CoherentAmp(xi, -q_before * q_before, q_before);
  const double after  = photon.CoherentAmp(xi, -q_after * q_after, q_after);
  REQUIRE(before > 0.0);
  REQUIRE(after < 0.0);
  REQUIRE(photon.Density(gra::nuclear::CoherenceType::Coherent, xi, -q_before * q_before, q_before).parallel ==
          Approx(before * before).epsilon(2.0e-14));
  REQUIRE(photon.Density(gra::nuclear::CoherenceType::Coherent, xi, -q_after * q_after, q_after).parallel ==
          Approx(after * after).epsilon(2.0e-14));
  REQUIRE(std::abs(photon.CoherentAmp(xi, q_after * q_after, q_after)) < 2.0e-14);
}

TEST_CASE("Nuclear photon source retains charge sign and helicity phase", "[gra::nuclear][photon][amplitude]") {
  const double               source          = std::sqrt(2.0);
  const std::complex<double> positive        = gra::nuclear::PhotonAmp(source, 0.25, 82.0, 0.3, 1);
  const std::complex<double> negative_charge = gra::nuclear::PhotonAmp(source, 0.25, -82.0, 0.3, 1);
  const std::complex<double> negative_source = gra::nuclear::PhotonAmp(-source, 0.25, 82.0, 0.3, 1);
  const std::complex<double> opposite        = gra::nuclear::PhotonAmp(source, 0.25, 82.0, 0.3, -1);
  REQUIRE(std::norm(positive) == Approx(4.0).epsilon(2.0e-14));
  REQUIRE(std::abs(negative_charge + positive) < 2.0e-14);
  REQUIRE(std::abs(negative_source + positive) < 2.0e-14);
  REQUIRE(std::abs(opposite - std::conj(positive)) < 2.0e-14);
  REQUIRE(std::abs(gra::nuclear::PhotonAmp(0.0, 0.25, 82.0, 0.3, 1)) < 2.0e-14);
}

// Reject an unresolved eigenstate variance and retain the requested physical moments
TEST_CASE("Nuclear fluctuation rules enforce moment convergence", "[gra::nuclear][glauber][ggcf]") {
  auto param         = GlauberControls(1.0e-3).fluctuation;
  param.cdf_max_iter = 1;
  REQUIRE_THROWS_WITH(gra::nuclear::CrossSectionRule(param), "CrossSectionRule: moment match did not converge");
  param.cdf_max_iter    = 10000;
  const auto   scale    = gra::nuclear::CrossSectionRule(param);
  const double mean     = gra::Sum(scale) / scale.size();
  const double variance = gra::SquaredNorm(scale) / scale.size() - mean * mean;
  REQUIRE(mean == Approx(1.0).margin(1.0e-14));
  REQUIRE(variance == Approx(param.omega).margin(1.0e-12));
}

// Keep failures in energy dependent target profiles recoverable during sampling
TEST_CASE("Photonuclear profile failures use amplitude bookkeeping", "[gra::nuclear][photo][failure]") {
  auto param                        = PhotoControls();
  param.fluctuation.cdf_max_iter    = 1;
  const auto                 oxygen = MakeOxygen();
  const gra::nuclear::MPhoto photo(oxygen, param, gra::nuclear::PhotoModel::Glauber);
  auto                       profile = gra::test::PhotoProfile(3.0, 0.1, 1.0e-3);
  REQUIRE_THROWS_AS(photo.Factors(profile, 0.1, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ),
                    gra::AmplitudeFailure);
  profile.omega     = 0.0;
  profile.sigma_eff = 1.0e6;
  REQUIRE_THROWS_AS(photo.Factors(profile, 0.1, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ),
                    gra::AmplitudeFailure);
  param.fluctuation.normal_shape_min = 0.0;
  REQUIRE_THROWS_AS(gra::nuclear::MPhoto(oxygen, param, gra::nuclear::PhotoModel::Glauber), std::invalid_argument);
}

TEST_CASE("GGCF survival reaches the optical limit continuously", "[gra::nuclear][glauber][ggcf]") {
  const auto                   oxygen = MakeOxygen();
  const gra::nuclear::MGlauber optical(oxygen, oxygen, GlauberControls(0.0));
  constexpr double omega = 1.0e-6;
  const gra::nuclear::MGlauber vanishing_fluctuation(oxygen, oxygen, GlauberControls(omega));
  const gra::nuclear::MGlauber fluctuating(oxygen, oxygen, GlauberControls(0.2));
  auto                         relaxed_param = GlauberControls(0.2);
  relaxed_param.fluctuation.cdf_tol          = 1.0e-8;
  const gra::nuclear::MGlauber relaxed(oxygen, oxygen, relaxed_param);
  for (const double impact : {0.0, 2.0, 4.0, 6.0, 10.0}) {
    REQUIRE(relaxed.GGCFAmp(impact) == Approx(fluctuating.GGCFAmp(impact)).margin(1.0e-7));
  }

  constexpr double b = 4.0;
  REQUIRE(optical.Overlap(0.0) > optical.Overlap(b));
  REQUIRE(optical.Overlap(b) > 0.0);
  REQUIRE(optical.GGCFAmp(b) == Approx(optical.OpticalAmp(b)).epsilon(2.0e-12));
  // Expand <exp(-tau f1 f2)> with independent unit-mean nucleon strengths
  const double tau = -std::log(optical.OpticalAmp(b));
  const double correction = 0.5 * tau * tau * (std::pow(1.0 + omega, 2) - 1.0);
  REQUIRE(vanishing_fluctuation.GGCFAmp(b) / optical.OpticalAmp(b) - 1.0 == Approx(correction).epsilon(0.01));
  REQUIRE(fluctuating.GGCFAmp(b) >= optical.OpticalAmp(b));
  RequireProbability(optical.OpticalProb(b));
  RequireProbability(fluctuating.GGCFProb(b));
  REQUIRE(optical.OpticalProb(19.9) > optical.OpticalProb(0.0));

  auto                         nucleus = std::make_shared<const gra::nuclear::MNucleus>(oxygen);
  const gra::nuclear::MGlauber proton_oxygen({gra::nuclear::HadronType::Proton, gra::nuclear::HadronType::Nucleus},
                                             {nullptr, nucleus}, GlauberControls(0.1));
  REQUIRE(proton_oxygen.IsProton(1));
  REQUIRE_FALSE(proton_oxygen.IsProton(2));
  REQUIRE(proton_oxygen.Nucleus(2)->A() == 16);
  REQUIRE(proton_oxygen.Overlap(0.0) > 0.0);
  RequireProbability(proton_oxygen.GGCFProb(b));
}

TEST_CASE("The configuration nucleon profile integrates to sigma_nn", "[gra::nuclear][glauber][config]") {
  auto                         param = GlauberControls(0.0);
  const gra::nuclear::MGlauber proton_proton({gra::nuclear::HadronType::Proton, gra::nuclear::HadronType::Proton},
                                             {nullptr, nullptr}, param);
  const double                 sigma_fm2 = NNCrossSection(proton_proton, param.profile.sigma);
  REQUIRE(sigma_fm2 == Approx(0.1 * param.profile.sigma).epsilon(2.0e-10));
}

TEST_CASE("Every fluctuating nucleon profile integrates to its sampled sigma",
          "[gra::nuclear][glauber][config][ggcf]") {
  auto param              = GlauberControls(0.2);
  param.fluctuation.nodes = 8;
  const gra::nuclear::MGlauber proton_proton({gra::nuclear::HadronType::Proton, gra::nuclear::HadronType::Proton},
                                             {nullptr, nullptr}, param);
  for (std::size_t upper = 0; upper < proton_proton.SigmaCount(); ++upper) {
    for (std::size_t lower = 0; lower < proton_proton.SigmaCount(); ++lower) {
      const double sigma    = proton_proton.Sigma(upper, lower);
      const double integral = NNCrossSection(proton_proton, sigma);
      REQUIRE(integral == Approx(0.1 * sigma).epsilon(3.0e-10));
      for (const double b : {0.0, 0.4, 1.1, 2.3}) {
        const double actual = proton_proton.ConfigPairAmp(b, 0.0, nullptr, nullptr, 0, 0, upper, lower);
        REQUIRE(actual == Approx(proton_proton.NNAmp(b, sigma)).margin(2.0e-13));
        const double exchanged = proton_proton.ConfigPairAmp(b, 0.0, nullptr, nullptr, 0, 0, lower, upper);
        REQUIRE(actual == Approx(exchanged).margin(2.0e-13));
      }
    }
  }
}

TEST_CASE("Configuration Glauber couples Gamma cross-section fluctuations", "[gra::nuclear][glauber][config][ggcf]") {
  const auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       config = ConfigControls(1);
  config.d_min      = 0.3;
  const auto bank   = SampleBankPtr(*oxygen, config, 90210);

  auto fluctuation_param              = GlauberControls(0.2);
  fluctuation_param.fluctuation.nodes = 16;
  const gra::nuclear::MGlauber fluctuating({gra::nuclear::HadronType::Proton, gra::nuclear::HadronType::Nucleus},
                                           {nullptr, oxygen}, fluctuation_param);
  REQUIRE(fluctuating.SampleCount(nullptr, bank.get()) == config.count);
  const std::size_t count = fluctuating.SigmaCount();
  REQUIRE(count == fluctuation_param.fluctuation.nodes);

  double mean   = 0.0;
  double second = 0.0;
  for (std::size_t sample = 0; sample < count; ++sample) {
    for (std::size_t other = 0; other < count; ++other) {
      const double sigma = fluctuating.Sigma(sample, other);
      mean += sigma;
      second += sigma * sigma;
    }
  }
  RequireProbability(gra::math::pow2(fluctuating.SampleAmp(4.0, nullptr, bank.get(), 0)));
  mean /= static_cast<double>(count * count);
  second /= static_cast<double>(count * count);
  const double relative_variance = second / (mean * mean) - 1.0;
  REQUIRE(mean == Approx(fluctuation_param.profile.sigma).epsilon(2.0e-14));
  REQUIRE(relative_variance == Approx(std::pow(1.0 + fluctuation_param.fluctuation.omega, 2) - 1.0).epsilon(2.0e-12));

  auto fixed_param              = fluctuation_param;
  fixed_param.profile.omega     = 0.0;
  fixed_param.fluctuation.omega = 0.0;
  const gra::nuclear::MGlauber fixed({gra::nuclear::HadronType::Proton, gra::nuclear::HadronType::Nucleus},
                                     {nullptr, oxygen}, fixed_param);
  auto                         small_param = fluctuation_param;
  small_param.profile.omega                = 1.0e-6;
  small_param.fluctuation.omega            = 1.0e-6;
  const gra::nuclear::MGlauber small({gra::nuclear::HadronType::Proton, gra::nuclear::HadronType::Nucleus},
                                     {nullptr, oxygen}, small_param);
  REQUIRE(small.ConfigAmp(4.0, nullptr, bank.get()) ==
          Approx(fixed.ConfigAmp(4.0, nullptr, bank.get())).epsilon(2.0e-5));
  // A finite-range profile need not be convex in its fluctuating cross section
  double mean_amp = 0.0;
  for (std::size_t upper = 0; upper < count; ++upper) {
    for (std::size_t lower = 0; lower < count; ++lower) {
      double amplitude = 1.0;
      for (const auto &nucleon : bank->At(0).Nucleons()) {
        amplitude *= fluctuating.NNAmp(std::hypot(4.0 + nucleon.x[0], nucleon.x[1]), fluctuating.Sigma(upper, lower));
      }
      mean_amp += amplitude;
    }
  }
  REQUIRE(fluctuating.ConfigAmp(4.0, nullptr, bank.get()) == Approx(mean_amp / (count * count)).epsilon(2.0e-13));
}

TEST_CASE("Glauber samples use Cartesian two-leg configuration indexing",
          "[gra::nuclear][glauber][config][cartesian]") {
  const auto oxygen       = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       first_param  = ConfigControls(2);
  first_param.d_min       = 0.3;
  auto second_param       = first_param;
  second_param.count      = 3;
  const auto first        = SampleBank(*oxygen, first_param, 20260826);
  const auto second       = SampleBank(*oxygen, second_param, 20260826 ^ 0x9e3779b9U);
  auto       param        = GlauberControls(0.0);
  param.fluctuation.nodes = 1;
  const gra::nuclear::MGlauber glauber({gra::nuclear::HadronType::Nucleus, gra::nuclear::HadronType::Nucleus},
                                       {oxygen, oxygen}, param);

  REQUIRE(glauber.SampleCount(&first, &second) == 6);
  for (std::size_t sample = 0; sample < 6; ++sample) {
    const std::size_t config1 = sample / second.Size();
    const std::size_t config2 = sample % second.Size();
    REQUIRE(glauber.SampleAmp(1.7, -0.4, &first, &second, sample) ==
            Approx(glauber.ConfigPairAmp(1.7, -0.4, &first, &second, config1, config2)).epsilon(2.0e-14));
  }
  REQUIRE_THROWS_AS(glauber.SampleAmp(1.7, -0.4, &first, &second, 6), std::out_of_range);
}

// Check spatially restricted amplitudes against every NN pair and beam exchange
TEST_CASE("Batched Glauber amplitudes equal scalar fixed-state products", "[gra::nuclear][glauber][config][batch]") {
  const auto oxygen       = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       config       = ConfigControls(1);
  config.d_min            = 0.3;
  const auto first        = SampleBank(*oxygen, config, 20260827);
  const auto second       = SampleBank(*oxygen, config, 20260827 ^ 0x9e3779b9U);
  auto       param        = GlauberControls(0.2);
  param.fluctuation.nodes = 8;
  const gra::nuclear::MGlauber glauber({gra::nuclear::HadronType::Nucleus, gra::nuclear::HadronType::Nucleus},
                                       {oxygen, oxygen}, param);
  const std::vector<std::array<double, 2>> impact = {{0.0, 0.0}, {1.7, -0.4}, {-4.2, 3.1}, {9.0, 8.0}, {30.0, 0.0}};
  const std::size_t                        sigma1 = 1;
  const std::size_t                        sigma2 = glauber.SigmaCount() - 2;
  const auto                               batch = glauber.ConfigPairAmp(impact, &first, &second, 0, 0, sigma1, sigma2);
  REQUIRE(batch.size() == impact.size());
  for (const auto &i : indices(impact)) {
    double direct = 1.0;
    for (const auto &upper : first.At(0).Nucleons()) {
      for (const auto &lower : second.At(0).Nucleons()) {
        const double radius = std::hypot(upper.x[0] - lower.x[0] - impact[i][0],
                                         upper.x[1] - lower.x[1] - impact[i][1]);
        direct *= glauber.NNAmp(radius, glauber.Sigma(sigma1, sigma2));
      }
    }
    REQUIRE(batch[i] == Approx(direct).epsilon(2.0e-13));
    REQUIRE(batch[i] == Approx(glauber.ConfigPairAmp(-impact[i][0], -impact[i][1], &second, &first, 0, 0, sigma2, sigma1))
                            .epsilon(2.0e-13));
    REQUIRE(batch[i] == Approx(glauber.ConfigPairAmp(impact[i][0], impact[i][1], &first, &second, 0, 0, sigma1, sigma2))
                            .epsilon(2.0e-13));
  }
  const auto averaged = glauber.ConfigPairAmp(impact, &first, &second, 0, 0);
  REQUIRE(averaged.size() == impact.size());
  for (const auto &i : indices(impact)) {
    REQUIRE(averaged[i] ==
            Approx(glauber.ConfigPairAmp(impact[i][0], impact[i][1], &first, &second, 0, 0)).epsilon(2.0e-13));
  }
  double fixed_mean = 0.0;
  for (std::size_t upper = 0; upper < glauber.SigmaCount(); ++upper) {
    for (std::size_t lower = 0; lower < glauber.SigmaCount(); ++lower) {
      fixed_mean += glauber.ConfigPairAmp(1.7, -0.4, &first, &second, 0, 0, upper, lower);
    }
  }
  fixed_mean /= static_cast<double>(glauber.SigmaCount() * glauber.SigmaCount());
  REQUIRE(fixed_mean == Approx(glauber.ConfigPairAmp(1.7, -0.4, &first, &second, 0, 0)).epsilon(2.0e-13));
}

TEST_CASE("Electromagnetic absorption retains its charge and impact dependence", "[gra::nuclear][breakup]") {
  const gra::nuclear::MBreakup breakup(BreakupControls(208, 82, 82.0, 1.0e7));
  REQUIRE(breakup.PhotoAbsorption(0.013373) > 644.0);
  REQUIRE(breakup.PhotoAbsorption(0.013373) < 648.0);
  REQUIRE(breakup.PhotoAbsorption(0.100) > 5.0);
  REQUIRE(breakup.PhotoAbsorption(0.350) > 80.0);
  constexpr double qd_energy = 50.0;
  const double     qd_pauli  = 8.3714e-2 - 9.8343e-3 * qd_energy + 4.1222e-4 * qd_energy * qd_energy -
                          3.4762e-6 * qd_energy * qd_energy * qd_energy + 9.3537e-9 * std::pow(qd_energy, 4);
  const double deuteron  = 61.2 * std::pow(qd_energy - 2.224, 1.5) / std::pow(qd_energy, 3);
  const double qd        = 6.5 * 126.0 * 82.0 / 208.0 * deuteron * qd_pauli;
  const double gdr_peak  = 2.0 * 1.33716 * (60.0 * 126.0 * 82.0 / 208.0) / (gra::math::PI * 3.938);
  const double gdr_ratio = qd_energy / 13.373;
  const double gdr_width = 3.938 / 13.373;
  const double gdr =
      gdr_peak * gdr_width * gdr_width / (std::pow(gdr_ratio - 1.0 / gdr_ratio, 2) + gdr_width * gdr_width);
  REQUIRE(breakup.PhotoAbsorption(0.050) == Approx(gdr + qd).epsilon(5.0e-13));
  const double     below_resonance = breakup.PhotoAbsorption(std::nextafter(0.100, 0.0));
  constexpr double area            = 23000.0;
  constexpr double width           = 110.0;
  constexpr double center          = 350.0;
  const double     pull            = (100.0 - center) / width;
  const double     resonance       = area * std::exp(-0.5 * pull * pull) / (width * std::sqrt(2.0 * gra::math::PI));
  REQUIRE(breakup.PhotoAbsorption(0.100) - below_resonance == Approx(resonance).epsilon(2.0e-12));
  REQUIRE(breakup.PhotoAbsorption(0.00599) == Approx(0.0).margin(1.0e-14));
  REQUIRE(breakup.PhotoAbsorption(0.006) > 0.0);
  REQUIRE(breakup.PhotoAbsorption(0.00736) > 0.0);
  REQUIRE(breakup.PhotoAbsorption(0.00738) > 0.0);
  gra::MRandom photon_random;
  photon_random.SetSeed(6091);
  bool sampled_below_threshold = false;
  for (unsigned int sample = 0; sample < 4096; ++sample) {
    const double omega = breakup.SampleEnergy(20.0, photon_random);
    REQUIRE(omega >= 0.006);
    sampled_below_threshold = sampled_below_threshold || omega < 0.0073678686;
  }
  REQUIRE(sampled_below_threshold);
  REQUIRE_THROWS_AS(breakup.PhotoAbsorption(-0.001), std::invalid_argument);
  REQUIRE_THROWS_AS(breakup.PhotoAbsorption(std::numeric_limits<double>::infinity()), std::invalid_argument);
  REQUIRE(std::isfinite(breakup.PhotoAbsorption(std::numeric_limits<double>::max())));

  const gra::nuclear::MBreakup half_charge(BreakupControls(208, 82, 41.0, 1.0e7));
  REQUIRE(half_charge.Mean(20.0) == Approx(0.25 * breakup.Mean(20.0)).epsilon(2.0e-12));

  // The additive continuum and its first derivative join continuously
  const auto controls = BreakupControls(208, 82, 82.0, 1.0e7);
  for (const double edge : {controls.photo.continuum.threshold, controls.photo.continuum.match, 8000.0}) {
    const double energy = edge * 1.0e-3;
    const double step   = energy * 1.0e-6;
    const double middle = breakup.PhotoAbsorption(energy);
    const double lower  = breakup.PhotoAbsorption(energy - step);
    const double upper  = breakup.PhotoAbsorption(energy + step);
    CAPTURE(energy);
    REQUIRE(lower == Approx(upper).epsilon(1.0e-5));
    REQUIRE((upper - 2.0 * middle + lower) / middle == Approx(0.0).margin(1.0e-10));
  }
  for (const double omega : {0.5, 1.0, 8.0, 80.0}) { REQUIRE(breakup.PhotoAbsorption(omega) > 0.0); }

  const gra::nuclear::MBreakup fine(BreakupControls(208, 82, 82.0, 1.0e7, 192));
  const gra::nuclear::MBreakup coarse(BreakupControls(208, 82, 82.0, 1.0e7, 32));
  auto                         tail_param = BreakupControls(208, 82, 82.0, 1.0e7);
  tail_param.response.tail_rel_tol *= 0.5;
  const gra::nuclear::MBreakup tail(tail_param);
  for (const double b : std::array<double, 4>{6.7, 20.0, 137.0, 2500.0}) {
    CAPTURE(b);
    REQUIRE(fine.Mean(b) == Approx(breakup.Mean(b)).epsilon(2.0e-10));
    REQUIRE(coarse.Mean(b) == Approx(fine.Mean(b)).epsilon(1.0e-7));
    REQUIRE(tail.Mean(b) == Approx(breakup.Mean(b)).epsilon(0.01));
  }
  REQUIRE(breakup.Mean(20.0) > breakup.Mean(137.0));
  REQUIRE(breakup.Mean(137.0) > breakup.Mean(2500.0));
}

// Check the continuous deposited energy against its analytic moments
TEST_CASE("TCM transfer conserves absorbed energy and its continuous moments", "[gra::nuclear][neutron][physics]") {
  const gra::nuclear::MBreakup breakup(BreakupControls(208, 82, 82.0, 1.0e7, 48));
  const gra::nuclear::MBreakup fine(BreakupControls(208, 82, 82.0, 1.0e7, 96));
  const double                 match = 1.0e-3 * breakup.Param().transfer.match;
  const double                 E0    = 1.0e-3 * breakup.Param().transfer.E0;
  for (const double omega : {0.001, 0.007, 0.010, 0.050, match, 10.0 * match}) {
    gra::MRandom random;
    gra::MRandom fine_random;
    random.SetSeed(53901);
    fine_random.SetSeed(53901);
    constexpr unsigned int samples  = 32768;
    double                 sum      = 0.0;
    double                 square   = 0.0;
    const double           transfer = std::min(omega, match);
    for (unsigned int i = 0; i < samples; ++i) {
      const double energy = breakup.SampleTransfer(omega, random);
      REQUIRE(energy >= 0.0);
      REQUIRE(energy <= transfer);
      REQUIRE(fine.SampleTransfer(omega, fine_random) == Approx(energy).margin(1.0e-15));
      sum += energy;
      square += energy * energy;
    }
    const double full   = std::exp(-transfer / E0);
    const double mean   = transfer * (full + (1.0 - full) / 2.0);
    const double second = transfer * transfer * (full + (1.0 - full) / 3.0);
    const double fourth = std::pow(transfer, 4) * (full + (1.0 - full) / 5.0);
    CAPTURE(omega, mean, sum / samples);
    CHECK(sum / samples == Approx(mean).margin(6.0 * std::sqrt((second - mean * mean) / samples)));
    CHECK(square / samples == Approx(second).margin(6.0 * std::sqrt((fourth - second * second) / samples)));
  }
}

// Check absorption convergence and the zero-excitation probability
TEST_CASE("EMD absorption tolerance preserves compound excitation", "[gra::nuclear][breakup][convergence]") {
  const auto numerics =
      nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("modeldata/TUNE0/NUMERICS.json")))
          .at("NUMERICS_NUCLEAR")
          .at("EMD");
  auto param             = BreakupControls(208, 82, 82.0, 1.0e7);
  param.response.nodes   = numerics.at("response").at("nodes").get<unsigned int>();
  param.response.rel_tol = numerics.at("response").at("rel_tol").get<double>();
  param.profile.nodes = numerics.at("profile").at("nodes").get<std::size_t>();
  param.profile.b_min = numerics.at("profile").at("b_min").get<double>();
  const gra::nuclear::MBreakup model(param);
  param.response.nodes *= 2;
  param.response.rel_tol *= 0.01;
  param.profile.nodes *= 2;
  const gra::nuclear::MBreakup reference(param);
  gra::MRandom                 random;
  random.SetSeed(73941);
  for (const double b : {14.0, 20.0, 100.0, 1000.0}) {
    CHECK(model.Mean(b) == Approx(reference.Mean(b)).epsilon(1.0e-3));
    constexpr unsigned int samples = 16384;
    unsigned int           empty   = 0;
    for (unsigned int i = 0; i < samples; ++i) {
      const double excitation = model.SampleExcitation(b, random);
      REQUIRE(std::isfinite(excitation));
      REQUIRE(excitation >= 0.0);
      empty += !(excitation > 0.0);
    }
    const double zero = std::exp(-reference.Mean(b));
    CHECK(static_cast<double>(empty) / samples ==
          Approx(zero).margin(6.0 * std::sqrt(zero * (1.0 - zero) / samples) + 1.0 / samples));
  }
}

TEST_CASE("Finite emitter charge densities regulate EMD at small impact", "[gra::nuclear][breakup][finite-size]") {
  auto                         lead        = std::make_shared<const gra::nuclear::MNucleus>(MakeLead());
  auto                         point_param = BreakupControls(208, 82, 82.0, 1.0e7);
  const gra::nuclear::MBreakup point(point_param);
  auto                         finite_param = point_param;
  finite_param.emitter                      = gra::nuclear::EmitterType::Nuclear;
  const gra::nuclear::MBreakup finite(finite_param, lead);

  REQUIRE(finite.Mean(0.0) == Approx(0.0).margin(1.0e-15));
  REQUIRE(finite.Mean(1.0) < point.Mean(1.0));
  REQUIRE(finite.Mean(30.0) == Approx(point.Mean(30.0)).epsilon(2.0e-13));

  auto proton_param    = point_param;
  proton_param.z_emit  = 1.0;
  proton_param.emitter = gra::nuclear::EmitterType::Proton;
  const gra::nuclear::MBreakup proton(proton_param);
  REQUIRE(proton.Mean(0.0) == Approx(0.0).margin(1.0e-15));
  REQUIRE(proton.Mean(0.1) < point.Mean(0.1) / (point_param.z_emit * point_param.z_emit));

  // A finite transverse charge density gives n(omega,b) proportional to b^2 near the origin
  constexpr double omega = 0.02;
  const double     b     = proton_param.profile.b_min;
  for (const double scale : {0.1, 0.01}) {
    REQUIRE(proton.PhotonDensity(omega, scale * b) / proton.PhotonDensity(omega, b) ==
            Approx(scale * scale).epsilon(1.0e-8));
    REQUIRE(finite.PhotonDensity(omega, scale * b) / finite.PhotonDensity(omega, b) ==
            Approx(scale * scale).epsilon(1.0e-8));
  }
}

// Check the shared E1 line strength and stable low and high energy limits
TEST_CASE("Dipole absorption preserves the TRK integral", "[gra::nuclear][breakup][decay]") {
  const auto param = BreakupControls(208, 82, 82.0, 1.0e7).photo.gdr;
  for (const auto &ion : std::array<std::array<unsigned int, 2>, 2>{{{208, 82}, {207, 82}}}) {
    const gra::nuclear::Dipole dipole(ion[0], ion[1], param);
    CHECK(dipole.Sigma(dipole.energy) == Approx(dipole.peak));
    const auto rule     = gra::math::GaussLegendreRule(256, 0.0, 1.0);
    double     integral = 0.0;
    for (const auto &i : indices(rule.first)) {
      const double u = rule.first[i];
      integral +=
          rule.second[i] * dipole.energy / gra::math::pow2(1.0 - u) * dipole.Sigma(dipole.energy * u / (1.0 - u));
    }
    double strength = param.systematics.strength;
    for (const auto &isotope : param.isotope) {
      if (isotope.a == ion[0] && isotope.z == ion[1]) { strength = isotope.strength; }
    }
    const double trk = 1.0e-3 * strength * param.systematics.trk * ion[1] * (ion[0] - ion[1]) / ion[0];
    CHECK(integral == Approx(trk).epsilon(1.0e-6));
    for (const double ratio : {1.0e-100, 1.0e-6, 1.0e6, 1.0e100}) {
      CHECK(std::isfinite(dipole.Sigma(dipole.energy * ratio)));
      CHECK(dipole.Sigma(dipole.energy * ratio) == Approx(dipole.Sigma(dipole.energy / ratio)).epsilon(1.0e-10));
    }
  }
}

// Require the measured photoabsorption fit to belong to the active target
TEST_CASE("Photonuclear continuum fits require matching target isotopes", "[gra::nuclear][breakup][isotope]") {
  {
    auto param = BreakupControls(208, 82, 82.0, 1.0e7);
    REQUIRE_NOTHROW(gra::nuclear::MBreakup(param));
    for (const auto &isotope : std::array<std::array<unsigned int, 2>, 2>{{{207, 82}, {208, 81}}}) {
      param.photo.continuum.a = isotope[0];
      param.photo.continuum.z = isotope[1];
      REQUIRE_THROWS_AS(gra::nuclear::MBreakup(param), std::invalid_argument);
    }
  }
}

TEST_CASE("Photonuclear coherent and incoherent target transitions are finite", "[gra::nuclear][photo][good-walker]") {
  const auto                 oxygen          = MakeOxygen();
  auto                       param           = PhotoControls();
  const auto                 profile         = gra::test::PhotoProfile(25.0, 0.1);
  const auto                 impulse_profile = gra::test::PhotoProfile(0.0, 0.1);
  const gra::nuclear::MPhoto photo(oxygen, param);
  const gra::nuclear::MPhoto impulse_mode(oxygen, param, gra::nuclear::PhotoModel::Impulse);
  REQUIRE(photo.Model() == gra::nuclear::PhotoModel::Glauber);
  REQUIRE(impulse_mode.Model() == gra::nuclear::PhotoModel::Impulse);

  const auto forward  = photo.Factors(profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  const auto finite_q = photo.Factors(profile, 0.16, -0.05, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  REQUIRE(std::isfinite(forward.coherent.real()));
  REQUIRE(std::isfinite(forward.coherent.imag()));
  REQUIRE(std::abs(forward.coherent) > 0.0);
  REQUIRE(forward.incoherent > 0.0);
  REQUIRE(std::isfinite(finite_q.incoherent));
  REQUIRE(finite_q.incoherent > 0.0);

  auto direct_param           = param;
  direct_param.table.qt_max   = 0.1;
  direct_param.table.qt_nodes = 33;
  const gra::nuclear::MPhoto direct(oxygen, direct_param);
  for (const auto direction : {gra::nuclear::PhotonDirection::PositiveZ, gra::nuclear::PhotonDirection::NegativeZ}) {
    for (const double qz : {-0.04, 0.0, 0.04}) {
      const auto cached    = photo.Factors(profile, 0.16, -0.05, qz, direction);
      const auto reference = direct.Factors(profile, 0.16, -0.05, qz, direction);
      REQUIRE(std::abs(cached.coherent - reference.coherent) / std::max(1.0, std::abs(reference.coherent)) < 1.0e-7);
      REQUIRE(std::abs(cached.incoherent - reference.incoherent) / std::max(1.0, reference.incoherent) < 1.0e-7);
    }
  }

  const gra::nuclear::MPhoto impulse(oxygen, param);
  const auto                 impulse_forward =
      impulse.Factors(impulse_profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  const auto selected_forward = impulse_mode.Factors(profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  REQUIRE(impulse_forward.coherent.real() == Approx(16.0).margin(2.0e-13));
  REQUIRE(std::abs(impulse_forward.coherent.imag()) < 2.0e-14);
  REQUIRE(impulse_forward.incoherent < 1.0e-12);
  REQUIRE(std::abs(selected_forward.coherent - impulse_forward.coherent) < 2.0e-13);
  REQUIRE(std::abs(selected_forward.incoherent - impulse_forward.incoherent) < 2.0e-13);
  REQUIRE(std::abs(forward.coherent) < std::abs(impulse_forward.coherent));

  const gra::M3Vec q     = {0.12, -0.04, 0.07};
  const double     qnorm = std::sqrt(gra::SquaredNorm(q));
  const double     form  = oxygen.MatterDensity().Form(qnorm);
  for (const auto direction : {gra::nuclear::PhotonDirection::PositiveZ, gra::nuclear::PhotonDirection::NegativeZ}) {
    const auto full     = impulse.Factors(impulse_profile, q[0], q[1], q[2], direction);
    const auto selected = impulse_mode.Factors(profile, q[0], q[1], q[2], direction);
    REQUIRE(full.coherent.real() == Approx(16.0 * form).epsilon(2.0e-12));
    REQUIRE(std::abs(full.coherent.imag()) < 2.0e-14);
    REQUIRE(full.incoherent == Approx(std::sqrt(16.0 * std::max(0.0, 1.0 - form * form))).epsilon(2.0e-12));
    REQUIRE(std::abs(selected.coherent - full.coherent) < 2.0e-13);
    REQUIRE(std::abs(selected.incoherent - full.incoherent) < 2.0e-13);
  }

  const auto positive = photo.Factors(profile, q[0], q[1], q[2], gra::nuclear::PhotonDirection::PositiveZ);
  const auto crossed  = photo.Factors(profile, q[0], q[1], -q[2], gra::nuclear::PhotonDirection::NegativeZ);
  REQUIRE(positive.coherent.real() == Approx(crossed.coherent.real()).epsilon(2.0e-10));
  REQUIRE(positive.coherent.imag() == Approx(crossed.coherent.imag()).epsilon(2.0e-10));
  REQUIRE(positive.incoherent == Approx(crossed.incoherent).epsilon(2.0e-10));

  const auto                 absorptive_profile = gra::test::PhotoProfile(25.0, 0.0);
  const gra::nuclear::MPhoto absorptive(oxygen, param);
  const auto plus  = absorptive.Factors(absorptive_profile, q[0], q[1], q[2], gra::nuclear::PhotonDirection::PositiveZ);
  const auto minus = absorptive.Factors(absorptive_profile, q[0], q[1], q[2], gra::nuclear::PhotonDirection::NegativeZ);
  REQUIRE(plus.coherent.real() == Approx(minus.coherent.real()).epsilon(2.0e-10));
  REQUIRE(plus.coherent.imag() == Approx(-minus.coherent.imag()).epsilon(2.0e-10));
  REQUIRE(plus.incoherent == Approx(minus.incoherent).epsilon(2.0e-10));

  auto config_param  = ConfigControls(8);
  config_param.d_min = 0.3;
  const auto bank    = SampleBank(oxygen, config_param, 424242);
  const auto bank_impulse =
      impulse.Factors(impulse_profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ, &bank);
  REQUIRE(bank_impulse.coherent.real() == Approx(impulse_forward.coherent.real()).epsilon(2.0e-14));
  REQUIRE(std::abs(bank_impulse.coherent.imag()) < 2.0e-14);
  REQUIRE(bank_impulse.incoherent < 2.0e-14);

  auto optical_config     = config_param;
  optical_config.count    = 512;
  optical_config.d_min    = 0.0;
  optical_config.sweeps   = 0;
  const auto optical_bank = SampleBank(oxygen, optical_config, 112358);
  const auto sampled =
      photo.ShadowStat(profile, optical_bank, q[0], q[1], q[2], gra::nuclear::PhotonDirection::PositiveZ);
  REQUIRE(std::abs(sampled.mean - positive.coherent) < 0.15 * std::max(1.0, std::abs(positive.coherent)));
  REQUIRE(std::abs(std::sqrt(sampled.variance) - positive.incoherent) < 0.35 * std::max(1.0, positive.incoherent));
  for (const auto direction : {gra::nuclear::PhotonDirection::PositiveZ, gra::nuclear::PhotonDirection::NegativeZ}) {
    const auto raw              = bank.At(0).MatterCurrent(q[0], q[1], q[2]);
    const auto impulse_current  = impulse.ShadowCurrent(impulse_profile, bank.At(0), q[0], q[1], q[2], direction);
    const auto selected_current = impulse_mode.ShadowCurrent(profile, bank.At(0), q[0], q[1], q[2], direction);
    REQUIRE(std::abs(impulse_current - raw) < 2.0e-14);
    REQUIRE(std::abs(selected_current - raw) < 2.0e-14);
    REQUIRE(std::abs(photo.ShadowCurrent(profile, bank.At(0), 0.0, 0.0, 0.0, direction)) <
            std::abs(impulse.ShadowCurrent(impulse_profile, bank.At(0), 0.0, 0.0, 0.0, direction)));

    std::complex<double> mean   = 0.0;
    double               second = 0.0;
    for (std::size_t sample = 0; sample < bank.Size(); ++sample) {
      const auto current = photo.ShadowCurrent(profile, bank.At(sample), q[0], q[1], q[2], direction);
      mean += current;
      second += std::norm(current);
    }
    mean /= static_cast<double>(bank.Size());
    second /= static_cast<double>(bank.Size());
    const double variance = gra::statistics::UnbiasedComplexVariance(second, mean, bank.Size());
    const auto   stat     = photo.ShadowStat(profile, bank, q[0], q[1], q[2], direction);
    REQUIRE(std::abs(stat.mean - mean) < 2.0e-14);
    REQUIRE(stat.second == Approx(std::norm(mean) + variance).epsilon(2.0e-14));
    REQUIRE(stat.variance == Approx(variance).epsilon(2.0e-14));
    const auto factors = photo.Factors(profile, q[0], q[1], q[2], direction, &bank);
    const auto smooth  = photo.Factors(profile, q[0], q[1], q[2], direction);
    REQUIRE(std::abs(factors.coherent - smooth.coherent) < 2.0e-14);
    REQUIRE(factors.incoherent == Approx(std::sqrt(stat.variance)).epsilon(2.0e-14));
  }
  REQUIRE_THROWS_AS(
      photo.ShadowCurrent(profile, bank.At(0), q[0], q[1], q[2], static_cast<gra::nuclear::PhotonDirection>(255)),
      std::invalid_argument);
  REQUIRE_THROWS_AS(photo.Factors(profile, q[0], q[1], q[2], static_cast<gra::nuclear::PhotonDirection>(255)),
                    std::invalid_argument);

  const std::complex<double> elementary(0.3, 1.2);
  const auto                 coherent =
      elementary * photo.Factors(profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ).coherent;
  const auto expected_coherent = elementary * forward.coherent;
  REQUIRE(coherent.real() == Approx(expected_coherent.real()).epsilon(1.0e-12));
  REQUIRE(coherent.imag() == Approx(expected_coherent.imag()).epsilon(1.0e-12));
  const auto incoherent =
      elementary * photo.Factors(profile, 0.16, -0.05, 0.0, gra::nuclear::PhotonDirection::PositiveZ).incoherent;
  REQUIRE(std::abs(incoherent) == Approx(std::abs(elementary) * finite_q.incoherent).epsilon(1.0e-12));
}

TEST_CASE("Photonuclear geometry constructs consistently across threads", "[gra::nuclear][photo][threading]") {
  const auto target = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param  = PhotoControls();
  param.b_nodes     = 41;
  param.z_nodes     = 43;
  const auto                                       profile = gra::test::PhotoProfile(25.0, 0.1);
  constexpr std::size_t                            count   = 4;
  std::array<gra::nuclear::PhotoTransition, count> transition;
  std::array<std::exception_ptr, count>            failure;
  std::vector<std::thread>                         worker;
  worker.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    worker.emplace_back([&, i]() {
      try {
        const gra::nuclear::MPhoto photo(target, param);
        transition[i] = photo.Factors(profile, 0.13, -0.07, 0.04, gra::nuclear::PhotonDirection::PositiveZ);
      } catch (...) { failure[i] = std::current_exception(); }
    });
  }
  for (auto &thread : worker) { thread.join(); }
  for (std::size_t i = 0; i < count; ++i) {
    REQUIRE(failure[i] == nullptr);
    REQUIRE(transition[i].coherent.real() == Approx(transition[0].coherent.real()).epsilon(2.0e-14));
    REQUIRE(transition[i].coherent.imag() == Approx(transition[0].coherent.imag()).epsilon(2.0e-14));
    REQUIRE(transition[i].incoherent == Approx(transition[0].incoherent).epsilon(2.0e-14));
  }
}

TEST_CASE("Leading twist gluon shadowing has stable effective cross sections", "[gra::nuclear][photo][shadow]") {
  auto                        photo = PhotoControls();
  const gra::nuclear::MShadow shadow(photo.shadow);
  const auto                  central = shadow.CrossSections(1.0e-3, 3.0);
  REQUIRE(central.sigma2 > 0.0);
  REQUIRE(central.sigma3 >= central.sigma2);
  REQUIRE(central.sigma3_in > 0.0);
  REQUIRE(central.sigma3_in < central.sigma3);
  REQUIRE(central.ratio == Approx(central.sigma2 / central.sigma3).epsilon(2.0e-14));
  const auto frozen   = shadow.CrossSections(1.0e-3, 2.4);
  const auto boundary = shadow.CrossSections(1.0e-3, photo.shadow.q2_min);
  REQUIRE(frozen.sigma2 == Approx(boundary.sigma2).epsilon(2.0e-14));
  REQUIRE(frozen.sigma3 == Approx(boundary.sigma3).epsilon(2.0e-14));

  // Integrate the forward DPDF in beta to cross-check flux powers and normalization
  gra::MLHAPDFStore store;
  const auto        pdf      = store.GetPDF(photo.shadow.pdf_set, photo.shadow.pdf_member);
  const auto        dpdf     = store.GetPDF(photo.shadow.dpdf_set, photo.shadow.dpdf_member);
  const double      x        = 1.0e-3;
  const auto        rule     = gra::math::GaussLegendreRule(512, x / photo.shadow.x_max, 1.0);
  double            integral = 0.0;
  for (const auto &i : indices(rule.first)) {
    const double beta  = rule.first[i];
    const double xp    = x / beta;
    const double t_min = -std::pow(gra::PDG::mp * xp, 2) / (1.0 - xp);
    const double slope = photo.shadow.flux_b - 2.0 * photo.shadow.alpha_prime * std::log(xp);
    integral += rule.second[i] * x / (beta * beta) * std::pow(xp, 1.0 - 2.0 * photo.shadow.alpha0) *
                std::exp(slope * t_min) * dpdf->xfxQ2(21, beta, 3.0);
  }
  const double eta      = std::tan(0.5 * gra::math::PI * (photo.shadow.alpha0 - 1.0));
  const double expected = 16.0 * gra::math::PI * integral * gra::PDG::GeV2mb /
                          ((1.0 + eta * eta) * photo.shadow.dpdf_diss * pdf->xfxQ2(21, x, 3.0));
  REQUIRE(central.sigma2 == Approx(expected).epsilon(0.003));

  // Intact-proton normalization rescales both moments and preserves their ratio
  auto rescaled_param = photo.shadow;
  rescaled_param.dpdf_diss *= 1.5;
  const gra::nuclear::MShadow rescaled(rescaled_param);
  for (const double x : {1.0e-5, 1.0e-3, 0.05}) {
    const auto original = shadow.CrossSections(x, 3.0);
    const auto value    = rescaled.CrossSections(x, 3.0);
    REQUIRE(1.5 * value.sigma2 == Approx(original.sigma2).epsilon(2.0e-13));
    REQUIRE(1.5 * value.sigma3 == Approx(original.sigma3).epsilon(2.0e-13));
    REQUIRE(value.ratio == Approx(original.ratio).epsilon(2.0e-13));
  }
  auto invalid_norm      = photo.shadow;
  invalid_norm.dpdf_diss = 0.5;
  REQUIRE_THROWS_AS(gra::nuclear::MShadow(invalid_norm), std::invalid_argument);

  auto order_mismatch     = photo.shadow;
  order_mismatch.dpdf_set = "GKG18_DPDF_FitB_LO";
  REQUIRE_THROWS_AS(gra::nuclear::MShadow(order_mismatch), std::invalid_argument);
  auto grid_mismatch   = photo.shadow;
  grid_mismatch.q2_min = 0.5;
  REQUIRE_THROWS_AS(gra::nuclear::MShadow(grid_mismatch), std::invalid_argument);
  for (const double x : {1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2}) {
    const auto value = shadow.CrossSections(x, 3.0);
    CAPTURE(x, value.sigma2, value.sigma3, value.ratio);
    REQUIRE(value.sigma2 >= 0.0);
    REQUIRE(value.sigma3 >= value.sigma2);
    REQUIRE(value.ratio == Approx(value.sigma2 / value.sigma3).epsilon(2.0e-14));
  }
  for (const double scale2 : {2.4, 2.7225, 3.0, 4.0}) {
    for (const double x :
         {2.5e-4, 4.0e-4, 6.0e-4, 8.0e-4, 1.0e-3, 1.3e-3, 1.6e-3, 2.0e-3, 5.0e-2, 8.0e-2, 9.0e-2, 9.8e-2}) {
      CAPTURE(x, scale2);
      const auto value = shadow.CrossSections(x, scale2);
      REQUIRE(value.sigma2 <= value.sigma3);
    }
  }
  const double log_q2_min = std::log10(photo.shadow.q2_min);
  for (std::size_t iq = 0; iq <= 64; ++iq) {
    const double log_scale2 =
        log_q2_min + (std::log10(photo.shadow.q2_max) - log_q2_min) * static_cast<double>(iq) / 64.0;
    const double scale2 = iq == 0 ? photo.shadow.q2_min : iq == 64 ? photo.shadow.q2_max : std::pow(10.0, log_scale2);
    for (std::size_t ix = 0; ix < 100; ++ix) {
      const double x = std::pow(10.0, -6.0 + 5.0 * static_cast<double>(ix) / 100.0);
      CAPTURE(x, scale2);
      const auto value = shadow.CrossSections(x, scale2);
      REQUIRE(value.sigma2 <= value.sigma3);
    }
  }
  for (const double x : {1.0e-6, 3.0e-6, 1.0e-5, 3.0e-5, 1.0e-4}) {
    const auto value = shadow.CrossSections(x, 3.0);
    REQUIRE(value.sigma3 == Approx(value.sigma2).epsilon(2.0e-14));
  }

  auto fine_param           = photo.shadow;
  fine_param.x_nodes        = 261;
  fine_param.scale_nodes    = 129;
  fine_param.integral_nodes = 96;
  const gra::nuclear::MShadow fine(fine_param);
  for (const auto &point :
       std::array<std::array<double, 2>, 4>{{{3.0e-4, 2.4}, {1.0e-3, 3.0}, {4.0e-3, 3.5}, {2.0e-2, 4.0}}}) {
    const auto value     = shadow.CrossSections(point[0], point[1]);
    const auto reference = fine.CrossSections(point[0], point[1]);
    CAPTURE(point[0], point[1], value.sigma2, reference.sigma2);
    REQUIRE(value.sigma2 == Approx(reference.sigma2).epsilon(0.003));
    REQUIRE(value.sigma3 == Approx(reference.sigma3).epsilon(0.006));
  }
}

TEST_CASE("Leading twist target moments preserve the optical normalization",
          "[gra::nuclear][photo][shadow][good-walker]") {
  const auto                 lead        = MakeLead();
  auto                       photo_param = PhotoControls();
  const auto                 profile     = gra::test::PhotoProfile(3.0, 0.1, 0.15);
  const gra::nuclear::MPhoto impulse(lead, photo_param, gra::nuclear::PhotoModel::Impulse);
  const gra::nuclear::MPhoto shadow(lead, photo_param, gra::nuclear::PhotoModel::LTA);
  const auto impulse_forward = impulse.Factors(profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  const auto shadow_forward  = shadow.Factors(profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  REQUIRE(std::abs(shadow_forward.coherent) < 0.9 * std::abs(impulse_forward.coherent));
  REQUIRE(shadow_forward.incoherent > 0.0);

  // Integrate the coherent density and its finite-difference response independently
  const gra::nuclear::MShadow leading(photo_param.shadow);
  const auto                  xs     = leading.CrossSections(profile.x, profile.scale2);
  const auto                  b_rule = gra::math::GaussLegendreRule(192, 0.0, lead.MatterDensity().Param().r_max);
  const auto                  z_rule =
      gra::math::GaussLegendreRule(192, -lead.MatterDensity().Param().r_max, lead.MatterDensity().Param().r_max);
  constexpr double mb_to_fm2         = 0.1;
  const double     a                 = static_cast<double>(lead.A());
  double           expected_mean     = 0.0;
  double           expected_second   = 0.0;
  double           expected_response = 0.0;
  for (const auto &ib : indices(b_rule.first)) {
    double thickness = 0.0;
    for (const auto &iz : indices(z_rule.first)) {
      thickness += z_rule.second[iz] * lead.MatterDensity().Rho(std::hypot(b_rule.first[ib], z_rule.first[iz]));
    }
    const double optical          = a * thickness;
    const double exponent         = 0.5 * mb_to_fm2 * xs.sigma3 * optical;
    const double rescatter        = exponent > 1.0e-12 ? -std::expm1(-exponent) / exponent : 1.0;
    const double first            = 1.0 - xs.ratio + xs.ratio * rescatter;
    const auto   coherent_density = [&](double t) {
      return (1.0 - xs.ratio) * t -
             2.0 * xs.ratio * std::expm1(-0.5 * mb_to_fm2 * xs.sigma3 * a * t) / (mb_to_fm2 * xs.sigma3 * a);
    };
    const double h        = 1.0e-4 * thickness;
    const double response = (coherent_density(thickness + h) - coherent_density(thickness - h)) / (2.0 * h);
    const double second   = response * response;
    const double area     = 2.0 * gra::math::PI * b_rule.second[ib] * b_rule.first[ib];
    expected_mean += a * area * thickness * first;
    expected_second += a * area * thickness * second;
    expected_response += area * thickness * response;
  }
  REQUIRE(shadow_forward.coherent.real() == Approx(expected_mean).epsilon(0.003));
  REQUIRE(std::abs(shadow_forward.coherent.imag()) < 2.0e-14);
  REQUIRE(shadow_forward.incoherent ==
          Approx(std::sqrt(expected_second - a * expected_response * expected_response)).epsilon(0.003));

  auto config     = ConfigControls(256);
  config.d_min    = 0.8;
  const auto bank = SampleBank(lead, config, 9801U);
  const auto sampled =
      shadow.Factors(profile, 0.035, -0.021, 0.002, gra::nuclear::PhotonDirection::PositiveZ, &bank, nullptr, true);
  const auto optical = shadow.Factors(profile, 0.035, -0.021, 0.002, gra::nuclear::PhotonDirection::PositiveZ);
  REQUIRE(sampled.current.has_value());
  const auto &current = *sampled.current;
  REQUIRE(current.sample.size() == bank.Size());
  const double coherent_error =
      std::abs(sampled.coherent - optical.coherent) / std::max(1.0, std::abs(optical.coherent));
  CAPTURE(coherent_error, sampled.incoherent, optical.incoherent);
  REQUIRE(coherent_error < 0.15);
  REQUIRE(std::abs(sampled.coherent - current.stat.mean) < 2.0e-14);
  REQUIRE(gra::math::pow2(sampled.incoherent) == Approx(current.stat.variance).epsilon(2.0e-13));
  REQUIRE(current.stat.second == Approx(std::norm(current.stat.mean) + current.stat.variance).epsilon(2.0e-13));
}

TEST_CASE("Photonuclear target quadrature is stable at production resolution",
          "[gra::nuclear][photo][shadow][convergence]") {
  const auto lead                       = MakeLead();
  auto       production_param           = PhotoControls();
  production_param.b_nodes              = 32;
  production_param.z_nodes              = 32;
  production_param.table.qt_max         = 2.5;
  production_param.table.qt_nodes       = 513;
  auto reference_param                  = production_param;
  production_param.fluctuation.cdf_tol  = 1.0e-8;
  production_param.table.series_abs_tol = 1.0e-6;
  reference_param.b_nodes               = 384;
  reference_param.z_nodes               = 384;
  // Force direct integration independently of the production momentum interpolation
  reference_param.table.qt_max   = 0.01;
  reference_param.table.qz_max   = 0.0001;
  reference_param.table.qt_nodes = reference_param.table.qz_nodes = 4;
  const auto profile                                              = gra::test::PhotoProfile(3.0, 0.1, 0.15);
  for (const auto model : {gra::nuclear::PhotoModel::LTA, gra::nuclear::PhotoModel::Glauber}) {
    const gra::nuclear::MPhoto production(lead, production_param, model);
    const gra::nuclear::MPhoto reference(lead, reference_param, model);
    for (const double qt : {0.0, 0.025, 0.055, 0.090, 0.3, 0.5, 1.0, 2.0, 2.47, 3.1}) {
      for (const double qz : {0.0015, 0.23, 0.31}) {
        const auto value = production.Factors(profile, qt, 0.0, qz, gra::nuclear::PhotonDirection::PositiveZ);
        const auto truth = reference.Factors(profile, qt, 0.0, qz, gra::nuclear::PhotonDirection::PositiveZ);
        CAPTURE(static_cast<int>(model), qt, qz, value.coherent, truth.coherent, value.incoherent, truth.incoherent);
        REQUIRE(std::abs(value.coherent - truth.coherent) < 1.0e-6 + 0.003 * std::abs(truth.coherent));
        REQUIRE(value.incoherent == Approx(truth.incoherent).epsilon(0.003));
      }
    }
    if (model == gra::nuclear::PhotoModel::LTA) {
      for (const double scale2 : {3.0, production_param.shadow.q2_max}) {
        auto strong   = profile;
        strong.x      = production_param.shadow.x_min;
        strong.scale2 = scale2;
        for (const double qt : {0.025, 0.1, 1.0}) {
          const auto value = production.Factors(strong, qt, 0.0, 0.23, gra::nuclear::PhotonDirection::PositiveZ);
          const auto truth = reference.Factors(strong, qt, 0.0, 0.23, gra::nuclear::PhotonDirection::PositiveZ);
          CAPTURE(qt, scale2, value.coherent, truth.coherent);
          REQUIRE(std::abs(value.coherent - truth.coherent) < 1.0e-6 + 0.003 * std::abs(truth.coherent));
          REQUIRE(value.incoherent == Approx(truth.incoherent).epsilon(0.003));
        }
      }
    }
  }
}

// Verify fixed-A closure as diffractive rescattering disappears
TEST_CASE("Leading twist tends to the finite nucleus impulse response", "[gra::nuclear][photo][shadow]") {
  const auto lead  = MakeLead();
  auto       param = PhotoControls();
  const gra::nuclear::MPhoto leading(lead, param, gra::nuclear::PhotoModel::LTA);
  const gra::nuclear::MPhoto impulse(lead, param, gra::nuclear::PhotoModel::Impulse);
  auto                       profile = gra::test::PhotoProfile();
  for (const double x : {param.shadow.x_max * (1.0 - 1.0e-10), param.shadow.x_max, 0.1144, 0.5}) {
    profile.x = x;
    for (const double qt : {0.0, 0.015, 0.1, 0.5}) {
      const auto lta = leading.Factors(profile, qt, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
      const auto ia  = impulse.Factors(profile, qt, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
      CAPTURE(qt, lta.incoherent, ia.incoherent);
      REQUIRE(std::abs(lta.coherent - ia.coherent) < 1.0e-5);
      REQUIRE(std::abs(lta.incoherent - ia.incoherent) < 1.0e-5);
    }
  }
}

// Compare independent nucleons to smooth closure without imposing optical variances
TEST_CASE("Leading twist sampled density response closes to the smooth variance", "[gra::nuclear][photo][shadow]") {
  const auto lead  = MakeLead();
  auto       param = PhotoControls();
  const gra::nuclear::MPhoto photo(lead, param, gra::nuclear::PhotoModel::LTA);
  const auto                 profile = gra::test::PhotoProfile();
  gra::MRandom               random;
  random.SetSeed(753019);
  const std::array<double, 3>                      qt = {0.0, 0.1, 0.5};
  std::array<std::vector<std::complex<double>>, 3> currents;
  for (unsigned int sample = 0; sample < 2048; ++sample) {
    std::vector<gra::nuclear::Nucleon> centers(lead.A());
    for (const auto &i : indices(centers)) {
      const double radius = lead.MatterDensity().Radius(random.U(0.0, 1.0));
      const double cosine = random.U(-1.0, 1.0), phi = random.U(0.0, 2.0 * gra::math::PI);
      const double transverse = radius * std::sqrt(1.0 - cosine * cosine);
      centers[i].x            = {transverse * std::cos(phi), transverse * std::sin(phi), radius * cosine};
      centers[i].type         = i < lead.Z() ? gra::nuclear::NucleonType::Proton : gra::nuclear::NucleonType::Neutron;
    }
    const gra::nuclear::MConfig config(lead.A(), lead.Z(), std::move(centers));
    for (const auto &i : indices(qt)) {
      currents[i].push_back(
          photo.ShadowCurrent(profile, config, qt[i], 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ));
    }
  }
  for (const auto &i : indices(qt)) {
    const auto smooth  = photo.Factors(profile, qt[i], 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
    const auto sampled = gra::nuclear::CurrentMoments(currents[i]);
    CAPTURE(qt[i], sampled.variance, smooth.incoherent * smooth.incoherent);
    REQUIRE(sampled.variance == Approx(smooth.incoherent * smooth.incoherent).epsilon(0.10));
  }
}

TEST_CASE("Photonuclear attenuation follows the incident photon direction", "[gra::nuclear][photo][direction]") {
  const auto                 oxygen  = MakeOxygen();
  auto                       param   = PhotoControls();
  const auto                 profile = gra::test::PhotoProfile(25.0, 0.0);
  const gra::nuclear::MPhoto photo(oxygen, param);
  const gra::nuclear::MPhoto leading(oxygen, param, gra::nuclear::PhotoModel::LTA);

  std::vector<gra::nuclear::Nucleon> nucleon(oxygen.A());
  for (const auto &i : indices(nucleon)) {
    const double fraction = static_cast<double>(i) / static_cast<double>(nucleon.size() - 1);
    nucleon[i].x    = {0.12 * static_cast<double>(i % 3), 0.09 * static_cast<double>(i % 5), -2.0 + 4.0 * fraction};
    nucleon[i].type = i < oxygen.Z() ? gra::nuclear::NucleonType::Proton : gra::nuclear::NucleonType::Neutron;
  }
  const gra::nuclear::MConfig front_loaded(oxygen.A(), oxygen.Z(), nucleon);
  for (const auto &produced : indices(nucleon)) {
    const auto downstream = front_loaded.DownstreamR2(produced, true);
    const auto upstream   = front_loaded.DownstreamR2(produced, false);
    REQUIRE(downstream.size() + upstream.size() == nucleon.size() - 1);
    for (const double radius2 : downstream) {
      REQUIRE(std::isfinite(radius2));
      REQUIRE(radius2 >= 0.0);
    }
    for (const double radius2 : upstream) {
      REQUIRE(std::isfinite(radius2));
      REQUIRE(radius2 >= 0.0);
    }
  }
  constexpr double qz              = 0.12;
  constexpr double angle           = 0.63;
  constexpr double qx              = 0.14;
  constexpr double qy              = -0.08;
  auto             rotated_nucleon = nucleon;
  for (auto &item : rotated_nucleon) {
    const double x = item.x[0];
    const double y = item.x[1];
    item.x[0]      = std::cos(angle) * x - std::sin(angle) * y;
    item.x[1]      = std::sin(angle) * x + std::cos(angle) * y;
  }
  const gra::nuclear::MConfig rotated(oxygen.A(), oxygen.Z(), rotated_nucleon);
  const double                rotated_qx = std::cos(angle) * qx - std::sin(angle) * qy;
  const double                rotated_qy = std::sin(angle) * qx + std::cos(angle) * qy;
  for (auto &item : nucleon) { item.x[2] = -item.x[2]; }
  const gra::nuclear::MConfig reflected(oxygen.A(), oxygen.Z(), nucleon);

  for (const auto *model : {&photo, &leading}) {
    const auto upper_target =
        model->ShadowCurrent(profile, front_loaded, 0.0, 0.0, qz, gra::nuclear::TargetPhotonDirection(1));
    const auto lower_target =
        model->ShadowCurrent(profile, front_loaded, 0.0, 0.0, qz, gra::nuclear::TargetPhotonDirection(2));
    for (const auto direction : {gra::nuclear::PhotonDirection::PositiveZ, gra::nuclear::PhotonDirection::NegativeZ}) {
      const auto original    = model->ShadowCurrent(profile, front_loaded, qx, qy, qz, direction);
      const auto transformed = model->ShadowCurrent(profile, rotated, rotated_qx, rotated_qy, qz, direction);
      REQUIRE(transformed.real() == Approx(original.real()).margin(2.0e-13));
      REQUIRE(transformed.imag() == Approx(original.imag()).margin(2.0e-13));
    }
    const auto reflected_lower =
        model->ShadowCurrent(profile, reflected, 0.0, 0.0, -qz, gra::nuclear::TargetPhotonDirection(2));
    if (model == &leading) {
      REQUIRE(std::abs(lower_target - upper_target) < 2.0e-12);
    } else {
      REQUIRE(std::abs(lower_target - upper_target) > 1.0e-6);
    }
    REQUIRE(std::abs(reflected_lower - upper_target) < 2.0e-12);
  }

  CHECK(gra::nuclear::TargetPhotonDirection(1) == gra::nuclear::PhotonDirection::NegativeZ);
  CHECK(gra::nuclear::TargetPhotonDirection(2) == gra::nuclear::PhotonDirection::PositiveZ);
  CHECK_THROWS_AS(gra::nuclear::TargetPhotonDirection(0), std::out_of_range);
}

TEST_CASE("UPC screening is inactive for lepton beams and controlled for pA and AA",
          "[gra::nuclear][upc][mixed-beam]") {
  const auto oxygen           = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param            = UPCControls();
  param.glauber               = GlauberControls(0.1);
  param.loop.r_max            = 1.0;
  param.loop.radial_intervals = 16;
  param.loop.azimuth_nodes    = 16;
  param.photo[0].b_nodes      = 32;
  param.photo[1].b_nodes      = 32;

  auto                     lepton_param = param;
  const gra::nuclear::MUPC electron_ion({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus},
                                        {nullptr, oxygen}, lepton_param);
  REQUIRE_FALSE(electron_ion.HasHadronicPair());
  REQUIRE(electron_ion.Glauber() == nullptr);
  REQUIRE(electron_ion.Nodes({}).empty());
  REQUIRE(electron_ion.SurvivalAmp(0.0) == Approx(1.0).margin(2.0e-14));

  const gra::nuclear::MUPC electron_proton({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Proton},
                                           {nullptr, nullptr}, lepton_param);
  REQUIRE_FALSE(electron_proton.HasHadronicPair());
  REQUIRE(electron_proton.Nodes({}).empty());
  REQUIRE(electron_proton.Glauber() == nullptr);

  const gra::nuclear::MUPC electron_electron({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Lepton},
                                             {nullptr, nullptr}, lepton_param);
  REQUIRE_FALSE(electron_electron.HasHadronicPair());
  REQUIRE(electron_electron.Nodes({}).empty());
  REQUIRE(electron_electron.Glauber() == nullptr);

  const gra::nuclear::MUPC proton_ion({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus},
                                      {nullptr, oxygen}, param);
  REQUIRE(proton_ion.HasHadronicPair());
  REQUIRE(proton_ion.HadronicConvolution());
  REQUIRE(proton_ion.Glauber() != nullptr);
  REQUIRE_FALSE(proton_ion.Nodes({}).empty());
  REQUIRE(proton_ion.ProfileTransform(0.0) > 0.0);
  RequireProbability(proton_ion.Glauber()->GGCFProb(4.0));

  gra::nuclear::UPCMode bare_mode;
  bare_mode.screening = false;
  const gra::nuclear::MUPC bare_proton_ion({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus},
                                           {nullptr, oxygen}, param, {}, bare_mode);
  REQUIRE(bare_proton_ion.HasHadronicPair());
  REQUIRE_FALSE(bare_proton_ion.HadronicConvolution());
  REQUIRE(bare_proton_ion.Nodes({}).empty());
  gra::nuclear::FinalState bare_final;
  HepMC3::GenEvent         bare_event(HepMC3::Units::GEV, HepMC3::Units::MM);
  gra::nuclear::AttachRecord(bare_proton_ion, bare_final, bare_event);
  REQUIRE(bare_event.heavy_ion() != nullptr);
  CHECK(bare_event.heavy_ion()->sigma_inel_NN == Approx(bare_proton_ion.Param().glauber.profile.sigma));

  const gra::nuclear::MUPC ion_ion(*oxygen, *oxygen, param);
  REQUIRE(ion_ion.HasHadronicPair());
  REQUIRE_FALSE(ion_ion.Nodes({}).empty());
  REQUIRE(ion_ion.ProfileTransform(0.0) > 0.0);
  REQUIRE(ion_ion.ProfileTransform(0.0) > proton_ion.ProfileTransform(0.0));
}

TEST_CASE("Smooth UPC profiles use event-centered loop quadrature",
          "[gra::nuclear][upc][optical][quadrature][rotation]") {
  const auto oxygen           = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param            = UPCControls();
  param.structure             = gra::nuclear::StructureType::Smooth;
  param.config                = {};
  param.survival              = gra::nuclear::SurvivalType::Optical;
  param.glauber               = GlauberControls(0.1);
  param.loop.r_max            = 0.3;
  param.loop.radial_intervals = 12;
  param.loop.azimuth_nodes    = 12;
  const gra::nuclear::MUPC upc({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus}, {nullptr, oxygen},
                               param);

  REQUIRE(upc.HadronicConvolution());
  REQUIRE(upc.Nodes({}).size() == 12 * 12);

  const std::array<gra::M3Vec, 2> transfer = {gra::M3Vec{0.08, -0.05, 0.02}, gra::M3Vec{-0.06, 0.04, -0.015}};
  constexpr double                angle    = 0.37;
  // Rotate one transverse momentum while preserving its longitudinal part
  const auto rotate = [](const gra::M3Vec &q) {
    return gra::M3Vec{std::cos(angle) * q[0] - std::sin(angle) * q[1], std::sin(angle) * q[0] + std::cos(angle) * q[1],
                      q[2]};
  };
  const auto node    = upc.Nodes(transfer);
  const auto rotated = upc.Nodes({rotate(transfer[0]), rotate(transfer[1])});
  REQUIRE(node.size() == 3 * upc.Nodes({}).size());
  REQUIRE(rotated.size() == node.size());

  bool                 focused     = false;
  std::complex<double> fixed_sum   = 0.0;
  std::complex<double> focused_sum = 0.0;
  for (const auto &item : upc.Nodes({})) { fixed_sum += item.weight; }
  for (const auto &i : indices(node)) {
    focused = focused || i >= upc.Nodes({}).size() || std::abs(node[i].kt - upc.Nodes({})[i].kt) > 1.0e-12 ||
              std::abs(node[i].phi - upc.Nodes({})[i].phi) > 1.0e-12;
    focused_sum += node[i].weight;
    REQUIRE(rotated[i].kt == Approx(node[i].kt).margin(1.0e-14));
    REQUIRE(rotated[i].kx == Approx(std::cos(angle) * node[i].kx - std::sin(angle) * node[i].ky).margin(2.0e-13));
    REQUIRE(rotated[i].ky == Approx(std::sin(angle) * node[i].kx + std::cos(angle) * node[i].ky).margin(2.0e-13));
    REQUIRE(rotated[i].weight.real() == Approx(node[i].weight.real()).margin(2.0e-13));
    REQUIRE(rotated[i].weight.imag() == Approx(node[i].weight.imag()).margin(2.0e-13));
  }
  REQUIRE(focused);
  auto reference_param                  = param;
  reference_param.loop.radial_intervals = 48;
  reference_param.loop.azimuth_nodes    = 48;
  const gra::nuclear::MUPC reference({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus},
                                     {nullptr, oxygen}, reference_param);
  std::complex<double>     reference_fixed   = 0.0;
  std::complex<double>     reference_focused = 0.0;
  std::complex<double>     focused_current   = 0.0;
  std::complex<double>     reference_current = 0.0;
  // Resolve both shifted coherent-current peaks in the same 2D integrand
  const auto current = [&](const gra::nuclear::LoopNode &item) {
    const double upper = gra::math::pow2(item.kx + transfer[0][0]) + gra::math::pow2(item.ky + transfer[0][1]) +
                         gra::math::pow2(transfer[0][2]);
    const double lower = gra::math::pow2(item.kx - transfer[1][0]) + gra::math::pow2(item.ky - transfer[1][1]) +
                         gra::math::pow2(transfer[1][2]);
    return 1.0 / (upper * lower);
  };
  for (const auto &item : node) { focused_current += item.weight * current(item); }
  for (const auto &item : FixedNodes(reference)) {
    reference_fixed += item.weight;
    reference_current += item.weight * current(item);
  }
  for (const auto &item : reference.Nodes(transfer)) { reference_focused += item.weight; }
  CAPTURE(fixed_sum, focused_sum, reference_fixed, reference_focused, focused_current, reference_current);
  REQUIRE(reference_focused.real() == Approx(reference_fixed.real()).epsilon(0.001));
  REQUIRE(focused_sum.real() == Approx(reference_fixed.real()).epsilon(0.03));
  REQUIRE(focused_sum.imag() == Approx(fixed_sum.imag()).margin(2.0e-13));
  REQUIRE(focused_current.real() == Approx(reference_current.real()).epsilon(0.03));
  REQUIRE(focused_current.imag() == Approx(reference_current.imag()).margin(2.0e-10));

  gra::nuclear::ScreenLayout layout;
  layout.type        = gra::nuclear::ScreenType::Scalar;
  layout.fixed.valid = true;
  gra::nuclear::ScreenPoint born;
  born.amplitude = {{2.0, -1.0}};
  born.transfer  = transfer;
  gra::nuclear::MUPCScreen screen(upc, layout, born);
  for (const auto &item : node) { screen.Add(item, born); }
  const auto                 result   = screen.Result();
  const std::complex<double> expected = (std::complex<double>(1.0, 0.0) + focused_sum) * born.amplitude[0];
  REQUIRE(result.amplitude.size() == 1);
  REQUIRE(result.amplitude[0].real() == Approx(expected.real()).epsilon(1.0e-13));
  REQUIRE(result.amplitude[0].imag() == Approx(expected.imag()).epsilon(1.0e-13));

  // Keep distinct MeV-width EPA poles in separate local channels
  const std::array<gra::M3Vec, 2> narrow_transfer       = {gra::M3Vec{0.020, 0.002, 0.001},
                                                           gra::M3Vec{-0.039, -0.001, -0.001}};
  const auto                      narrow_node           = upc.Nodes(narrow_transfer);
  const auto                      narrow_reference_node = reference.Nodes(narrow_transfer);
  REQUIRE(narrow_node.size() == 3 * upc.Nodes({}).size());
  REQUIRE(narrow_reference_node.size() == 3 * reference.Nodes({}).size());

  // |J_1T J_2T| = q'_1T q'_2T / [(q'^2_1T + q^2_1z)(q'^2_2T + q^2_2z)]
  const auto current_norm = [&](const gra::nuclear::LoopNode &item) {
    const double upper_x = item.kx + narrow_transfer[0][0];
    const double upper_y = item.ky + narrow_transfer[0][1];
    const double lower_x = narrow_transfer[1][0] - item.kx;
    const double lower_y = narrow_transfer[1][1] - item.ky;
    const double upper2  = upper_x * upper_x + upper_y * upper_y;
    const double lower2  = lower_x * lower_x + lower_y * lower_y;
    return std::sqrt(upper2 * lower2) /
           ((upper2 + gra::math::pow2(narrow_transfer[0][2])) * (lower2 + gra::math::pow2(narrow_transfer[1][2])));
  };
  std::complex<double> narrow_current           = 0.0;
  std::complex<double> narrow_reference_current = 0.0;
  for (const auto &item : narrow_node) { narrow_current += item.weight * current_norm(item); }
  for (const auto &item : narrow_reference_node) { narrow_reference_current += item.weight * current_norm(item); }
  CAPTURE(narrow_current, narrow_reference_current);
  REQUIRE(narrow_current.real() == Approx(narrow_reference_current.real()).epsilon(0.04));
  REQUIRE(narrow_current.imag() == Approx(narrow_reference_current.imag()).margin(2.0e-10));
}

TEST_CASE("UPC momentum screening matches an independent impact-space integral",
          "[gra::nuclear][upc][screen][fourier][physics]") {
  const auto oxygen                 = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param                  = UPCControls();
  param.structure                   = gra::nuclear::StructureType::Smooth;
  param.config                      = {};
  param.survival                    = gra::nuclear::SurvivalType::Optical;
  param.glauber                     = GlauberControls(0.0);
  param.glauber.b_nodes             = 128;
  param.convolution.smooth_b_nodes  = 128;
  param.convolution.sample_b_nodes  = 128;
  param.convolution.smooth_kt_nodes = 128;
  param.convolution.sample_kt_nodes = 128;
  param.loop.r_max                  = 0.4;
  param.loop.radial_intervals       = 64;
  param.loop.azimuth_nodes          = 16;
  const gra::nuclear::MUPC upc(*oxygen, *oxygen, param);

  // A(k) = exp(-a k^2) transforms to
  // a(b) = exp[-b^2/(4 a hbarc^2)]/(4 pi a)
  constexpr double           slope = 250.0;
  const std::complex<double> born  = std::polar(1.0, 0.71);
  gra::nuclear::ScreenLayout layout;
  layout.type        = gra::nuclear::ScreenType::Scalar;
  layout.fixed.valid = true;
  gra::nuclear::ScreenPoint point;
  point.amplitude = {born};
  gra::nuclear::MUPCScreen screen(upc, layout, point);
  for (const auto &node : upc.Nodes({})) {
    point.amplitude[0] = born * std::exp(-slope * node.kt * node.kt);
    screen.Add(node, point);
  }
  const auto result = screen.Result();

  // Fourier inversion gives int d^2b S(b) a(b)/hbarc^2 directly
  const auto [b_node, b_weight] = gra::math::GaussLegendreRule(256, 0.0, param.convolution.b_max);
  constexpr double hbarc        = gra::PDG::GeV2fm;
  const double     width2       = 4.0 * slope * hbarc * hbarc;
  double           direct       = 0.0;
  for (const auto &i : indices(b_node)) {
    direct += b_weight[i] * b_node[i] * std::exp(-b_node[i] * b_node[i] / width2) * upc.SurvivalAmp(b_node[i]) /
              (2.0 * slope * hbarc * hbarc);
  }
  direct += std::exp(-param.convolution.b_max * param.convolution.b_max / width2);

  REQUIRE(result.amplitude.size() == 1);
  REQUIRE(std::abs(result.amplitude[0] - direct * born) < 0.015);
  const std::complex<double> loop       = result.amplitude[0] - born;
  const double               decomposed = std::norm(born) + 2.0 * std::real(std::conj(born) * loop) + std::norm(loop);
  REQUIRE(result.helicity_norm[0] == Approx(decomposed).epsilon(2.0e-14));
  REQUIRE(std::real(std::conj(born) * loop) < 0.0);
}

TEST_CASE("Shifted Gaussian EPA currents match their exact impact-space convolution",
          "[gra::nuclear][upc][screen][fourier][rotation][current]") {
  const auto oxygen                 = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param                  = UPCControls();
  param.structure                   = gra::nuclear::StructureType::Smooth;
  param.config                      = {};
  param.survival                    = gra::nuclear::SurvivalType::Optical;
  param.glauber                     = GlauberControls(0.0);
  param.glauber.b_nodes             = 192;
  param.convolution.smooth_b_nodes  = 192;
  param.convolution.sample_b_nodes  = 192;
  param.convolution.smooth_kt_nodes = 256;
  param.convolution.sample_kt_nodes = 256;
  param.loop.r_min                  = 0.0;
  param.loop.r_max                  = 0.5;
  param.loop.radial_intervals       = 12;
  param.loop.azimuth_nodes          = 12;
  const gra::nuclear::MUPC upc(*oxygen, *oxygen, param);
  auto                     reference_param = param;
  reference_param.loop.radial_intervals    = 48;
  reference_param.loop.azimuth_nodes       = 48;
  const gra::nuclear::MUPC reference(*oxygen, *oxygen, reference_param);

  const std::array<gra::M3Vec, 2>           transfer = {gra::M3Vec{0.02000, 0.00200, 0.00010},
                                                        gra::M3Vec{-0.02035, -0.00165, -0.00100}};
  constexpr std::array<double, 2>           slope    = {60.0, 85.0};
  const auto                                screened = GaussianCurrentScreen(reference, transfer, slope);
  const auto                                plus     = GaussianCurrentClosure(upc, transfer, slope, 1);
  const auto                                minus    = GaussianCurrentClosure(upc, transfer, slope, -1);
  const std::array<std::complex<double>, 4> expected = {plus[0], plus[1], minus[0], minus[1]};
  REQUIRE(screened.size() == expected.size());
  for (const auto &i : indices(expected)) {
    CAPTURE(i, screened[i], expected[i]);
    REQUIRE(std::abs(screened[i] - expected[i]) < 0.004 * std::max(0.01, std::abs(expected[i])));
  }

  // The 12 x 12 rule resolves two nearly coincident sub-MeV EPA poles
  const auto epa_low       = GaussianEPAScreen(upc, transfer, slope);
  const auto epa_reference = GaussianEPAScreen(reference, transfer, slope);
  REQUIRE(epa_low.size() == epa_reference.size());
  for (const auto &i : indices(epa_low)) {
    CAPTURE(i, epa_low[i], epa_reference[i]);
    REQUIRE(std::abs(epa_low[i] - epa_reference[i]) < 0.04 * std::max(1.0, std::abs(epa_reference[i])));
  }

  // A current pole at qT = 0 remains distinct from the soft hadronic channel
  const std::array<gra::M3Vec, 2> origin_transfer = {gra::M3Vec{0.0, 0.0, 0.00010}, transfer[1]};
  REQUIRE(upc.Nodes(origin_transfer).size() == 3 * upc.Nodes({}).size());

  constexpr double angle = 0.63;
  // Rotate one transfer while preserving its longitudinal current scale
  const auto rotate = [](const gra::M3Vec &q) {
    return gra::M3Vec{std::cos(angle) * q[0] - std::sin(angle) * q[1], std::sin(angle) * q[0] + std::cos(angle) * q[1],
                      q[2]};
  };
  const std::array<gra::M3Vec, 2> rotated_transfer = {rotate(transfer[0]), rotate(transfer[1])};
  const auto                      rotated          = GaussianCurrentScreen(reference, rotated_transfer, slope);
  for (const auto &i : indices(screened)) {
    const int  helicity = i < 2 ? 1 : -1;
    const auto phase    = std::polar(1.0, static_cast<double>(helicity) * angle);
    CAPTURE(i, screened[i], rotated[i], phase);
    REQUIRE(std::abs(rotated[i] - phase * screened[i]) < 2.0e-11 * std::max(0.01, std::abs(screened[i])));
  }

  const auto                   epa_rotated = GaussianEPAScreen(upc, rotated_transfer, slope);
  constexpr std::array<int, 4> harmonic    = {2, 0, 0, -2};
  for (const auto &i : indices(epa_low)) {
    const auto phase = std::polar(1.0, static_cast<double>(harmonic[i]) * angle);
    CAPTURE(i, epa_low[i], epa_rotated[i], phase);
    REQUIRE(std::abs(epa_rotated[i] - phase * epa_low[i]) < 2.0e-11 * std::max(1.0, std::abs(epa_low[i])));
  }
}

TEST_CASE("Configuration screening retains correlated charge-current samples",
          "[gra::nuclear][upc][good-walker][config]") {
  const auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       config = ConfigControls(4);
  config.d_min      = 0.4;
  const auto bank   = SampleBankPtr(*oxygen, config, 8128);

  auto param                       = UPCControls();
  param.survival                   = gra::nuclear::SurvivalType::MCGGCF;
  param.config                     = config;
  param.glauber                    = GlauberControls(0.1);
  param.glauber.b_nodes            = 32;
  param.convolution.smooth_b_nodes = 32;
  param.convolution.sample_b_nodes = 32;
  param.loop.r_min                 = 0.01;
  param.loop.r_max                 = 0.2;
  param.loop.radial_intervals      = 8;
  param.loop.azimuth_nodes         = 4;
  param.photo[1].b_nodes           = 32;
  const gra::nuclear::MUPC coherent({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus},
                                    {nullptr, oxygen}, param, {nullptr, bank});
  REQUIRE(coherent.HasSamples());
  const std::size_t sample_count = coherent.SampleCount();
  REQUIRE(sample_count == config.count);
  const auto coherent_node = coherent.Nodes({});
  REQUIRE_FALSE(coherent_node.empty());
  bool nonzero_weight   = false;
  bool sample_variation = false;
  for (const auto &node : coherent_node) {
    REQUIRE(node.sample_weight.size() == sample_count);
    std::complex<double> mean_weight = 0.0;
    for (const auto &sample : indices(node.sample_weight)) {
      const auto weight = node.sample_weight[sample];
      REQUIRE(std::isfinite(weight.real()));
      REQUIRE(std::isfinite(weight.imag()));
      mean_weight += weight;
      nonzero_weight = nonzero_weight || std::abs(weight) > 0.0;
      if (sample > 0) { sample_variation = sample_variation || std::abs(weight - node.sample_weight[0]) > 1.0e-14; }
    }
    mean_weight /= static_cast<double>(sample_count);
    const double weight_scale = std::max(1.0, std::abs(node.weight));
    REQUIRE(std::abs(mean_weight - node.weight) < 2.0e-11 * weight_scale);
  }
  REQUIRE(nonzero_weight);
  REQUIRE(sample_variation);
  REQUIRE(coherent_node.size() >= param.loop.azimuth_nodes);
  bool azimuth_variation = false;
  for (std::size_t phi = 1; phi < param.loop.azimuth_nodes; ++phi) {
    for (std::size_t sample = 0; sample < sample_count; ++sample) {
      azimuth_variation = azimuth_variation ||
                          std::abs(coherent_node[phi * gra::math::PolarNodeCount(param.loop)].sample_weight[sample] -
                                   coherent_node[0].sample_weight[sample]) > 1.0e-14;
    }
  }
  REQUIRE(azimuth_variation);

  const gra::M3Vec q        = {0.11, -0.04, 0.07};
  const auto       selected = coherent.EmissionRatios(2, gra::nuclear::CoherenceType::Coherent, q);
  const auto       resolved = coherent.EmissionComponents(2, q);
  REQUIRE(selected.size() == resolved[0].size());
  for (const auto &sample : indices(selected)) { REQUIRE(std::abs(selected[sample] - resolved[0][sample]) < 1.0e-15); }

  const auto           coherent_ratio = coherent.EmissionRatios(2, gra::nuclear::CoherenceType::Coherent, q);
  std::complex<double> coherent_mean  = 0.0;
  for (std::size_t sample = 0; sample < sample_count; ++sample) { coherent_mean += coherent_ratio[sample]; }
  coherent_mean /= static_cast<double>(sample_count);
  REQUIRE(coherent_mean.real() == Approx(1.0).epsilon(2.0e-12));
  REQUIRE(coherent_mean.imag() == Approx(0.0).margin(2.0e-12));

  param.emission[1] = gra::nuclear::CoherenceType::Incoherent;
  const gra::nuclear::MUPC incoherent({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus},
                                      {nullptr, oxygen}, param, {nullptr, bank});
  const auto               incoherent_ratio =
      incoherent.EmissionRatios(2, gra::nuclear::CoherenceType::Incoherent, {0.11, -0.04, 0.07});
  std::complex<double> fluctuation_mean   = 0.0;
  double               fluctuation_second = 0.0;
  for (std::size_t sample = 0; sample < sample_count; ++sample) {
    const auto ratio = incoherent_ratio[sample];
    fluctuation_mean += ratio;
    fluctuation_second += std::norm(ratio);
  }
  fluctuation_mean /= static_cast<double>(sample_count);
  fluctuation_second /= static_cast<double>(sample_count);
  REQUIRE(std::abs(fluctuation_mean) < 2.0e-12);
  REQUIRE(fluctuation_second == Approx(1.0).epsilon(2.0e-12));
}

// Check ensemble convergence without imposing seed-dependent relative precision
TEST_CASE("Configuration survival ensemble converges at the nuclear edge", "[gra::nuclear][upc][convergence][config]") {
  const auto                         oxygen        = MakeOxygen();
  auto                               glauber_param = GlauberControls(0.15);
  glauber_param.fluctuation.nodes                  = 8;
  const gra::nuclear::MGlauber       glauber(oxygen, oxygen, glauber_param);
  constexpr std::array<std::size_t, 3> counts = {8, 32, 64};
  std::array<std::vector<double>, counts.size()> coherent, incoherent;
  for (const auto &i : indices(counts)) {
    for (std::uint32_t replica = 0; replica < 32; ++replica) {
      const std::uint32_t seed = 1009U + 7919U * replica;
      CAPTURE(counts[i], seed);
      const auto moment = SurvivalSectors(glauber, oxygen, counts[i], seed, 6.0);
      coherent[i].push_back(moment.sector[0][0]);
      incoherent[i].push_back(moment.sector[0][1] + moment.sector[1][0] + moment.sector[1][1]);
      for (const auto &upper : moment.sector) {
        for (const double sector : upper) { RequireProbability(sector); }
      }
      REQUIRE(moment.probability == Approx(moment.sector[0][0]).epsilon(2.0e-12));
      REQUIRE(gra::Sum(moment.sector[0]) + gra::Sum(moment.sector[1]) == Approx(moment.mean_second).epsilon(2.0e-12));
    }
  }

  // Four times more independent configurations should halve the leading RMS error
  // Require only sqrt(2) improvement to allow finite-ensemble and nonlinear corrections
  for (const auto *sector : {&coherent, &incoherent}) {
    const double error8   = EnsembleError((*sector)[0], sector->back());
    const double error32  = EnsembleError((*sector)[1], sector->back());
    const double spread8  = SeedSpread((*sector)[0]);
    const double spread32 = SeedSpread((*sector)[1]);
    CAPTURE(error8, error32, spread8, spread32);
    REQUIRE(error32 < error8 / std::sqrt(2.0));
    REQUIRE(spread32 < spread8 / std::sqrt(2.0));
  }

  const auto first  = SurvivalSectors(glauber, oxygen, 16, 3011U, 6.0);
  const auto replay = SurvivalSectors(glauber, oxygen, 16, 3011U, 6.0);
  REQUIRE(first.probability == Approx(replay.probability).margin(1.0e-15));
  REQUIRE(first.mean_second == Approx(replay.mean_second).margin(1.0e-15));
  for (const auto &upper : indices(first.sector)) {
    for (const auto &lower : indices(first.sector[upper])) {
      REQUIRE(first.sector[upper][lower] == Approx(replay.sector[upper][lower]).margin(1.0e-15));
    }
  }
}

// Compare equilibration with paired seed uncertainties for correlated nucleon pairs
TEST_CASE("Hard-core geometry is stable from eight to thirty-two sweeps", "[gra::nuclear][upc][convergence][config]") {
  const auto lead = MakeLead();
  constexpr std::array<double GeometryMoment::*, 3> field = {
      &GeometryMoment::radius2, &GeometryMoment::pair12, &GeometryMoment::pair20};
  std::array<gra::statistics::RunningMoments, field.size()> difference, reference;
  for (std::uint32_t replica = 0; replica < 32; ++replica) {
    const std::uint32_t seed = 1009U + 7919U * replica;
    auto param   = ConfigControls(32);
    param.d_min  = 0.8;
    param.sweeps = 8;
    const auto short_run = ConfigGeometry(SampleBank(lead, param, seed));
    param.sweeps = 32;
    const auto long_run = ConfigGeometry(SampleBank(lead, param, seed));
    for (const auto &i : indices(field)) {
      difference[i].Add(short_run.*field[i] - long_run.*field[i]);
      reference[i].Add(long_run.*field[i]);
    }
  }
  for (const auto &i : indices(field)) {
    const auto &delta = difference[i];
    const double n = delta.Count();
    const double error = std::sqrt(delta.M2() / (n * (n - 1.0)));
    CAPTURE(i, delta.Mean(), reference[i].Mean(), error);
    // Retain the 0.6 percent equilibration target and allow five standard errors
    REQUIRE(std::abs(delta.Mean()) < 0.006 * std::abs(reference[i].Mean()) + 5.0 * error);
  }
}

TEST_CASE("Nuclear fusion quadrature is stable in momentum and azimuth",
          "[gra::nuclear][upc][convergence][quadrature]") {
  // Resolve the Fourier transform before varying the convolution quadrature
  constexpr unsigned int transform_nodes = 1024;

  constexpr std::uint32_t seed          = 8128U;
  constexpr std::size_t   count         = 8;
  const auto              low_k         = ScreenMoment(ConfigUPC(count, seed, 32, 24, transform_nodes));
  const auto              medium        = ScreenMoment(ConfigUPC(count, seed, 64, 24, transform_nodes));
  const auto              k_reference   = ScreenMoment(ConfigUPC(count, seed, 128, 24, transform_nodes));
  const auto              low_phi       = ScreenMoment(ConfigUPC(count, seed, 64, 12, transform_nodes));
  const auto              phi_reference = ScreenMoment(ConfigUPC(count, seed, 64, 96, transform_nodes));
  const auto              reference     = ScreenMoment(ConfigUPC(count, seed, 256, 96, transform_nodes));
  for (const auto &moment : {low_k, medium, k_reference, low_phi, phi_reference, reference}) {
    RequireScreenMoment(moment);
    REQUIRE(moment.survival == Approx(medium.survival).epsilon(2.0e-13));
  }

  const double low_k_error      = ScreenDistance(low_k, k_reference);
  const double medium_k_error   = ScreenDistance(medium, k_reference);
  const double low_phi_error    = ScreenDistance(low_phi, phi_reference);
  const double medium_phi_error = ScreenDistance(medium, phi_reference);
  const double reference_error  = ScreenDistance(medium, reference);
  CAPTURE(low_k_error, medium_k_error, low_phi_error, medium_phi_error, reference_error);
  REQUIRE(medium_k_error < low_k_error);
  REQUIRE(medium_k_error < 0.002);
  REQUIRE(medium_phi_error < low_phi_error);
  REQUIRE(medium_phi_error < 0.004);
  REQUIRE(reference_error < 0.005);
}

TEST_CASE("AA configurations use the complete Cartesian ensemble", "[gra::nuclear][upc][config][cartesian]") {
  constexpr std::size_t            count = 5;
  const auto                       upc   = ConfigUPC(count, 8128U, 8, 4);
  const std::array<std::size_t, 2> shape = {count, count};
  REQUIRE(upc.SampleShape() == shape);
  REQUIRE(upc.SampleCount() == count * count);
  for (std::size_t sample = 0; sample < upc.SampleCount(); ++sample) {
    REQUIRE(upc.SampleIndex(sample, 1) == sample / count);
    REQUIRE(upc.SampleIndex(sample, 2) == sample % count);
  }
  for (const auto &node : upc.Nodes({})) { REQUIRE(node.sample_weight.size() == count * count); }
  const std::array<gra::M3Vec, 2> zero_transfer{};
  const auto                      zero_nodes = upc.Nodes(zero_transfer);
  REQUIRE(zero_nodes.size() == 8 * 4);
  const std::array<gra::M3Vec, 2> outside_transfer = {gra::M3Vec{2.5, 0.0, 0.0}, gra::M3Vec{-2.5, 0.0, 0.0}};
  const auto                      outside_nodes    = upc.Nodes(outside_transfer);
  REQUIRE(outside_nodes.size() == 8 * 4);
  for (const auto &node : outside_nodes) {
    REQUIRE(std::isfinite(node.weight.real()));
    REQUIRE(std::isfinite(node.weight.imag()));
  }
  const std::array<gra::M3Vec, 2> transfer = {gra::M3Vec{0.08, -0.05, 0.02}, gra::M3Vec{-0.06, 0.04, -0.015}};
  constexpr double                angle    = 0.37;
  // Rotate one transverse momentum while preserving its longitudinal part
  const auto rotate = [](const gra::M3Vec &q) {
    return gra::M3Vec{std::cos(angle) * q[0] - std::sin(angle) * q[1], std::sin(angle) * q[0] + std::cos(angle) * q[1],
                      q[2]};
  };
  const auto nodes         = upc.Nodes(transfer);
  const auto rotated_nodes = upc.Nodes({rotate(transfer[0]), rotate(transfer[1])});
  REQUIRE(rotated_nodes.size() == nodes.size());
  for (const auto &i : indices(nodes)) {
    REQUIRE(rotated_nodes[i].kt == Approx(nodes[i].kt).margin(1.0e-14));
    REQUIRE(rotated_nodes[i].kx ==
            Approx(std::cos(angle) * nodes[i].kx - std::sin(angle) * nodes[i].ky).margin(2.0e-13));
    REQUIRE(rotated_nodes[i].ky ==
            Approx(std::sin(angle) * nodes[i].kx + std::cos(angle) * nodes[i].ky).margin(2.0e-13));
  }
  RequireScreenMoment(ScreenMoment(upc));
}

TEST_CASE("Lepton-ion configuration currents survive without hadronic profiles",
          "[gra::nuclear][upc][good-walker][config][mixed-beam]") {
  const auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       config = ConfigControls(8);
  config.d_min      = 0.35;
  const auto bank   = SampleBankPtr(*oxygen, config, 271828);

  auto event_param     = UPCControls();
  event_param.survival = gra::nuclear::SurvivalType::Optical;
  const gra::nuclear::MUPC model({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus}, {nullptr, oxygen},
                                 event_param);
  REQUIRE_FALSE(model.HasSamples());
  gra::MRandom event_random;
  event_random.SetSeed(271828);
  REQUIRE(model.Sample(event_random)->HasSamples());

  for (const std::size_t ion_leg : {0U, 1U}) {
    CAPTURE(ion_leg);
    auto param                   = UPCControls();
    param.survival               = gra::nuclear::SurvivalType::Optical;
    param.config                 = config;
    param.emission[ion_leg]      = gra::nuclear::CoherenceType::Inclusive;
    param.target[ion_leg]        = gra::nuclear::CoherenceType::Inclusive;
    param.photo[ion_leg].b_nodes = 32;
    param.photo[ion_leg].z_nodes = 32;

    std::array<gra::nuclear::BeamType, 2> type = {gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Lepton};
    std::array<std::shared_ptr<const gra::nuclear::MNucleus>, 2>    nucleus{};
    std::array<std::shared_ptr<const gra::nuclear::MConfigBank>, 2> banks{};
    type[ion_leg]    = gra::nuclear::BeamType::Nucleus;
    nucleus[ion_leg] = oxygen;
    banks[ion_leg]   = bank;

    const gra::nuclear::MUPC upc(type, nucleus, param, banks);
    REQUIRE_FALSE(upc.HasHadronicPair());
    REQUIRE(upc.Glauber() == nullptr);
    REQUIRE(upc.HasSamples());
    REQUIRE(upc.SampleCount() == config.count);
    REQUIRE(upc.Nodes({}).empty());

    const int        leg = static_cast<int>(ion_leg + 1);
    const gra::M3Vec q   = {0.11, -0.04, 0.07};
    RequireCurrentMoments(upc.EmissionComponents(leg, q), config.count);
    RequireCurrentMoments({upc.TargetRatios(leg, gra::nuclear::CoherenceType::Coherent, q),
                           upc.TargetRatios(leg, gra::nuclear::CoherenceType::Incoherent, q)},
                          config.count);
    const auto  direction = gra::nuclear::TargetPhotonDirection(leg);
    const auto *photo     = upc.Photo(leg);
    REQUIRE(photo != nullptr);
    const auto profile    = gra::test::PhotoProfile(18.0, 0.1);
    const auto transition = photo->Factors(profile, q[0], q[1], q[2], direction, bank.get(), nullptr, true);
    REQUIRE(transition.current.has_value());
    const auto &current     = *transition.current;
    const auto  target      = upc.TargetCurrentRatios(leg, gra::nuclear::CoherenceType::Coherent, current);
    const auto  fluctuation = upc.TargetCurrentRatios(leg, gra::nuclear::CoherenceType::Incoherent, current);
    const auto  moment      = gra::statistics::ComplexMoments(current.sample);
    for (const auto &sample : indices(current.sample)) {
      REQUIRE(std::abs(target[sample] - std::complex<double>(1.0, 0.0)) < 2.0e-13);
      const auto reconstructed = moment.mean * target[sample] + std::sqrt(moment.variance) * fluctuation[sample];
      REQUIRE(std::abs(reconstructed - current.sample[sample]) < 2.0e-13);
    }
    // Common complex rescaling preserves the normalized current phase, including tiny amplitudes
    const auto phase = std::polar(1.0, 0.73);
    for (const auto sector : {gra::nuclear::CoherenceType::Inclusive, gra::nuclear::CoherenceType::Incoherent}) {
      const auto reference = upc.TargetCurrentRatios(leg, sector, current);
      for (const double scale : {1.0e-20, 1.0, 1.0e20}) {
        auto scaled = current;
        for (auto &value : scaled.sample) { value *= scale * phase; }
        const auto ratio = upc.TargetCurrentRatios(leg, sector, scaled);
        for (const auto &i : indices(ratio)) {
          CHECK(std::abs(ratio[i] - phase * reference[i]) < 2.0e-12);
        }
      }
      auto zero = current;
      std::fill(zero.sample.begin(), zero.sample.end(), std::complex<double>{});
      for (const auto &value : upc.TargetCurrentRatios(leg, sector, zero)) {
        CHECK(std::abs(value) < std::numeric_limits<double>::min());
      }
    }
    const int  lepton_leg = static_cast<int>(2 - ion_leg);
    const auto lepton     = upc.EmissionRatios(lepton_leg, gra::nuclear::CoherenceType::Coherent, {0.11, -0.04, 0.07});
    REQUIRE(lepton.size() == config.count);
    for (const auto &ratio : lepton) { REQUIRE(std::abs(ratio - std::complex<double>(1.0, 0.0)) < 2.0e-14); }
  }
}

TEST_CASE("Optical screening retains sampled nuclear current ensembles", "[gra::nuclear][upc][good-walker][optical]") {
  const auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       config = ConfigControls(4);
  config.d_min      = 0.35;
  const auto bank   = SampleBankPtr(*oxygen, config, 271829);

  auto param      = UPCControls();
  param.survival  = gra::nuclear::SurvivalType::Optical;
  param.config    = config;
  param.target[1] = gra::nuclear::CoherenceType::Incoherent;
  param.glauber   = GlauberControls(0.0);
  const gra::nuclear::MUPC upc({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus}, {nullptr, oxygen},
                               param, {nullptr, bank});

  REQUIRE(upc.HasSamples());
  REQUIRE(upc.SampleShape() == std::array<std::size_t, 2>{1, 4});
  REQUIRE(upc.SampleCount() == 4);
}

TEST_CASE("Ion-ion configuration currents use the full Cartesian bank",
          "[gra::nuclear][upc][good-walker][config][cartesian]") {
  const auto oxygen       = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       first_config = ConfigControls(5);
  first_config.d_min      = 0.8;
  auto second_config      = first_config;
  second_config.count     = GENERATE(3U, 5U);
  const auto first_bank   = SampleBankPtr(*oxygen, first_config, 1618033);
  const auto second_bank  = SampleBankPtr(*oxygen, second_config, 1618033 ^ 0x9e3779b9U);

  auto param                       = UPCControls();
  param.survival                   = gra::nuclear::SurvivalType::MCGGCF;
  param.config                     = first_config;
  param.glauber                    = GlauberControls(0.0);
  param.glauber.b_nodes            = 32;
  param.convolution.smooth_b_nodes = 32;
  param.convolution.sample_b_nodes = 32;
  param.glauber.fluctuation.nodes  = 1;
  param.loop.r_min                 = 0.01;
  param.loop.r_max                 = 0.12;
  param.loop.radial_intervals      = 8;
  param.loop.azimuth_nodes         = 4;
  param.photo[0].b_nodes           = 32;
  param.photo[1].b_nodes           = 32;
  const gra::nuclear::MUPC complete({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus},
                                    {oxygen, oxygen}, param, {first_bank, second_bank});

  const std::array<std::size_t, 2> expected_shape = {first_config.count, second_config.count};
  REQUIRE(complete.SampleShape() == expected_shape);
  REQUIRE(complete.SampleCount() == first_config.count * second_config.count);
  std::array<std::vector<std::size_t>, 2> degree = {std::vector<std::size_t>(first_config.count, 0),
                                                    std::vector<std::size_t>(second_config.count, 0)};
  std::vector<bool>                       edge(complete.SampleCount(), false);
  for (std::size_t sample = 0; sample < complete.SampleCount(); ++sample) {
    const std::size_t upper = complete.SampleIndex(sample, 1);
    const std::size_t lower = complete.SampleIndex(sample, 2);
    REQUIRE(upper < first_config.count);
    REQUIRE(lower < second_config.count);
    REQUIRE_FALSE(edge[upper * second_config.count + lower]);
    edge[upper * second_config.count + lower] = true;
    ++degree[0][upper];
    ++degree[1][lower];
  }
  for (const auto &leg : indices(degree)) {
    for (const std::size_t value : degree[leg]) { REQUIRE(value == expected_shape[1 - leg]); }
  }

  const auto *glauber = complete.Glauber();
  REQUIRE(glauber != nullptr);
  std::vector<std::complex<double>> kernel(complete.SampleCount(), 0.0);
  double                            survival = 0.0;
  for (const auto &sample : indices(kernel)) {
    kernel[sample] = glauber->ConfigPairAmp(4.0, 0.0, first_bank.get(), second_bank.get(),
                                            complete.SampleIndex(sample, 1), complete.SampleIndex(sample, 2));
    survival += kernel[sample].real();
  }
  survival /= static_cast<double>(kernel.size());
  REQUIRE(complete.SurvivalAmp(4.0) == Approx(survival).epsilon(2.0e-14));
  const auto sector      = gra::nuclear::NuclearGoodWalkerProject(kernel, expected_shape);
  double     mean_second = 0.0;
  for (const auto &value : kernel) { mean_second += std::norm(value); }
  mean_second /= static_cast<double>(kernel.size());
  for (const auto &upper : sector) {
    for (const double sector : upper) { REQUIRE(sector >= 0.0); }
  }
  REQUIRE(gra::Sum(sector[0]) + gra::Sum(sector[1]) == Approx(mean_second).epsilon(2.0e-14));
  const auto complete_node = complete.Nodes({});
  REQUIRE_FALSE(complete_node.empty());
  for (const auto &node : complete_node) { REQUIRE(node.sample_weight.size() == complete.SampleCount()); }
  bool anisotropic = false;
  for (std::size_t phi = 1; phi < param.loop.azimuth_nodes; ++phi) {
    for (std::size_t sample = 0; sample < complete.SampleCount(); ++sample) {
      anisotropic =
          anisotropic || std::abs(complete_node[phi * gra::math::PolarNodeCount(param.loop)].sample_weight[sample] -
                                  complete_node[0].sample_weight[sample]) > 1.0e-14;
    }
  }
  REQUIRE(anisotropic);
}

TEST_CASE("Nuclear screening projects resolved photon-source sectors", "[gra::nuclear][screen][good-walker]") {
  const auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       config = ConfigControls(8);
  config.d_min      = 0.35;
  const auto bank   = SampleBankPtr(*oxygen, config, 314159);

  auto param             = UPCControls();
  param.survival         = gra::nuclear::SurvivalType::Optical;
  param.config           = config;
  param.emission[1]      = gra::nuclear::CoherenceType::Inclusive;
  param.photo[1].b_nodes = 32;
  const gra::nuclear::MUPC upc({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus}, {nullptr, oxygen},
                               param, {nullptr, bank});

  gra::nuclear::ScreenLayout layout;
  layout.type             = gra::nuclear::ScreenType::Fusion;
  layout.fusion.rows      = 1;
  layout.fusion.sector[0] = {gra::nuclear::CoherenceType::Coherent};
  layout.fusion.sector[1] = {gra::nuclear::CoherenceType::Coherent, gra::nuclear::CoherenceType::Incoherent};
  gra::nuclear::ScreenPoint born;
  born.amplitude = {{2.0, -1.0}, {2.0, -1.0}};
  born.transfer  = {{{0.0, 0.0, 0.0}, {0.11, -0.04, 0.07}}};

  const gra::nuclear::MUPCScreen screen(upc, layout, born);
  const auto                     result = screen.Result();
  REQUIRE(result.amplitude.size() == 2);
  REQUIRE(result.helicity_norm.size() == 2);
  // Nuclear projection returns raw |A_h|^2 before process-specific spin averaging
  for (const auto &component : result.helicity_norm) { REQUIRE(component == Approx(5.0).epsilon(2.0e-12)); }
  REQUIRE(std::abs(result.amplitude[0] - born.amplitude[0]) < 2.0e-12);
  REQUIRE(std::abs(result.amplitude[1]) < 2.0e-14);

  // Require identical Born and loop accumulation for helicity and color projections
  born.flow = {{{1.0, 2.0}, {1.0, 2.0}}};
  for (const bool sampled : {false, true}) {
    gra::nuclear::MUPCScreen accumulated(upc, layout, born);
    gra::nuclear::LoopNode   node;
    node.weight = {-0.3, 0.4};
    if (sampled) { node.sample_weight.assign(upc.SampleCount(), node.weight); }
    accumulated.Add(node, born);
    const auto scaled = accumulated.Result();
    const auto factor = 1.0 + node.weight;
    for (const double norm : scaled.helicity_norm) {
      REQUIRE(norm == Approx(5.0 * std::norm(factor)).epsilon(2.0e-12));
    }
    REQUIRE(scaled.color_norm.size() == 1);
    REQUIRE(scaled.color_norm[0] == Approx(10.0 * std::norm(factor)).epsilon(2.0e-12));
    REQUIRE(std::abs(scaled.amplitude[0] - factor * born.amplitude[0]) < 2.0e-12);
    REQUIRE(std::abs(scaled.amplitude[1]) < 2.0e-14);
  }
}

// Compare the screened charge currents with projection of the complete Cartesian ensemble
TEST_CASE("Nuclear fusion screens complete currents before selecting final sectors",
          "[gra::nuclear][screen][good-walker][screen-fix]") {
  using namespace gra::nuclear;
  const auto                                oxygen = std::make_shared<const MNucleus>(MakeOxygen());
  const auto                                config = ConfigControls(4);
  const auto                                first  = SampleBankPtr(*oxygen, config, 314159);
  const auto                                second = SampleBankPtr(*oxygen, config, 271828);
  const std::array<std::complex<double>, 2> hard   = {{{1.0, 0.3}, {-0.4, 0.7}}};
  const std::array<std::complex<double>, 2> pol1   = {{{0.6, -0.2}, {0.1, 0.9}}};
  const std::array<std::complex<double>, 2> pol2   = {{{0.3, 0.5}, {-0.8, 0.2}}};
  const double               spin_norm = gra::SquaredNorm(hard) * gra::SquaredNorm(pol1) * gra::SquaredNorm(pol2);
  const std::complex<double> color{0.4, -0.7};

  for (const auto type :
       {std::array{BeamType::Proton, BeamType::Nucleus}, std::array{BeamType::Nucleus, BeamType::Proton},
        std::array{BeamType::Nucleus, BeamType::Nucleus}}) {
    const std::array nucleus = {type[0] == BeamType::Nucleus ? oxygen : nullptr,
                                type[1] == BeamType::Nucleus ? oxygen : nullptr};
    const std::array bank = {nucleus[0] != nullptr ? first : nullptr, nucleus[1] != nullptr ? second : nullptr};
    for (const auto selected : {CoherenceType::Inclusive, CoherenceType::Coherent, CoherenceType::Incoherent}) {
      auto param     = UPCControls();
      param.config   = config;
      param.survival = SurvivalType::MCGGCF;
      gra::test::Glauber(param.glauber, 0.1);
      for (const auto &leg : indices(type)) {
        param.emission[leg] = nucleus[leg] != nullptr ? selected : CoherenceType::Coherent;
      }
      const MUPC   upc(type, nucleus, param, bank);
      ScreenLayout layout;
      layout.type        = ScreenType::Fusion;
      layout.fusion.rows = pol1.size();
      for (const auto &leg : indices(type)) {
        layout.fusion.sector[leg] = CoherenceSectors(upc.SourceSector(leg + 1, param.emission[leg]));
      }

      // Resolve each measured current into its mean and fluctuation coefficients
      const auto point = [&](const std::array<gra::M3Vec, 2> &q) {
        std::array<std::vector<std::complex<double>>, 2> source;
        for (const auto &leg : indices(source)) {
          const auto stat = gra::statistics::ComplexMoments(ChargeBank(bank[leg].get(), q[leg]));
          source[leg]     = {stat.mean};
          if (bank[leg] != nullptr) { source[leg].push_back(std::sqrt(stat.variance)); }
        }
        ScreenPoint output;
        output.transfer  = q;
        output.amplitude = gra::KroneckerProduct(hard, gra::KroneckerProduct(gra::KroneckerProduct(source[0], pol1),
                                                                             gra::KroneckerProduct(source[1], pol2)));
        output.flow      = {output.amplitude};
        gra::Scale(output.flow[0], color);
        return output;
      };

      for (const bool zero_transfer : {false, true}) {
        std::array<gra::M3Vec, 2> q = {{{0.07, 0.03, -0.02}, {0.11, -0.04, 0.07}}};
        if (zero_transfer) { q = {}; }
        const auto born  = point(q);
        const auto nodes = upc.Nodes(q);
        for (const bool sampled : {false, true}) {
          CAPTURE(type, selected, zero_transfer, sampled);
          MUPCScreen screen(upc, layout, born);
          auto       direct = gra::KroneckerProduct(ChargeBank(bank[0].get(), q[0]), ChargeBank(bank[1].get(), q[1]));
          for (auto node : nodes) {
            if (!sampled) { node.sample_weight.clear(); }
            auto shifted = q;
            shifted[0][0] += node.kx;
            shifted[0][1] += node.ky;
            shifted[1][0] -= node.kx;
            shifted[1][1] -= node.ky;
            screen.Add(node, point(shifted));
            const auto current =
                gra::KroneckerProduct(ChargeBank(bank[0].get(), shifted[0]), ChargeBank(bank[1].get(), shifted[1]));
            for (const auto &i : indices(direct)) {
              direct[i] += (sampled ? node.sample_weight[i] : node.weight) * current[i];
            }
          }
          const auto result   = screen.Result();
          const auto weight   = FinalSectorWeights(layout, result.helicity_norm);
          const auto expected = NuclearGoodWalkerProject(direct, upc.SampleShape());
          for (const auto &pair : indices(weight)) {
            const auto upper    = pair / 2 == 0 ? CoherenceType::Coherent : CoherenceType::Incoherent;
            const auto lower    = pair % 2 == 0 ? CoherenceType::Coherent : CoherenceType::Incoherent;
            const bool accepted = (param.emission[0] == CoherenceType::Inclusive || param.emission[0] == upper) &&
                                  (param.emission[1] == CoherenceType::Inclusive || param.emission[1] == lower);
            const double reference = accepted ? spin_norm * expected[pair / 2][pair % 2] : 0.0;
            CHECK(weight[pair] == Approx(reference).epsilon(2.0e-11).margin(1.0e-24));
            CHECK(result.color_sector[0][pair] ==
                  Approx(std::norm(color) * reference).epsilon(2.0e-11).margin(1.0e-24));
          }
          CHECK(result.color_norm[0] == Approx(std::norm(color) * gra::Sum(weight)).epsilon(2.0e-11));
        }
      }
    }
  }
}

// Compare photonuclear screening with the directly shadowed configuration currents
TEST_CASE("Photonuclear projection follows convolution of both complete directions",
          "[gra::nuclear][photo][screen][good-walker][screen-fix]") {
  using namespace gra::nuclear;
  const auto       oxygen = std::make_shared<const MNucleus>(MakeOxygen());
  const auto       config = ConfigControls(4);
  const std::array bank   = {SampleBankPtr(*oxygen, config, 3123), SampleBankPtr(*oxygen, config, 5782)};
  const std::array<std::complex<double>, 2> hard = {{{1.0, 0.3}, {-0.6, 0.2}}};
  const std::array profile = {gra::test::PhotoProfile(5.0, 0.04), gra::test::PhotoProfile(25.0, 0.12)};
  for (const auto selected : {CoherenceType::Inclusive, CoherenceType::Coherent, CoherenceType::Incoherent}) {
    auto param     = UPCControls();
    param.config   = config;
    param.survival = SurvivalType::MCGGCF;
    param.emission = {CoherenceType::Coherent, CoherenceType::Inclusive};
    param.target   = {selected, selected};
    gra::test::Glauber(param.glauber, 0.1);
    const MUPC   upc({BeamType::Nucleus, BeamType::Nucleus}, {oxygen, oxygen}, param, bank);
    ScreenLayout layout;
    layout.type  = ScreenType::Photo;
    layout.photo = PhotoChannels(upc);

    // Evaluate both mean and fluctuation target factors with their actual shadowed currents
    const auto point = [&](const std::array<gra::M3Vec, 2> &q) {
      ScreenPoint output;
      output.transfer = q;
      std::array<std::array<std::complex<double>, 2>, 2> emission, target;
      for (const auto &leg : indices(bank)) {
        const auto stat           = gra::statistics::ComplexMoments(ChargeBank(bank[leg].get(), q[leg]));
        emission[leg]             = {stat.mean, std::sqrt(stat.variance)};
        auto value                = upc.Photo(leg + 1)->Factors(profile[leg], q[leg][0], q[leg][1], q[leg][2],
                                                                TargetPhotonDirection(leg + 1), bank[leg].get(), nullptr, true);
        target[leg]               = {value.coherent, value.incoherent};
        output.photo_current[leg] = std::move(value.current);
      }
      for (const auto &channel : layout.photo) {
        const std::size_t emitter = channel.direction == PhotoDirection::Upper ? 0 : 1;
        const std::size_t a       = channel.emission == CoherenceType::Coherent ? 0 : 1;
        const std::size_t b       = channel.target == CoherenceType::Coherent ? 0 : 1;
        output.amplitude.push_back(hard[emitter] * emission[emitter][a] * target[1 - emitter][b]);
      }
      return output;
    };
    const std::array<gra::M3Vec, 2> q    = {{{0.07, 0.03, -0.02}, {0.11, -0.04, 0.07}}};
    const auto                      born = point(q);
    for (const bool sampled : {false, true}) {
      CAPTURE(selected, sampled);
      MUPCScreen                                       screen(upc, layout, born);
      std::array<std::vector<std::complex<double>>, 2> direct;
      for (auto &amplitude : direct) { amplitude.assign(upc.SampleCount(), 0.0); }
      // Sum the complete currents without using the screening source ratios
      const auto accumulate = [&](const ScreenPoint &value, const LoopNode &node) {
        for (const auto &emitter : indices(direct)) {
          const auto  charge  = ChargeBank(bank[emitter].get(), value.transfer[emitter]);
          const auto &current = *value.photo_current[1 - emitter];
          auto        target  = current.sample;
          // Preserve the unbiased target variance used by MPhoto::Factors
          const double correction = std::sqrt(static_cast<double>(target.size()) / (target.size() - 1));
          for (auto &sample : target) { sample = current.stat.mean + correction * (sample - current.stat.mean); }
          for (const auto &sample : indices(direct[emitter])) {
            const auto weight = node.sample_weight.empty() ? node.weight : node.sample_weight[sample];
            direct[emitter][sample] += weight * hard[emitter] * charge[upc.SampleIndex(sample, emitter + 1)] *
                                       target[upc.SampleIndex(sample, 2 - emitter)];
          }
        }
      };
      LoopNode initial;
      initial.weight = std::complex<double>(1.0, 0.0);
      accumulate(born, initial);
      for (auto node : upc.Nodes(q)) {
        if (!sampled) { node.sample_weight.clear(); }
        auto shifted = q;
        shifted[0][0] += node.kx;
        shifted[0][1] += node.ky;
        shifted[1][0] -= node.kx;
        shifted[1][1] -= node.ky;
        const auto value = point(shifted);
        screen.Add(node, value);
        accumulate(value, node);
      }
      const auto result            = screen.Result();
      const auto weight            = FinalSectorWeights(layout, result.helicity_norm);
      const auto selected_channels = PhotoChannels(param);
      for (const auto &pair : indices(weight)) {
        std::vector<std::complex<double>> complete(upc.SampleCount(), 0.0);
        for (const auto &channel : selected_channels) {
          if (channel.Pair() == pair) {
            gra::AddScaled(complete, direct[channel.direction == PhotoDirection::Upper ? 0 : 1], 1.0);
          }
        }
        const auto expected = NuclearGoodWalkerProject(complete, upc.SampleShape());
        CHECK(weight[pair] == Approx(expected[pair / 2][pair % 2]).epsilon(2.0e-11).margin(1.0e-24));
      }
    }
  }
}

TEST_CASE("Cartesian Good-Walker projection separates both nuclear legs",
          "[gra::nuclear][screen][good-walker][cartesian]") {
  using gra::nuclear::NuclearGoodWalkerProject;
  const std::array<std::size_t, 2>        shape          = {2, 2};
  const std::vector<std::complex<double>> amplitude      = {{4.0, 7.0}, {-2.0, -1.0}, {4.0, -5.0}, {-2.0, 3.0}};
  const std::vector<std::complex<double>> upper_permuted = {{4.0, -5.0}, {-2.0, 3.0}, {4.0, 7.0}, {-2.0, -1.0}};
  const std::vector<std::complex<double>> lower_permuted = {{-2.0, -1.0}, {4.0, 7.0}, {-2.0, 3.0}, {4.0, -5.0}};

  // Require all four orthogonal sectors and closure for one leg ordering
  const auto reference = NuclearGoodWalkerProject(amplitude, shape);
  REQUIRE(reference[0][0] == Approx(2.0).epsilon(2.0e-14));
  REQUIRE(reference[0][1] == Approx(9.0).epsilon(2.0e-14));
  REQUIRE(reference[1][0] == Approx(4.0).epsilon(2.0e-14));
  REQUIRE(reference[1][1] == Approx(16.0).epsilon(2.0e-14));
  REQUIRE(gra::Sum(reference[0]) + gra::Sum(reference[1]) == Approx(31.0).epsilon(2.0e-14));

  // Require invariance under independent upper and lower bank permutations
  for (const auto &permuted : {upper_permuted, lower_permuted}) {
    const auto result = NuclearGoodWalkerProject(permuted, shape);
    for (const auto &upper : indices(reference)) {
      for (const auto &lower : indices(reference[upper])) {
        REQUIRE(result[upper][lower] == Approx(reference[upper][lower]).epsilon(2.0e-14));
      }
    }
  }

  // Require the exact factorized finite bank moments
  const std::array<std::complex<double>, 3> x = {std::complex<double>{1.0, 0.2}, {2.0, -0.4}, {-0.5, 0.7}};
  const std::array<std::complex<double>, 4> y = {std::complex<double>{0.8, -0.1}, {-1.1, 0.3}, {0.4, 0.6}, {1.3, -0.2}};
  std::vector<std::complex<double>>         factorized;
  std::complex<double>                      mean_x   = 0.0;
  std::complex<double>                      mean_y   = 0.0;
  double                                    second_x = 0.0;
  double                                    second_y = 0.0;
  for (const auto &value : x) {
    mean_x += value;
    second_x += std::norm(value);
  }
  for (const auto &value : y) {
    mean_y += value;
    second_y += std::norm(value);
  }
  mean_x /= static_cast<double>(x.size());
  mean_y /= static_cast<double>(y.size());
  second_x /= static_cast<double>(x.size());
  second_y /= static_cast<double>(y.size());
  for (const auto &upper : x) {
    for (const auto &lower : y) { factorized.push_back(upper * lower); }
  }
  const double coherent_x        = std::norm(mean_x);
  const double coherent_y        = std::norm(mean_y);
  const double variance_x        = second_x - coherent_x;
  const double variance_y        = second_y - coherent_y;
  const auto   factorized_result = NuclearGoodWalkerProject(factorized, {x.size(), y.size()});
  REQUIRE(factorized_result[0][0] == Approx(coherent_x * coherent_y).epsilon(2.0e-14));
  REQUIRE(factorized_result[0][1] == Approx(coherent_x * variance_y).epsilon(2.0e-14));
  REQUIRE(factorized_result[1][0] == Approx(variance_x * coherent_y).epsilon(2.0e-14));
  REQUIRE(factorized_result[1][1] == Approx(variance_x * variance_y).epsilon(2.0e-14));
}

TEST_CASE("Photonuclear directions interfere only within one final sector", "[gra::nuclear][screen][photo]") {
  auto param         = UPCControls();
  param.emission[0]  = gra::nuclear::CoherenceType::Inclusive;
  param.target[1]    = gra::nuclear::CoherenceType::Inclusive;
  const auto channel = gra::nuclear::PhotoChannels(param);
  REQUIRE(channel.size() == 5);

  const std::vector<std::complex<double>> amplitude = {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}};
  const auto                              combined  = gra::nuclear::CombinePhotoChannels(amplitude, channel);
  REQUIRE(combined.size() == amplitude.size());
  REQUIRE(gra::SquaredNorm(combined) == Approx(65.0));
  REQUIRE(std::count_if(combined.begin(), combined.end(), [](const auto &item) { return std::abs(item) > 0.0; }) == 4);
}

TEST_CASE("UPC sampling requires a resolvable Good-Walker variance", "[gra::nuclear][upc][good-walker]") {
  const auto oxygen   = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param    = UPCControls();
  param.config.count  = 1;
  param.current_count = 0;
  param.emission      = {gra::nuclear::CoherenceType::Incoherent, gra::nuclear::CoherenceType::Coherent};
  gra::nuclear::UPCMode mode;
  mode.target    = {false, false};
  mode.screening = false;
  REQUIRE_THROWS_WITH(gra::nuclear::MUPC({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus},
                                         {oxygen, oxygen}, param, {}, mode),
                      Catch::Matchers::Contains("Good-Walker sampling requires at least two nuclear configurations"));
}

TEST_CASE("Focused UPC screening accepts only its implemented quadrature rules", "[gra::nuclear][upc][quadrature]") {
  const auto oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param  = UPCControls();
  param.loop.r_min  = 1.0e-6;

  REQUIRE_NOTHROW(
      gra::nuclear::MUPC({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus}, {oxygen, oxygen}, param));

  for (const std::string radial : {"1/3", "3/8", "Boole"}) {
    DYNAMIC_SECTION("radial rule " << radial) {
      auto unsupported                   = param;
      unsupported.loop.radial_integrator = radial;
      REQUIRE_THROWS_AS(gra::nuclear::MUPC({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus},
                                           {oxygen, oxygen}, unsupported),
                        std::invalid_argument);
    }
  }

  param.loop.azimuth_integrator = "1/3";
  REQUIRE_THROWS_AS(
      gra::nuclear::MUPC({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus}, {oxygen, oxygen}, param),
      std::invalid_argument);
}

TEST_CASE("Nuclear final-sector weights preserve orthogonal C and I pairs", "[gra::nuclear][screen][record]") {
  auto param     = UPCControls();
  param.emission = {gra::nuclear::CoherenceType::Inclusive, gra::nuclear::CoherenceType::Inclusive};
  param.target   = {gra::nuclear::CoherenceType::Inclusive, gra::nuclear::CoherenceType::Inclusive};

  gra::nuclear::ScreenLayout photo;
  photo.type  = gra::nuclear::ScreenType::Photo;
  photo.photo = gra::nuclear::PhotoChannels(param);
  std::vector<double>        component(photo.photo.size(), 0.0);
  gra::nuclear::FinalWeights expected{};
  for (const auto &i : indices(component)) {
    component[i]        = static_cast<double>(i + 1);
    const auto &channel = photo.photo[i];
    expected[channel.Pair()] += component[i];
  }
  CHECK(gra::nuclear::FinalSectorWeights(photo, component) == expected);

  gra::nuclear::ScreenLayout fusion;
  fusion.type          = gra::nuclear::ScreenType::Fusion;
  fusion.fusion.rows   = 2;
  fusion.fusion.sector = {std::vector<gra::nuclear::CoherenceType>{gra::nuclear::CoherenceType::Coherent,
                                                                   gra::nuclear::CoherenceType::Incoherent},
                          std::vector<gra::nuclear::CoherenceType>{gra::nuclear::CoherenceType::Coherent,
                                                                   gra::nuclear::CoherenceType::Incoherent}};
  const std::vector<double> fusion_component(16, 1.0);
  const auto                fusion_weight = gra::nuclear::FinalSectorWeights(fusion, fusion_component);
  for (const double weight : fusion_weight) { CHECK(weight == Approx(4.0)); }

  const gra::nuclear::FinalWeights weight = {1.0, 2.0, 3.0, 4.0};
  const auto                       cc     = gra::nuclear::SampleFinalState(weight, 0.00);
  const auto                       ci     = gra::nuclear::SampleFinalState(weight, 0.15);
  const auto                       ic     = gra::nuclear::SampleFinalState(weight, 0.45);
  const auto                       ii     = gra::nuclear::SampleFinalState(weight, 0.85);
  REQUIRE(cc.valid);
  REQUIRE(ci.valid);
  REQUIRE(ic.valid);
  REQUIRE(ii.valid);
  CHECK((cc.leg == std::array<gra::nuclear::CoherenceType, 2>{gra::nuclear::CoherenceType::Coherent,
                                                              gra::nuclear::CoherenceType::Coherent}));
  CHECK((ci.leg == std::array<gra::nuclear::CoherenceType, 2>{gra::nuclear::CoherenceType::Coherent,
                                                              gra::nuclear::CoherenceType::Incoherent}));
  CHECK((ic.leg == std::array<gra::nuclear::CoherenceType, 2>{gra::nuclear::CoherenceType::Incoherent,
                                                              gra::nuclear::CoherenceType::Coherent}));
  CHECK((ii.leg == std::array<gra::nuclear::CoherenceType, 2>{gra::nuclear::CoherenceType::Incoherent,
                                                              gra::nuclear::CoherenceType::Incoherent}));
}

TEST_CASE("Configuration photon densities use the full rest transfer", "[gra::nuclear][photon][config]") {
  using gra::nuclear::CoherenceType;
  const auto                  nucleus = MakeOxygen();
  const auto                  bank    = SampleBank(nucleus, ConfigControls(16), 12345);
  const gra::nuclear::MPhoton photon(nucleus);
  const double                xi = 1.0e-4;
  const double                t  = -0.04;
  const double                pt = 0.2;
  const gra::M3Vec            reference{0.2, 0.0, 0.0};
  const auto                  base     = photon.Density(CoherenceType::Incoherent, xi, t, pt, reference, &bank);
  const double                variance = bank.ChargeStat(reference[0], reference[1], reference[2]).variance;
  REQUIRE(base.Trace() > 0.0);
  for (const auto &q : std::array<gra::M3Vec, 3>{{{0.0, 0.2, 0.0}, {0.0, 0.0, 0.2}, {0.12, -0.16, 0.0}}}) {
    const auto   density = photon.Density(CoherenceType::Incoherent, xi, t, pt, q, &bank);
    const double ratio   = bank.ChargeStat(q[0], q[1], q[2]).variance / variance;
    REQUIRE(density.parallel == Approx(base.parallel * ratio).epsilon(2.0e-13));
    REQUIRE(density.perpendicular == Approx(base.perpendicular * ratio).epsilon(2.0e-13));
    const auto reversed = photon.Density(CoherenceType::Incoherent, xi, t, pt, {-q[0], -q[1], -q[2]}, &bank);
    REQUIRE(reversed.Trace() == Approx(density.Trace()).epsilon(2.0e-13));
    const auto inclusive = photon.Density(CoherenceType::Inclusive, xi, t, pt, q, &bank);
    const auto coherent  = photon.Density(CoherenceType::Coherent, xi, t, pt);
    REQUIRE(inclusive.Trace() == Approx(coherent.Trace() + density.Trace()).epsilon(2.0e-13));
  }
}

TEST_CASE("Coherent photoproduction has no finite-bank variance floor", "[gra::nuclear][photo][config]") {
  const auto                 nucleus = MakeOxygen();
  const gra::nuclear::MPhoto photo(nucleus, PhotoControls(), gra::nuclear::PhotoModel::Impulse);
  const auto                 profile = gra::test::PhotoProfile(0.0, 0.0, 0.0);
  REQUIRE_THROWS_AS(photo.Factors(profile, 0.2, 0.0, 0.0, static_cast<gra::nuclear::PhotonDirection>(255)),
                    std::invalid_argument);
  for (const std::size_t count : {4U, 16U}) {
    for (const unsigned int seed : {42U, 314159U, 12345U}) {
      const auto bank = SampleBank(nucleus, ConfigControls(count), seed);
      for (const double q : {0.0, 0.2, 0.8}) {
        const auto direction = gra::nuclear::PhotonDirection::PositiveZ;
        const auto smooth    = photo.Factors(profile, q, 0.0, 0.0, direction);
        const auto factor    = photo.Factors(profile, q, 0.0, 0.0, direction, &bank, nullptr, true);
        REQUIRE(std::abs(factor.coherent - smooth.coherent) < 2.0e-13);
        REQUIRE(factor.current.has_value());
        const auto stat = gra::statistics::ComplexMoments(factor.current->sample);
        REQUIRE(std::abs(stat.mean - smooth.coherent) < 2.0e-13);
        const auto raw = photo.ShadowStat(profile, bank, q, 0.0, 0.0, direction);
        REQUIRE(factor.incoherent * factor.incoherent == Approx(raw.variance).margin(2.0e-12));
        const double unbiased = gra::statistics::UnbiasedComplexVariance(stat.second, stat.mean, count);
        REQUIRE(unbiased == Approx(raw.variance).margin(2.0e-12));
      }
    }
  }
}

TEST_CASE("UPC cache follows effective convolution controls", "[gra::nuclear][cache]") {
  const auto nucleus      = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param        = UPCControls();
  param.structure         = gra::nuclear::StructureType::Smooth;
  param.config.count      = 0;
  param.survival          = gra::nuclear::SurvivalType::Optical;
  param.table_fingerprint = "nuclear-effective-convolution-regression";
  gra::nuclear::UPCMode mode;
  mode.target                                      = {false, false};
  const std::array<gra::nuclear::BeamType, 2> type = {gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus};
  const std::array<std::shared_ptr<const gra::nuclear::MNucleus>, 2> beam = {nullptr, nucleus};
  const gra::nuclear::MUPC                                           original(type, beam, param, {}, mode);
  param.loop.r_max *= 0.5;
  param.glauber.profile.sigma *= 1.1;
  for (auto &b : param.glauber.profile.b_node) { b *= std::sqrt(1.1); }
  const gra::nuclear::MUPC cached(type, beam, param, {}, mode);
  param.table_fingerprint.clear();
  const gra::nuclear::MUPC fresh(type, beam, param, {}, mode);
  REQUIRE(cached.Nodes({}).size() == fresh.Nodes({}).size());
  REQUIRE(cached.Nodes({}).back().kt < original.Nodes({}).back().kt);
  for (const auto &i : indices(cached.Nodes({}))) {
    REQUIRE(cached.Nodes({})[i].kt == Approx(fresh.Nodes({})[i].kt).epsilon(2.0e-13));
    REQUIRE(std::abs(cached.Nodes({})[i].weight - fresh.Nodes({})[i].weight) < 2.0e-12);
  }
  for (const double q : {0.0, 0.03, 0.08}) {
    REQUIRE(cached.ProfileTransform(q) == Approx(fresh.ProfileTransform(q)).margin(2.0e-10));
  }
}

TEST_CASE("Nuclear screening rejects malformed nodes before accumulation", "[gra::nuclear][screen][validation]") {
  const auto            nucleus = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto                  param   = UPCControls();
  gra::nuclear::UPCMode mode;
  mode.screening = false;
  mode.target    = {false, false};
  const gra::nuclear::MUPC  upc({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus}, {nullptr, nucleus},
                                param, {}, mode);
  gra::nuclear::ScreenPoint born;
  born.amplitude = {{1.0, 0.2}, {-0.3, 0.5}};
  born.flow      = {born.amplitude};
  gra::nuclear::MUPCScreen screen(upc, {}, born);
  const auto               before = screen.Result();
  gra::nuclear::LoopNode   node;
  node.weight = 0.1;
  for (const unsigned int failure : {0U, 1U, 2U, 3U, 4U, 5U, 6U}) {
    auto point   = born;
    auto invalid = node;
    if (failure == 0) { point.amplitude.clear(); }
    if (failure == 1) { point.amplitude.push_back(0.0); }
    if (failure == 2) { point.flow[0].pop_back(); }
    if (failure == 3) { point.flow.clear(); }
    if (failure == 4) { point.transfer[1][0] = std::numeric_limits<double>::quiet_NaN(); }
    if (failure == 5) { invalid.weight = std::numeric_limits<double>::infinity(); }
    if (failure == 6) { invalid.sample_weight.assign(upc.SampleCount() + 1, 1.0); }
    REQUIRE_THROWS_AS(screen.Add(invalid, point), gra::AmplitudeFailure);
    const auto after = screen.Result();
    for (const auto &h : indices(before.amplitude)) {
      REQUIRE(std::abs(after.amplitude[h] - before.amplitude[h]) < 1.0e-15);
      REQUIRE(after.helicity_norm[h] == Approx(before.helicity_norm[h]).margin(1.0e-15));
    }
  }
  screen.Add(node, born);
  const auto valid = screen.Result();
  REQUIRE(gra::AllFinite(valid.helicity_norm));
}

TEST_CASE("Bare photonuclear interference retains charge-matter covariance", "[gra::nuclear][photo][screen]") {
  using gra::nuclear::CoherenceType;
  const auto nucleus           = MakeOxygen();
  auto       param             = UPCControls();
  param.photo_model            = gra::nuclear::PhotoModel::Impulse;
  param.emission               = {CoherenceType::Inclusive, CoherenceType::Inclusive};
  param.target                 = param.emission;
  const auto            first  = SampleBankPtr(nucleus, param.config, 12345);
  const auto            second = SampleBankPtr(nucleus, param.config, 67890);
  gra::nuclear::UPCMode mode;
  mode.screening = false;
  const gra::nuclear::MUPC   upc(nucleus, nucleus, param, first, second, mode);
  gra::nuclear::ScreenLayout layout;
  layout.type  = gra::nuclear::ScreenType::Photo;
  layout.photo = gra::nuclear::PhotoChannels(upc);
  gra::nuclear::ScreenPoint born;
  born.transfer = {gra::M3Vec{0.2, 0.1, 0.03}, gra::M3Vec{-0.1, 0.15, -0.02}};
  born.amplitude.assign(layout.photo.size(), 0.0);
  for (const auto &i : indices(layout.photo)) {
    if (layout.photo[i].Pair() == gra::nuclear::FinalPairIndex(CoherenceType::Incoherent, CoherenceType::Coherent)) {
      born.amplitude[i] = 1.0;
    }
  }
  std::vector<std::complex<double>> charge(first->Size());
  std::vector<std::complex<double>> matter(first->Size());
  const auto                       &q = born.transfer[0];
  for (const auto &i : indices(charge)) {
    charge[i] = first->At(i).ChargeCurrent(q[0], q[1], q[2]);
    matter[i] = first->At(i).MatterCurrent(q[0], q[1], q[2]);
  }
  const auto           c          = gra::statistics::ComplexMoments(charge);
  const auto           m          = gra::statistics::ComplexMoments(matter);
  std::complex<double> covariance = 0.0;
  for (const auto &i : indices(charge)) { covariance += (charge[i] - c.mean) * std::conj(matter[i] - m.mean); }
  covariance /= static_cast<double>(charge.size());
  const double                   expected = 2.0 + 2.0 * covariance.real() / std::sqrt(c.variance * m.variance);
  const gra::nuclear::MUPCScreen screen(upc, layout, born);
  REQUIRE(gra::Sum(screen.Result().helicity_norm) == Approx(expected).margin(2.0e-12));
  REQUIRE(std::abs(expected - 4.0) > 1.0e-4);
  REQUIRE_THROWS_AS(gra::nuclear::CombinePhotoChannels(born.amplitude, layout.photo), gra::AmplitudeFailure);
  auto smooth         = param;
  smooth.structure    = gra::nuclear::StructureType::Smooth;
  smooth.config.count = 0;
  REQUIRE_THROWS_AS(gra::nuclear::MUPC(nucleus, nucleus, smooth, nullptr, nullptr, mode), std::invalid_argument);
  for (auto &amplitude : born.amplitude) { amplitude *= std::polar(1.0, 0.73); }
  for (auto &transfer : born.transfer) {
    for (auto &component : transfer) { component = -component; }
  }
  const gra::nuclear::MUPCScreen reversed(upc, layout, born);
  REQUIRE(gra::Sum(reversed.Result().helicity_norm) == Approx(expected).margin(2.0e-12));
}

TEST_CASE("Scalar photonuclear screening accepts optical direction terms", "[gra::nuclear][photo][screen]") {
  const auto nucleus = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param   = UPCControls();
  param.structure    = gra::nuclear::StructureType::Smooth;
  param.config.count = 0;
  param.target[1]    = gra::nuclear::CoherenceType::Inclusive;
  gra::nuclear::UPCMode mode;
  mode.screening = false;
  const gra::nuclear::MUPC   upc({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus}, {nullptr, nucleus},
                                 param, {}, mode);
  gra::nuclear::ScreenLayout layout;
  layout.photo = gra::nuclear::PhotoChannels(upc);
  gra::nuclear::ScreenPoint born;
  born.amplitude.assign(layout.photo.size(), {0.7, -0.2});
  gra::nuclear::PhotoTerm term;
  term.direction   = gra::nuclear::PhotoDirection::Upper;
  term.amplitude   = born.amplitude;
  born.photo_terms = {term};
  gra::nuclear::MUPCScreen screen(upc, layout, born);
  gra::nuclear::LoopNode   node;
  node.weight  = 0.25;
  auto shifted = born;
  gra::Scale(shifted.amplitude, std::complex<double>(0.4, 0.3));
  shifted.photo_terms[0].amplitude = shifted.amplitude;
  screen.Add(node, shifted);
  const auto result = screen.Result();
  for (const auto &h : indices(born.amplitude)) {
    const auto expected = born.amplitude[h] + node.weight * shifted.amplitude[h];
    REQUIRE(std::abs(result.amplitude[h] - expected) < 2.0e-13);
    REQUIRE(result.helicity_norm[h] == Approx(std::norm(expected)).margin(2.0e-13));
  }
  shifted.photo_terms[0].amplitude[0] += 0.2;
  REQUIRE_THROWS_AS(screen.Add(node, shifted), gra::AmplitudeFailure);
}

TEST_CASE("Identical beam nuclei retain independent target response controls", "[gra::nuclear][photo][cache]") {
  const auto nucleus      = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param        = UPCControls();
  param.structure         = gra::nuclear::StructureType::Smooth;
  param.config.count               = 0;
  param.table_fingerprint          = "independent-photonuclear-leg-regression";
  param.photo[0].table.qt_nodes = 4;
  param.photo[1].table.qt_nodes = 65;
  gra::nuclear::UPCMode mode;
  mode.screening = false;
  const gra::nuclear::MUPC   upc({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus}, {nucleus, nucleus},
                                 param, {}, mode);
  const auto                 profile   = gra::test::PhotoProfile(25.0, 0.1);
  const auto                 direction = gra::nuclear::PhotonDirection::PositiveZ;
  const gra::nuclear::MPhoto direct(nucleus, param.photo[1]);
  const auto                 expected = direct.Factors(profile, 0.06, 0.02, 0.01, direction);
  const auto                 upper    = upc.Photo(1)->Factors(profile, 0.06, 0.02, 0.01, direction);
  const auto                 lower    = upc.Photo(2)->Factors(profile, 0.06, 0.02, 0.01, direction);
  REQUIRE(std::abs(lower.coherent - expected.coherent) < 2.0e-12);
  REQUIRE(lower.incoherent == Approx(expected.incoherent).margin(2.0e-12));
  REQUIRE(std::abs(upper.coherent - lower.coherent) > 1.0e-4);
}

TEST_CASE("Photon direction sums reject non-finite amplitudes", "[gra::nuclear][photo][validation]") {
  const auto channels = gra::nuclear::PhotoChannels(UPCControls());
  REQUIRE(channels.size() == 2);
  const double maximum = std::numeric_limits<double>::max();
  REQUIRE_THROWS_AS(gra::nuclear::CombinePhotoChannels({maximum, maximum}, channels), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(gra::nuclear::CombinePhotoChannels({std::numeric_limits<double>::quiet_NaN(), 0.0}, channels),
                    gra::AmplitudeFailure);
}

// Check orthogonal hard neutron sums without additional photon absorption
TEST_CASE("Coherent photon density follows the hard interference", "[gra::nuclear][photo][screen]") {
  const auto nucleus   = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto       param     = UPCControls();
  param.reaction       = gra::nuclear::ReactionType::Internal;
  param.additional_emd = true;
  param.structure      = gra::nuclear::StructureType::Smooth;
  param.config.count   = 0;
  gra::nuclear::UPCMode mode;
  mode.screening = false;
  const gra::nuclear::MUPC   upc({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus}, {nucleus, nucleus},
                                 param, {}, mode);
  gra::nuclear::ScreenLayout layout;
  layout.photo = gra::nuclear::PhotoChannels(upc);
  REQUIRE(layout.photo.size() == 2);
  for (const double phase : {0.0, 0.7, 2.4}) {
    gra::nuclear::ScreenPoint born;
    born.amplitude    = {{0.7, -0.2}, std::polar(0.4, phase)};
    const auto result = gra::nuclear::MUPCScreen(upc, layout, born).Result();
    CHECK(result.photo[0] == Approx(std::norm(born.amplitude[0])).margin(2e-13));
    CHECK(result.photo[1] == Approx(std::norm(born.amplitude[1])).margin(2e-13));
    CHECK(result.photo[2] == Approx(2.0 * std::real(std::conj(born.amplitude[0]) * born.amplitude[1])).margin(2e-13));
    CHECK(gra::Sum(result.photo) == Approx(gra::Sum(result.helicity_norm)).margin(2e-13));
    std::swap(born.amplitude[0], born.amplitude[1]);
    gra::Scale(born.amplitude, std::polar(1.0, 0.9));
    const auto reversed = gra::nuclear::MUPCScreen(upc, layout, born).Result();
    REQUIRE(gra::Sum(result.helicity_norm) == Approx(std::norm(born.amplitude[0] + born.amplitude[1])).margin(2e-13));
    CHECK(reversed.photo[0] == Approx(result.photo[1]).margin(2e-13));
    CHECK(reversed.photo[1] == Approx(result.photo[0]).margin(2e-13));
    CHECK(reversed.photo[2] == Approx(result.photo[2]).margin(2e-13));
    REQUIRE(gra::Sum(reversed.helicity_norm) == Approx(gra::Sum(result.helicity_norm)).margin(2e-13));
  }
}

// Check optical attenuation against its forward thickness integral without a depth table
TEST_CASE("Optical photonuclear current obeys the forward thickness identity", "[gra::nuclear][photo]") {
  const auto   nucleus = MakeOxygen();
  const auto   profile = gra::test::PhotoProfile(20.0, 0.0, 0.0);
  const double a = nucleus.A(), opacity = 0.05 * profile.sigma_eff * a;
  const auto   rule  = gra::math::GaussLegendreRule(256, 0.0, nucleus.MatterDensity().Param().r_max);
  const double sigma = 0.1 * profile.sigma_eff;
  const double sigma_in =
      sigma - sigma * sigma / (16.0 * gra::math::PI * profile.slope * gra::PDG::GeV2fm * gra::PDG::GeV2fm);
  double mean = 0.0, local = 0.0, response = 0.0;
  for (const auto &i : indices(rule.first)) {
    const double t    = nucleus.MatterDensity().Thick(rule.first[i]);
    const double area = 2.0 * gra::math::PI * rule.first[i] * rule.second[i];
    mean += a * area * (-std::expm1(-opacity * t)) / opacity;
    local += area * t * std::exp(-a * sigma_in * t);
    response += area * t * std::exp(-opacity * t);
  }
  const double          variance = a * (local - response * response);
  std::array<double, 2> error{1.0, 1.0};
  for (const unsigned int nodes : {128U, 256U}) {
    auto param    = PhotoControls();
    param.b_nodes = param.z_nodes = nodes;
    param.table.qt_nodes = param.table.qz_nodes = param.table.series_terms = 4;
    const gra::nuclear::MPhoto  model(nucleus, param);
    const auto                  value = model.Factors(profile, 0.0, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
    const std::array<double, 2> current = {std::abs(value.coherent.real() / mean - 1.0),
                                           std::abs(value.incoherent * value.incoherent / variance - 1.0)};
    for (const auto &i : indices(error)) { REQUIRE(current[i] < error[i]); }
    error = current;
    REQUIRE(std::abs(value.coherent.imag()) < 1e-13);
  }
  for (const double value : error) { REQUIRE(value < 1e-7); }
}

// Reject missing nuclear identities before configuration samplers can dereference them
TEST_CASE("UPC validates beam identities before preparing nuclear models", "[gra::nuclear][upc][validation]") {
  auto                  param = UPCControls();
  gra::nuclear::UPCMode mode;
  mode.screening = false;
  mode.emission = mode.target                      = {false, false};
  const std::array<gra::nuclear::BeamType, 2> type = {gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus};
  REQUIRE_THROWS_AS(gra::nuclear::MUPC(type, {}, param, {}, mode), std::invalid_argument);
}

// Check the published long-coherence attenuation independently of the optical moment cache
TEST_CASE("Smooth Glauber incoherent attenuation uses the full target thickness", "[gra::nuclear][photo][physics]") {
  const auto                 nucleus = MakeLead();
  const gra::nuclear::MPhoto photo(nucleus, PhotoControls());
  const auto                 rule = gra::math::GaussLegendreRule(256, 0.0, nucleus.MatterDensity().Param().r_max);
  for (const double sigma : {3.0, 15.0}) {
    const auto   profile   = gra::test::PhotoProfile(sigma, 0.0, 0.0);
    const double inelastic = sigma - sigma * sigma / (16.0 * gra::math::PI * profile.slope * gra::PDG::GeV2mb);
    double       expected  = 0.0;
    for (const auto &i : indices(rule.first)) {
      const double t = nucleus.A() * nucleus.MatterDensity().Thick(rule.first[i]);
      expected += 2.0 * gra::math::PI * rule.first[i] * rule.second[i] * t * std::exp(-0.1 * inelastic * t);
    }
    const auto result = photo.Factors(profile, 0.3, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
    CHECK(result.incoherent * result.incoherent == Approx(expected).epsilon(0.002));
  }
}

// Check the correlated production and rescattering eigenstates against independent thickness integrals
// [REFERENCE: Frankfurt et al., Phys. Lett. B752 (2016) 51, Eq. (10)]
TEST_CASE("Photonuclear eigenstates weight production and attenuation together", "[gra::nuclear][photo][physics]") {
  using namespace gra::nuclear;
  const auto nucleus = MakeOxygen();
  auto       param   = PhotoControls();
  param.b_nodes = param.z_nodes = 128;
  param.table.qt_nodes = param.table.qz_nodes = 4;
  const auto radial = gra::math::GaussLegendreRule(256, 0.0, nucleus.MatterDensity().Param().r_max);
  for (const double omega : {0.0, 0.3}) {
    param.fluctuation.omega = omega;
    const auto scales       = CrossSectionRule(param.fluctuation);
    for (const double eta : {0.0, 0.2}) {
      const auto                 profile = gra::test::PhotoProfile(20.0, eta, omega);
      const double               a = nucleus.A(), sigma = 0.1 * profile.sigma_eff;
      const std::complex<double> c       = 0.5 * a * sigma * std::complex<double>(1.0, -eta);
      const double               elastic = a * sigma * sigma * (1.0 + eta * eta) /
                             (8.0 * gra::math::PI * profile.slope * gra::PDG::GeV2fm * gra::PDG::GeV2fm);
      std::complex<double> mean = 0.0, response = 0.0;
      double               local = 0.0;
      for (const auto &i : indices(radial.first)) {
        const double t    = nucleus.MatterDensity().Thick(radial.first[i]);
        const double area = 2.0 * gra::math::PI * radial.first[i] * radial.second[i];
        for (const double s : scales) {
          mean += a * area * (1.0 - std::exp(-c * s * t)) / (c * static_cast<double>(scales.size()));
          response += area * t * s * std::exp(-c * s * t) / static_cast<double>(scales.size());
          for (const double u : scales) {
            local += area * t * s * u * std::exp((-c * s - std::conj(c) * u + elastic * s * u / (s + u)) * t).real() /
                     static_cast<double>(scales.size() * scales.size());
          }
        }
      }
      for (const unsigned int terms : {4U, 128U}) {
        param.table.series_terms = terms;
        const MPhoto model(nucleus, param, PhotoModel::Glauber);
        const auto   value = model.Factors(profile, 0.0, 0.0, 0.0, PhotonDirection::PositiveZ);
        CAPTURE(omega, eta, terms, value.coherent, mean);
        REQUIRE(std::abs(value.coherent - mean) < 2.0e-4 * std::abs(mean));
        REQUIRE(value.incoherent * value.incoherent == Approx(a * (local - std::norm(response))).epsilon(0.002));
      }
      if (eta < 0.1 && omega > 0.1) {
        const MPhoto model(nucleus, param, PhotoModel::Glauber);
        auto         fixed = profile;
        fixed.omega        = 0.0;
        REQUIRE(std::abs(mean) < std::abs(model.Factors(fixed, 0.0, 0.0, 0.0, PhotonDirection::PositiveZ).coherent));
      }
    }
  }
}

// Check complex configuration currents against a coherent sum of fixed cross-section eigenstates
TEST_CASE("Photonuclear configuration currents retain production eigenvalues", "[gra::nuclear][photo][physics]") {
  using namespace gra::nuclear;
  const auto nucleus      = MakeOxygen();
  auto       param        = PhotoControls();
  param.fluctuation.omega = 0.3;
  const auto   scales     = CrossSectionRule(param.fluctuation);
  const MPhoto model(nucleus, param, PhotoModel::Glauber);
  const auto   profile = gra::test::PhotoProfile(20.0, 0.2, param.fluctuation.omega);
  const auto   bank    = SampleBank(nucleus, ConfigControls(4), 731U);
  for (const auto direction : {PhotonDirection::PositiveZ, PhotonDirection::NegativeZ}) {
    for (std::size_t sample = 0; sample < bank.Size(); ++sample) {
      std::complex<double> expected = 0.0;
      for (const double s : scales) {
        auto eigenstate  = profile;
        eigenstate.omega = 0.0;
        eigenstate.sigma_eff *= s;
        eigenstate.slope *= s;
        expected += s * model.ShadowCurrent(eigenstate, bank.At(sample), 0.07, -0.03, 0.02, direction);
      }
      expected /= static_cast<double>(scales.size());
      const auto current = model.ShadowCurrent(profile, bank.At(sample), 0.07, -0.03, 0.02, direction);
      REQUIRE(std::abs(current - expected) < 1.0e-10 * std::max(1.0, std::abs(expected)));
    }
  }
}

// Check the even form-factor slope and the resulting forward incoherent photon density
TEST_CASE("Nuclear form interpolation preserves the forward charge variance", "[gra::nuclear][numerics]") {
  auto param              = MakeLead().Param();
  param.charge.form_q_max = 8.0;
  const gra::nuclear::MNucleus nucleus(param);
  const auto&                  density = nucleus.ChargeDensity();
  const gra::nuclear::MPhoton  photon(nucleus);
  const double                 slope = gra::math::pow2(density.Rms() / gra::PDG::GeV2fm) / 6.0;
  for (const double q : {1.0e-4, 1.0e-5, 1.0e-7}) {
    CAPTURE(q);
    REQUIRE((1.0 - density.Form(q)) / (q * q) ==
            Approx(slope).epsilon(3.0e-6).margin(std::numeric_limits<double>::epsilon() / (q * q)));
    const double xi         = q / (2.0 * nucleus.Mass());
    const auto   incoherent = photon.Density(gra::nuclear::CoherenceType::Incoherent, xi, -q * q, q);
    REQUIRE(incoherent.parallel > 0.0);
    REQUIRE(incoherent.perpendicular > 0.0);
  }
}

// Check enclosed charge against an independent integral of the longitudinal density projection
TEST_CASE("Nuclear cylinder charge agrees with the transverse thickness", "[gra::nuclear][numerics]") {
  const auto  nucleus = MakeLead();
  const auto& density = nucleus.ChargeDensity();
  for (const double b : {1.0e-5, 0.001, 0.01, 0.1, 1.0, 3.0, 5.0, 6.6, 8.0}) {
    const auto rule     = gra::math::GaussLegendreRule(192, 0.0, b);
    double     expected = 0.0;
    for (const auto& i : indices(rule.first)) {
      expected += 2.0 * gra::math::PI * rule.first[i] * rule.second[i] * density.Thick(rule.first[i]);
    }
    CAPTURE(b);
    REQUIRE(density.Cylinder(b) == Approx(expected).epsilon(1.0e-9));
  }
}

// Check the analytic x K1(x) limit without evaluating a divergent Bessel function
TEST_CASE("Nuclear EPA kernel remains finite at vanishing photon energy", "[gra::nuclear][numerics]") {
  for (const double x : {0.0, std::numeric_limits<double>::denorm_min(), 1.0e-300, 1.0e-20, 1.0e-10}) {
    CAPTURE(x);
    REQUIRE(gra::nuclear::PhotonKernel(x, 100.0) == Approx(1.0).margin(2.0e-14));
  }
}

// Check the continuous impulse limit when the leading twist rescattering cross section vanishes
TEST_CASE("Leading twist reaches its upper x limit without roundoff rejection", "[gra::nuclear][numerics]") {
  const auto                  param = PhotoControls();
  const gra::nuclear::MShadow shadow(param.shadow);
  const double                x  = std::nextafter(param.shadow.x_max, 0.0);
  const auto                  xs = shadow.CrossSections(x, 3.0);
  REQUIRE(xs.sigma2 >= 0.0);
  REQUIRE(xs.sigma3 > 0.0);
  REQUIRE(xs.sigma3_in > 0.0);
  REQUIRE(xs.sigma3_in <= xs.sigma3);
  REQUIRE(xs.sigma3 < 1.0e-12);
  const auto                 nucleus = MakeOxygen();
  const gra::nuclear::MPhoto lta(nucleus, param, gra::nuclear::PhotoModel::LTA);
  const gra::nuclear::MPhoto impulse(nucleus, param, gra::nuclear::PhotoModel::Impulse);
  auto                       profile = gra::test::PhotoProfile();
  profile.x                          = x;
  const auto expected = impulse.Factors(profile, 0.04, -0.03, 0.02, gra::nuclear::PhotonDirection::PositiveZ);
  const auto value    = lta.Factors(profile, 0.04, -0.03, 0.02, gra::nuclear::PhotonDirection::PositiveZ);
  REQUIRE(std::abs(value.coherent - expected.coherent) < 1.0e-12);
  REQUIRE(value.incoherent == Approx(expected.incoherent).margin(1.0e-12));
}

// Reject unresolved focused screening rules before any event-dependent annulus can be sampled
TEST_CASE("UPC rejects empty and insufficient focused quadratures at initialization", "[gra::nuclear][numerics]") {
  const auto            oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  gra::nuclear::UPCMode mode;
  mode.target = {false, false};
  for (const auto counts : {std::array<unsigned int, 2>{0, 4}, {1, 4}, {8, 0}}) {
    auto param                  = UPCControls();
    param.loop.radial_intervals = counts[0];
    param.loop.azimuth_nodes    = counts[1];
    CAPTURE(counts[0], counts[1]);
    REQUIRE_THROWS_AS(gra::nuclear::MUPC({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus},
                                         {nullptr, oxygen}, param, {}, mode),
                      std::invalid_argument);
  }
}

// Validate sampled hotspot inputs while the immutable UPC model is constructed
TEST_CASE("UPC rejects invalid hotspot controls before event sampling", "[gra::nuclear][numerics]") {
  const auto            oxygen = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  gra::nuclear::UPCMode mode;
  mode.target      = {false, false};
  mode.screening   = false;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (const auto hotspot : {gra::nuclear::HotSpotParam{0, 1.0, 1.0, 0.2},
                             {3, -1.0, 1.0, 0.2},
                             {3, 1.0, -1.0, 0.2},
                             {3, 1.0, 1.0, -0.2},
                             {3, 1.0, 1.0, nan}}) {
    auto param             = UPCControls();
    param.structure        = gra::nuclear::StructureType::Hotspot;
    param.photo[1].hotspot = hotspot;
    REQUIRE_THROWS_AS(gra::nuclear::MUPC({gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Nucleus},
                                         {nullptr, oxygen}, param, {}, mode),
                      std::invalid_argument);
  }
}

// Check neutral isovector production against fixed proton and neutron number sums
TEST_CASE("Photonuclear isospin preserves coherent signs and correlated fluctuations", "[gra::nuclear][photo][isospin]") {
  using namespace gra::nuclear;
  const auto oxygen = MakeOxygen();
  const auto bank = SampleBank(oxygen, ConfigControls(32), 51317);
  const gra::M3Vec q = {0.12, -0.07, 0.035};
  auto scalar = gra::test::PhotoProfile(20.0, 0.12);
  auto vector = scalar;
  vector.isospin = PhotoIsospin::Isovector;
  for (const auto model : {PhotoModel::Impulse, PhotoModel::Glauber}) {
    const MPhoto photo(oxygen, PhotoControls(), model);
    const auto a = photo.Factors(scalar, q[0], q[1], q[2], PhotonDirection::PositiveZ, &bank, nullptr, true);
    const auto b = photo.Factors(vector, q[0], q[1], q[2], PhotonDirection::PositiveZ, &bank, nullptr, true);
    REQUIRE(std::abs(b.coherent) < 1.0e-11);
    REQUIRE(b.incoherent > 0.0);
    REQUIRE(a.current.has_value());
    REQUIRE(b.current.has_value());
    const std::complex<double> g0{0.7, 0.2}, g1{-0.3, 0.8};
    std::vector<std::complex<double>> mixed, direct;
    for (const auto i : indices(a.current->sample)) {
      const auto& config = bank.At(i);
      const auto j0 = photo.ShadowCurrent(scalar, config, q[0], q[1], q[2], PhotonDirection::PositiveZ);
      const auto j1 = photo.ShadowCurrent(vector, config, q[0], q[1], q[2], PhotonDirection::PositiveZ);
      mixed.push_back(g0 * a.current->sample[i] + g1 * b.current->sample[i]);
      direct.push_back(g0 * j0 + g1 * j1);
      if (model == PhotoModel::Impulse) {
        CHECK(std::abs(j1 - (2.0 * config.ChargeCurrent(q[0], q[1], q[2]) - config.MatterCurrent(q[0], q[1], q[2]))) < 1.0e-12);
      }
      auto swapped = config.Nucleons();
      auto reflected = swapped;
      auto rotated = swapped;
      for (const auto k : indices(swapped)) {
        swapped[k].type = swapped[k].type == NucleonType::Proton ? NucleonType::Neutron : NucleonType::Proton;
        reflected[k].x[2] *= -1.0;
        rotated[k].x = {-rotated[k].x[1], rotated[k].x[0], rotated[k].x[2]};
      }
      const auto opposite = photo.ShadowCurrent(vector, MConfig(16, 8, swapped), q[0], q[1], q[2], PhotonDirection::PositiveZ);
      const auto backward = photo.ShadowCurrent(vector, MConfig(16, 8, reflected), q[0], q[1], -q[2], PhotonDirection::NegativeZ);
      const auto turned = photo.ShadowCurrent(vector, MConfig(16, 8, rotated), -q[1], q[0], q[2], PhotonDirection::PositiveZ);
      CHECK(std::abs(opposite + j1) < 1.0e-12);
      CHECK(std::abs(backward - j1) < 1.0e-12);
      CHECK(std::abs(turned - j1) < 1.0e-12);
    }
    const auto actual = CurrentMoments(mixed), expected = CurrentMoments(direct);
    CHECK(actual.variance == Approx(expected.variance).epsilon(1.0e-12));
    CHECK(std::abs(actual.mean - (g0 * a.coherent + g1 * b.coherent)) < 1.0e-12);
  }
}

// Separate a neutron excess from a difference in the proton and neutron radii
TEST_CASE("Photonuclear isovector densities retain neutron skin and optical attenuation", "[gra::nuclear][photo][isospin]") {
  using namespace gra::nuclear;
  auto scalar = gra::test::PhotoProfile(20.0, 0.12);
  auto vector = scalar;
  vector.isospin = PhotoIsospin::Isovector;
  auto isotope = MakeOxygen().Param();
  isotope.pdg = EncodeNuclearPDG(16, 6);
  const MNucleus target(isotope);
  for (const auto model : {PhotoModel::Impulse, PhotoModel::Glauber}) {
    auto controls = PhotoControls();
    const MPhoto photo(target, controls, model);
    for (const double qt : {0.0, 0.08, controls.table.qt_max + 0.01}) {
      for (const auto direction : {PhotonDirection::PositiveZ, PhotonDirection::NegativeZ}) {
        const auto a = photo.Factors(scalar, qt, 0.0, 0.03, direction, nullptr, nullptr, false, CoherenceType::Coherent);
        const auto b = photo.Factors(vector, qt, 0.0, 0.03, direction, nullptr, nullptr, false, CoherenceType::Coherent);
        CHECK(std::abs(b.coherent + 0.25 * a.coherent) < 1.0e-10);
      }
    }
  }
  const auto skin = MakeNeutronSkin();
  const MPhoto impulse(skin, PhotoControls(), PhotoModel::Impulse);
  for (const double q : {0.0, 0.05, 0.12}) {
    const auto amplitude = impulse.Factors(vector, q, 0.0, 0.0, PhotonDirection::PositiveZ);
    const double expected = 2.0 * skin.Z() * skin.ChargeDensity().Form(q) - skin.A() * skin.MatterDensity().Form(q);
    CHECK(std::abs(amplitude.coherent - expected) < 1.0e-11);
  }
}

// Check the fixed-cross-section limit and covariance of explicit photonuclear fluctuations
TEST_CASE("Photonuclear fluctuations have an explicit zero-fluctuation limit", "[gra::nuclear][ggcf][photo][physics]") {
  using namespace gra::nuclear;
  const auto nucleus = MakeOxygen();
  auto param = PhotoControls();
  const MPhoto ordinary(nucleus, param, PhotoModel::Glauber);
  param.fluctuation.omega = 0.2;
  const MPhoto fluctuating(nucleus, param, PhotoModel::Glauber);
  auto profile = gra::test::PhotoProfile(20.0, 0.1, 0.0);
  const auto fixed = ordinary.Factors(profile, 0.02, -0.03, 0.01, PhotonDirection::PositiveZ).coherent;
  const auto limit = fluctuating.Factors(profile, 0.02, -0.03, 0.01, PhotonDirection::PositiveZ).coherent;
  profile.omega = 0.2;
  const auto ggcf = fluctuating.Factors(profile, 0.02, -0.03, 0.01, PhotonDirection::PositiveZ).coherent;
  CHECK(std::abs(fixed - limit) < 1.0e-10 * std::abs(fixed));
  CHECK(std::abs(ggcf - fixed) > 1.0e-3 * std::abs(fixed));
  for (const auto *model : {&ordinary, &fluctuating}) {
    for (const double omega : {0.0, 0.2}) {
      profile.omega = omega;
      const auto reference = model->Factors(profile, 0.02, -0.03, 0.01, PhotonDirection::PositiveZ).coherent;
      const auto rotated = model->Factors(profile, -0.03, -0.02, 0.01, PhotonDirection::PositiveZ).coherent;
      const auto reversed = model->Factors(profile, 0.02, -0.03, -0.01, PhotonDirection::NegativeZ).coherent;
      CHECK(std::abs(rotated - reference) < 1.0e-10 * std::abs(reference));
      CHECK(std::abs(reversed - reference) < 1.0e-10 * std::abs(reference));
    }
  }
}

// Compare the complete complex convolution to an analytic two-source Gaussian amplitude
TEST_CASE("Excitation operators preserve complex photon amplitudes and rotations", "[gra::nuclear][EMD][amplitude]") {
  const auto nucleus = std::make_shared<const gra::nuclear::MNucleus>(MakeOxygen());
  auto param = UPCControls();
  param.structure = gra::nuclear::StructureType::Smooth;
  param.config.count = 0;
  param.loop.radial_intervals = 96;
  param.loop.azimuth_nodes = 64;
  gra::nuclear::UPCMode mode;
  mode.screening = false;
  const gra::nuclear::MUPC model({gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus},
                                 {nucleus, nucleus}, param, {}, mode);
  const double beta = 0.005, scale = 4.0 * beta * gra::PDG::GeV2fm * gra::PDG::GeV2fm;
  auto channel = std::make_shared<gra::nuclear::ExcitationChannel>();
  channel->born = GENERATE(0.0, 0.4);
  channel->radius = 20.0 / std::sqrt(beta);
  for (unsigned int i = 0; i <= 2048; ++i) {
    const double k = 1.0e-10 * std::pow(param.loop.r_max / 1.0e-10, static_cast<double>(i) / 2048);
    channel->momentum.push_back(k);
    channel->bare.push_back(-4.0 * gra::math::PI / scale * std::exp(-k * k / scale));
  }
  channel->screened = channel->bare;
  const auto upc = model.WithExcitation(channel, false);
  const std::array<std::complex<double>, 2> coefficient = {{{0.7, 0.4}, {-0.3, 0.8}}};
  const std::array<double, 2> width = {50.0, 120.0};
  for (const bool photo : {false, true}) {
    gra::nuclear::ScreenLayout layout;
    if (photo) { layout.photo = gra::nuclear::PhotoChannels(*upc); }
    for (const double angle : {0.0, 0.7, 2.1}) {
      const double c = std::cos(angle), s = std::sin(angle);
      const std::array<double, 2> q1 = {0.05 * c, 0.05 * s}, q2 = {-0.03 * c - 0.02 * s, -0.03 * s + 0.02 * c};
      const auto point = [&](double kx, double ky) {
        gra::nuclear::ScreenPoint p;
        p.transfer = {{{q1[0] - kx, q1[1] - ky, 0.01}, {q2[0] + kx, q2[1] + ky, 0.015}}};
        const double norm = gra::math::pow2(q1[0] - kx) + gra::math::pow2(q1[1] - ky) +
                            gra::math::pow2(q2[0] + kx) + gra::math::pow2(q2[1] + ky);
        for (const auto& h : indices(width)) { p.amplitude.push_back(coefficient[h] * std::exp(-width[h] * norm)); }
        return p;
      };
      const auto born = point(0.0, 0.0);
      gra::nuclear::MUPCScreen screen(*upc, layout, born);
      for (const auto& node : upc->Nodes(born.transfer)) { screen.Add(node, point(node.kx, node.ky)); }
      const auto result = screen.Result();
      std::array<std::complex<double>, 2> expected;
      for (const auto& h : indices(expected)) {
        const double a = width[h], denominator = 2.0 * a + 1.0 / scale;
        expected[h] = coefficient[h] / (scale * denominator) *
                      std::exp(-a * (0.05 * 0.05 + 0.03 * 0.03 + 0.02 * 0.02) +
                               a * a * (0.08 * 0.08 + 0.02 * 0.02) / denominator) +
                      channel->born * born.amplitude[h];
      }
      if (photo) {
        REQUIRE(result.amplitude.size() == 2);
        CHECK(std::abs(result.amplitude[1]) < 1.0e-14);
        CHECK(std::abs(result.amplitude[0] - expected[0] - expected[1]) < 2.0e-4);
      } else {
        REQUIRE(result.amplitude.size() == expected.size());
        for (const auto& h : indices(expected)) { CHECK(std::abs(result.amplitude[h] - expected[h]) < 2.0e-4); }
      }
    }
  }
  // Attaching a channel cannot mutate the shared inclusive model
  CHECK(model.BornWeight() == Approx(1.0));
  CHECK_FALSE(model.Convolution());
  CHECK(upc->Convolution());
}

// Check the EMD times configuration-survival identity before the Good-Walker projection
TEST_CASE("EMD configuration kernels reproduce direct impact integrals", "[gra::nuclear][EMD][mc_ggcf][amplitude]") {
  const auto model = ConfigUPC(2, 817241, 64, 48, 512, 128, 53);
  const double beta = 0.01, slope = 1000.0, hbarc = gra::PDG::GeV2fm;
  const std::complex<double> phase(0.6, -0.8);
  auto channel = std::make_shared<gra::nuclear::ExcitationChannel>();
  channel->born = GENERATE(0.0, 0.4);
  channel->radius = 200.0;
  for (unsigned int i = 0; i <= 4096; ++i) {
    const double k = 1.0e-10 * std::pow(model.Param().loop.r_max / 1.0e-10, static_cast<double>(i) / 4096);
    channel->momentum.push_back(k);
    channel->bare.push_back(-gra::math::PI / (beta * hbarc * hbarc) * std::exp(-k * k / (4.0 * beta * hbarc * hbarc)));
    const double b = 1.0e-5 * std::pow(channel->radius / 1.0e-5, static_cast<double>(i) / 4096);
    channel->impact.push_back(b);
    channel->amplitude.push_back(channel->born + std::exp(-beta * b * b));
  }
  const auto event = model.WithExcitation(channel, true);
  const auto [radius, weight] = gra::math::GaussLegendreRule(192, 0.0, 100.0);
  std::vector<std::array<double, 2>> impact;
  const unsigned int angles = 52;
  for (const auto b : radius) {
    for (unsigned int j = 0; j < angles; ++j) {
      const double phi = 2.0 * gra::math::PI * (j + 0.5) / static_cast<double>(angles);
      impact.push_back({b * std::cos(phi), b * std::sin(phi)});
    }
  }
  for (const auto& q : std::array<gra::M3Vec, 3>{{{0.0, 0.0, 0.0}, {0.02, 0.01, 0.0}, {-0.02, -0.01, 0.0}}}) {
    const double ratio = -0.4;
    const std::array<gra::M3Vec, 2> transfer = {q, {ratio * q[0], ratio * q[1], 0.0}};
    const auto hard = [&](double kx, double ky) {
      const double norm = gra::math::pow2(transfer[0][0] - kx) + gra::math::pow2(transfer[0][1] - ky) +
                          gra::math::pow2(transfer[1][0] + kx) + gra::math::pow2(transfer[1][1] + ky);
      return phase * std::exp(-slope * norm);
    };
    std::vector<std::complex<double>> amplitude(event->SampleCount(), channel->born * hard(0.0, 0.0));
    for (const auto& k : event->Nodes(transfer)) { gra::AddScaled(amplitude, k.sample_weight, hard(k.kx, k.ky)); }
    for (const auto& sample : indices(amplitude)) {
      const auto survival = model.Glauber()->ConfigPairAmp(impact, model.Bank(1), model.Bank(2),
                                                            model.SampleIndex(sample, 1), model.SampleIndex(sample, 2));
      std::complex<double> integral = 0.0;
      for (const auto& i : indices(radius)) {
        for (unsigned int j = 0; j < angles; ++j) {
          const auto& b = impact[i * angles + j];
          const double angle = (1.0 - ratio) * (q[0] * b[0] + q[1] * b[1]) / (2.0 * hbarc);
          integral += weight[i] * radius[i] * survival[i * angles + j] *
                      (channel->born + std::exp(-beta * radius[i] * radius[i])) * std::polar(1.0, angle) *
                      std::exp(-radius[i] * radius[i] / (8.0 * slope * hbarc * hbarc)) / static_cast<double>(angles);
        }
      }
      const auto expected = phase * integral * std::exp(-0.5 * slope * gra::math::pow2(1.0 + ratio) * (q[0] * q[0] + q[1] * q[1])) /
                            (4.0 * slope * hbarc * hbarc);
      CAPTURE(q[0], q[1], sample, amplitude[sample], expected);
      CHECK(std::abs(amplitude[sample] - expected) < 0.005 * std::abs(expected));
    }
    const auto sectors = gra::nuclear::NuclearGoodWalkerProject(amplitude, event->SampleShape());
    CHECK(gra::Sum(sectors[0]) + gra::Sum(sectors[1]) == Approx(gra::SquaredNorm(amplitude) / amplitude.size()).epsilon(1.0e-12));
  }
}

// Exercise the neutron grammar and all comparisons on actual integer multiplicities
TEST_CASE("Neutron selection grammar and count comparisons", "[nuclear][neutron][steering]") {
  using gra::nuclear::ParseNeutronSelection;
  using gra::nuclear::NeutronName;
  for (std::size_t n = 0; n < 5; ++n) {
    CHECK(ParseNeutronSelection(" * ").Accept(n));
    CHECK(ParseNeutronSelection("n == 2").Accept(n) == (n == 2));
    CHECK(ParseNeutronSelection("n!=2").Accept(n) == (n != 2));
    CHECK(ParseNeutronSelection("n > 2").Accept(n) == (n > 2));
    CHECK(ParseNeutronSelection("n>=2").Accept(n) == (n >= 2));
    CHECK(ParseNeutronSelection("n < 2").Accept(n) == (n < 2));
    CHECK(ParseNeutronSelection("n<=2").Accept(n) == (n <= 2));
  }
  for (const auto text : {"*", "n == 0", "n != 1", "n > 0", "n >= 2", "n < 3", "n <= 4"}) {
    CHECK(NeutronName(ParseNeutronSelection(text)) == text);
  }
  CHECK(NeutronName(ParseNeutronSelection(" n>=002 ")) == "n >= 2");
  for (const auto text : {"n=0", "n==0", "n ==0", "n = 0", " n=0 "}) {
    const auto selection = ParseNeutronSelection(text);
    CHECK(selection.Accept(0));
    CHECK_FALSE(selection.Accept(1));
    CHECK(NeutronName(selection) == "n == 0");
  }
  for (const auto text : {"any", "ANY", "0n", "xn", "Xn", "", "N > 0", "n === 1", "n > -1", "n == +1",
                          "n == 1.0", "n == 1 2", "n > = 1", "n > 0 && n < 2", "n == 999999999999999999999999999"}) {
    CHECK_THROWS_AS(ParseNeutronSelection(text), std::invalid_argument);
  }
}
