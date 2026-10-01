// Forward baryon excitation and cylinder fragmentation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <iterator>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Regge/MFragment.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "json.hpp"
using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

using gra::PDG::mp;
using gra::PDG::mpi;

using gra::PDG::PDG_gamma;
using gra::PDG::PDG_n;
using gra::PDG::PDG_p;
using gra::PDG::PDG_pi0;
using gra::PDG::PDG_pim;
using gra::PDG::PDG_pip;

namespace gra {
namespace {

// Reject unknown fields in one fragmentation steering block
void ValidateFields(const nlohmann::json &block,
                    const std::set<std::string> &allowed,
                    const std::string &context) {
  if (!block.is_object()) {
    throw std::invalid_argument(context + " must be an object");
  }
  for (const auto &[field, value] : block.items()) {
    (void)value;
    if (!allowed.contains(field)) {
      throw std::invalid_argument(context + " has unknown field " + field);
    }
  }
}

// Compute true when one probability is finite and lies in [0,1]
bool ValidProbability(const double value) {
  return std::isfinite(value) && value >= 0.0 && value <= 1.0;
}

// Read one strictly positive unsigned numerical control
unsigned int ReadCount(const nlohmann::json &value,
                       const std::string &context) {
  unsigned long long count = 0;
  if (value.is_number_unsigned()) {
    count = value.get<unsigned long long>();
  } else if (value.is_number_integer()) {
    const long long signed_count = value.get<long long>();
    if (signed_count <= 0) {
      throw std::invalid_argument(context + " must be positive");
    }
    count = static_cast<unsigned long long>(signed_count);
  } else {
    throw std::invalid_argument(context + " must be an integer");
  }
  if (count == 0 || count > std::numeric_limits<unsigned int>::max()) {
    throw std::invalid_argument(context + " is outside unsigned range");
  }
  return static_cast<unsigned int>(count);
}

// Validate and normalize one vector of relative probabilities
template <typename Container>
void NormalizeWeights(Container &weights, const std::string &context) {
  if (!gra::AllFinite(weights) ||
      std::any_of(weights.begin(), weights.end(),
                  [](const double value) { return value < 0.0; })) {
    throw std::invalid_argument(context +
                                " requires finite nonnegative weights");
  }
  weights = gra::NormalizedSum(weights);
}

} // namespace

// Compute the class prescriptions with the photon-target override resolved
DissociationModel MNstarParam::Model(const std::string &istate) const {
  const auto entry = model.find(istate);
  return entry == model.end() ? DissociationModel{} : entry->second;
}

// Read one immutable PARAM_NSTAR block
void MNstarParam::Configure(const nlohmann::json &block,
                            const std::string &source) {
  try {
    ValidateFields(block,
                   {"MODEL", "single_side_prob", "mass_margin", "FEWBODY",
                    "CYLINDER", "STRING"},
                   "PARAM_NSTAR");

    const std::map<std::string, DissociationType> types = {
        {"soft", DissociationType::Soft}, {"hera", DissociationType::Hera},
        {"structure", DissociationType::Structure}, {"triple_regge", DissociationType::TripleRegge}};
    const auto &models = block.at("MODEL");
    if (!models.is_object()) { throw std::invalid_argument("PARAM_NSTAR.MODEL must be a process map"); }
    model.clear();
    for (const auto &[key, value] : models.items()) {
      const auto read = [&](const nlohmann::json &entry) {
        const auto type = types.find(entry.get<std::string>());
        if (type == types.end()) { throw std::invalid_argument("PARAM_NSTAR.MODEL: unknown prescription for " + key); }
        return type->second;
      };
      DissociationModel rule;
      if (value.is_string()) {
        rule.hadron = rule.photo = read(value);
      } else {
        if (!value.is_object() || (!value.contains("[*]") && !(key == "ygg" && value.contains("[22,P]")))) {
          throw std::invalid_argument("PARAM_NSTAR.MODEL." + key + " requires a prescription or [*] (or [22,P] for ygg)");
        }
        if (value.contains("[*]")) { rule.hadron = rule.photo = read(value.at("[*]")); }
        for (const auto &[selector, entry] : value.items()) {
          if (selector == "[*]") { continue; }
          if (selector == "[22]") {
            if (key != "MP" && key != "XP" && key != "GP" && key != "TP" && key != "ygg" && key != "yy") {
              throw std::invalid_argument("PARAM_NSTAR.MODEL: [22] requires a photon-emitting class");
            }
            if (read(entry) != DissociationType::Structure) {
              throw std::invalid_argument("PARAM_NSTAR.MODEL: [22] describes photon emission and requires structure");
            }
          } else if (selector == "[22,P]") {
            if (key != "MP" && key != "XP" && key != "GP" && key != "TP" && key != "ygg") {
              throw std::invalid_argument("PARAM_NSTAR.MODEL: [22,P] requires a photoproduction class");
            }
            rule.photo = read(entry);
          } else {
            throw std::invalid_argument("PARAM_NSTAR.MODEL: unknown exchange selector " + selector);
          }
        }
      }
      model.emplace(key, rule);
    }

    single_side_prob = block.at("single_side_prob").get<double>();
    mass_margin = block.at("mass_margin").get<double>();
    if (!std::isfinite(single_side_prob) || !(single_side_prob > 0.0 && single_side_prob < 1.0)) {
      throw std::invalid_argument("PARAM_NSTAR.single_side_prob must be in (0,1) to sample both beams");
    }
    if (!std::isfinite(mass_margin) || mass_margin < 0.0) {
      throw std::invalid_argument("PARAM_NSTAR has invalid mass controls");
    }
    const auto &few = block.at("FEWBODY");
    ValidateFields(few, {"body_br", "two_body_br", "three_body_br", "unweight"},
                   "PARAM_NSTAR.FEWBODY");
    fewbody.body_br = few.at("body_br").get<decltype(fewbody.body_br)>();
    fewbody.two_body_br =
        few.at("two_body_br").get<decltype(fewbody.two_body_br)>();
    fewbody.three_body_br =
        few.at("three_body_br").get<decltype(fewbody.three_body_br)>();
    fewbody.unweight = few.at("unweight").get<bool>();
    NormalizeWeights(fewbody.body_br, "PARAM_NSTAR.FEWBODY.body_br");
    NormalizeWeights(fewbody.two_body_br, "PARAM_NSTAR.FEWBODY.two_body_br");
    NormalizeWeights(fewbody.three_body_br,
                     "PARAM_NSTAR.FEWBODY.three_body_br");
    const auto &cyl = block.at("CYLINDER");
    ValidateFields(cyl,
                   {"charged_prob",
                    "pair_prob",
                    "particle_prob",
                    "species_ratio",
                    "pt_distribution",
                    "mult_a",
                    "mult_b",
                    "q2_ref",
                    "mass_unit",
                    "q_power",
                    "T",
                    "impact_power",
                    "max_pt",
                    "exp_lambda",
                    "outer_trials",
                    "mult_trials",
                    "pick_trials",
                    "tube_trials",
                    "pt_trials",
                    "pt_bins",
                    "alpha_iter",
                    "alpha_tol"},
                   "PARAM_NSTAR.CYLINDER");
    cylinder.charged_prob = cyl.at("charged_prob").get<double>();
    cylinder.pair_prob = cyl.at("pair_prob").get<double>();
    cylinder.particle_prob = cyl.at("particle_prob").get<double>();
    cylinder.species_ratio =
        cyl.at("species_ratio").get<decltype(cylinder.species_ratio)>();
    cylinder.pt_distribution = cyl.at("pt_distribution").get<std::string>();
    cylinder.mult_a = cyl.at("mult_a").get<double>();
    cylinder.mult_b = cyl.at("mult_b").get<double>();
    cylinder.q2_ref = cyl.at("q2_ref").get<double>();
    cylinder.mass_unit = cyl.at("mass_unit").get<double>();
    cylinder.q_power = cyl.at("q_power").get<double>();
    cylinder.T = cyl.at("T").get<double>();
    cylinder.impact_power = cyl.at("impact_power").get<double>();
    cylinder.max_pt = cyl.at("max_pt").get<double>();
    cylinder.exp_lambda = cyl.at("exp_lambda").get<double>();
    cylinder.outer_trials =
        ReadCount(cyl.at("outer_trials"), "PARAM_NSTAR.CYLINDER.outer_trials");
    cylinder.mult_trials =
        ReadCount(cyl.at("mult_trials"), "PARAM_NSTAR.CYLINDER.mult_trials");
    cylinder.pick_trials =
        ReadCount(cyl.at("pick_trials"), "PARAM_NSTAR.CYLINDER.pick_trials");
    cylinder.tube_trials =
        ReadCount(cyl.at("tube_trials"), "PARAM_NSTAR.CYLINDER.tube_trials");
    cylinder.pt_trials =
        ReadCount(cyl.at("pt_trials"), "PARAM_NSTAR.CYLINDER.pt_trials");
    cylinder.pt_bins =
        ReadCount(cyl.at("pt_bins"), "PARAM_NSTAR.CYLINDER.pt_bins");
    cylinder.alpha_iter =
        ReadCount(cyl.at("alpha_iter"), "PARAM_NSTAR.CYLINDER.alpha_iter");
    cylinder.alpha_tol = cyl.at("alpha_tol").get<double>();
    NormalizeWeights(cylinder.species_ratio,
                     "PARAM_NSTAR.CYLINDER.species_ratio");
    const double q_min = std::pow(2.0, cylinder.q_power);
    if (!(cylinder.charged_prob > 0.0) || cylinder.charged_prob > 1.0 ||
        !ValidProbability(cylinder.pair_prob) ||
        !ValidProbability(cylinder.particle_prob) ||
        (cylinder.pt_distribution != "powexp" &&
         cylinder.pt_distribution != "exp") ||
        !std::isfinite(cylinder.mult_a) || cylinder.mult_a < 0.0 ||
        !std::isfinite(cylinder.mult_b) || cylinder.mult_b <= 0.0 ||
        !std::isfinite(cylinder.q2_ref) || cylinder.q2_ref <= 0.0 ||
        !std::isfinite(cylinder.mass_unit) || cylinder.mass_unit <= 0.0 ||
        !std::isfinite(cylinder.q_power) || cylinder.q_power < 0.0 ||
        !std::isfinite(cylinder.T) || cylinder.T <= 0.0 ||
        !std::isfinite(cylinder.impact_power) || cylinder.impact_power < 0.0 ||
        !std::isfinite(cylinder.max_pt) || cylinder.max_pt <= 0.0 ||
        !std::isfinite(cylinder.exp_lambda) || cylinder.exp_lambda <= 0.0 ||
        cylinder.outer_trials == 0 || cylinder.mult_trials == 0 ||
        cylinder.pick_trials == 0 || cylinder.tube_trials == 0 ||
        cylinder.pt_trials == 0 || cylinder.pt_bins == 0 ||
        (cylinder.pt_distribution == "powexp" &&
         (!std::isfinite(q_min) || !(q_min > 1.0) || cylinder.pt_bins < 2)) ||
        cylinder.alpha_iter == 0 || !std::isfinite(cylinder.alpha_tol) ||
        cylinder.alpha_tol <= 0.0) {
      throw std::invalid_argument(
          "PARAM_NSTAR.CYLINDER contains invalid values");
    }

    const auto &str = block.at("STRING");
    ValidateFields(
        str,
        {"first_valence_prob", "scalar_diquark_mass", "vector_diquark_mass"},
        "PARAM_NSTAR.STRING");
    string.first_valence_prob = str.at("first_valence_prob").get<double>();
    string.scalar_diquark_mass = str.at("scalar_diquark_mass").get<double>();
    string.vector_diquark_mass = str.at("vector_diquark_mass").get<double>();
    if (!ValidProbability(string.first_valence_prob) ||
        !std::isfinite(string.scalar_diquark_mass) ||
        string.scalar_diquark_mass <= 0.0 ||
        !std::isfinite(string.vector_diquark_mass) ||
        string.vector_diquark_mass <= 0.0) {
      throw std::invalid_argument("PARAM_NSTAR.STRING contains invalid values");
    }
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "MNstarParam::Configure: Error reading " + source + ": " +
        error.what());
  }
}

// Construct one immutable forward excitation parameter block
MNstarParamPtr ReadNstarParam(const MModelTune &tune) {
  auto param = std::make_shared<MNstarParam>();
  param->Configure(tune.General("PARAM_NSTAR"), tune.GeneralFile());
  return param;
}

// Compute the run owned immutable forward excitation parameter block
MNstarParamPtr GetNstarParam(MModelCache &cache) {
  return cache.Get<MNstarParam>("nstar", [&cache] {
    return ReadNstarParam(cache.Tune());
  });
}

// Compute decay status [RESERVATION!]
//
//
void MFragment::GetDecayStatus(const std::vector<int> &pdgcode,
                               std::vector<bool> &isstable) {
  isstable.resize(pdgcode.size(), true);
  for (const auto &i : indices(isstable)) {
    isstable[i] = pdgcode[i] != 111 && std::abs(pdgcode[i]) != 311;
  }
}

// Sample from the Hagedorn transverse-momentum parametrization
// density is proportional to pT mT [1+(q-1)mT/T]^{q/(1-q)}
//
// This function is useful to generate distribution with
// exponential at low pt, powerlaw at high pt
// Does not generate jet like topologies,
// so this can be used as a strict null model (H0) for soft like topologies
//
// [slow function, performance should be upgraded]
//
// Parameters <=>
//  p0 = T/(q-1)
//   n = 1/(q-1)
//
bool MFragment::ExpPowRND(double q, double T, double maxpt,
                          const std::vector<double> &mass,
                          std::vector<double> &x, MRandom &rng,
                          const MNstarParam &param) {
  const unsigned int MAXTRIAL = param.cylinder.pt_trials;
  const unsigned int Nbins = param.cylinder.pt_bins;

  x.assign(mass.size(), 0.0);
  if (MAXTRIAL == 0 || Nbins == 0 || !std::isfinite(q) || q <= 1.0 ||
      !std::isfinite(T) || T <= 0.0 || !std::isfinite(maxpt) || maxpt <= 0.0 ||
      !std::all_of(mass.begin(), mass.end(), [](const double value) {
        return std::isfinite(value) && value >= 0.0;
      })) {
    return false;
  }
  const double ptstep = maxpt / Nbins;

  // Calculate pt-bin values
  std::vector<double> ptval(Nbins);
  for (std::size_t i = 0; i < Nbins; ++i) {
    ptval[i] = i * ptstep;
  }

  // Random integer from [0,NBins-1]
  std::uniform_int_distribution<int> random_bin(0, Nbins - 1);

  // Calculate for pion, kaon, proton masses
  const std::vector<double> fixmass = {PDG::mpi0, PDG::mpi, PDG::mK, PDG::mp};
  std::vector<std::vector<double>> dsdpt(fixmass.size(),
                                         std::vector<double>(Nbins, 0.0));
  std::vector<double> maxval(fixmass.size(), 0.0);

  // -------------------------------------------------------------------------------
  // Very slow

  // Pre-evaluate dsigma/dpt |_y=0 (mid rapidity)
  for (const auto &k : indices(fixmass)) {
    const double m2 = pow2(fixmass[k]);
    for (std::size_t i = 0; i < Nbins; ++i) {
      const double mtval = msqrt(pow2(ptval[i]) + m2);

      // PDF
      dsdpt[k][i] = ptval[i] * mtval *
                    std::pow(1.0 + (q - 1.0) * mtval / T, q / (1.0 - q));
      // Save maximum
      maxval[k] = (dsdpt[k][i] > maxval[k]) ? dsdpt[k][i] : maxval[k];
    }
  }
  // -------------------------------------------------------------------------------

  // For each particle, generate pt
  for (const auto &p : indices(mass)) {
    // Find the best distribution
    unsigned int best = 0;
    double mindist = 1e32;
    for (const auto &k : indices(fixmass)) {
      const double dist = std::abs(fixmass[k] - mass[p]);
      if (dist < mindist) {
        mindist = dist;
        best = k;
      }
    }

    // Acceptance-Rejection, very slow
    bool accepted = false;
    for (unsigned int trials = 0; trials < MAXTRIAL; ++trials) {
      const int BIN = random_bin(rng.rng);
      if (rng.U(0.0, 1.0) * maxval[best] < dsdpt[best][BIN]) {
        x[p] = ptval[BIN];
        accepted = true;
        break;
      }
    }
    if (!accepted) {
      x.assign(mass.size(), 0.0);
      return false;
    }
  }
  return true;
}

// N-body fragmentation with tube (cylinder) phase space
// (approximate target distribution is flat over rapidity, pt from a given
// distribution)
//
// [REFERENCE: Jadach, Computer Physics Communications, 9 (1975) 297-304]
// [REFERENCE: UA5, http://cds.cern.ch/record/173907/files/198701320.pdf]
//
// This function contains a mixture of dynamics and kinematics, i.e.,
// is not a pure phase space and is suitable for soft fragmentation studies
//
// Compute: 1.0 for a valid fragmentation and -1.0 for a kinematically impossible

double MFragment::TubeFragment(const M4Vec &mother, double M0,
                               const std::vector<double> &m,
                               std::vector<M4Vec> &p, double q, double T,
                               double lambda, double maxpt, MRandom &rng,
                               const MNstarParam &param) {
  const std::string &pt_distribution = param.cylinder.pt_distribution;
  if (m.size() < 2) {
    throw std::invalid_argument(
        "MFragment::TubeFragment: at least two particles are required");
  }
  if (pt_distribution != "powexp" && pt_distribution != "exp") {
    throw std::invalid_argument(
        "MFragment::TubeFragment: Unknown pt-distribution parameter = " +
        pt_distribution);
  }
  const bool finite_mother =
      std::isfinite(mother.E()) && std::isfinite(mother.Px()) &&
      std::isfinite(mother.Py()) && std::isfinite(mother.Pz());
  if (!finite_mother || !std::isfinite(M0) || M0 <= 0.0 ||
      !std::isfinite(maxpt) || maxpt <= 0.0 ||
      (pt_distribution == "powexp" &&
       (!std::isfinite(q) || q <= 1.0 || !std::isfinite(T) || T <= 0.0)) ||
      (pt_distribution == "exp" && (!std::isfinite(lambda) || lambda <= 0.0))) {
    return -1.0;
  }
  double threshold = 0.0;
  for (const double mass : m) {
    if (!std::isfinite(mass) || mass < 0.0) {
      return -1.0;
    }
    threshold += mass;
  }
  if (!(threshold < M0)) {
    return -1.0;
  }

  // Number of re-trials in a case of failing kinematics
  const unsigned int MAXTRIAL = param.cylinder.tube_trials;
  unsigned int trials = 0;

  const unsigned int N = m.size();
  std::valarray<double> m2(N);
  for (const auto &i : indices(m)) {
    m2[i] = pow2(m[i]);
  }
  // transverse momentum px, py
  std::valarray<double> px(N);
  std::valarray<double> py(N);
  std::valarray<double> pt2(N);

  // rapidity and transverse mass
  std::valarray<double> y(N);
  std::valarray<double> mt(N);

  while (true) {
    // Random sample px,py and initial rapidity
    std::vector<double> pt(N);

    if (pt_distribution == "powexp") {
      if (!ExpPowRND(q, T, maxpt, m, pt, rng, param)) {
        ++trials;
        if (trials > MAXTRIAL) {
          return -1.0;
        }
        continue;
      }
    } else if (pt_distribution == "exp") {
      for (const auto &i : indices(pt)) {
        pt2[i] = rng.ExpRandom(lambda);
        pt[i] = msqrt(pt2[i]);
      }
    }

    for (const auto &i : indices(pt)) {
      // Sample angle, px, py
      const double phi = rng.U(0, 2.0 * gra::math::PI);
      px[i] = pt[i] * std::cos(phi);
      py[i] = pt[i] * std::sin(phi);

      // Rapidity
      y[i] = rng.U(0.0, 1.0);
    }
    // ---------------------------------------------------------------
    // Span [0,1] without correlating particle labels with rapidity order
    const double rapidity_min = y.min();
    const double rapidity_span = y.max() - rapidity_min;
    if (!(rapidity_span > std::numeric_limits<double>::epsilon())) {
      ++trials;
      if (trials > MAXTRIAL) {
        return -1.0;
      }
      continue;
    }
    y = (y - rapidity_min) / rapidity_span;

    // Set transverse momentum with zero sum
    px = px - px.sum() / N;
    py = py - py.sum() / N;
    pt2 = px * px + py * py;

    // Transverse mass for all particles
    mt = sqrt(m2 + pt2);

    // ---------------------------------------------------------------
    // Solve rapidity scale factor 'alpha'
    double alpha = 0.0;
    if (!SolveAlpha(alpha, M0, m, mt, y, param)) {
      ++trials;
      if (trials > MAXTRIAL) {
        return -1.0;
      } else {
        continue;
      }
    }
    // Scale all rapidities
    y = alpha * y;

    // Set longitudinal momentum and energy using 'alpha'
    // get boost factor to +- 0
    const double Q = (mt * exp(y)).sum();
    y = y + std::log(M0 / Q);
    const std::valarray<double> pz = mt * sinh(y);
    const std::valarray<double> E = mt * cosh(y);

    // ---------------------------------------------------------------
    // Set all particles 4-momenta and boost to the original frame
    p.resize(N);
    const int sign = 1; // positive
    for (const auto &i : indices(p)) {
      p[i] = M4Vec(px[i], py[i], pz[i], E[i]);
      gra::kinematics::LorentzBoost(mother, M0, p[i], sign);
    }

    // Check EM-conservation
    M4Vec p_sum(0, 0, 0, 0);
    for (const auto &n : p) {
      p_sum += n;
    }
    bool valid = gra::math::CheckEMC(p_sum - mother);

    // Check rapidity (floating points can fail after boost in forward)
    for (const auto &i : indices(p)) {
      if (std::isnan(p[i].Rap()) || std::isinf(p[i].Rap())) {
        valid = false;
        break;
      }
    }
    if (!valid) {
      ++trials;
      if (trials > MAXTRIAL) {
        return -1.0;
      } else {
        continue;
      }
    }
    break; // Fragmentation complete
  }
  return 1.0;
}

// ---------------------------------------------------------------
// Solve rapidity scale factor via Newton iteration
// (alternative strategies are viable, too, this does not converge
// always with non-gaussian pt-distributions)
//
// Apply Newton iteration to find a root of f(a)
// a_1 = a_0 - f(a_0)/f'(a_0)
//   iteratively via
// a_{n+1} = a_n - f(a_n) / f'(a_n)
//
bool MFragment::SolveAlpha(double &alpha, double M0,
                           const std::vector<double> &m,
                           const std::valarray<double> &mt,
                           const std::valarray<double> &y,
                           const MNstarParam &param) {
  const double STOP_EPS = param.cylinder.alpha_tol;
  const unsigned int MAXITER = param.cylinder.alpha_iter;

  // Validate the Newton system before indexing its boundary values
  const std::size_t N = m.size();
  if (N < 2 || mt.size() != N || y.size() != N || !std::isfinite(M0) ||
      M0 <= 0.0) {
    return false;
  }
  for (const auto &i : indices(m)) {
    if (!std::isfinite(m[i]) || m[i] < 0.0 || !std::isfinite(mt[i]) ||
        mt[i] <= 0.0 || !std::isfinite(y[i])) {
      return false;
    }
  }

  // Starting value, retaining the mass-based estimate for massive boundaries
  const double C = std::log(M0 * M0);
  const auto [low, high] = std::minmax_element(std::begin(y), std::end(y));
  const auto first = std::distance(std::begin(y), low), last = std::distance(std::begin(y), high);
  double boundary_scale = m[first] * m[last];
  if (!(boundary_scale > 0.0)) {
    boundary_scale = mt[first] * mt[last];
  }
  if (!(boundary_scale > 0.0) || !std::isfinite(boundary_scale)) {
    return false;
  }
  alpha = C - std::log(boundary_scale);
  std::vector<double> E = {0, 0, 0, 0};

  unsigned int iter = 0;
  while (true) {
    const std::valarray<double> x = exp(alpha * y);

    E[0] = (mt * x).sum();
    E[1] = (mt / x).sum();
    E[2] = (mt * y * x).sum(); // Derivative d/dy
    E[3] = (mt * y / x).sum(); // Derivative d/dy

    // Iterate the solution, d/dx ln(x) = 1/x
    const double product = E[0] * E[1];
    const double derivative = E[0] * E[3] - E[1] * E[2];
    if (!(product > 0.0) || !std::isfinite(product) ||
        !std::isfinite(derivative) ||
        std::abs(derivative) <= std::numeric_limits<double>::epsilon()) {
      return false;
    }
    const double DY = product * (C - std::log(product)) / derivative;
    alpha -= DY;

    if (!std::isfinite(DY) || !std::isfinite(alpha)) {
      return false;
    }
    if (std::abs(DY) <= STOP_EPS * std::max(1.0, std::abs(alpha))) {
      return true;
    }
    if ((++iter) > MAXITER) {
      return false;
    }
  }
}

// N* decay table [set manually according to experimental data]
//
// Q is the proton / antiproton charge (1,-1)
// M0 is the N* mass
//
void MFragment::NstarDecayTable(int Q, double M0, std::vector<int> &pdgcodes,
                                MRandom &rng, const MNstarParam &param) {
  pdgcodes.clear();
  if ((Q != -1 && Q != 1) || !std::isfinite(M0) || M0 <= PDG::mp + PDG::mpi) {
    return;
  }
  int decaymode = 0;

  // Only 2-body decay possible, mass below 3-body threshold
  if (M0 < (PDG::mp + 2 * PDG::mpi)) {
    decaymode = 0;

    // 2- or 3-body decay possible
  } else {
    // 2->body / 3-body branching ratios from PDG
    std::discrete_distribution<> d(param.fewbody.body_br.begin(),
                                   param.fewbody.body_br.end());
    decaymode = d(rng.rng); // Draw random
  }

  // ----------------------------------------------------------------------
  // N*(J^P = 1/2^+) Decay Parameters from PDG
  // Only major decays implemented

  std::vector<int> decayids;

  // 2-Body channel
  if (decaymode == 0) {
    // Subchannels
    std::discrete_distribution<> subd(param.fewbody.two_body_br.begin(),
                                      param.fewbody.two_body_br.end());
    const int channel = subd(rng.rng);

    // PDG-ID of subchannels
    const std::vector<std::vector<int>> ID = {
        {Q * PDG_n, Q * PDG_pip}, // (anti)neutron & pi(-)+
        {Q * PDG_p, PDG_pi0}};    // (anti)proton  & pi0

    // Choose the decay channel
    pdgcodes = ID[channel];
  }

  // 3-body channel
  if (decaymode == 1) {
    // Subchannels
    std::discrete_distribution<> subd(param.fewbody.three_body_br.begin(),
                                      param.fewbody.three_body_br.end());
    const int channel = subd(rng.rng);

    // PDG-ID of subchannels
    const std::vector<std::vector<int>> ID = {{Q * PDG_n, Q * PDG_pip, PDG_pi0},
                                              {Q * PDG_p, PDG_pip, PDG_pim},
                                              {Q * PDG_p, PDG_pi0, PDG_pi0}};

    // Perhaps to add
    //{Q*PDG_delta0, Q*PDG_pip, PDG_pi0},    // delta0 & pi+ & pi0
    //{Q*PDG_deltap, PDG_pip, PDG_pim},      // delta+ & pi+ & pi-
    //{Q*PDG_deltap, PDG_pi0, PDG_pi0}};     // delta+ & pi0 & pi0

    // Draw the decay channel
    pdgcodes = ID[channel];
  }
}

// Simple statistical toy particle pick-up, nothing more
//
bool MFragment::PickParticles(double M, unsigned int N, int B, int S, int Q,
                              std::vector<double> &mass,
                              std::vector<int> &pdgcode, const MPDG &PDG,
                              MRandom &rng, const MNstarParam &param) {
  if (!std::isfinite(M) || M <= 0.0 || N < 2) {
    return false;
  }

  // At least |B| nucleons are required, and every other hadron is at least as heavy as the lightest meson
  const auto baryons = static_cast<unsigned int>(std::abs(B));
  if (baryons > N) { return false; }
  const double mN = std::min(PDG.FindByPDG(2212).mass, PDG.FindByPDG(2112).mass);
  const double mM = std::min({PDG.FindByPDG(111).mass, PDG.FindByPDG(211).mass,
                             PDG.FindByPDG(311).mass, PDG.FindByPDG(321).mass});
  const long double threshold = static_cast<long double>(baryons) * mN +
                                static_cast<long double>(N - baryons) * mM;
  if (static_cast<long double>(M) <= threshold) { return false; }

  const unsigned int MAXTRIAL = param.cylinder.pick_trials;

  std::uniform_int_distribution<int> random_species(0, 2);

  // Pion, Kaon, Proton ratios [MEASURED]
  const auto &ratio3 = param.cylinder.species_ratio;

  const std::vector<int> charged_pdg = {211, 321,
                                        2212}; // pi+-, K+-, (anti)proton
  const std::vector<int> charged_B = {0, 0, 1};
  const std::vector<int> charged_S = {0, 1, 0};

  const std::vector<int> neutral_pdg = {111, 311, 2112}; // pi0, K0, neutron
  const std::vector<int> neutral_B = {0, 0, 1};
  const std::vector<int> neutral_S = {0, 1, 0};

  const double QProb = param.cylinder.charged_prob;
  const double DOUBLE = param.cylinder.pair_prob;
  // -----------------------------------------------------------------

  // Now resize!
  mass.resize(N, 0);
  pdgcode.resize(N, 0);
  std::vector<int> Qcharges(N, 0);
  std::vector<int> Bcharges(N, 0);
  std::vector<int> Scharges(N, 0);
  unsigned int trials = 0;

  while (true) {
    unsigned int i = 0;

    do { // Hadron picking

      // Charged
      if (rng.U(0, 1) < QProb) {
        // Charged pair
        if (rng.U(0, 1) < DOUBLE && i + 1 < N) {
          // Sample particle flavour
          while (true) {
            const int bin = random_species(rng.rng);

            if (rng.U(0, 1) < ratio3[bin]) {
              Qcharges[i] = 1;
              Bcharges[i] = charged_B[bin];
              Scharges[i] = charged_S[bin];
              pdgcode[i] = charged_pdg[bin];
              ++i; // Next slot

              Qcharges[i] = -1;
              Bcharges[i] = -charged_B[bin];
              Scharges[i] = -charged_S[bin];
              pdgcode[i] = -charged_pdg[bin];
              ++i;
              break;
            }
          }

          // Single
        } else {
          while (true) {
            const int bin = random_species(rng.rng);

            if (rng.U(0, 1) < ratio3[bin]) {
              const int sign =
                  rng.U(0, 1) < param.cylinder.particle_prob ? 1 : -1;

              Qcharges[i] = sign;
              Bcharges[i] = sign * charged_B[bin];
              Scharges[i] = sign * charged_S[bin];
              pdgcode[i] = sign * charged_pdg[bin];
              ++i;
              break;
            }
          }
        }

        // Neutral
      } else {
        while (true) {
          // Sample particle flavour
          const int bin = random_species(rng.rng);
          if (rng.U(0, 1) < ratio3[bin]) {
            const bool self_conjugate =
                neutral_B[bin] == 0 && neutral_S[bin] == 0;
            const int sign =
                self_conjugate || rng.U(0, 1) < param.cylinder.particle_prob
                    ? 1
                    : -1;
            Qcharges[i] = 0;
            Bcharges[i] = sign * neutral_B[bin];
            Scharges[i] = sign * neutral_S[bin];
            pdgcode[i] = sign * neutral_pdg[bin];
            ++i;
            break;
          }
        }
      }
    } while (i < N); // Picking loop

    // Sum charges
    const int Q_sum = gra::Sum(Qcharges);
    const int B_sum = gra::Sum(Bcharges);
    const int S_sum = gra::Sum(Scharges);

    // Get corresponding masses
    for (std::size_t i = 0; i < N; ++i) {
      const auto particle = PDG.PDG_table.find(pdgcode[i]);
      if (particle == PDG.PDG_table.end()) {
        return false;
      }
      mass[i] = particle->second.mass;
    }
    const double M_sum = gra::Sum(mass);

    // Check charge and mass threshold
    bool Q_check = (Q_sum == Q) ? true : false;
    bool B_check = (B_sum == B) ? true : false;
    bool S_check = (S_sum == S) ? true : false;
    bool M_check = (M_sum < M) ? true : false;

    if (Q_check && B_check && S_check && M_check) {
      break; // we are ok!
    } else {
      ++trials;
      if (trials > MAXTRIAL) {
        return false;
      }
    }
  } // Out loop

  return true;
}

} // namespace gra
