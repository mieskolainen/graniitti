// Analytic GP central Regge spin amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityNorm.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;

namespace gra {

namespace gpom {

using detail::SourceFor;
using detail::ResonanceScale;
using detail::ContractFusion;
using detail::ConProd;

namespace {

// Freeze only the GP angular spin, keeping S+alpha1+alpha2+2 > 1 and 2alpha+1 > 0
// The initialized floor at most one preserves every physical spin-one and spin-two pole
// [REFERENCE: NIST DLMF 16.2.2, https://dlmf.nist.gov/16.2.E2]
double SpinAlpha(double alpha, const regge::Param& param) { return std::max(alpha, param.gp_alpha_min); }

}  // namespace

// Build per-channel analytic resonance production matrices
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const regge::Param& param, const PARAM_RES& res,
                              AmpCache* amp_cache) {
  AmpCache            local;
  auto&               cache    = amp_cache ? *amp_cache : local;
  const auto          nmu      = static_cast<std::size_t>(res.p.spinX2 + 1);
  const auto          nbeam    = spin::SpinHalfTransitions(lts.process.FORWARD_NOFLIP).size();
  const double        momentum = lts.process.DERIVATIVE_FACTOR ? lts.q1_in_X.P3mod() : 0.0;
  std::vector<HelAmp> out;
  out.reserve(res.production.size());
  for (const auto& production : res.production) {
    const auto& hel  = production.hel;
    const auto& tree = production.tree;
    if (!lts.process.SPINGEN) {
      out.push_back(
          spin::Blind(nbeam * nbeam, nmu, ResonanceScale(hel, tree, momentum, lts.process.DERIVATIVE_FACTOR)));
      continue;
    }
    auto& source = SourceFor(lts, param, {lts.process.MMAX, tree[0].p.pdg, tree[1].p.pdg}, cache);
    out.push_back(ContractFusion(hel, source, momentum, lts.process.DERIVATIVE_FACTOR));
  }
  return out;
}

// Build per-channel analytic continuum production matrices
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, const regge::Param& param, AmpCache* amp_cache) {
  AmpCache             local;
  auto&                cache = amp_cache ? *amp_cache : local;
  std::vector<HelPair> channels;
  channels.reserve(lts.process.CONT_PRODUCTION.size());
  bool photon = false;
  for (const auto& exchange : lts.process.CONT_PRODUCTION) {
    photon = photon || std::find(exchange.cbegin(), exchange.cend(), PDG::PDG_gamma) != exchange.cend();
  }
  std::array<M4Vec, 2> final_in_X{};
  if (lts.process.SPINGEN && photon) {
    final_in_X = {
        kinematics::BoostToRestFrame(lts.decaytree[0].p4, lts.pfinal[0], "MReggeGP::Continuum upper final state"),
        kinematics::BoostToRestFrame(lts.decaytree[1].p4, lts.pfinal[0], "MReggeGP::Continuum lower final state")};
  }
  for (const auto& i : indices(lts.process.CONT_PRODUCTION)) {
    channels.push_back(
        ConProd(lts, param, lts.process.CONTINUUM_GP[i], lts.process.CONT_PRODUCTION[i], final_in_X, cache));
  }
  return channels;
}

// Compute the photon transition at the requested transfer using the production kinematics
std::complex<double> PhotoCoupling(const LORENTZSCALAR &lts, const regge::Param &param,
                                  const RES_PRODUCTION &production, double t) {
  const bool upper_gamma = production.tree[0].p.pdg == PDG::PDG_gamma;
  const int exchange = production.tree[upper_gamma ? 1 : 0].p.pdg;
  const auto trajectory = param.exchanges.at(regge::TrajectoryIndex(param, exchange)).soft_exchange;
  const double alpha = SpinAlpha(param.soft_model->Alpha(trajectory, t), param);
  const double photon_t = upper_gamma ? lts.t1 : lts.t2;
  const double mass2 = lts.pfinal[0].M2();
  const double momentum = kinematics::SqrtKallenLambda(mass2, photon_t, t) / (2.0 * std::sqrt(mass2));
  return gpom::PhotoCoupling(production, alpha, momentum, lts.process.DERIVATIVE_FACTOR);
}

// Compute the transverse helicity-conserving GP vertex at the requested trajectory spin
// The target m=0 vertex has no helicity barrier or nonsense factor
std::complex<double> PhotoCoupling(const RES_PRODUCTION& production, double alpha, double momentum,
                                   bool derivative_factor) {
  const auto&   hel = production.hel;
  const bool upper_photon = production.tree[0].p.pdg == PDG::PDG_gamma;
  auto source = detail::PoleSource(hel, upper_photon, !upper_photon);
  auto& basis = source.basis;
  basis.alpha1 = upper_photon ? 1.0 : alpha;
  basis.alpha2 = upper_photon ? alpha : 1.0;
  const auto           central   = Fusion(hel, source, momentum, derivative_factor);
  const int            JX        = static_cast<int>(std::llround(hel.J));
  std::complex<double> amplitude = 0.0;
  for (const int photon : {-1, 1}) {
    const int m1 = basis.upper_photon ? photon : 0;
    const int m2 = basis.lower_photon ? photon : 0;
    const int mu = m1 - m2;
    if (std::abs(mu) > JX) { continue; }
    const std::size_t row =
        static_cast<std::size_t>(m1 + basis.MMAX) * basis.nm + static_cast<std::size_t>(m2 + basis.MMAX);
    amplitude += 0.5 * central[row][static_cast<std::size_t>(JX - mu)];
  }
  return amplitude;
}

// Compute the local GP exchange projection range
std::vector<double> LocalBasis(const ReggeContinuumPole& pole, const std::size_t leg) {
  const auto& basis = pole.gp_vertex[leg];
  return basis.Jz_values;
}

// Evaluate one analytic GP local final-pair kernel
HelAmp PairKernel(const ReggeContinuumPole& pole, const regge::Param& param, const MDecayBranch& first, const MDecayBranch& second,
                  const double alpha_upper, const double alpha_lower) {
  const auto upper      = Crossed(pole.gp_vertex[0], SpinAlpha(alpha_upper, param), M4Vec{}, false);
  const auto lower      = Crossed(pole.gp_vertex[1], SpinAlpha(alpha_lower, param), M4Vec{}, true);
  auto       subchannel = spin::Subchannel(upper, lower, first.p, second.p, false, true);
  subchannel.Reshape(upper.helicity.Jz_values.size(), lower.helicity.Jz_values.size());
  return subchannel;
}

// Supporting GP spin mathematics and analytic vertices

// Compute the analytic trajectory m-basis matrix column index
// Analytic exchange helicity m runs from -MMAX to +MMAX
std::size_t AnalyticMIndex(int m, int MMAX, const std::string& context) {
  if (MMAX < 0) { throw std::invalid_argument(context + ": NUMERICS_REGGE.MMAX must be non-negative"); }
  if (MMAX > (std::numeric_limits<int>::max() - 1) / 2) {
    throw std::invalid_argument(context + ": NUMERICS_REGGE.MMAX is too large");
  }
  if (m < -MMAX || m > MMAX) {
    throw std::invalid_argument(context + ": analytic exchange helicity m=" + std::to_string(m) +
                                " exceeds MMAX=" + std::to_string(MMAX));
  }
  return static_cast<std::size_t>(m + MMAX);
}

// Precompute the pole-normalized central-resonance orbital SU(2) factors
void InitResonanceLS(HELMatrix& hc, const int JX, const int j1x2, const int j2x2) {
  if (!hc.UsesLSCouplings() || JX < 0 || hc.alpha_ls.Empty() || !std::isfinite(hc.analytic_Lambda) ||
      hc.analytic_Lambda <= 0.0) {
    throw std::invalid_argument("MReggeGP::InitResonanceLS requires an analytic LS resonance");
  }
  hc.J  = static_cast<double>(JX);
  hc.s1 = 0.5 * static_cast<double>(j1x2);
  hc.s2 = 0.5 * static_cast<double>(j2x2);

  GPOrbitalBasis basis;
  basis.nmu = static_cast<std::size_t>(2 * JX + 1);
  basis.terms.reserve(hc.alpha_ls.Size());
  basis.su2.reserve(hc.alpha_ls.Size() * basis.nmu);
  for (const spin::LSTerm& term : hc.alpha_ls) {
    const auto position = std::lower_bound(basis.spins.begin(), basis.spins.end(), term.two_s);
    if (position == basis.spins.end() || *position != term.two_s) { basis.spins.insert(position, term.two_s); }
  }

  for (const spin::LSTerm& term : hc.alpha_ls) {
    if (term.l > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        term.two_s > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::invalid_argument("MReggeGP::InitResonanceLS: LS quantum number is too large");
    }
    const auto        position   = std::lower_bound(basis.spins.begin(), basis.spins.end(), term.two_s);
    const std::size_t spin_index = static_cast<std::size_t>(position - basis.spins.begin());
    basis.terms.push_back({term.l, term.two_s, spin_index, basis.su2.size()});
    const double pole = spin::RawPoleLSNormalization(j1x2, j2x2, 2 * JX, term.l, static_cast<int>(term.two_s));
    const double norm = pole * msqrt((2.0 * static_cast<double>(term.l) + 1.0) / (2.0 * static_cast<double>(JX) + 1.0));
    const double S    = 0.5 * static_cast<double>(term.two_s);
    for (int mu = -JX; mu <= JX; ++mu) {
      const double su2 = norm * wigner::CG(static_cast<double>(term.l), S, 0.0, static_cast<double>(mu),
                                           static_cast<double>(JX), static_cast<double>(mu));
      if (!std::isfinite(su2)) {
        throw std::invalid_argument("MReggeGP::InitResonanceLS: non-finite orbital coefficient");
      }
      basis.su2.push_back(su2);
    }
  }
  hc.gp_orbital = std::move(basis);
}

namespace {

// Build pole normalizations for one sparse crossed LS coefficient block
GPOrbitalBasis CrossedLSBasis(const HELMatrix& hc, const spin::LSCoefficients& terms, const int JX) {
  const double j1x2 = 2.0 * hc.s1;
  const double j2x2 = 2.0 * hc.s2;
  if (!std::isfinite(j1x2) || !std::isfinite(j2x2) || std::abs(j1x2 - std::round(j1x2)) > 1.0e-12 ||
      std::abs(j2x2 - std::round(j2x2)) > 1.0e-12) {
    throw std::invalid_argument("MReggeGP::InitCrossedLS requires physical leg spins");
  }
  GPOrbitalBasis basis;
  basis.terms.reserve(terms.Size());
  basis.su2.reserve(terms.Size());
  for (const spin::LSTerm& term : terms) {
    basis.terms.push_back({term.l, term.two_s, 0, basis.su2.size()});
    basis.su2.push_back(spin::RawPoleLSNormalization(static_cast<int>(std::llround(j1x2)),
                                                     static_cast<int>(std::llround(j2x2)), 2 * JX, term.l,
                                                     static_cast<int>(term.two_s)));
  }
  return basis;
}

}  // namespace

// Precompute the pole normalization of continuous-J crossed LS vertices
void InitCrossedLS(HELMatrix& hc, const int JX) {
  if (!hc.UsesReggeDomain() || !hc.UsesLSCouplings() || JX < 0) {
    throw std::invalid_argument("MReggeGP::InitCrossedLS requires an analytic crossed LS residue");
  }
  hc.J = static_cast<double>(JX);
  if (hc.exchange_basis == ExchangeBasisType::ReggeHelicity) {
    const std::size_t nm = static_cast<std::size_t>(2 * hc.analytic_MMAX + 1);
    if (hc.m_ls.size() != nm || std::none_of(hc.m_ls.cbegin(), hc.m_ls.cend(),
                                             [](const spin::LSCoefficients& terms) { return !terms.Empty(); })) {
      throw std::invalid_argument("MReggeGP::InitCrossedLS requires nonempty per-m LS residues");
    }
    hc.m_orbital.assign(nm, {});
    for (const auto& m : indices(hc.m_ls)) {
      if (!hc.m_ls[m].Empty()) { hc.m_orbital[m] = CrossedLSBasis(hc, hc.m_ls[m], JX); }
    }
    hc.gp_orbital = {};
    return;
  }
  if (hc.exchange_basis != ExchangeBasisType::HelicityTransport || hc.analytic_MMAX != 1 || hc.alpha_ls.Empty()) {
    throw std::invalid_argument("MReggeGP::InitCrossedLS requires a physical photon or Regge m basis");
  }
  hc.gp_orbital = CrossedLSBasis(hc, hc.alpha_ls, JX);
  hc.m_orbital.clear();
}

// Compute one proton-side Regge-helicity vertex factor for an integer m label
// The unit-vertex option keeps azimuthal m phases without qT^|m|
// weights
std::complex<double> Residue(int m, double delta_lambda, double qt, double phi, double s0,
                             bool second_exchange_daughter, bool use_exchange_helicity_barrier) {
  return spin::HelicityFactor(m, delta_lambda, qt, phi, math::msqrt(s0), PDG::mp, second_exchange_daughter,
                              use_exchange_helicity_barrier);
}

// Compute the canonical pole spin for a photon or mapped analytic trajectory
int PoleSpinX2(const regge::Param& param, const MPDG& pdg_table, const int pdg) {
  if (pdg == PDG::PDG_gamma) { return 2; }
  return regge::PoleRepresentative(param, pdg_table, pdg).spinX2;
}

// Compute the anchor-free falling factorial for one exchange projection
// N_m(alpha) = product_{k=0}^{|m|-1}(alpha-k)
double NonsenseZero(double alpha_t, int m) {
  if (!std::isfinite(alpha_t)) { throw AmplitudeFailure("MReggeGP::NonsenseZero: non-finite alpha"); }
  const int helicity = std::abs(m);
  double    factor   = 1.0;
  for (int k = 0; k < helicity; ++k) { factor *= alpha_t - static_cast<double>(k); }
  return factor;
}

namespace {

// Store one coupled exchange-spin contribution
struct PreparedSpin {
  std::size_t spin_index     = 0;
  std::size_t orbital_offset = 0;
};

// Store the factorized LS sum for one resonance
struct PreparedLS {
  std::size_t               nmu = 0;
  std::vector<PreparedSpin> spins;
  HelVec                    orbital;
};

// Build an analytic Regge vertex as R = f r^T
spin::ForwardSourceFactors ReggeFactors(const std::vector<std::pair<double, double>>& rows,
                                        const std::vector<int>& m_values, const std::vector<double>& nonsense,
                                        const M4Vec& q, double s0, bool second_exchange_daughter,
                                        bool use_exchange_helicity_barrier, const std::string& context) {
  if (m_values.size() != nonsense.size()) {
    throw std::invalid_argument(context + ": analytic residue basis dimension mismatch");
  }

  const double           qt                        = q.Pt();
  const double           phi                       = q.Phi();
  const double           collider_transfer_azimuth = (second_exchange_daughter ? q : -q).Phi();
  spin::ForwardSourceFactors residue;
  residue.exchange_helicity.assign(m_values.size(), 0.0);
  for (const auto& col : indices(m_values)) {
    residue.exchange_helicity[col] =
        gpom::Residue(m_values[col], 0.0, qt, phi, s0, second_exchange_daughter, use_exchange_helicity_barrier) *
        nonsense[col];
  }

  residue.beam_transition.assign(rows.size(), 0.0);
  for (const auto& row : indices(rows)) {
    residue.beam_transition[row] = spin::SpinHalfForwardHelicitySectionFactor(
        rows[row].first, rows[row].second, collider_transfer_azimuth, second_exchange_daughter);
  }
  return residue;
}

}  // namespace

// Build the nonphoton analytic vertex as the outer product R = f r^T
// The SOFT Good Walker source supplies the proton flip magnitude
HelAmp ResidueMatrix(const std::vector<std::pair<double, double>>& rows, const std::vector<int>& m_values,
                     const std::vector<double>& nonsense, const M4Vec& q, double s0, bool second_exchange_daughter,
                     bool use_exchange_helicity_barrier, const std::string& context) {
  const auto residue =
      ReggeFactors(rows, m_values, nonsense, q, s0, second_exchange_daughter, use_exchange_helicity_barrier, context);
  return OuterProduct(residue.beam_transition, residue.exchange_helicity);
}

// Build the common analytic exchange helicity basis for GP[RES] and GP[CON]
ReggeBasis BuildBasis(const LORENTZSCALAR& lts, const regge::Param& param, int MMAX, int up_pdg, int dn_pdg,
                      const std::string& context) {
  (void)AnalyticMIndex(0, MMAX, context);

  ReggeBasis basis;
  basis.MMAX                          = MMAX;
  basis.nm                            = static_cast<std::size_t>(2 * MMAX + 1);
  basis.upper_photon                  = (up_pdg == PDG::PDG_gamma);
  basis.lower_photon                  = (dn_pdg == PDG::PDG_gamma);
  basis.use_exchange_helicity_barrier = (lts.process.FORWARD_VERTEX != ForwardVertexMode::UnitResidue);
  basis.upper_rows                    = spin::SpinHalfTransitions(lts.process.FORWARD_NOFLIP);
  basis.lower_rows                    = spin::SpinHalfTransitions(lts.process.FORWARD_NOFLIP);
  basis.m_values.reserve(basis.nm);
  for (int m = -MMAX; m <= MMAX; ++m) { basis.m_values.push_back(m); }

  basis.alpha1 = basis.upper_photon ? 1.0 : SpinAlpha(regge::Alpha(param, up_pdg, lts.t1), param);
  basis.alpha2 = basis.lower_photon ? 1.0 : SpinAlpha(regge::Alpha(param, dn_pdg, lts.t2), param);
  // Compensate the CGRegge Gamma branch under identical-leg exchange
  if (up_pdg == dn_pdg) {
    basis.exchange_section = std::exp(-2.0 * math::zi * math::PI * (basis.alpha1 - basis.alpha2));
  }
  basis.upper_nonsense.assign(basis.nm, 1.0);
  basis.lower_nonsense.assign(basis.nm, 1.0);

  for (int m = -MMAX; m <= MMAX; ++m) {
    const std::size_t col = static_cast<std::size_t>(m + MMAX);
    if (!basis.upper_photon) { basis.upper_nonsense[col] = NonsenseZero(basis.alpha1, m); }
    if (!basis.lower_photon) { basis.lower_nonsense[col] = NonsenseZero(basis.alpha2, m); }
  }

  return basis;
}

// Build upper and lower forward vertices with their respective Jacob-Wick sections
ForwardSource BuildSource(const LORENTZSCALAR& lts, const regge::Param& param, const SourceKey& key,
                          const char* context) {
  ForwardSource source;
  source.key = key;
  source.basis = BuildBasis(lts, param, key.MMAX, key.up_pdg, key.dn_pdg, context);
  const auto& basis = source.basis;
  const auto leg = [&](bool lower, HelAmp& matrix, HelVec& row, HelVec& column) {
    const auto& rows = lower ? basis.lower_rows : basis.upper_rows;
    if (lower ? basis.lower_photon : basis.upper_photon) {
      matrix = qed::PhotonSourceMatrixTransitions(lts, lower ? 2 : 1, rows, basis.m_values, lts.process.PHOTON_VERTEX, lower);
    } else {
      auto residue = ReggeFactors(rows, basis.m_values, lower ? basis.lower_nonsense : basis.upper_nonsense,
                                   lower ? lts.q2 : lts.q1, param.s0, lower, basis.use_exchange_helicity_barrier, context);
      row = std::move(residue.beam_transition);
      column = std::move(residue.exchange_helicity);
      if (basis.upper_photon || basis.lower_photon) { matrix = OuterProduct(row, column); }
    }
  };
  leg(false, source.upper_residue, source.upper_row_factor, source.upper_column_factor);
  leg(true, source.lower_residue, source.lower_row_factor, source.lower_column_factor);
  return source;
}

namespace {

// Compute the discrete reflection phase of one coupled physical-pole spin
int CoupledReflectionPhase(const long long j1x2, const long long j2x2, const long long sx2, const char* context) {
  const long long exponent_x2 = j1x2 + j2x2 - sx2;
  if (exponent_x2 % 2 != 0) { throw std::logic_error(std::string(context) + ": non-integral reflection phase"); }
  return (exponent_x2 / 2) % 2 == 0 ? 1 : -1;
}

// Locate or append a spin block with stable indices for this amplitude
std::size_t CacheReggeSpin(ForwardSource& source, std::size_t two_s, int reflection_phase) {
  for (const auto i : indices(source.regge_cg)) {
    const auto& block = source.regge_cg[i];
    if (block.two_s != two_s) { continue; }
    if (block.reflection_phase != reflection_phase) {
      throw std::logic_error("CacheReggeSpin: inconsistent pole-spin reflection");
    }
    return i;
  }
  source.regge_cg.push_back(
      {two_s, reflection_phase, std::vector<std::optional<std::complex<double>>>(source.basis.nm * source.basis.nm)});
  return source.regge_cg.size() - 1;
}

// Compute one pole-reflection projected Gamma-continued exchange-spin coefficient
// C = zeta_ex[CG(m1,-m2)+eta CG(-m1,m2)]/2
std::complex<double> ReggeSpin(ForwardSource& source, std::size_t spin_index, std::size_t pair_index, int m1, int m2) {
  ReggeCGBlock& block = source.regge_cg[spin_index];
  if (pair_index >= block.coefficient.size()) {
    throw std::logic_error("MReggeGP::ReggeSpin: pair index outside cache");
  }
  std::optional<std::complex<double>>& cached = block.coefficient[pair_index];
  if (cached.has_value()) { return *cached; }

  const std::size_t reflected_i1 = AnalyticMIndex(-m1, source.basis.MMAX, "MReggeGP::ReggeSpin reflected upper spin");
  const std::size_t reflected_i2 = AnalyticMIndex(-m2, source.basis.MMAX, "MReggeGP::ReggeSpin reflected lower spin");
  const std::size_t reflected_pair_index                = reflected_i1 * source.basis.nm + reflected_i2;
  std::optional<std::complex<double>>& reflected_cached = block.coefficient[reflected_pair_index];

  const std::size_t two_s  = block.two_s;
  const long long   mu     = static_cast<long long>(m1) - m2;
  const std::size_t abs_mu = static_cast<std::size_t>(std::abs(mu));
  if (abs_mu > two_s / 2 || (reflected_pair_index == pair_index && block.reflection_phase < 0)) {
    cached           = 0.0;
    reflected_cached = 0.0;
    return *cached;
  }

  const double               j1 = source.basis.upper_photon ? 1.0 : source.basis.alpha1;
  const double               j2 = source.basis.lower_photon ? 1.0 : source.basis.alpha2;
  const double               S  = 0.5 * static_cast<double>(two_s);
  const std::complex<double> direct =
      wigner::CGRegge(j1, j2, static_cast<double>(m1), -static_cast<double>(m2), S, static_cast<double>(mu));
  // The zero-helicity pair is its own reflection
  const std::complex<double> reflected =
      reflected_pair_index == pair_index
          ? direct
          : wigner::CGRegge(j1, j2, -static_cast<double>(m1), static_cast<double>(m2), S, -static_cast<double>(mu));
  cached = source.basis.exchange_section * 0.5 * (direct + static_cast<double>(block.reflection_phase) * reflected);
  if (reflected_pair_index != pair_index) { reflected_cached = static_cast<double>(block.reflection_phase) * *cached; }
  return *cached;
}

// Factor the LS sum into one orbital coefficient per exchange spin and mu
PreparedLS PrepareResonanceLS(const HELMatrix& hel, double momentum, bool derivative_factor, ForwardSource& source) {
  const GPOrbitalBasis& basis = hel.gp_orbital;

  const int  j1x2 = spin::SpinLabelX2(hel.s1, "MReggeGP::PrepareResonanceLS upper pole spin");
  const int  j2x2 = spin::SpinLabelX2(hel.s2, "MReggeGP::PrepareResonanceLS lower pole spin");
  PreparedLS prepared;
  prepared.nmu = basis.nmu;
  prepared.spins.reserve(basis.spins.size());
  prepared.orbital.assign(basis.spins.size() * prepared.nmu, 0.0);
  for (const auto i : indices(basis.spins)) {
    const auto two_s = basis.spins[i];
    const int  phase = CoupledReflectionPhase(j1x2, j2x2, two_s, "PrepareResonanceLS");
    prepared.spins.push_back({CacheReggeSpin(source, two_s, phase), i * prepared.nmu});
  }

  const double scaled_momentum = derivative_factor ? momentum / hel.analytic_Lambda : 1.0;
  if (!std::isfinite(scaled_momentum)) {
    throw AmplitudeFailure("MReggeGP::PrepareResonanceLS: non-finite derivative factor");
  }
  std::size_t term_index = 0;
  for (const spin::LSTerm& term : hel.alpha_ls) {
    const GPOrbitalTerm& orbital = basis.terms[term_index++];

    const double               derivative  = math::IntegerPower(scaled_momentum, static_cast<unsigned int>(orbital.l));
    const std::complex<double> coefficient = term.coefficient * derivative;
    const std::size_t          destination = orbital.spin_index * prepared.nmu;
    for (std::size_t mu = 0; mu < prepared.nmu; ++mu) {
      prepared.orbital[destination + mu] += coefficient * basis.su2[orbital.su2_offset + mu];
    }
  }
  return prepared;
}

// Compute one Gamma continued resonance LS vertex
std::complex<double> ResonanceLS(const PreparedLS& ls, std::size_t mu_index, int m1, int m2, std::size_t pair_index,
                                 ForwardSource& source) {
  std::complex<double> out = 0.0;
  for (const PreparedSpin& spin : ls.spins) {
    const auto orbital = ls.orbital[spin.orbital_offset + mu_index];
    if (math::IsZero(orbital)) { continue; }
    out += orbital * ReggeSpin(source, spin.spin_index, pair_index, m1, m2);
  }
  return out;
}

// Traverse the allowed fusion terms with one produced spin column per exchange pair
template <typename Action>
void ForEachFusion(const HELMatrix& hel, ForwardSource& source, double momentum, bool derivative_factor,
                   Action&& action) {
  const int        JX    = static_cast<int>(std::llround(hel.J));
  const auto&      basis = source.basis;
  const PreparedLS ls =
      hel.UsesLSCouplings() ? PrepareResonanceLS(hel, momentum, derivative_factor, source) : PreparedLS{};
  if (hel.UsesLSCouplings()) {
    for (int m1 = -basis.MMAX; m1 <= basis.MMAX; ++m1) {
      if (basis.upper_photon && std::abs(m1) != 1) { continue; }
      const std::size_t i1 = static_cast<std::size_t>(m1 + basis.MMAX);
      for (int m2 = -basis.MMAX; m2 <= basis.MMAX; ++m2) {
        if (basis.lower_photon && std::abs(m2) != 1) { continue; }
        const int mu = m1 - m2;
        if (std::abs(mu) > JX) { continue; }
        const std::size_t i2         = static_cast<std::size_t>(m2 + basis.MMAX);
        const std::size_t pair_index = i1 * basis.nm + i2;
        const std::size_t mu_index   = static_cast<std::size_t>(mu + JX);
        const double metric = spin::JacobWickSecondLegReversalPhase(static_cast<double>(JX), -static_cast<double>(mu));
        action(i1, i2, static_cast<std::size_t>(JX - mu),
               metric * ResonanceLS(ls, mu_index, m1, m2, pair_index, source));
      }
    }
  } else {
    for (const auto& [i1, i2] : hel.T_active) {
      const int    m1     = static_cast<int>(i1) - basis.MMAX;
      const int    m2     = static_cast<int>(i2) - basis.MMAX;
      const int    mu     = m1 - m2;
      const double metric = spin::JacobWickSecondLegReversalPhase(static_cast<double>(JX), -static_cast<double>(mu));
      action(i1, i2, static_cast<std::size_t>(JX - mu), metric * hel.T[i1][i2]);
    }
  }
}

}  // namespace

namespace {

// Evaluate one continuous-J crossed LS coefficient block
std::complex<double> CrossedLS(const HELMatrix& hel, const spin::LSCoefficients& terms, const GPOrbitalBasis& basis,
                               double alpha, double lambda1, double lambda2) {
  const double lambda = lambda1 - lambda2;
  if (std::fpclassify(2.0 * alpha + 1.0) == FP_ZERO) { throw AmplitudeFailure("MReggeGP::Continuum singular Regge spin norm"); }
  const int            Jx2        = spin::SpinLabelX2(hel.J, "MReggeGP::CrossedLS trajectory pole spin");
  std::complex<double> out        = 0.0;
  std::size_t          term_index = 0;
  for (const spin::LSTerm& term : terms) {
    const GPOrbitalTerm&       prepared = basis.terms[term_index++];
    const double               S        = 0.5 * static_cast<double>(term.two_s);
    const double               spin     = wigner::CG(hel.s1, hel.s2, lambda1, -lambda2, S, lambda);
    const std::complex<double> direct_orbital =
        wigner::CGRegge(static_cast<double>(term.l), S, 0.0, lambda, alpha, lambda);
    const std::complex<double> reflected_orbital =
        wigner::CGRegge(static_cast<double>(term.l), S, 0.0, -lambda, alpha, -lambda);
    const int                  reflection_phase = CoupledReflectionPhase(2LL * static_cast<long long>(term.l),
                                                                         static_cast<long long>(term.two_s), Jx2, "MReggeGP::CrossedLS");
    const std::complex<double> orbital =
        0.5 * (direct_orbital + static_cast<double>(reflection_phase) * reflected_orbital);
    const std::complex<double> norm =
        std::sqrt(std::complex<double>((2.0 * static_cast<double>(term.l) + 1.0) / (2.0 * alpha + 1.0), 0.0));
    out += term.coefficient * basis.su2[prepared.su2_offset] * norm * spin * orbital;
  }
  return out;
}

}  // namespace

// Convert one analytic photon input cache to a physical pole helicity tensor
HELMatrix PhotonPoleHelicity(const HELMatrix& input) {
  const std::size_t n1  = spin::SpinStateCount(input.s1, "MReggeGP::PhotonPoleHelicity first leg");
  const std::size_t n2  = spin::SpinStateCount(input.s2, "MReggeGP::PhotonPoleHelicity second leg");
  HELMatrix         out = input;
  out.domain            = HelicityDomain::PhysicalPole;
  out.coupling_basis    = CouplingBasis::Helicity;
  out.analytic_MMAX     = -1;
  out.J                 = 1.0;
  out.Jz_values         = {-1.0, 1.0};
  out.T                 = HelAmp(n1, n2, 0.0);
  out.T_set             = MMatrix<bool>(n1, n2, false);
  out.T_active.clear();
  for (std::size_t row = 0; row < input.lambda_values.size_row(); ++row) {
    const double         lambda1 = input.lambda_values[row][0];
    const double         lambda2 = input.lambda_values[row][1];
    std::complex<double> reduced = 0.0;
    bool                 active  = input.UsesLSCouplings();
    if (input.UsesLSCouplings()) {
      reduced = CrossedLS(input, input.alpha_ls, input.gp_orbital, 1.0, lambda1, lambda2);
    } else {
      for (const int m : {-1, 1}) {
        const std::size_t col = static_cast<std::size_t>(m + 1);
        if (!input.T_set[row][col]) { continue; }
        if (!active) {
          reduced = input.T[row][col];
          active  = true;
        }
      }
    }
    const std::size_t i1 = input.lambda_idx[row][0];
    const std::size_t i2 = input.lambda_idx[row][1];
    out.T[i1][i2]        = reduced;
    out.T_set[i1][i2]    = active;
  }
  return out;
}

namespace {

// Evaluate one physical photon subvertex with the common pole-frame sewing
spin::EvaluatedPoleSubvertex PhotonSubvertex(const HELMatrix& input, const M4Vec& final_in_X,
                                             const M4Vec& parent_dir_in_X, const bool second_exchange_daughter) {
  auto helicity = PhotonPoleHelicity(input);
  return spin::EvaluatePoleSubvertex(std::move(helicity), final_in_X, parent_dir_in_X, second_exchange_daughter);
}

// Embed helicity columns into the common exchange projection range
spin::EvaluatedPoleSubvertex EmbedSubvertex(spin::EvaluatedPoleSubvertex input, const int MMAX) {
  if (MMAX < 0 || input.frame.size_col() != input.helicity.Jz_values.size()) {
    throw std::invalid_argument("MReggeGP::EmbedSubvertex requires matching helicity columns");
  }
  const std::size_t nm = static_cast<std::size_t>(2 * MMAX + 1);
  HelAmp            frame(input.frame.size_row(), nm, 0.0);
  for (const auto& column : indices(input.helicity.Jz_values)) {
    const int         m      = static_cast<int>(std::llround(input.helicity.Jz_values[column]));
    const std::size_t target = AnalyticMIndex(m, MMAX, "MReggeGP::EmbedSubvertex");
    for (std::size_t row = 0; row < input.frame.size_row(); ++row) { frame[row][target] = input.frame[row][column]; }
  }
  input.frame                   = std::move(frame);
  input.helicity.exchange_basis = ExchangeBasisType::ReggeHelicity;
  input.helicity.analytic_MMAX  = MMAX;
  input.helicity.Jz_values.clear();
  input.helicity.Jz_values.reserve(nm);
  for (int m = -MMAX; m <= MMAX; ++m) { input.helicity.Jz_values.push_back(static_cast<double>(m)); }
  return input;
}

}  // namespace

// Evaluate one angle-free crossed Regge-helicity vertex
spin::EvaluatedPoleSubvertex Crossed(const HELMatrix& input, const double alpha, const M4Vec& axis,
                                     const bool second_exchange_daughter) {
  spin::EvaluatedPoleSubvertex out;
  out.helicity = input;
  if (input.UsesLSCouplings()) {
    const std::size_t nm = static_cast<std::size_t>(2 * input.analytic_MMAX + 1);
    out.helicity.T       = HelAmp(input.lambda_values.size_row(), nm, 0.0);
    out.helicity.T_set   = MMatrix<bool>(input.lambda_values.size_row(), nm, false);
    out.helicity.T_active.clear();
    for (std::size_t row = 0; row < input.lambda_values.size_row(); ++row) {
      const double lambda1 = input.lambda_values[row][0];
      const double lambda2 = input.lambda_values[row][1];
      for (const auto& m : indices(input.m_ls)) {
        if (input.m_ls[m].Empty()) { continue; }
        out.helicity.T[row][m]     = CrossedLS(input, input.m_ls[m], input.m_orbital[m], alpha, lambda1, lambda2);
        out.helicity.T_set[row][m] = true;
      }
    }
  }
  out.frame = spin::ReggeCrossedFrame(out.helicity, axis, second_exchange_daughter);
  return out;
}

namespace {

// Compute the common crossed pole density, summing independent Regge m vertices
double ContinuumScale(const HELMatrix& hel, bool photon) {
  if (photon) {
    const auto pole = PhotonPoleHelicity(hel);
    HelAmp     reduced(pole.lambda_values.size_row(), 1, 0.0);
    for (std::size_t row = 0; row < reduced.size_row(); ++row) {
      reduced[row][0] = pole.T[pole.lambda_idx[row][0]][pole.lambda_idx[row][1]];
    }
    return std::sqrt(spin::ForwardHelicityDensity(reduced, pole.lambda_values, spin::ForwardLegType::Hadron,
                                                  spin::ForwardLegType::Hadron));
  }
  const auto pole = Crossed(hel, hel.J, M4Vec{}, false);
  return std::sqrt(spin::ForwardHelicityDensity(pole.frame, pole.helicity.lambda_values, spin::ForwardLegType::Hadron,
                                                spin::ForwardLegType::Hadron));
}

}  // namespace

// Evaluate a photon or Regge vertex in the common exchange projection range
spin::EvaluatedPoleSubvertex CrossedSubvertex(const HELMatrix& vertex, int exchange_pdg, const M4Vec& final_in_X,
                                              const M4Vec& axis, bool lower, double alpha, int MMAX) {
  if (exchange_pdg == PDG::PDG_gamma) { return EmbedSubvertex(PhotonSubvertex(vertex, final_in_X, axis, lower), MMAX); }
  auto out = Crossed(vertex, alpha, axis, lower);
  return vertex.analytic_MMAX == MMAX ? std::move(out) : EmbedSubvertex(std::move(out), MMAX);
}

// Compute the full pole tensor for physical reference densities and photon couplings
HelAmp Fusion(const HELMatrix& hel, ForwardSource& source, double momentum, bool derivative_factor) {
  const std::size_t nmu = static_cast<std::size_t>(2 * std::llround(hel.J) + 1);
  HelAmp            central(source.basis.nm * source.basis.nm, nmu, 0.0);
  ForEachFusion(hel, source, momentum, derivative_factor,
                [&](std::size_t i1, std::size_t i2, std::size_t mu, std::complex<double> value) {
                  central[i1 * source.basis.nm + i2][mu] = value;
                });
  return central;
}

// Share a source only within the supplied amplitude evaluation
ForwardSource& detail::SourceFor(const LORENTZSCALAR& lts, const regge::Param& param, SourceKey key, AmpCache& cache) {
  for (auto& source : cache.sources) {
    if (source.key == key) { return source; }
  }
  cache.sources.push_back(BuildSource(lts, param, key, "GP source"));
  return cache.sources.back();
}

// Construct the physical-pole exchange basis used by normalization and photon transitions
ForwardSource detail::PoleSource(const HELMatrix& hel, bool upper_photon, bool lower_photon) {
  ForwardSource source;
  auto& basis = source.basis;
  basis.MMAX = hel.analytic_MMAX;
  basis.nm = static_cast<std::size_t>(2 * basis.MMAX + 1);
  basis.alpha1 = hel.s1;
  basis.alpha2 = hel.s2;
  basis.upper_photon = upper_photon;
  basis.lower_photon = lower_photon;
  return source;
}

// Compute the common physical pole density when production spin is disabled
double detail::ResonanceScale(const HELMatrix& hel, const std::vector<MDecayBranch>& tree, double momentum,
                      bool derivative_factor) {
  auto source = detail::PoleSource(hel, tree[0].p.pdg == PDG::PDG_gamma, tree[1].p.pdg == PDG::PDG_gamma);
  const auto& basis = source.basis;
  const auto      central = Fusion(hel, source, momentum, derivative_factor);
  MMatrix<double> helicities(basis.nm * basis.nm, 2, 0.0);
  for (int m1 = -basis.MMAX; m1 <= basis.MMAX; ++m1) {
    for (int m2 = -basis.MMAX; m2 <= basis.MMAX; ++m2) {
      const std::size_t row =
          static_cast<std::size_t>(m1 + basis.MMAX) * basis.nm + static_cast<std::size_t>(m2 + basis.MMAX);
      helicities[row][0] = static_cast<double>(m1);
      helicities[row][1] = static_cast<double>(m2);
    }
  }
  return std::sqrt(spin::ForwardHelicityDensity(
      central, helicities, basis.upper_photon ? spin::ForwardLegType::RealPhoton : spin::ForwardLegType::Hadron,
      basis.lower_photon ? spin::ForwardLegType::RealPhoton : spin::ForwardLegType::Hadron));
}

// Contract fusion terms directly without allocating an exchange-pair tensor
HelAmp detail::ContractFusion(const HELMatrix& hel, ForwardSource& source, double momentum, bool derivative_factor) {
  const std::size_t nmu         = static_cast<std::size_t>(2 * std::llround(hel.J) + 1);
  const auto        destination = spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(
             source.basis.upper_rows.size(), source.basis.lower_rows.size());
  if (!source.basis.upper_photon && !source.basis.lower_photon) {
    HelVec spin(nmu, 0.0);
    ForEachFusion(hel, source, momentum, derivative_factor,
                  [&](std::size_t i1, std::size_t i2, std::size_t mu, std::complex<double> value) {
                    spin[mu] += source.upper_column_factor[i1] * source.lower_column_factor[i2] * value;
                  });
    return OuterProduct(MappedKroneckerProduct(source.upper_row_factor, source.lower_row_factor, destination), spin);
  }
  const auto upper = source.upper_residue.Transpose();
  const auto lower = source.lower_residue.Transpose();
  HelAmp     out(destination.size(), nmu, 0.0);
  ForEachFusion(hel, source, momentum, derivative_factor,
                [&](std::size_t i1, std::size_t i2, std::size_t mu, std::complex<double> value) {
                  out.AddColumnKroneckerProduct(mu, upper.Row(i1), lower.Row(i2), value);
                });
  return out.PermuteRows(destination);
}

// Build analytic t/u continuum production matrices from crossed vertices
HelPair detail::ConProd(const LORENTZSCALAR& lts, const regge::Param& param, const std::vector<HELMatrix>& vertex,
                const std::vector<int>& exchange_channel, const std::array<M4Vec, 2>& final_in_X, AmpCache& cache) {
  const std::size_t n_final = spin::FinalStateHelicityCount(lts.decaytree[0].p, "MReggeGP::Continuum left") *
                              spin::FinalStateHelicityCount(lts.decaytree[1].p, "MReggeGP::Continuum right");
  const std::size_t proton_rows  = spin::SpinHalfTransitions(lts.process.FORWARD_NOFLIP).size();
  const bool        upper_photon = exchange_channel[0] == PDG::PDG_gamma;
  const bool        lower_photon = exchange_channel[1] == PDG::PDG_gamma;
  if (!lts.process.SPINGEN) {
    const double t_scale = ContinuumScale(vertex[0], upper_photon) * ContinuumScale(vertex[1], lower_photon);
    const double u_scale = ContinuumScale(vertex[2], upper_photon) * ContinuumScale(vertex[3], lower_photon);
    return {spin::Blind(proton_rows * proton_rows, n_final, t_scale),
            spin::Blind(proton_rows * proton_rows, n_final, u_scale)};
  }

  int common_mmax = 0;
  for (const auto& hel : vertex) { common_mmax = std::max(common_mmax, hel.analytic_MMAX); }
  const auto&  source = SourceFor(lts, param, {common_mmax, exchange_channel[0], exchange_channel[1]}, cache);
  const double alpha1 = source.basis.alpha1, alpha2 = source.basis.alpha2;
  // Match nonphoton m sections to the collider axes of their forward sources
  // Physical photon spin transport retains its central rest frame axes
  const M4Vec& upper_axis = upper_photon ? lts.q1_in_X : lts.q1;
  M4Vec        lower_axis = lower_photon ? lts.q2_in_X : lts.q2;
  lower_axis.Flip3();
  const std::array<spin::EvaluatedPoleSubvertex, 4> sub = {
      CrossedSubvertex(vertex[0], exchange_channel[0], final_in_X[0], upper_axis, false, alpha1, common_mmax),
      CrossedSubvertex(vertex[1], exchange_channel[1], final_in_X[1], lower_axis, true, alpha2, common_mmax),
      CrossedSubvertex(vertex[2], exchange_channel[0], final_in_X[1], upper_axis, false, alpha1, common_mmax),
      CrossedSubvertex(vertex[3], exchange_channel[1], final_in_X[0], lower_axis, true, alpha2, common_mmax)};
  // Sew the t/u vertices after projecting onto the current forward sources
  const auto project = [&](const auto& upper, const auto& lower) {
    return std::pair{
        spin::ProjectedSubchannel(sub[0], sub[1], upper, lower, lts.decaytree[0].p, lts.decaytree[1].p, false, true),
        spin::ProjectedSubchannel(sub[2], sub[3], upper, lower, lts.decaytree[0].p, lts.decaytree[1].p, true, true)};
  };
  if (upper_photon || lower_photon) { return project(source.upper_residue, source.lower_residue); }
  const auto subchannels = project(source.upper_column_factor, source.lower_column_factor);
  const auto destination = spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(
      source.upper_row_factor.size(), source.lower_row_factor.size());
  const auto beam = MappedKroneckerProduct(source.upper_row_factor, source.lower_row_factor, destination);
  return {OuterProduct(beam, subchannels.first), OuterProduct(beam, subchannels.second)};
}

}  // namespace gpom

}  // namespace gra
