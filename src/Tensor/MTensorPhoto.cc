// Charged meson pair production: pi+ pi- and K+ K- in the Tensor Pomeron model
// 
// Coherent gamma-Pomeron and gamma-Reggeon Drell-Soding amplitudes from both
// beam directions. Optional gamma-gamma fusion includes scalar-QED
// t- and u-channel meson exchange and the contact term
// 
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#include "Graniitti/Tensor/MTensorPhoto.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Spin/MHelicityBasis.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra {
namespace {

using Current = MTensorPhoto::Current;

// Compute one checked immutable Tensor parameter block
const MTensorPomeronParam& RequirePhotoParameters(const MTensorPomeronParamPtr& parameters) {
  if (parameters == nullptr || !parameters->initialized) {
    throw std::invalid_argument("MTensorPhoto: immutable Tensor parameters are missing");
  }
  return *parameters;
}

// Store one oriented gamma-star proton event and unequal Regge subenergies
struct PhotoKinematics {
  int                   pdg = 0;
  M4Vec                 q;
  M4Vec                 p;
  M4Vec                 prime;
  M4Vec                 kplus;
  M4Vec                 kminus;
  M4Vec                 psum;
  M4Vec                 d;
  M4Vec                 rplus;
  M4Vec                 rminus;
  double                t      = 0.0;
  double                nubar2 = 0.0;
  double                rp_p   = 0.0;
  double                rm_p   = 0.0;
  double                d_p    = 0.0;
  double                kappa  = 0.0;
  std::array<double, 2> log_nu{};
  bool                  valid = false;
};

// Compute the direct charged-meson momenta in fixed charge order
std::pair<M4Vec, M4Vec> ChargedMesons(const LORENTZSCALAR& lts) {
  const auto positive = lts.decaytree[0].p.pdg > 0 ? 0 : 1;
  return {lts.decaytree[positive].p4, lts.decaytree[1 - positive].p4};
}

// Resolve one photon direction and its positive unequal subenergies
PhotoKinematics PhotoEvent(const LORENTZSCALAR& lts, const bool photon_from_upper, const ForwardLegState& target) {
  PhotoKinematics out;
  out.pdg                    = std::abs(lts.decaytree[0].p.pdg);
  const auto [kplus, kminus] = ChargedMesons(lts);
  out.kplus                  = kplus;
  out.kminus                 = kminus;
  out.q                      = photon_from_upper ? lts.q1 : lts.q2;
  out.p                      = target.incoming;
  out.prime                  = target.outgoing;
  out.psum                   = out.prime + out.p;
  out.d                      = out.kminus - out.kplus;
  out.rplus                  = out.d + out.q;
  out.rminus                 = out.d - out.q;
  out.t                      = (out.prime - out.p).M2();

  // Use the equal-proton-mass identity without subtracting large invariants
  // [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:2508.06334v2, Eq. (A.6)]
  const auto   p   = out.psum.Contravariant<long double>();
  const double nu1 = static_cast<double>(0.5L * gra::MinkowskiProduct(p, out.kplus.Contravariant<long double>()));
  const double nu2 = static_cast<double>(0.5L * gra::MinkowskiProduct(p, out.kminus.Contravariant<long double>()));
  
  // Preserve small contractions between nearly collinear high energy momenta
  out.rp_p         = static_cast<double>(gra::MinkowskiProduct(p, out.rplus.Contravariant<long double>()));
  out.rm_p         = static_cast<double>(gra::MinkowskiProduct(p, out.rminus.Contravariant<long double>()));
  out.d_p          = static_cast<double>(gra::MinkowskiProduct(p, out.d.Contravariant<long double>()));
  const double sum = math::pow2(nu1) + math::pow2(nu2);
  if (!(nu1 > 0.0) || !(nu2 > 0.0) || !(sum > 0.0) || !std::isfinite(sum) || !std::isfinite(out.t)) { return out; }
  out.nubar2             = 0.5 * sum;
  out.kappa              = (math::pow2(nu1) - math::pow2(nu2)) / sum;
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
  if (std::abs(out.kappa) > 1.0 + tolerance) { return out; }
  out.kappa = std::clamp(out.kappa, -1.0, 1.0);

  // Preserve both positive subenergies even when kappa rounds to unity
  if (std::abs(out.kappa) < 0.5) {
    out.log_nu = {std::log1p(out.kappa), std::log1p(-out.kappa)};
  } else {
    const double nubar = std::sqrt(out.nubar2);
    out.log_nu         = {2.0 * std::log(nu1 / nubar), 2.0 * std::log(nu2 / nubar)};
  }
  out.valid = true;
  return out;
}

// Compute the monopole meson form factor used in the paper amplitudes
// F_M(q2) = m0^2/(m0^2-q2)
double MesonFF(const double q2, const double m0_2) {
  const double denominator = m0_2 - q2;
  if (!(denominator > 0.0) || !std::isfinite(denominator)) { return 0.0; }
  return m0_2 / denominator;
}

// Evaluate the unequal-subenergy continuation stably
// g = [(1-kappa)^(-lambda)-1]/(lambda kappa)
double ReggeG(const double lambda, const double kappa, const double log_nu) {
  constexpr double small = 1.0e-6;
  if (std::abs(kappa) < small && std::abs(lambda * kappa) < small) {
    return 1.0 + 0.5 * (lambda + 1.0) * kappa + (lambda + 1.0) * (lambda + 2.0) * math::pow2(kappa) / 6.0;
  }
  const double exponent = -lambda * log_nu;
  if (std::abs(exponent) < small) {
    return -log_nu / kappa * (1.0 + exponent * (0.5 + exponent * (1.0 / 6.0 + exponent / 24.0)));
  }
  return std::expm1(exponent) / (lambda * kappa);
}

// Compute the lower-index target Dirac current and scalar bilinear
std::pair<Current, std::complex<double>> TargetBilinears(const MTensorPomeron& tensor, const PhotoKinematics& event,
                                                         const ForwardLegState& state, const bool noflip,
                                                         const std::size_t initial_helicity,
                                                         const std::size_t final_helicity) {
  Current current{};
  if (noflip) {
    if (initial_helicity != final_helicity) { return {current, 0.0}; }
    for (const auto& mu : indices(current)) { current[mu] = event.psum % mu; }
    return {current, 2.0 * PDG::mp};
  }

  const auto   incoming   = tensor.SpinorStates(event.p, "u");
  const auto   outgoing   = tensor.SpinorStates(event.prime, "ubar");
  const auto&  u          = incoming[initial_helicity];
  const auto&  ubar       = outgoing[final_helicity];
  const double lambda_in  = static_cast<double>(spin::BinaryHelicityLabelX2(initial_helicity)) / 2.0;
  const double lambda_out = static_cast<double>(spin::BinaryHelicityLabelX2(final_helicity)) / 2.0;
  const M4Vec  transfer =
      state.leg == ForwardBeamLeg::Upper ? state.outgoing - state.incoming : state.incoming - state.outgoing;
  const std::complex<double> phase = MDirac::ElasticSpinHalfColliderCurrentPhase(
      MDirac::FermionKind::Particle, state.Index(), lambda_in, lambda_out, transfer.Phi());
  for (const auto& mu : indices(current)) { current[mu] = phase * tensor.gamma_lo[mu].BilinearForm(ubar, u); }
  return {current, phase * tensor.I4.BilinearForm(ubar, u)};
}

// Contract one contravariant vector with a lower-index complex current
std::complex<double> Contract(const M4Vec& vector, const Current& current) {
  const std::array<std::complex<long double>, 4> extended = {current[0], current[1], current[2], current[3]};
  return static_cast<std::complex<double>>(gra::BilinearProduct(vector.Contravariant<long double>(), extended));
}

// Compute the dressed photon-meson external-leg current
Current MesonPoleCurrent(const M4Vec& meson, const M4Vec& q, const double m0_2) {
  Current                    out{};
  const double               q2 = q.M2();
  const double               fm = MesonFF(q2, m0_2);
  const std::complex<double> denominator(-2.0 * (meson * q) + q2, 1.0e-15);
  if (!(std::abs(denominator) > 0.0)) { return out; }
  const M4Vec numerator = meson * 2.0 - q;
  for (const auto& mu : indices(out)) { out[mu] = fm * (numerator % mu) / denominator + (q % mu) / (m0_2 - q2); }
  return out;
}

// Compute the common rank-two meson-proton Regge amplitude factor at 2 nubar
std::complex<double> RankTwoResidue(const MTensorPomeron& tensor, const MTensorPomeronParam& parameter,
                                    const PhotoKinematics& event, const int exchange_pdg) {
  const auto& meson  = parameter.exchange.FindVertex(exchange_pdg, event.pdg);
  const auto& proton = parameter.exchange.FindVertex(exchange_pdg, PDG::PDG_p);
  if (meson.active_g_tensor.empty() || proton.active_g_tensor.empty()) { return 0.0; }
  const double gm    = parameter.exchange.PseudoscalarCoupling(exchange_pdg, meson.g_tensor.front());
  const double gp    = parameter.exchange.BaryonCoupling(exchange_pdg, proton.g_tensor.front());
  const double nubar = std::sqrt(event.nubar2);
  return 6.0 * gm * gp * regge::TransferFF(event.t, proton.ff_transfer) *
         tensor.TensorPropagatorFactor(exchange_pdg, 2.0 * nubar, event.t);
}

// Compute the common vector-Reggeon meson-proton amplitude factor at 2 nubar
std::complex<double> RankOneResidue(const MTensorPomeronParam& parameter, const PhotoKinematics& event,
                                    const int exchange_pdg) {
  const auto& meson  = parameter.exchange.FindVertex(exchange_pdg, event.pdg);
  const auto& proton = parameter.exchange.FindVertex(exchange_pdg, PDG::PDG_p);
  if (meson.active_g_tensor.empty() || proton.active_g_tensor.empty()) { return 0.0; }
  const double               gm = parameter.exchange.VectorPseudoscalarCoupling(exchange_pdg, meson.g_tensor.front());
  const double               gp = parameter.exchange.VectorBaryonCoupling(exchange_pdg, proton.g_tensor.front());
  const double               nubar      = std::sqrt(event.nubar2);
  const std::complex<double> propagator = parameter.exchange.PropagatorFactor(exchange_pdg, 2.0 * nubar, event.t);
  return -math::zi * gm * gp * regge::TransferFF(event.t, proton.ff_transfer) * propagator;
}

// Compute the full strong meson factors F_i and (F_1-F_2)/(t_pi1-t_pi2)
// The original F_M cancels in F_M G_i = F_i, including the contact term
// [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:2609.01285, Eqs. (2.2), (2.7), (2.15)-(2.16)]
std::array<double, 3> OffshellFactors(const PhotoKinematics& event, const MTensorPomeronParam& parameter) {
  const auto& [transfer, offshell] = parameter.photo.offshell.at(event.pdg);
  // Supply virtuality minus the pion pole directly to preserve collinear precision
  auto factors = regge::MassFFPair(event.q.M2() - 2.0 * (event.q * event.kplus),
                                   event.q.M2() - 2.0 * (event.q * event.kminus), 0.0, offshell);
  const double form = regge::TransferFF(event.t, transfer);
  for (auto& factor : factors) { factor *= form; }
  return factors;
}

// Add one gauge-restored C-even rank-two Drell-Soding exchange
// [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:2508.06334v2, Eqs. (2.21)-(2.23)]
void AddRankTwo(Current& out, const MTensorPomeron& tensor, const MTensorPomeronParam& parameter,
                const PhotoKinematics& event, const Current& target, const std::complex<double> scalar,
                const int exchange_pdg, const std::array<double, 3>& f) {
  const auto&                exchange = parameter.exchange.FindExchange(exchange_pdg);
  const double               alpha    = 1.0 + exchange.delta + exchange.ap * event.t;
  const double               lambda   = 0.5 * (2.0 - alpha);
  const double               gp       = ReggeG(lambda, event.kappa, event.log_nu[1]);
  const double               gm       = ReggeG(lambda, -event.kappa, event.log_nu[0]);
  const double               wa       = std::exp(-lambda * event.log_nu[1]);
  const double               wb       = std::exp(-lambda * event.log_nu[0]);
  const std::complex<double> residue  = RankTwoResidue(tensor, parameter, event, exchange_pdg);
  if (std::fpclassify(std::abs(residue)) == FP_ZERO) { return; }

  const double               e        = qed::e_QED();
  const Current              ha       = MesonPoleCurrent(event.kplus, event.q, parameter.photo.m0_2);
  const Current              hb       = MesonPoleCurrent(event.kminus, event.q, parameter.photo.m0_2);
  const std::complex<double> jrplus   = Contract(event.rplus, target);
  const std::complex<double> jrminus  = Contract(event.rminus, target);
  const std::complex<double> strong_a = 2.0 * event.rp_p * jrplus - event.rplus.M2() * PDG::mp * scalar;
  const std::complex<double> strong_b = 2.0 * event.rm_p * jrminus - event.rminus.M2() * PDG::mp * scalar;
  const std::complex<double> idot     = math::zi * residue;
  const double               pdiff    = -event.d_p;
  const double               dpsum    = event.d_p;
  const std::complex<double> dj       = Contract(event.d, target);
  const double               skew     = pdiff / (16.0 * event.nubar2);

  for (const auto& mu : indices(out)) {
    out[mu] += e * idot * (f[0] * ha[mu] * wa * strong_a - f[1] * hb[mu] * wb * strong_b
                           + (event.d % mu) * f[2] * (wa * strong_a + wb * strong_b));
    const std::complex<double> vector_contact =
        2.0 * dpsum * target[mu] + 2.0 * (event.psum % mu) * dj +
        (event.psum % mu) * (2.0 - alpha) * skew * (gp * jrplus * event.rp_p + gm * jrminus * event.rm_p);
    const double scalar_contact = -2.0 * (event.d % mu) - (event.psum % mu) * 0.5 * (2.0 - alpha) * skew *
                                                              (gp * event.rplus.M2() + gm * event.rminus.M2());
    out[mu] += (f[0] + f[1]) * e * idot * (vector_contact + scalar_contact * PDG::mp * scalar);
  }
}

// Add one gauge-restored C-odd vector-Reggeon Drell-Soding exchange
// [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:2508.06334v2, Eqs. (2.47)-(2.49)]
void AddRankOne(Current& out, const MTensorPomeronParam& parameter,
                const PhotoKinematics& event, const Current& target, const int exchange_pdg, const std::array<double, 3>& f) {
  const auto&                exchange = parameter.exchange.FindExchange(exchange_pdg);
  const double               alpha    = 1.0 + exchange.delta + exchange.ap * event.t;
  const double               lambda   = 0.5 * (1.0 - alpha);
  const double               gp       = ReggeG(lambda, event.kappa, event.log_nu[1]);
  const double               gm       = ReggeG(lambda, -event.kappa, event.log_nu[0]);
  const double               wa       = std::exp(-lambda * event.log_nu[1]);
  const double               wb       = std::exp(-lambda * event.log_nu[0]);
  const std::complex<double> residue  = RankOneResidue(parameter, event, exchange_pdg);
  if (std::fpclassify(std::abs(residue)) == FP_ZERO) { return; }

  const double               e       = qed::e_QED();
  const Current              ha      = MesonPoleCurrent(event.kplus, event.q, parameter.photo.m0_2);
  const Current              hb      = MesonPoleCurrent(event.kminus, event.q, parameter.photo.m0_2);
  const std::complex<double> jrplus  = Contract(event.rplus, target);
  const std::complex<double> jrminus = Contract(event.rminus, target);
  const double               pdiff   = -event.d_p;
  const double               skew    = pdiff / (16.0 * event.nubar2);
  for (const auto& mu : indices(out)) {
    out[mu] += e * residue * (f[0] * ha[mu] * wa * jrplus - f[1] * hb[mu] * wb * jrminus
                             + (event.d % mu) * f[2] * (wa * jrplus + wb * jrminus));
    out[mu] +=
        0.5 * (f[0] + f[1]) * e * residue * (2.0 * target[mu] + (event.psum % mu) * (1.0 - alpha) * skew * (gp * jrplus + gm * jrminus));
  }
}

// Add the form factor dressed scalar-QED meson continuum
// [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:2508.06334v2, Eqs. (2.42)-(2.45)]
void AddPhotonExchange(Current& out, const MTensorPomeron& tensor, const MTensorPomeronParam& parameter,
                       const PhotoKinematics& event, const ForwardLegState& state, const std::size_t initial_helicity,
                       const std::size_t final_helicity) {
  if (parameter.FORWARD_NOFLIP && initial_helicity != final_helicity) { return; }
  
  // Use the same electromagnetic current on both proton legs, also in diagonal rows
  const double lambda_in  = static_cast<double>(spin::BinaryHelicityLabelX2(initial_helicity)) / 2.0;
  const double lambda_out = static_cast<double>(spin::BinaryHelicityLabelX2(final_helicity)) / 2.0;
  const double f2         = parameter.photo.proton_pauli ? tensor.F2_(event.t) : 0.0;
  const auto   current =
      tensor.DiracPauliCurrent(event.p, event.prime, event.prime - event.p, MDirac::FermionKind::Particle, lambda_in,
                               lambda_out, PDG::mp, tensor.F1_(event.t), f2);
  const M4Vec transfer =
      state.leg == ForwardBeamLeg::Upper ? state.outgoing - state.incoming : state.incoming - state.outgoing;
  const auto phase = MDirac::ElasticSpinHalfColliderCurrentPhase(MDirac::FermionKind::Particle, state.Index(),
                                                                 lambda_in, lambda_out, transfer.Phi());
  Current    target_lower{};
  for (const auto& mu : indices(target_lower)) { target_lower[mu] = phase * current[mu]; }

  Current target_upper{};
  for (const auto& nu : indices(target_lower)) { target_upper[nu] = nu == 0 ? target_lower[nu] : -target_lower[nu]; }
  const double ht       = (event.q - event.kplus).M2();
  const double hu       = (event.q - event.kminus).M2();
  const double hs       = (event.kplus + event.kminus).M2();
  const double mass2    = math::pow2(event.kplus.M());
  const double lambda2  = math::pow2(parameter.photo.lambda.at(event.pdg));
  const double ft       = std::exp((ht - mass2) / lambda2);
  const double fu       = std::exp((hu - mass2) / lambda2);
  const double fs       = std::exp(-(hs - 4.0 * mass2) / lambda2);
  const double offshell = (math::pow2(ft) + math::pow2(fu)) / (1.0 + math::pow2(fs));
  if (std::abs(event.t) < 1.0e-18 || std::abs(ht - mass2) < 1.0e-18 || std::abs(hu - mass2) < 1.0e-18) { return; }

  const M4Vec          vt = event.q - event.kplus + event.kminus;
  const M4Vec          vu = event.q + event.kplus - event.kminus;
  std::complex<double> jt = 0.0;
  std::complex<double> ju = 0.0;
  for (const auto& nu : indices(target_upper)) {
    jt += (vt % nu) * target_upper[nu];
    ju += (vu % nu) * target_upper[nu];
  }
  const double factor = math::pow3(qed::e_QED()) * MesonFF(event.q.M2(), parameter.photo.m0_2) *
                        MesonFF(event.t, parameter.photo.m0_2) * offshell / event.t;
  for (const auto& mu : indices(out)) {
    out[mu] += factor * (((event.q - event.kplus * 2.0) % mu) * jt / (ht - mass2) +
                         ((event.q - event.kminus * 2.0) % mu) * ju / (hu - mass2) - 2.0 * target_lower[mu]);
  }
}

}  // namespace

// Contract a conserved target current with physical photon sources and nuclear transitions
nuclear::PhotoTerm TensorPhotoTerm(const MTensorPomeron& tensor, const MTensorPomeronParam& param,
                                   const LORENTZSCALAR& lts, const ForwardLegState& target,
                                   const TensorPhotoCurrent& current, const flux::PhotoTargetProfile& profile,
                                   bool continuum) {
  nuclear::PhotoTerm out;
  const bool         upper_gamma = target.leg == ForwardBeamLeg::Lower;
  out.direction                  = upper_gamma ? nuclear::PhotoDirection::Upper : nuclear::PhotoDirection::Lower;
  const auto photon  = ResolveForwardLegState(lts, upper_gamma ? ForwardBeamLeg::Upper : ForwardBeamLeg::Lower);
  const auto factors = flux::ResolvePhotoTarget(target, profile, 1.0, 1.0);
  out.photo_current[static_cast<std::size_t>(target.Index() - 1)] = factors.current;
  const auto        source_count                                  = flux::PhotoSourceCount(photon);
  const std::size_t eigen_count                                   = lts.upc_model ? 2 : 1;
  const auto        vertex =
      continuum && !param.photo.proton_pauli ? TensorPhotonVertex::Dirac : TensorPhotonVertex::DiracPauli;
  std::array<MDirac::Spinor, 2> incoming{}, outgoing{};
  if (!photon.IsNuclear()) {
    incoming = tensor.SpinorStates(photon.incoming, "u");
    outgoing = tensor.SpinorStates(photon.outgoing, "ubar");
  }
  std::array<std::vector<FTensor::Tensor1<std::complex<double>, 4>>, 4> sources;
  for (const auto h : spin::BinaryHelicityIndices()) {
    for (const auto prime : spin::BinaryHelicityIndices()) {
      sources[spin::BinaryPairHelicityIndex(h, prime)] =
          tensor.iG_yForwardSources(lts, photon, outgoing[prime], incoming[h], h, prime, vertex);
    }
  }
  for (const auto ha : spin::BinaryHelicityIndices()) {
    for (const auto hb : spin::BinaryHelicityIndices()) {
      for (const auto h1 : spin::BinaryHelicityIndices()) {
        for (const auto h2 : spin::BinaryHelicityIndices()) {
          if (param.FORWARD_NOFLIP && (ha != h1 || hb != h2)) { continue; }
          const auto photon_h = spin::BinaryPairHelicityIndex(upper_gamma ? ha : hb, upper_gamma ? h1 : h2);
          const auto target_h = spin::BinaryPairHelicityIndex(upper_gamma ? hb : ha, upper_gamma ? h2 : h1);
          for (std::size_t eigen = 0; eigen < eigen_count; ++eigen) {
            for (std::size_t source = 0; source < source_count; ++source) {
              const auto           index     = eigen * source_count + source;
              std::complex<double> amplitude = 0.0;
              if (index < sources[photon_h].size()) {
                for (std::size_t mu = 0; mu < 4; ++mu) {
                  amplitude += sources[photon_h][index](mu) * current[target_h](mu);
                }
              }
              for (const auto factor : factors.factor) { out.amplitude.push_back(amplitude * factor); }
            }
          }
        }
      }
    }
  }
  return out;
}

// Bind shared Tensor vertices and immutable photoproduction steering
MTensorPhoto::MTensorPhoto(const MTensorPomeron& tensor, MTensorPomeronParamPtr parameters)
    : tensor(tensor), parameter_handle(std::move(parameters)), parameter(RequirePhotoParameters(parameter_handle)) {}

// Compute one gamma-star proton Drell-Soding current for fixed target spins
MTensorPhoto::Current MTensorPhoto::DSCurrent(const LORENTZSCALAR& lts, const bool photon_from_upper,
                                              const std::size_t initial_helicity,
                                              const std::size_t final_helicity) const {
  Current    out{};
  const auto physical = ResolveForwardLegState(lts, photon_from_upper ? ForwardBeamLeg::Lower : ForwardBeamLeg::Upper);
  const auto target_state = flux::ElementaryPhotoTarget(lts, physical);
  if (std::abs(target_state.emitter.pdg) != PDG::PDG_p) { return out; }
  const PhotoKinematics event = PhotoEvent(lts, photon_from_upper, target_state);
  const double          q2    = -event.q.M2();
  if (!std::isfinite(q2)) { throw AmplitudeFailure("MTensorPhoto: non-finite photon virtuality"); }
  if (q2 < 0.0 || q2 > parameter.photo.q2_max) { return out; }
  if (!event.valid) { throw AmplitudeFailure("MTensorPhoto: invalid Regge subenergies"); }

  const auto [target, scalar] = TargetBilinears(
      tensor, event, target_state, parameter.FORWARD_NOFLIP || physical.IsNuclear(), initial_helicity, final_helicity);
  const auto f = OffshellFactors(event, parameter);
  for (const int exchange_pdg : parameter.photo.exchanges) {
    const auto& exchange = parameter.exchange.FindExchange(exchange_pdg);
    if (exchange.rank == 2) {
      AddRankTwo(out, tensor, parameter, event, target, scalar, exchange_pdg, f);
    } else {
      AddRankOne(out, parameter, event, target, exchange_pdg, f);
    }
  }

  const bool both_protons = std::abs(lts.beam1.pdg) == PDG::PDG_p && std::abs(lts.beam2.pdg) == PDG::PDG_p;
  
  // Count gamma-gamma once and apply the same virtuality cut to both photons
  const double target_q2 = -event.t;
  if (parameter.photo.photon_exchange && target_q2 >= 0.0 && target_q2 <= parameter.photo.q2_max &&
      (!both_protons || photon_from_upper)) {
    AddPhotonExchange(out, tensor, parameter, event, target_state, initial_helicity, final_helicity);
  }
  return out;
}

// Build conserved meson currents in both directions with the common beam contraction
std::vector<nuclear::PhotoTerm> MTensorPhoto::Terms(const LORENTZSCALAR& lts) const {
  std::vector<nuclear::PhotoTerm> terms;
  for (const auto target_leg : {ForwardBeamLeg::Lower, ForwardBeamLeg::Upper}) {
    const auto target = ResolveForwardLegState(lts, target_leg);
    if (!flux::SupportsPhotoTarget(target)) { continue; }
    TensorPhotoCurrent current;
    for (const auto h : spin::BinaryHelicityIndices()) {
      for (const auto prime : spin::BinaryHelicityIndices()) {
        const auto value = DSCurrent(lts, target_leg == ForwardBeamLeg::Lower, h, prime);
        auto&      row   = current[spin::BinaryPairHelicityIndex(h, prime)];
        for (const auto& mu : indices(value)) { row(mu) = value[mu]; }
      }
    }
    terms.push_back(TensorPhotoTerm(tensor, parameter, lts, target, current, {}, true));
  }
  return terms;
}

// Evaluate both physical photon directions in the collider process
double MTensorPhoto::Amp2(LORENTZSCALAR& lts) const {
  lts.proton_good_walker.reset();
  flux::SumPhotoTerms(lts, Terms(lts));
  return gra::SquaredNorm(lts.hamp) / 4.0;
}

}  // namespace gra
