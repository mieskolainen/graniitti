// Spin polarization functions for scattering
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHELICITYSCATTER_H
#define MHELICITYSCATTER_H

#include <array>
#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// Own headers
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Spin/MHelicityBasis.h"
#include "Graniitti/Spin/MPoleLS.h"

namespace gra {
namespace spin {

// Select the forward spin source and exchange-helicity barrier
struct ForwardSpec {
  ForwardVertexMode mode  = ForwardVertexMode::HelicityResidue;
  double            s0    = 1.0;
  ExchangeBasisType basis = ExchangeBasisType::HelicityTransport;
};

// Store one exact rank-one forward source in beam and exchange helicity space
struct ForwardSourceFactors {
  HelVec beam_transition;
  HelVec exchange_helicity;
};

// Build the exchange-helicity part of one fixed-spin forward vertex
HelVec ExchangeHelicityResidue(const std::vector<double>& projections, const M4Vec& transfer,
                               bool second_exchange_daughter, ForwardSpec forward);

// Build the spherical metric or the identity for a single exchange vertex
HelAmp ExchangeMetric(ExchangeBasisType basis, const std::vector<double>& projections, double pole_spin);

// Build an exchange vertex without the proton spin factors
HelVec ExchangeBoundary(const std::vector<double>& projections, const M4Vec& transfer, bool second_exchange_daughter,
                        ForwardSpec forward);

// Canonical hard proton-pair rows use (i1,i2,f1,f2) with negative helicity 0
struct CanonicalProtonPairSpinLayout {
  // Store a stack allocated permutation from Kronecker rows to hard rows
  struct KroneckerRowMap {
    std::array<std::size_t, 16> destination{};
    std::size_t                 count = 0;

    // Compute the active permutation dimension
    std::size_t size() const { return count; }

    // Compute one destination row
    std::size_t operator[](std::size_t row) const { return destination.at(row); }
  };

  // Compute one full hard row from the two initial and two final helicities
  static constexpr std::size_t HardRow(std::size_t i1, std::size_t i2, std::size_t f1, std::size_t f2) {
    return PairHelicityTransitionIndex(BinaryPairHelicityIndex(i1, i2), BinaryPairHelicityIndex(f1, f2));
  }

  // Compute one compact no-flip hard row with f1=i1 and f2=i2
  static constexpr std::size_t CompactNoFlipHardRow(std::size_t i1, std::size_t i2) {
    return BinaryPairHelicityIndex(i1, i2);
  }

  // Convert producer leg rows with negative helicity first to one canonical row
  static std::size_t FromNegativeFirstLegRows(std::size_t upper_row, std::size_t lower_row, std::size_t upper_rows,
                                              std::size_t lower_rows);

  // Compute the full or compact destination rows with the lower leg fastest
  static KroneckerRowMap KroneckerDestinationRows(std::size_t upper_rows, std::size_t lower_rows);
};

// Apply spin steering without changing the setup-time vertex normalization
HelAmp Steer(const HelAmp& central, const HelVec& spin_steering);

// Build an incoherent spin-independent production operator
HelAmp Blind(std::size_t initial_states, std::size_t spin_states, std::complex<double> scale = 1.0);

// Build a spin-independent operator with exact forward transition sections
HelAmp Blind(const HelVec& upper_transition, const HelVec& lower_transition, std::size_t spin_states,
             std::complex<double> scale = 1.0);

// Select physical rows from one forward helicity source
std::vector<std::size_t> Rows(const MDecayBranch& branch, bool drop_flip);

// Build selected rows of one forward photon or exchanged-state source
HelAmp Forward(const LORENTZSCALAR& lts, const MDecayBranch& branch, const M4Vec& incoming_beam,
               const M4Vec& outgoing_system, bool second_exchange_daughter, const std::vector<std::size_t>& rows,
               const std::string& photon_vertex, ForwardSpec forward);

// Compute exact fixed-spin forward factors or no value for photon sources
std::optional<ForwardSourceFactors> ForwardFactors(const MDecayBranch& branch, const M4Vec& incoming_beam,
                                                   const M4Vec& outgoing_system, bool second_exchange_daughter,
                                                   const std::vector<std::size_t>& rows, ForwardSpec forward);

// Transform a pair into the requested production spin frame
void ProductionFrame(std::vector<M4Vec>& p, const LORENTZSCALAR& lts, const std::string& frame,
                     const M4Vec& helicity_dir);

// Compute the produced-spin rotation from the central CM axes to the selected frame
HelAmp ProductionRotation(const LORENTZSCALAR& lts, const std::string& frame, double J);

// Sew two evaluated vertices with an optional internal pole numerator
HelAmp Subchannel(const EvaluatedPoleSubvertex& upper, const EvaluatedPoleSubvertex& lower, const MParticle& left,
                  const MParticle& right, bool swap, bool skip_unphysical = false,
                  const InternalHelicityMetric* metric = nullptr);

// Project exchange helicities before sewing the internal pole
HelVec ProjectedSubchannel(const EvaluatedPoleSubvertex& upper, const EvaluatedPoleSubvertex& lower,
                           const HelVec& upper_source, const HelVec& lower_source, const MParticle& left,
                           const MParticle& right, bool swap, bool skip_unphysical = false,
                           const InternalHelicityMetric* metric = nullptr);

// Project forward sources before sewing into canonical proton rows
HelAmp ProjectedSubchannel(const EvaluatedPoleSubvertex& upper, const EvaluatedPoleSubvertex& lower,
                           const HelAmp& upper_source, const HelAmp& lower_source, const MParticle& left,
                           const MParticle& right, bool swap, bool skip_unphysical = false,
                           const InternalHelicityMetric* metric = nullptr);

// Contract two forward sources into canonical proton-pair hard rows
HelAmp Contract(const HelAmp& upper, const HelAmp& lower, const HelAmp& central, std::complex<double> scale = 1.0);

// Contract two exact rank-one forward sources without materializing them
HelAmp Contract(const ForwardSourceFactors& upper, const ForwardSourceFactors& lower, const HelAmp& central,
                std::complex<double> scale = 1.0);

// Expand projected continuum subchannels into canonical proton-pair hard rows
HelPair Contract(const ForwardSourceFactors& upper, const ForwardSourceFactors& lower, const HelVec& sub_t,
                 const HelVec& sub_u);

}  // namespace spin
}  // namespace gra

#endif  // MHELICITYSCATTER_H
