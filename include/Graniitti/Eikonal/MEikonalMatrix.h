// Coupled-channel matrix eikonal state
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEIKONAL_MATRIX_H
#define MEIKONAL_MATRIX_H

// C++
#include <array>
#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Eikonal/MEikonalHelicity.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Regge/MSoftModel.h"

namespace gra {

struct MEikonalNumerics;

// Immutable coupled Good-Walker and helicity eikonal tables
class MEikonalMatrix {
public:
  using Complex = std::complex<double>;
  using Matrix = MMatrix<Complex>;

  // Named defining entries of the reference-plane pair spin operator
  struct PairHelicityBank {
    Matrix phi1;
    Matrix phi2;
    Matrix phi3;
    Matrix phi4;
    Matrix phi5;
    Matrix phi5_first;

    // Compute one named pair-space helicity matrix
    const Matrix &Component(ElasticHelicityComponent component) const;

    // Compute one mutable named pair-space helicity matrix during construction
    Matrix &Component(ElasticHelicityComponent component);
  };

  // Store every negative-first proton spin transition in Good Walker pair space
  using PairSpinBank = std::array<Matrix, 16>;

  // Construct all coupled-channel tables for one initialized eikonal object
  static std::shared_ptr<const MEikonalMatrix>
  Build(const SoftModelPtr &model, double s,
        const std::vector<MParticle> &initialstate,
        const MEikonalNumerics &numerics, bool build_amplitude);

  // Access the immutable SOFT model used to construct these tables
  const SoftModelPtr &SoftModelHandle() const noexcept { return soft_model; }

  // Compute the Good-Walker channel count
  std::size_t ChannelCount() const { return tables.N; }

  // Access the immutable impact parameter interpolation nodes
  const std::vector<double> &ImpactParameterNodes() const noexcept {
    return tables.b_node;
  }

  // Access the immutable momentum transfer interpolation nodes
  const std::vector<double> &MomentumTransferNodes() const noexcept {
    return tables.q2_node;
  }

  // Compute separated crossing-even and crossing-odd diagonal profiles
  std::pair<Complex, Complex> DensityProfiles(double bt, std::size_t i,
                                              std::size_t k) const;

  // Compute the spin-averaged diagonal eigenchannel amplitude
  Complex ScreeningAmplitude(double kt2, std::size_t i, std::size_t k) const;

  // Compute physical helicity amplitudes for one Good-Walker final state
  ElasticHelicityAmplitudes HelicityAmplitudes(double kt2, std::size_t f1,
                                               std::size_t f2) const;

  // Compute the exact physical spin matrix for one Good Walker transition
  ProtonHelicityMatrix HelicityMatrix(double kt2, std::size_t f1,
                                      std::size_t f2, std::size_t i1 = 0,
                                      std::size_t i2 = 0,
                                      double azimuth = 0.0) const;

  // Compute exclusive amplitudes built from the selected screening exchanges
  ElasticHelicityAmplitudes
  ScreeningHelicityAmplitudes(double kt2, std::size_t f1, std::size_t f2,
                              std::size_t i1 = 0, std::size_t i2 = 0) const;

  // Compute the exact selected-exchange screening spin matrix
  ProtonHelicityMatrix ScreeningHelicityMatrix(double kt2, std::size_t f1,
                                               std::size_t f2,
                                               std::size_t i1 = 0,
                                               std::size_t i2 = 0,
                                               double azimuth = 0.0) const;

  // Compute the six named defining entries of the radial screening operator
  PairHelicityBank PairScreeningHelicityBank(double kt2) const;

  // Compute every radial pair-space screening spin transition
  PairSpinBank PairScreeningSpinBank(double kt2) const;

  // Compute physical impact-parameter helicity amplitudes for one final state
  ElasticHelicityAmplitudes ImpactHelicityAmplitudes(double bt, std::size_t f1,
                                                     std::size_t f2) const;

  // Compute the exact impact-parameter spin matrix for one physical transition
  ProtonHelicityMatrix ImpactHelicityMatrix(double bt, std::size_t f1,
                                            std::size_t f2, std::size_t i1 = 0,
                                            std::size_t i2 = 0,
                                            double azimuth = 0.0) const;

  // Compute one physical impact-parameter helicity S-matrix transition
  Matrix PhysicalImpactSMatrix(double bt, std::size_t f1 = 0,
                               std::size_t f2 = 0, double azimuth = 0.0) const;

  // Compute separately screened crossing helicity amplitudes
  ElasticCrossingHelicityAmplitudes
  CrossingHelicityAmplitudes(double kt2, std::size_t f1, std::size_t f2) const;

  // Pomeron-only effective opacity spectrum defined from the actual S matrix:
  //   S_P^dagger S_P |a> = exp(-lambda_a) |a>.
  struct CutOpacitySpectrum {
    // A zero survival singular value gives positive infinite opacity
    std::vector<double> eigenvalues;
    std::vector<double> incoming_weights;
  };

  // Compute the Pomeron-only pre-unitarization chi operator at one impact
  // parameter
  Matrix CutOpacity(double bt) const;

  // Compute the q-eikonal-consistent effective cut-opacity spectrum
  CutOpacitySpectrum CutSpectrum(double bt) const;

  // Compute unpolarized incoming proton states in the cut-opacity basis
  std::vector<std::vector<double>> IncomingCutStates() const;

  // Compute the largest eikonal S-matrix singular value
  double MaxSingularValue() const { return tables.max_singular_value; }

  // Compute the numerical tolerance used for S-matrix contractivity
  double UnitarityTolerance() const { return unitarity_tolerance; }

  // Compute the versioned physics and numerical table fingerprint
  const std::string &RuntimeFingerprint() const noexcept {
    return runtime_fingerprint;
  }

private:
  // Hold all mutable implementation tables until const publication
  struct Tables {
    std::size_t N = 0;
    std::size_t D = 0;
    std::vector<double> U;
    std::vector<double> b_node;
    std::vector<double> q2_node;
    std::vector<Matrix> chi_even_b;
    std::vector<Matrix> chi_odd_b;
    std::vector<Matrix> chi_cut_b;
    std::array<std::vector<Matrix>, 16> amplitude_spin_b;
    std::array<std::vector<Matrix>, 16> amplitude_spin_q;
    std::array<Matrix, 16> amplitude_spin_zero;
    std::array<std::vector<Matrix>, 16> screening_spin_b;
    std::array<std::vector<Matrix>, 16> screening_spin_q;
    std::array<Matrix, 16> screening_spin_zero;
    std::array<std::vector<Matrix>, 16> crossing_even_spin_q;
    std::array<Matrix, 16> crossing_even_spin_zero;
    std::array<std::vector<Matrix>, 16> crossing_odd_spin_q;
    std::array<Matrix, 16> crossing_odd_spin_zero;
    double max_singular_value = 0.0;

    // Compute one Good-Walker rotation element
    double Uat(std::size_t row, std::size_t col) const {
      return U.at(row * N + col);
    }
  };

  // Construct mutable table storage bound to one immutable SOFT model
  explicit MEikonalMatrix(SoftModelPtr model) : soft_model(std::move(model)) {}

  // Private matrix and interpolation table storage
  Tables tables;

  // Immutable SOFT model shared with the owning eikonal object
  SoftModelPtr soft_model;

  // Numerical tolerance retained for cut-spectrum validation
  double unitarity_tolerance = 0.0;

  // Versioned key of the immutable strong matrix table
  std::string runtime_fingerprint;
};

} // namespace gra

#endif
