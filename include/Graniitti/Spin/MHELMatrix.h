// Decay helicity coupling structures [HEADER ONLY file]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHELMATRIX_H
#define MHELMATRIX_H

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Regge/MFormFactor.h"
#include "Graniitti/Spin/MWigner.h"

namespace gra {

// Store helicity amplitudes as a matrix, vector or ordered t/u pair
using HelAmp  = MMatrix<std::complex<double>>;
using HelVec  = std::vector<std::complex<double>>;
using HelPair = std::pair<HelAmp, HelAmp>;

}  // namespace gra

namespace gra::spin {

// Identify one orbital and total-spin coupling
struct LSCoupling {
  std::size_t l = 0;
  std::size_t two_s = 0;
};

// Store one nonzero or explicitly selected LS coefficient
struct LSTerm {
  std::size_t l = 0;
  std::size_t two_s = 0;
  std::complex<double> coefficient = 0.0;
};

// Store unique LS coefficients in ascending L and 2S order
class LSCoefficients {
public:
  using Container = std::vector<LSTerm>;
  using iterator = Container::iterator;
  using const_iterator = Container::const_iterator;

  // Check whether no LS coefficient is stored
  bool Empty() const { return terms.empty(); }

  // Compute the number of stored LS coefficients
  std::size_t Size() const { return terms.size(); }

  // Remove every LS coefficient
  void Clear() { terms.clear(); }

  // Compute the first LS coefficient
  iterator begin() { return terms.begin(); }

  // Compute one past the last LS coefficient
  iterator end() { return terms.end(); }

  // Compute the first constant LS coefficient
  const_iterator begin() const { return terms.begin(); }

  // Compute one past the last constant LS coefficient
  const_iterator end() const { return terms.end(); }

  // Compute the first constant LS coefficient for range access
  const_iterator cbegin() const { return terms.cbegin(); }

  // Compute one past the last constant LS coefficient for range access
  const_iterator cend() const { return terms.cend(); }

  // Insert one coefficient while rejecting a duplicate LS index
  bool Insert(std::size_t l, std::size_t two_s,
              std::complex<double> coefficient) {
    const auto position = LowerBound(l, two_s);
    if (position != terms.end() && position->l == l &&
        position->two_s == two_s) {
      return false;
    }
    terms.insert(position, {l, two_s, coefficient});
    return true;
  }

  // Set one coefficient while retaining unique sorted LS indices
  void Set(std::size_t l, std::size_t two_s, std::complex<double> coefficient) {
    const auto position = LowerBound(l, two_s);
    if (position != terms.end() && position->l == l &&
        position->two_s == two_s) {
      position->coefficient = coefficient;
      return;
    }
    terms.insert(position, {l, two_s, coefficient});
  }

  // Compute whether one LS index is stored
  bool Contains(std::size_t l, std::size_t two_s) const {
    return Find(l, two_s) != nullptr;
  }

  // Compute one mutable LS term or null when it is absent
  LSTerm *Find(std::size_t l, std::size_t two_s) {
    const auto position = LowerBound(l, two_s);
    return position != terms.end() && position->l == l &&
                   position->two_s == two_s
               ? &*position
               : nullptr;
  }

  // Compute one constant LS term or null when it is absent
  const LSTerm *Find(std::size_t l, std::size_t two_s) const {
    const auto position = LowerBound(l, two_s);
    return position != terms.end() && position->l == l &&
                   position->two_s == two_s
               ? &*position
               : nullptr;
  }

  // Compute one mutable coefficient and reject an absent LS index
  std::complex<double> &At(std::size_t l, std::size_t two_s) {
    LSTerm *term = Find(l, two_s);
    if (term == nullptr) {
      throw std::out_of_range("LSCoefficients::At: LS index is absent");
    }
    return term->coefficient;
  }

  // Compute one constant coefficient and reject an absent LS index
  const std::complex<double> &At(std::size_t l, std::size_t two_s) const {
    const LSTerm *term = Find(l, two_s);
    if (term == nullptr) {
      throw std::out_of_range("LSCoefficients::At: LS index is absent");
    }
    return term->coefficient;
  }

  // Scale every stored LS coefficient
  void Scale(std::complex<double> scale) {
    for (LSTerm &term : terms) {
      term.coefficient *= scale;
    }
  }

  // Compute the squared norm of all stored LS coefficients
  // ||alpha||^2 = sum_(LS) |alpha_LS|^2
  double Norm2() const {
    double norm2 = 0.0;
    for (const LSTerm &term : terms) {
      norm2 += std::norm(term.coefficient);
    }
    return norm2;
  }

  // Compute whether every stored LS coefficient is finite
  bool IsFinite() const {
    return std::all_of(terms.cbegin(), terms.cend(), [](const LSTerm &term) {
      return std::isfinite(term.coefficient.real()) &&
             std::isfinite(term.coefficient.imag());
    });
  }

  // Remove coefficients at or below one magnitude threshold
  void RemoveBelow(double threshold) {
    if (!std::isfinite(threshold) || threshold < 0.0) {
      throw std::invalid_argument(
          "LSCoefficients::RemoveBelow: threshold must be nonnegative");
    }
    std::erase_if(terms, [threshold](const LSTerm &term) {
      return std::abs(term.coefficient) <= threshold;
    });
  }

private:
  // Compute the ordered insertion position of one LS index
  iterator LowerBound(std::size_t l, std::size_t two_s) {
    return std::lower_bound(terms.begin(), terms.end(), std::pair{l, two_s},
                            [](const LSTerm &term, const auto &index) {
                              return std::pair{term.l, term.two_s} < index;
                            });
  }

  // Compute the ordered constant position of one LS index
  const_iterator LowerBound(std::size_t l, std::size_t two_s) const {
    return std::lower_bound(terms.cbegin(), terms.cend(), std::pair{l, two_s},
                            [](const LSTerm &term, const auto &index) {
                              return std::pair{term.l, term.two_s} < index;
                            });
  }

  Container terms;
};

} // namespace gra::spin

namespace gra {

// One normalized Jacob-Wick LS contribution to a helicity matrix
struct LSHelicityComponent {
  std::size_t l = 0;
  std::size_t two_s = 0;
  std::complex<double> alpha = 0.0;
  MMatrix<std::complex<double>> matrix;
};

// Index one finite GP orbital SU(2) row
struct GPOrbitalTerm {
  std::size_t l = 0;
  std::size_t two_s = 0;
  std::size_t spin_index = 0;
  std::size_t su2_offset = 0;
};

// Store the finite GP orbital SU(2) factors prepared during initialization
struct GPOrbitalBasis {
  std::size_t nmu = 0;
  std::vector<std::size_t> spins;
  std::vector<GPOrbitalTerm> terms;
  std::vector<double> su2;
};

// Store one finite Jacob-Wick helicity row map
struct JWRotationBasis {
  bool ready = false;
  std::vector<double> difference;
  std::vector<std::size_t> row;
  std::shared_ptr<const wigner::Rotation> rotation;
};

// Select physical helicity transport, angle-free Regge m or one reduced column
enum class ExchangeBasisType {
  HelicityTransport,
  ReggeHelicity,
  ReducedRegge
};

// Select a physical pole or analytic Regge helicity domain
enum class HelicityDomain { PhysicalPole, ReggeTrajectory };

// Select unset, LS or reduced helicity coupling coordinates
enum class CouplingBasis { Unspecified, LS, Helicity };

// Decay helicity amplitude information
struct HELMatrix {
  // Clebsch-Gordan amplitude matrix (calculated by SU(2) decomposition
  // routines)
  gra::MMatrix<std::complex<double>> T;

  // Explicitly supplied reduced helicity entries for direct-helicity input
  gra::MMatrix<bool> T_set;

  // Compact active direct-helicity matrix indices prepared during setup
  std::vector<std::pair<std::size_t, std::size_t>> T_active;

  // Coordinate domain carried by the helicity operator
  HelicityDomain domain = HelicityDomain::PhysicalPole;

  // Coupling coordinates used to construct the operator
  CouplingBasis coupling_basis = CouplingBasis::Unspecified;

  // Compute whether reduced helicity coordinates define the coupling
  bool UsesHelicityCouplings() const noexcept {
    return coupling_basis == CouplingBasis::Helicity;
  }

  // Compute whether the operator carries analytic Regge coordinates
  bool UsesReggeDomain() const noexcept {
    return domain == HelicityDomain::ReggeTrajectory;
  }

  // Compute whether LS coordinates define the coupling
  bool UsesLSCouplings() const noexcept {
    return coupling_basis == CouplingBasis::LS;
  }

  // Exchange-spin representation used by the production contraction
  ExchangeBasisType exchange_basis = ExchangeBasisType::HelicityTransport;

  // Momentum scale used by dimensionless analytic LS couplings
  double analytic_Lambda = 1.0;

  // MMAX used by the analytic exchange m=-MMAX,...,+MMAX column basis
  int analytic_MMAX = -1;

  // Resonance spin state projections
  std::vector<double> Jz_values;

  // Init values with -1 for initialization checks
  double J = -1;
  double s1 = -1;
  double s2 = -1;

  // Final state helicity projections
  MMatrix<double> lambda_values;
  MMatrix<std::size_t> lambda_idx;

  // Finite Jacob-Wick row map prepared with the helicity projections
  JWRotationBasis jw_rotation;

  // Decay alpha_ls or cached production g_ls in ascending L and 2S order
  spin::LSCoefficients alpha_ls;

  // Sparse crossed LS coefficients aligned with analytic exchange m columns
  std::vector<spin::LSCoefficients> m_ls;

  // Finite GP orbital factors aligned with alpha_ls
  GPOrbitalBasis gp_orbital;

  // Crossed GP orbital factors aligned with m_ls
  std::vector<GPOrbitalBasis> m_orbital;

  // Pole-normalized LS components used by optional mass scaling
  std::vector<LSHelicityComponent> ls_components;

  // Parity conservation
  bool P_symmetry = true;
  bool C_symmetry = false;

  // --------------------------------------------------------------------
  // Tensor Pomeron couplings
  std::vector<double> g_decay_TP; // Decay couplings
  regge::FFParam ff_decay;
  // --------------------------------------------------------------------

  // Branching Ratio to a particular decay
  double BR = 1.0;
  bool BR_set = false;
  double zeta = 0.0;                  // Free phase parameter [rad]
  std::complex<double> g_decay = 1.0; // Equivalent decay coupling
};

} // namespace gra

#endif
