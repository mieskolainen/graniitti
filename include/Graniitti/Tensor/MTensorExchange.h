// Tensor Pomeron and tensor Reggeon exchange steering
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSOREXCHANGE_H
#define MTENSOREXCHANGE_H

// C++
#include <array>
#include <complex>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

// Own
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Regge/MReggeParam.h"
#include "Graniitti/Regge/MSoftModel.h"

// Libraries
#include "json.hpp"

namespace gra {

// Select one effective exchange field of the tensor-pomeron framework
enum class TensorExchangeType {
  Pomeron,
  TensorReggeon,
  Odderon,
  VectorReggeon
};

// Store one soft exchange trajectory selected by its fixed-spin PDG alias
struct TensorExchangeParam {
  int pdg = 0;
  TensorExchangeType type = TensorExchangeType::Pomeron;
  int rank = 2;
  double delta = 0.0;
  double ap = 0.0;
  int eta = 1;
  double scale = 1.0;
};

// Store one exchange-to-particle-pair vertex from CON_TP.json
struct TensorVertexParam {
  int exchange_pdg = 0;
  std::array<int, 2> pair_pdgs = {0, 0};
  std::vector<double> g_tensor;
  std::vector<std::size_t> active_g_tensor;
  regge::FFParam ff_transfer;
  regge::FFParam ff_offshell;
  bool threshold = false;
};

// Store one final-state selection of unordered soft exchange pairs
struct TensorContinuumChannel {
  std::array<int, 2> final_pdgs = {0, 0};
  bool fallback = false;
  std::vector<std::array<int, 2>> exchange_pairs;
};

// Immutable exchange trajectories, vertices and channel selection
class MTensorExchangeModel {
public:
  // Read Tensor exchange and continuum steering
  void Configure(const nlohmann::json &general, const nlohmann::json &vertex_table,
                 const gra::MPDG &pdg_table, double coupling_min);

  // Compute one configured soft exchange by fixed-spin PDG alias
  const TensorExchangeParam &FindExchange(int pdg) const;

  // Compute the SOFT exchange matched to one Tensor exchange alias
  SoftExchangeId SoftId(int pdg, const SoftModel &model) const;

  // Compute one configured exchange-to-hadron continuum vertex
  const TensorVertexParam &FindVertex(int exchange_pdg, int hadron_pdg) const;

  // Compute one configured exchange-to-particle-pair vertex
  const TensorVertexParam &FindVertex(int exchange_pdg, int first_pdg,
                                      int second_pdg) const;

  // Compute the selected unordered exchange pairs for one final-state pair
  const std::vector<std::array<int, 2>> &
  FindContinuumPairs(int first_pdg, int second_pdg) const;

  // Compute common internal lines implied by two ordinary pair vertices
  std::vector<int> FindTransfers(int first_exchange, int second_exchange,
                                 int first_hadron, int second_hadron) const;

  // Compute active common internal lines implied by two pair vertices
  std::vector<int> FindActiveTransfers(int first_exchange, int second_exchange,
                                       int first_hadron,
                                       int second_hadron) const;

  // Expand one unordered exchange pair into its physical beam orderings
  static std::vector<std::array<int, 2>>
  OrderedPairs(const std::array<int, 2> &pair);

  // Compute the scalar Regge factor of one soft exchange propagator
  std::complex<double> PropagatorFactor(int exchange_pdg, double s,
                                        double t) const;

  // Convert a card coupling to the common rank-two baryon vertex
  double BaryonCoupling(int exchange_pdg, double coupling) const;

  // Convert a card coupling to the common rank-two pseudoscalar vertex
  double PseudoscalarCoupling(int exchange_pdg, double coupling) const;

  // Convert a card coupling to the common vector baryon vertex
  double VectorBaryonCoupling(int exchange_pdg, double coupling) const;

  // Convert a card coupling to the common vector pseudoscalar vertex
  double VectorPseudoscalarCoupling(int exchange_pdg, double coupling) const;

  // Compute all continuum hadrons with the requested spinX2
  std::vector<int> HadronPdgs(int spinX2) const;

  // Check whether exchange and continuum parsing is complete
  bool Initialized() const noexcept { return initialized; }

private:
  std::map<int, TensorExchangeParam> exchanges;
  std::map<std::array<int, 3>, TensorVertexParam> vertices;
  std::vector<TensorContinuumChannel> continuum;
  std::map<int, int> hadron_spins;
  bool initialized = false;
};

} // namespace gra

#endif
