// Serial four-particle and six-particle Regge ladders
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEMULTI_H
#define MREGGEMULTI_H

#include <complex>
#include <array>
#include <span>
#include <string>
#include <vector>

#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MReggeParam.h"
#include "Graniitti/Spin/MHelicityScatter.h"

namespace gra::regge {

// Store one allowed local pair vertex resolved during process construction
struct ContinuumVertexPlan {
  const VertexParam *param            = nullptr;
  int         upper_alias      = 0;
  int         lower_alias      = 0;
  std::size_t upper_trajectory = 0;
  std::size_t lower_trajectory = 0;
  ReggeContinuumPole             pole;
  std::array<std::vector<double>, 2> basis;
};

// Store one ordered final-state pair resolved during process construction
struct ContinuumPairPlan {
  std::size_t                       first       = 0;
  std::size_t                       second      = 0;
  const PairParam                  *param       = nullptr;
  double                            first_mass2 = 0.0;
  std::vector<ContinuumVertexPlan> vertex;
};

// Store one connected trajectory chain and its exact propagator cache slots
struct ContinuumChainPlan {
  std::array<std::size_t, 3> vertex        = {};
  std::array<std::size_t, 2> internal_slot = {};
  std::size_t                boundary      = 0;
};

// Store one common pair of forward Regge beam factors
struct ContinuumBoundaryPlan {
  int upper_exchange = 0;
  int first_particle = 0;
  int lower_exchange = 0;
  int last_particle  = 0;
};

// Store one ordered Feynman graph and its exact pair and meson cache slots
struct ContinuumDiagramPlan {
  std::array<int, 6>              order       = {};
  std::array<std::size_t, 3>      pair        = {};
  std::array<std::size_t, 3>      meson_slot  = {};
  std::vector<ContinuumChainPlan> chain;
};

// Store one immutable four- or six-body continuum evaluation plan
struct ContinuumPlan {
  std::size_t                               central_count = 0;
  std::array<std::size_t, 3>                meson_slots   = {};
  std::array<std::size_t, 2>                internal_slots = {};
  std::vector<std::vector<int>>             permutations;
  std::vector<Topology>                     topologies;
  std::vector<ContinuumPairPlan>             pair;
  std::vector<ContinuumDiagramPlan>          diagram;
  std::vector<ContinuumBoundaryPlan>         boundary;
};

// Build the immutable four-body continuum graph and cache-slot plan
ContinuumPlan BuildContinuum4Plan(const LORENTZSCALAR &lts, const Param &param, ReggeProductionModel model);

// Build the immutable six-body continuum graph and cache-slot plan
ContinuumPlan BuildContinuum6Plan(const LORENTZSCALAR &lts, const Param &param, ReggeProductionModel model);

}  // namespace gra::regge

namespace gra::ladder {

// Compute one prepared local continuum pole by its exact structural key
const ReggeContinuumPole &Pole(const LORENTZSCALAR &lts, const ReggeContinuumPoleKey &key, const std::string &context);

// Compute the model-specific local continuum contraction basis
std::vector<double> Basis(const ReggeContinuumPole &pole, ReggeProductionModel model, std::size_t leg);

// Evaluate one model-specific local final-pair kernel
MMatrix<std::complex<double>> Kernel(const ReggeContinuumPole &pole, const regge::Param &param, ReggeProductionModel model, const MDecayBranch &first, const MDecayBranch &second, const M4Vec &upper_exchange, const M4Vec &lower_exchange, double alpha_upper,
                                     double alpha_lower);

// Validate that adjacent kernels use the same retained spin basis
void CheckBasis(const std::vector<double> &upper, const std::vector<double> &lower, const std::string &context);

// Contract a serial kernel chain between its two boundary vertices
std::complex<double> Contract(const std::vector<std::complex<double>> &upper, std::span<const MMatrix<std::complex<double>>> kernels, std::span<const MMatrix<std::complex<double>>> metrics, const std::vector<std::complex<double>> &lower);

}  // namespace gra::ladder

#endif
