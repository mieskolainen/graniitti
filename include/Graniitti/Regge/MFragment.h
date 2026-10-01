// Forward baryon excitation and toy fragmentation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFRAGMENTATION_H
#define MFRAGMENTATION_H

// C++
#include <array>
#include <map>
#include <memory>
#include <string>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/MModelCache.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Sampling/MRandom.h"

// Libraries
#include "json.hpp"

namespace gra {

// Beam remnant fragmentation selected by process steering
enum class BeamFragType { None, FewBody, Cylinder, Diquark };

// Forward excitation prescriptions selected before event generation
enum class DissociationType { None, Soft, Hera, Structure, TripleRegge };

// Select the hadronic default and the target transition in gamma-Pomeron production
struct DissociationModel {
  DissociationType hadron = DissociationType::None;
  DissociationType photo = DissociationType::None;
};

// Immutable forward excitation and fragmentation parameters
struct MNstarParam {
  struct FewBody {
    std::array<double, 2> body_br = {};
    std::array<double, 2> two_body_br = {};
    std::array<double, 3> three_body_br = {};
    bool unweight = false;
  };

  struct Cylinder {
    double charged_prob = 0.0;
    double pair_prob = 0.0;
    double particle_prob = 0.0;
    std::array<double, 3> species_ratio = {};
    std::string pt_distribution;
    double mult_a = 0.0;
    double mult_b = 0.0;
    double q2_ref = 0.0;
    double mass_unit = 0.0;
    double q_power = 0.0;
    double T = 0.0;
    double impact_power = 0.0;
    double max_pt = 0.0;
    double exp_lambda = 0.0;
    unsigned int outer_trials = 0;
    unsigned int mult_trials = 0;
    unsigned int pick_trials = 0;
    unsigned int tube_trials = 0;
    unsigned int pt_trials = 0;
    unsigned int pt_bins = 0;
    unsigned int alpha_iter = 0;
    double alpha_tol = 0.0;
  };

  struct String {
    double first_valence_prob = 0.0;
    double scalar_diquark_mass = 0.0;
    double vector_diquark_mass = 0.0;
  };

  std::map<std::string, DissociationModel> model;
  double single_side_prob = 0.0;
  double mass_margin = 0.0;
  FewBody fewbody;
  Cylinder cylinder;
  String string;

  // Compute the class prescriptions with the photon-target override resolved
  DissociationModel Model(const std::string &istate) const;

  // Read one immutable PARAM_NSTAR block
  void Configure(const nlohmann::json &block, const std::string &source);
};

using MNstarParamPtr = std::shared_ptr<const MNstarParam>;

// Construct one immutable forward excitation parameter block
MNstarParamPtr ReadNstarParam(const MModelTune &tune);

// Compute the run owned immutable forward excitation parameter block
MNstarParamPtr GetNstarParam(MModelCache &cache);

// Generate simple forward fragmentation final states
class MFragment {
public:
  // Generate N-body cylinder phase space
  static double TubeFragment(const M4Vec &mother, double M0,
                             const std::vector<double> &m,
                             std::vector<M4Vec> &p, double q, double T,
                             double lambda, double maxpt, MRandom &rng,
                             const MNstarParam &param);

  // Solve the cylinder rapidity scale
  static bool SolveAlpha(double &alpha, double M0, const std::vector<double> &m,
                         const std::valarray<double> &mt,
                         const std::valarray<double> &y,
                         const MNstarParam &param);

  // Sample the power-law exponential transverse-momentum model
  static bool ExpPowRND(double q, double T, double maxpt,
                        const std::vector<double> &mass, std::vector<double> &x,
                        MRandom &rng, const MNstarParam &param);

  // Compute decay status flags for simple forward daughters
  static void GetDecayStatus(const std::vector<int> &pdgcode,
                             std::vector<bool> &isstable);

  // Select one low-mass N-star decay channel
  static void NstarDecayTable(int Q, double m0, std::vector<int> &pdgcode,
                              MRandom &rng, const MNstarParam &param);

  // Select one cylinder particle set with conserved quantum numbers
  static bool PickParticles(double M, unsigned int N, int B, int S, int Q,
                            std::vector<double> &mass,
                            std::vector<int> &pdgcode, const MPDG &PDG,
                            MRandom &rng, const MNstarParam &param);
};

} // namespace gra

#endif
