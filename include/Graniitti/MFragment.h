// Toy fragmentation class
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFRAGMENTATION_H
#define MFRAGMENTATION_H

// C++
#include <complex>
#include <map>
#include <random>
#include <vector>

// Own
#include "Graniitti/MPDG.h"
#include "Graniitti/MRandom.h"

namespace gra {
class MFragment {
 public:
  MFragment() {}
  ~MFragment() {}

  static double TubeFragment(const M4Vec &mother, double M0, const std::vector<double> &m,
                             std::vector<M4Vec> &p, double q, double T, double lambda, double maxpt,
                             MRandom &rng, const std::string &pt_distribution);

  static bool SolveAlpha(double &alpha, double M0, const std::vector<double> &m,
                         const std::valarray<double> &mt, const std::valarray<double> &y);

  static void ExpPowRND(double q, double T, double maxpt, const std::vector<double> &mass,
                        std::vector<double> &x, MRandom &rng);

  static void GetDecayStatus(const std::vector<int> &pdgcode, std::vector<bool> &isstable);

  static void GetForwardMass(double &mass1, double &mass2, bool &excite1, bool &excite2,
                             unsigned int excite, MRandom &random);
  static void GetSingleForwardMass(double &mass, MRandom &random);
  static void NstarDecayTable(int Q, double m0, std::vector<int> &pdgcode, MRandom &rng);
  
  static bool PickParticles(double M, unsigned int N, int B, int S, int Q,
                            std::vector<double> &mass, std::vector<int> &pdgcode, const MPDG &PDG,
                            MRandom &rng);
};

}  // namespace gra

#endif
