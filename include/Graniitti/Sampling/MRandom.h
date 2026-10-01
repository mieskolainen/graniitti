// Random number class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MRANDOM_H
#define MRANDOM_H

// C++
#include <complex>
#include <random>
#include <vector>

namespace gra {
class MRandom {
 public:
  // Calling constructors of member functions
  MRandom() : flat(0, 1), gaussian(0, 1) {
    
    std::seed_seq seedseq({0}); // Default
    rng.seed(seedseq);
  }
  ~MRandom() = default;
  
  // Set random number engine seed
  void SetSeed(uint32_t seed) {

    // Seed via seed sequence
    std::seed_seq seedseq({seed});
    rng.seed(seedseq);
    flat.reset();
    gaussian.reset();

    // Save it
    RNDSEED = seed;
  }

  // Compute current random seed
  uint32_t GetSeed() const { return RNDSEED; }

  // Random sampling functions
  double U(double a, double b);
  double G(double mu, double sigma);

  double ExpRandom(double lambda);
  double ExpBoundedRandom(double a, double b, double lambda);
  double PowerBoundedRandom(double a, double b, double alpha);
  int    NBDRandom(double avgN, double k, int maxvalue);
  int    LogRandom(double p, int maxvalue);
  
  void   DirRandom(const std::vector<double> &alpha, std::vector<double> &y);
  // Sample a Poisson count, optionally conditioned on a positive count
  int    PoissonRandom(double lambda, bool positive = false);
  double RelativisticBWRandom(double m0, double Gamma, double LIMIT = 5.0, double M_MIN = 0.0);
  double CauchyRandom(double m0, double Gamma, double LIMIT = 5.0, double M_MIN = 0.0);

  // Integrate a finite-width relativistic Breit-Wigner over mass squared
  static double RelativisticBWMass2Integral(double m0, double Gamma, double m2min, double m2max);
  
  double ExpPdf(double x, double lambda);
  double ExpBoundedPdf(double x, double a, double b, double lambda);
  double PowerBoundedPdf(double x, double a, double b, double alpha);
  double NBDPdf(int n, double avgN, double k);
  double LogPdf(int k, double p);
  
  // For generators, check: https://nullprogram.com/blog/2017/09/21/

  // 64-bit Mersenne Twister by Matsumoto and Nishimura, fast, basic
  // For some (possible) sources problems, see:
  // [REFERENCE: Harase, https://arxiv.org/abs/1708.06018]
  // std::mt19937_64 rng;

  // 48-bit RANLUX (a bit slower)
  // [REFERENCE: Luscher, https://arxiv.org/abs/hep-lat/9309020]
  std::ranlux48 rng;

  std::ranlux48 get_generator() const {
    return rng;
  }
  
  // Distribution engines
  std::uniform_real_distribution<double> flat;
  std::normal_distribution<double>       gaussian;
  uint32_t                               RNDSEED = 0;  // Random seed set

 private:
  // Nothing
};

}  // namespace gra

#endif
