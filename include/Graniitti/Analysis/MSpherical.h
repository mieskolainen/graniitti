// Functional methods for Spherical Harmonic Expansions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef ANALYSIS_MSPHERICAL_H
#define ANALYSIS_MSPHERICAL_H

// C++
#include <complex>
#include <iostream>
#include <vector>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MStatistics.h"

using gra::aux::indices;

namespace gra {
namespace spherical {
// Metadata
struct Meta {
  std::string NAME;     // Input name ID
  std::string LEGEND;   // Legend string
  std::string YAXIS;    // Yaxis string
  std::string MODE;     // MC or DATA
  bool FASTSIM = false; // Fast simulation on
  std::string FRAME;    // Lorentz frame
  double SCALE = 1.0;   // Scale/normalization value

  std::vector<std::string> TITLES; // Phase space titles

  // So we can use this in std::map<>
  bool operator<(const Meta &rhs) const { return NAME < rhs.NAME; }

  void Print() const {
    std::cout << "NAME:     " << NAME << std::endl;
    std::cout << "LEGEND:   " << LEGEND << std::endl;
    std::cout << "YAXIS:    " << YAXIS << std::endl;
    std::cout << "MODE:     " << MODE << std::endl;
    std::cout << "FRAME:    " << FRAME << std::endl;
    std::cout << "FASTSIM:  " << (FASTSIM ? "true" : "false") << std::endl;
    std::cout << "SCALE:    " << SCALE << std::endl;

    std::cout << std::endl;
    for (const auto &i : indices(TITLES)) {
      printf("TITLES[%lu] = %s \n", i, TITLES[i].c_str());
    }
  }
};

// Microevent structure
struct Omega {
  // Decay daughter
  // rest frame variables
  double costheta = 0.0;
  double phi = 0.0;

  // Invariant / system lab frame variables
  double M = 0.0;
  double Pt = 0.0;
  double Y = 0.0;

  double weight = 1.0;

  bool fiducial = false;
  bool selected = false;
};

// Container
struct Data {
  // Metadata
  spherical::Meta META;

  // Events
  std::vector<spherical::Omega> EVENTS;
};

// Detector data for one hypercell, e.g., in (M,Pt,Y)
struct SH_DET {
  // Moment mixing matrix
  MMatrix<double> MIXlm;

  // Delete-group jackknife replicas of the moment mixing matrix
  std::vector<MMatrix<double>> MIXlm_jackknife;

  // Efficiency decomposition coefficients
  std::vector<double> E_lm;
  std::vector<double> E_lm_error;

  double generated_weight = 0.0;
  bool valid = false;
};

// Data for one hypercell, e.g., in (M,Pt,Y)
struct SH {
  // Directly (algebraic) observed moments
  std::vector<double> t_lm_MPP;
  std::vector<double> t_lm_MPP_error;
  MMatrix<double> t_lm_MPP_covariance;

  // Extended Maximum Likelihood fitted moments
  std::vector<double> t_lm_EML;
  std::vector<double> t_lm_EML_error;
  MMatrix<double> t_lm_EML_covariance;

  bool MPP_valid = false;
  bool EML_valid = false;
};

struct MixingEstimate {
  MMatrix<double> value;
  std::vector<MMatrix<double>> jackknife;
  double generated_weight = 0.0;
};

struct MomentEstimate {
  std::vector<double> value;
  MMatrix<double> covariance;
  double sum_weight = 0.0;
  double sum_weight2 = 0.0;
  std::size_t entries = 0;
};

struct WeightSummary {
  double sum_weight = 0.0;
  double sum_weight2 = 0.0;
  std::size_t entries = 0;

  // Add one weight while retaining scale-normalized ESS statistics
  // Accumulate sum w and sum w^2
  void Add(double weight) {
    weights.Add(weight);
    sum_weight = static_cast<double>(weights.Sum());
    sum_weight2 = static_cast<double>(weights.SquareSum());
    entries = weights.Count();
  }

  // Compute the scale-invariant effective number of entries
  // Compute N_eff=(sum_i w_i)^2/sum_i w_i^2
  double EffectiveEntries() const {
    return static_cast<double>(weights.EffectiveSampleSize());
  }

private:
  statistics::ScaledWeightSums weights;
};

MixingEstimate GetGMixing(const std::vector<Omega> &events,
                          const std::vector<std::size_t> &ind, int LMAX,
                          const std::string &mode);

std::pair<std::vector<double>, std::vector<double>>
GetELM(const std::vector<Omega> &MC, const std::vector<std::size_t> &ind,
       int LMAX, const std::string &mode);

MomentEstimate SphericalMoments(const std::vector<Omega> &input,
                                const std::vector<std::size_t> &ind, int LMAX,
                                const std::string &mode);

MMatrix<double> YLM(const std::vector<Omega> &events, int LMAX);

double HarmDotProd(const std::vector<double> &G, const std::vector<double> &x,
                   const std::vector<bool> &ACTIVE, int LMAX);

double HarmDotProdError(const std::vector<double> &G,
                        const MMatrix<double> &covariance,
                        const std::vector<bool> &ACTIVE, int LMAX);

void PrintOutMoments(const std::vector<double> &x,
                     const std::vector<double> &x_error,
                     const std::vector<bool> &ACTIVE, int LMAX);

std::vector<std::size_t>
GetIndices(const std::vector<Omega> &events, const std::vector<double> &M,
           const std::vector<double> &Pt, const std::vector<double> &Y,
           const std::vector<bool> &include_upper = {false, false, false});

WeightSummary SummarizeWeights(const std::vector<Omega> &events,
                               const std::vector<std::size_t> &ind,
                               const std::string &mode);

MMatrix<double> ActiveResponsePseudoInverse(const MMatrix<double> &response,
                                            const std::vector<bool> &active,
                                            double svd_regularization);

int LinearInd(int l, int m);
double CalcError(double f2, double f, double N);
void PrintMatrix(FILE *fp, const std::vector<std::vector<double>> &A);

MMatrix<double> Y_real_synthesize(const std::vector<double> &c_lm,
                                  const std::vector<bool> &ACTIVE,
                                  std::size_t N, std::vector<double> &costheta,
                                  std::vector<double> &phi,
                                  bool normalized = false);

std::vector<double> ErrorProp(const MMatrix<double> &A,
                              const std::vector<double> &x);

MMatrix<double> CovarianceProp(const MMatrix<double> &A,
                               const MMatrix<double> &covariance);

std::vector<double> CovarianceErrors(const MMatrix<double> &covariance);

MMatrix<double> ResponseJackknifeCovariance(
    const std::vector<MMatrix<double>> &inverse_response,
    const std::vector<MMatrix<double>> &forward_response,
    const std::vector<double> &observed, const std::vector<bool> &active,
    double svd_regularization);

} // namespace spherical
} // namespace gra

#endif
