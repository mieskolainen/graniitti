// Load complex fixed-helicity MadLoop amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_MADLOOP_H
#define AMP_MG5_MADLOOP_H

#include <complex>
#include <cstddef>
#include <mutex>
#include <string>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Model.h"

namespace gra::mg5 {

// Store one complete fixed-helicity MadLoop evaluation
struct MadLoopResult {
  std::vector<std::complex<double>> amplitude;
  std::vector<std::complex<double>> pole;
  std::vector<double> squared;
  std::vector<int> return_code;
};

// Own one isolated or explicitly serialized generated MadLoop runtime
class MadLoop {
public:
  // Load and initialize one generated MadLoop runtime directory
  MadLoop(const std::string &directory, const std::string &name,
          std::size_t nexternal, std::size_t nhelicity, std::size_t ncolor);

  // Release the generated runtime and serialization handle
  ~MadLoop();

  MadLoop(const MadLoop &) = delete;
  MadLoop &operator=(const MadLoop &) = delete;
  MadLoop(MadLoop &&) = delete;
  MadLoop &operator=(MadLoop &&) = delete;

  // Compute an empty string when the runtime files are readable
  static std::string Check(const std::string &directory);

  // Compute whether the generated runtime was loaded and initialized
  bool Ready() const noexcept;

  // Compute the load or initialization diagnostic
  std::string Error() const;

  // Compute the runtime parameter-card strong coupling
  double AlphaS() const;

  // Compute the runtime parameter-card electromagnetic coupling
  double AlphaQED() const;

  // Initialize all external UFO parameters from their generated SLHA indices
  void InitParameters(const SLHAReader &card);

  // Compute the physical model poles for phase space
  ParticleMap Particles() const;

  // Evaluate every complex color and helicity amplitude
  bool Evaluate(const std::vector<double> &momentum, MadLoopResult &result);

private:
  using InitFunction = void (*)(const char *, int, const char *, int, int *);
  using EvalFunction = void (*)(const double *, double *, double *, double *,
                                double *, double *, int *, int *);
  using ModelFunction = void (*)(double *, double *, double *, int);
  using SizeFunction = void (*)(int *);

  // Read or update evaluated model parameters under the runtime lock
  void Model(int update);

  // Load symbols and initialize immutable model and resource data
  void Load(const std::string &directory);

  // Validate finite amplitudes and cancellation of loop poles
  bool Validate(const MadLoopResult &result);

  // Acquire the process-wide fallback lock for shared Fortran state
  bool Lock();

  // Release the process-wide fallback lock
  void Unlock();

  std::string name_;
  std::size_t nexternal_ = 0;
  std::size_t nhelicity_ = 0;
  std::size_t ncolor_ = 0;
  void *handle_ = nullptr;
  InitFunction init_ = nullptr;
  EvalFunction eval_ = nullptr;
  ModelFunction model_ = nullptr;
  std::vector<std::pair<std::string, std::vector<int>>> parameters_;
  std::vector<double> values_;
  std::vector<double> defaults_;
  ParticleMap particles_;
  std::vector<int> pdgs_;
  int lock_fd_ = -1;
  bool isolated_ = false;
  bool ready_ = false;
  double alpha_s_ = 0.0;
  double alpha_qed_ = 0.0;
  bool scale_qcd_ = false;
  bool scale_qed_ = false;
  std::string error_;
  mutable std::mutex mutex_;
};

} // namespace gra::mg5

#endif
