// Hard diffractive Pomeron PDF access
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHARDPOMERONPDF_H
#define MHARDPOMERONPDF_H

// C++
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// Own
#include "Graniitti/PDF/MLHAPDF.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Tech/MFixedStore.h"

namespace gra {

// Steering parameters for factorized hard-Pomeron densities
struct MHardPomeronPDFParam {
  std::string         DPDF_SET       = "GKG18_DPDF_FitB_LO";
  int                 DPDF_MEMBER    = 0;
  double              alpha0         = 1.0988;
  double              alpha_prime    = 0.0;
  double              B_flux         = 7.0;
  double              norm           = 1.0;
  std::vector<int>    parton_flavours = {21, 1, 2, 3, 4, 5};
  std::vector<double> xi_range       = {1e-4, 0.1};
  std::vector<double> beta_range     = {1e-5, 1.0};
  std::vector<double> t_range        = {-1.0, 0.0};
  double              Q2_min         = 1.0;
  double              remnant_mass   = 6.0;

  // Validate model-card ranges before using LHAPDF
  void Validate() const;
};

// Factorized f^D_i/p(xi,beta,t,Q2) = f_IP/p(xi,t) f_i/IP(beta,Q2)
class MHardPomeronPDF {
 public:
  // Construct with the default model-card parameters
  MHardPomeronPDF();

  // Construct from a GRANIITTI GENERAL.json model card
  explicit MHardPomeronPDF(const std::string &modelfile);

  // Construct from one already captured immutable model snapshot
  explicit MHardPomeronPDF(const SoftModelPtr &soft_model);

  // Construct from one model snapshot and a run owned LHAPDF store
  MHardPomeronPDF(const SoftModelPtr &soft_model, MLHAPDFStore &pdf_store);

  // Read parameters and replace the PDF only after successful initialization
  void ReadParameters(const std::string &modelfile);

  // Read only the configured hard-process parton flavours without loading LHAPDF
  static std::vector<int> ReadPartonFlavours(const std::string &modelfile);

  // Compute alpha_P(t)
  double AlphaP(double t) const;

  // Compute the Pomeron flux f_IP/p(xi,t)
  double Flux(double xi, double t) const;

  // Compute f_i/IP(beta,Q2), with LHAPDF xfx divided by beta
  double PartonDensity(int pid, double beta, double Q2) const;

  // Compute the factorized diffractive density
  double DiffractiveDensity(int pid, double xi, double beta, double t, double Q2) const;

  // Compute x_hard = xi beta for event-record bookkeeping
  double HardX(double xi, double beta) const;

  // Compute max(hard_scale2, Q2_min) for DPDF evaluation
  double FactorizationQ2(double hard_scale2) const;

  // Compute alpha_s from the same DPDF member used for diffractive densities
  double AlphaS(double Q2) const;

  // Compute whether one proton PDF uses the same alpha_s evolution
  bool MatchAlphaS(const LHAPDF::PDF &pdf, double Q2) const;

  // Validate PDF evolution across the requested hard-scale interval before sampling
  void ValidateAlphaS(const LHAPDF::PDF &pdf, double Q2_min, double Q2_max) const;

  // Map one unit random variable to the configured xi range
  double MapXi(double unit) const;

  // Map one unit random variable to the configured beta range
  double MapBeta(double unit) const;

  // Map one unit random variable to the configured t range
  double MapT(double unit) const;

  // Compute the exponential t slope of the flux at fixed xi
  double FluxTSlope(double xi) const;

  // Access the configured xi support
  const std::vector<double> &XiRange() const;

  // Access the configured beta support
  const std::vector<double> &BetaRange() const;

  // Access the configured t support
  const std::vector<double> &TRange() const;

  // Compute the integration volume for one diffractive proton leg
  double DomainVolume() const;

  // Compute the configured hard-process parton flavours
  const std::vector<int> &PartonFlavours() const;

  // Compute the ordinary proton remnant mass in GeV
  double RemnantMass() const;

  // Print the configured DPDF domain and set metadata
  void PrintSummary(std::ostream &os) const;

 private:
  // Read and validate parameters from one exact GENERAL JSON text
  void ConfigureFromJson(const std::string &source_file,
                                  const std::string &json_text);

  // Initialize the shared read-only LHAPDF member
  void InitPDF();

  // Initialize the PDF member through one run owned store
  void InitPDF(MLHAPDFStore &pdf_store);

  // Test whether a value is inside an inclusive two-element range
  bool InRange(double x, const std::vector<double> &range) const;

  MHardPomeronPDFParam               param;
  std::shared_ptr<const LHAPDF::PDF> dpdf = nullptr;
};

// Process-wide cache for read-only hard-Pomeron PDF wrappers
class MHardPomeronPDFStore {
 public:
  using HardPomeronPtr = std::shared_ptr<const MHardPomeronPDF>;

  // Construct a standalone wrapper store with its own LHAPDF store
  MHardPomeronPDFStore();

  // Construct a run owned wrapper store sharing one LHAPDF store
  explicit MHardPomeronPDFStore(MLHAPDFStore &pdf_store);

  // Compute one shared read-only hard-Pomeron PDF wrapper
  HardPomeronPtr GetHardPomeronPDF(const std::string &modelfile);

  // Compute a wrapper from one already captured immutable model snapshot
  HardPomeronPtr GetHardPomeronPDF(const SoftModelPtr &soft_model);

 private:
  struct Key {
    std::string modelfile;
    std::string model_json;

    // Provide strict weak ordering for cache lookup
    bool operator<(const Key &other) const;
  };

  // Build the full cache key from one immutable model snapshot
  Key MakeKey(const SoftModelPtr &soft_model) const;

  // Construct one initialized wrapper from an immutable model snapshot
  HardPomeronPtr LoadHardPomeronPDF(const SoftModelPtr &soft_model) const;

  MFixedStore<Key, MHardPomeronPDF> store;
  std::unique_ptr<MLHAPDFStore> owned_pdf;
  MLHAPDFStore *pdf_store = nullptr;
};

}  // namespace gra

#endif
