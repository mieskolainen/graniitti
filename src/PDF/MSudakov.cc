// Shuvaev PDF and Sudakov suppression class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cerrno>
#include <compare>
#include <cstdio>
#include <cstring>
#include <complex>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <system_error>
#include <tuple>
#include <thread>
#include <utility>
#include <vector>

// C file processing
#include <fcntl.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/PDF/MSudakov.h"
#include "Graniitti/Tech/MTimer.h"
#include "Graniitti/Tech/MJsonZip.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;

namespace gra {

using math::msqrt;
using math::pow2;
using math::zi;

namespace {

// Compute the explicit cache schema for numerical formula changes
std::string SudakovCacheVersion() { return "MSudakov-cache-v1-jsonzip"; }

// Compute the PDF-library and set-data identity affecting all tabulated values
std::string SudakovPDFIdentity(const std::string &pdfset, int dataversion) {
  return pdfset + "|LHAPDF=" + LHAPDF::version() + "|data=" + std::to_string(dataversion);
}

// Hold an advisory interprocess lock for one cache filename
class SudakovCacheFileLock {
 public:
  // Acquire an exclusive lock shared by all writers of one cache file
  explicit SudakovCacheFileLock(const std::string &filename) {
    const std::string lockfile = filename + ".lock";
    fd = ::open(lockfile.c_str(), O_CREAT | O_RDWR, 0664);
    if (fd < 0) {
      // Capture errno in an owning exception without using strerror's shared buffer
      const int error = errno;
      throw std::system_error(error, std::generic_category(), "SudakovCacheFileLock: cannot open " + lockfile);
    }

    while (::flock(fd, LOCK_EX) != 0) {
      if (errno == EINTR) { continue; }
      // Capture the error before close can change errno and keep its text local to this exception
      const int code = errno;
      const std::system_error error(code, std::generic_category(), "SudakovCacheFileLock: cannot lock " + lockfile);
      ::close(fd);
      fd = -1;
      throw error;
    }
  }

  // Release the advisory cache-file lock
  ~SudakovCacheFileLock() {
    if (fd >= 0) {
      ::flock(fd, LOCK_UN);
      ::close(fd);
    }
  }

  SudakovCacheFileLock(const SudakovCacheFileLock &) = delete;
  SudakovCacheFileLock &operator=(const SudakovCacheFileLock &) = delete;

 private:
  int fd = -1;
};

// Compute a process-and-thread-unique temporary filename beside the final cache file
std::string SudakovTemporaryFilename(const std::string &filename) {
  std::ostringstream suffix;
  suffix << filename << ".tmp." << static_cast<long long>(::getpid()) << "."
         << std::this_thread::get_id();
  return suffix.str();
}

// Remove only the owned temporary cache file unless it was atomically published
class SudakovTemporaryFileGuard {
 public:
  // Start tracking one not-yet-published temporary cache file
  explicit SudakovTemporaryFileGuard(std::string filename_in) : filename(std::move(filename_in)) {}

  // Remove a failed or abandoned temporary cache file
  ~SudakovTemporaryFileGuard() {
    if (!published) { ::unlink(filename.c_str()); }
  }

  // Mark the temporary file as successfully published
  void Release() { published = true; }

  SudakovTemporaryFileGuard(const SudakovTemporaryFileGuard &) = delete;
  SudakovTemporaryFileGuard &operator=(const SudakovTemporaryFileGuard &) = delete;

 private:
  std::string filename;
  bool        published = false;
};

// Compute a scale-aware floating-point comparison tolerance
double CacheCoordinateTolerance(double expected) {
  return 128.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(expected));
}

// Check a physical value against one interpolation axis with boundary roundoff tolerance
bool IsInsideGridDomain(double value, double minimum, double maximum, bool logarithmic) {
  if (!std::isfinite(value) || (logarithmic && !(value > 0.0))) { return false; }
  const double coordinate = logarithmic ? std::log(value) : value;
  const double lower      = logarithmic ? std::log(minimum) : minimum;
  const double upper      = logarithmic ? std::log(maximum) : maximum;
  const double tolerance =
      std::max(CacheCoordinateTolerance(lower), CacheCoordinateTolerance(upper));
  return coordinate >= lower - tolerance && coordinate <= upper + tolerance;
}

// Integrate each fixed-flavour radiation interval without sampling across a threshold
// [REFERENCE: Coughlin, Forshaw, arxiv.org/abs/0912.3280v2, Eq. (74)]
template <typename F>
double SudakovExponent(double q2, double mu2, double charm2, double bottom2, unsigned int intervals, F&& integrand) {
  std::vector<double> bounds = {q2, mu2};
  for (const double threshold : {charm2, bottom2}) {
    if (threshold > q2 && threshold < mu2) { bounds.push_back(threshold); }
  }
  std::sort(bounds.begin(), bounds.end());
  const double range = std::log(mu2 / q2);
  double       sum   = 0.0;
  for (std::size_t b = 1; b < bounds.size(); ++b) {
    const double width = std::log(bounds[b] / bounds[b - 1]);
    const auto   count = 2U * std::max(1U, static_cast<unsigned int>(std::ceil(intervals * width / (2.0 * range))));
    const double step  = width / count;
    const double scale = std::sqrt(bounds[b] * bounds[b - 1]);
    std::vector<double> values(count + 1);
    for (const auto& i : indices(values)) {
      const double kt2 = std::clamp(bounds[b - 1] * std::exp(i * step), bounds[b - 1], bounds[b]);
      values[i] = integrand(kt2, scale);
    }
    sum += math::CS13Integral(values, step);
  }
  return sum;
}

// Reconstruct the veto and its derivatives from K = -log(T)/log^2(mu^2/Q^2)
std::pair<double, double> SudakovFromRadiator(double q2, double mu, const std::pair<double, double>& radiator,
                                              double curvature = 0.0, double* second = nullptr) {
  if (second != nullptr) { *second = 0.0; }
  if (q2 >= math::pow2(mu)) { return {1.0, 0.0}; }
  const double L     = std::log(math::pow2(mu) / q2);
  const auto [K, dK] = radiator;
  const double T     = std::exp(-L * L * K);
  const double slope = 2.0 * L * K / q2 - L * L * dK;
  if (second != nullptr) {
    *second = T * (slope * slope - 2.0 * (1.0 + L) * K / (q2 * q2) + 4.0 * L * dK / q2 - L * L * curvature);
  }
  return {T, T * slope};
}

// Enforce the empty radiation interval without a discontinuous table splice
std::pair<double, double> SudakovVeto(const IArray2D& grid, double q2, double mu, double* second = nullptr) {
  if (second != nullptr) { *second = 0.0; }
  if (q2 >= math::pow2(mu)) { return {1.0, 0.0}; }
  double     curvature = 0.0;
  const auto radiator  = grid.Interpolate2D(q2, mu, second != nullptr ? &curvature : nullptr);
  return SudakovFromRadiator(q2, mu, radiator, curvature, second);
}

// Check the rectangular storage shape required by one interpolation array
bool HasValidArrayShape(const IArray2D &arr) {
  if (arr.N[0] < 2 || arr.N[1] < 2 || arr.F.size() != arr.N[0] + 1 ||
      arr.first_axis_physical.size() != arr.N[0] + 1) {
    return false;
  }
  for (const auto &row : arr.F) {
    if (row.size() != arr.N[1] + 1) { return false; }
  }
  return true;
}

// Describe the exact interpolation array using the common JSON cache conventions
nlohmann::json ArrayMetadata(const IArray2D &arr) {
  auto cache = arr.Metadata();
  cache.update({{"version", 1}, {"type", "IArray2D"}, {"key", arr.GetHashString()},
                {"shape", {arr.N[0] + 1, arr.N[1] + 1, 4}}});
  return cache;
}

}  // namespace

// Constructor
MSudakov::MSudakov() {}

// Destructor
MSudakov::~MSudakov() {}

// Construct a standalone Sudakov store with one private PDF store
MSudakovStore::MSudakovStore()
    : owned_pdf(std::make_unique<MLHAPDFStore>()),
      pdf_store(owned_pdf.get()) {}

// Construct a Sudakov store sharing the run owned PDF store
MSudakovStore::MSudakovStore(MLHAPDFStore &pdf_store_in)
    : pdf_store(&pdf_store_in) {}

// Provide strict weak ordering for cache lookup
bool MSudakovStore::Key::operator<(const Key &other) const {
  if (soft_model.owner_before(other.soft_model)) { return true; }
  if (other.soft_model.owner_before(soft_model)) { return false; }
  return std::tie(sqrts, pdf_identity) <
         std::tie(other.sqrts, other.pdf_identity);
}

// Compute one shared read-only Sudakov/Shuvaev table set
MSudakovStore::SudakovPtr MSudakovStore::GetSudakov(double sqrts, const std::string &pdfset,
                                                    const SoftModelPtr &soft_model) {
  if (!std::isfinite(sqrts) || !(sqrts > 0.0)) {
    throw std::invalid_argument("MSudakovStore::GetSudakov: invalid sqrt(s) " +
                                std::to_string(sqrts));
  }
  if (pdfset.empty() || pdfset == "null") {
    throw std::invalid_argument("MSudakovStore::GetSudakov: invalid PDF set name '" + pdfset + "'");
  }

  if (!soft_model) {
    throw std::invalid_argument(
        "MSudakovStore::GetSudakov: missing SOFT model");
  }
  const auto pdf = pdf_store->GetPDF(pdfset, 0);
  const std::string pdf_identity =
      SudakovPDFIdentity(pdfset, pdf->dataversion());
  const Key key = MakeKey(sqrts, pdf_identity, soft_model);
  return store.GetOrLoad(key, [this, sqrts, &pdfset, &soft_model, &pdf] {
    try {
      SudakovPtr sudakov = LoadSudakov(sqrts, pdfset, soft_model, pdf);
      std::cout << "MSudakovStore::LoadSudakov SUCCESS\n";
      return sudakov;
    } catch (const std::exception &e) {
      std::cerr << "MSudakovStore::LoadSudakov FAILED: "
                << e.what() << '\n';
      throw;
    }
  });
}

// Build the full cache key from runtime energy and current numerics
MSudakovStore::Key MSudakovStore::MakeKey(
    double sqrts, const std::string &pdf_identity,
    const SoftModelPtr &soft_model) const {
  return {soft_model, sqrts, pdf_identity};
}

// Construct one fully initialized Sudakov object
MSudakovStore::SudakovPtr MSudakovStore::LoadSudakov(
    double sqrts, const std::string &pdfset,
    const SoftModelPtr &soft_model, const MLHAPDFStore::PDFPtr &pdf) const {
  if (!soft_model->HasNumericsSnapshot()) {
    throw std::invalid_argument(
        "MSudakovStore::LoadSudakov: missing immutable NUMERICS snapshot");
  }
  std::shared_ptr<MSudakov> sudakov = std::make_shared<MSudakov>();
  sudakov->Init(sqrts, pdfset, soft_model,
                soft_model->NumericsSourceFile(),
                soft_model->NumericsSourceJson(), true, pdf);
  return sudakov;
}

// Init
void MSudakov::Init(double _sqrts, const std::string &PDFSET,
                    const SoftModelPtr &soft_model_snapshot,
                    const std::string &numerics_source,
                    const std::string &numerics_json, bool init_arrays,
                    std::shared_ptr<const LHAPDF::PDF> pdf) {
  if (!soft_model_snapshot) {
    throw std::invalid_argument("MSudakov::Init: missing SOFT model");
  }
  soft_model = soft_model_snapshot;
  Model.ConfigureFromJson(soft_model->SourceFile(),
                                   soft_model->SourceJson());
  Numerics.ConfigureFromJson(numerics_source, numerics_json);

  InitLHAPDF(PDFSET, std::move(pdf));
  const double q_match = Model.Q0;
  if (q_match < std::sqrt(pdf_q2_MIN)) {
    throw std::invalid_argument("MSudakov::Init: Q0 must be at least the LHAPDF QMin");
  }
  Numerics.q2_MIN      = pow2(q_match);
  Numerics.mu_MIN      = q_match;
  if (!(Numerics.q2_MIN < Numerics.q2_MAX)) {
    throw std::invalid_argument(
        "MSudakov::Init: physical matching scale must be below q2_MAX");
  }

  if (!std::isfinite(_sqrts) || _sqrts <= Numerics.mu_MIN) {
    throw std::invalid_argument("MSudakov::Init: sqrt(s) must be finite and larger than mu_MIN");
  }

  // Setup energy dependent boundaries
  Numerics.sqrts = _sqrts;
  Numerics.mu_MAX = _sqrts;

  // Use the lower perturbative scale to set the minimum x grid value
  Numerics.x_MIN = math::pow2(Numerics.mu_MIN / Numerics.sqrts);

  if (init_arrays) { InitArrays(); }
}

// Initialize interpolation arrays and write missing cache files
void MSudakov::InitArrays() {
  const std::string pdf_identity = SudakovPDFIdentity(PDFSETNAME, PdfPtr->dataversion());

  // Init Sudakov variables
  {
    std::cout << "Initializing <Sudakov> array:" << std::endl;

    // q2,mu
    veto.sqrts = Numerics.sqrts;  // FIRST THIS
    veto.pdf = PDFSETNAME;
    veto.quantity = "sudakov_radiator";
    veto.Set(0, "q2", Numerics.q2_MIN, Numerics.q2_MAX, Numerics.SUDA_N[0],
             Numerics.SUDA_log_ON[0]);
    veto.Set(1, "mu", Numerics.mu_MIN, Numerics.mu_MAX, Numerics.SUDA_N[1],
             Numerics.SUDA_log_ON[1]);
    veto.InitArray();  // Initialize (call last!)

    const unsigned long hash = gra::aux::djb2hash(
        veto.GetHashString() + std::to_string(Numerics.config_hash) + SudakovCacheVersion() + "|fixed-flavour-v1|" +
        pdf_identity);
    const std::string   filename = gra::aux::GetBasePath(2) + "/sudakov/SUDA_" +
                                 gra::aux::ToString(veto.sqrts, 0) + "_" + PDFSETNAME + "_" +
                                 std::to_string(hash) + ".json";

    // Pointer to member function: ReturnType
    // (ClassType::*)(ParameterTypes...)
    std::pair<double, double> (MSudakov::*f)(double, double) = &MSudakov::SudakovRadiator;
    LoadOrBuildArray(veto, filename, f);
  }

  // Init Shuvaev PDF variables
  {
    std::cout << "Initializing <Shuvaev> array:" << std::endl;

    // q2,x
    spdf.sqrts = Numerics.sqrts;  // FIRST THIS
    spdf.pdf = PDFSETNAME;
    spdf.quantity = "shuvaev";
    spdf.Set(0, "q2", Numerics.q2_MIN, Numerics.q2_MAX, Numerics.SHUV_N[0],
             Numerics.SHUV_log_ON[0]);
    spdf.Set(1, "x", Numerics.x_MIN, Numerics.x_MAX, Numerics.SHUV_N[1], Numerics.SHUV_log_ON[1]);
    spdf.InitArray();  // Initialize (call last!)

    const unsigned long hash = gra::aux::djb2hash(
        spdf.GetHashString() + std::to_string(Numerics.config_hash) + SudakovCacheVersion() +
        pdf_identity);
    const std::string   filename = gra::aux::GetBasePath(2) + "/sudakov/SHUV_" +
                                 gra::aux::ToString(spdf.sqrts, 0) + "_" + PDFSETNAME + "_" +
                                 std::to_string(hash) + ".json";

    // Pointer to member function: ReturnType
    // (ClassType::*)(ParameterTypes...)
    std::pair<double, double> (MSudakov::*f)(double, double) = &MSudakov::Shuvaev_H;
    LoadOrBuildArray(spdf, filename, f);
  }
  // --------------------------------------------------
  // Finally
  initialized = true;
  std::cout << std::endl;

  if (Numerics.DEBUG) { TestPDF(); }
}

// Compute gluon xg(x,Q2) directly from LHAPDF
double MSudakov::xg_xQ2_raw(double x, double q2) const {
  const int pid = 21;  // gluon
  try {
    return PdfPtr->xfxQ2(pid, x, q2);
  } catch (const std::exception &error) {
    throw std::invalid_argument("MSudakov::xg_xQ2_raw: Problem with x = " + std::to_string(x) +
                                ", q2 = " + std::to_string(q2) + ": " + error.what());
  }
}

// Compute gluon xg(x,Q2) inside the selected PDF domain
double MSudakov::xg_xQ2(double x, double q2) const {
  if (q2 < pdf_q2_MIN * (1.0 - 1e-12)) {
    throw std::out_of_range("MSudakov::xg_xQ2: Q2 lies below the LHAPDF starting scale");
  }
  return xg_xQ2_raw(x, q2);
}

// PDF access
//
// [REFERENCE: LHAPDF6, arxiv.org/abs/1412.7420]
void MSudakov::InitLHAPDF(const std::string &pdfname,
                          std::shared_ptr<const LHAPDF::PDF> pdf) {
  PDFSETNAME = pdfname;
  if (pdf == nullptr) {
    MLHAPDFStore pdf_store;
    pdf = pdf_store.GetPDF(pdfname, 0);
  }
  PdfPtr = std::move(pdf);

  const double q_min = PdfPtr->info().get_entry_as<double>("QMin", 0.0);
  if (!(q_min > 0.0) || !std::isfinite(q_min)) {
    throw std::invalid_argument("MSudakov::InitLHAPDF: PDF metadata must provide positive QMin");
  }
  pdf_q2_MIN = math::pow2(q_min);

  charm_mass  = PdfPtr->info().get_entry_as<double>("MCharm", 0.0);
  bottom_mass = PdfPtr->info().get_entry_as<double>("MBottom", 0.0);
  if (!(charm_mass > 0.0) || !(bottom_mass > charm_mass) ||
      !std::isfinite(charm_mass) || !std::isfinite(bottom_mass)) {
    throw std::invalid_argument(
        "MSudakov::InitLHAPDF: PDF metadata must provide ordered MCharm and MBottom");
  }
}

// PDF test routine
void MSudakov::TestPDF() const {
  const double MINLOGX = std::log10(Numerics.x_MIN);
  const double MAXLOGX = std::log10(Numerics.x_MAX);
  const int    NX      = 5;  // Number of points - 1
  const double stepX   = (MAXLOGX - MINLOGX) / NX;

  const double MINLOGQ2 = std::log10(Numerics.q2_MIN);
  const double MAXLOGQ2 = std::log10(Numerics.q2_MAX);
  const int    NQ2      = 5;  // Number of points - 1
  const double stepQ2   = (MAXLOGQ2 - MINLOGQ2) / NQ2;

  const double MINLOGMU = std::log10(Numerics.mu_MIN);
  const double MAXLOGMU = std::log10(Numerics.mu_MAX);
  const int    NMU      = 5;  // Number of points - 1
  const double stepMU   = (MAXLOGMU - MINLOGMU) / NMU;

  // Test loop
  for (std::size_t i = 0; i < NMU + 1; ++i) {
    const double log10mu = MINLOGMU + i * stepMU;
    const double mu      = std::pow(10, log10mu);

    printf("[mu = %0.1f GeV] : alpha_s(Q = mu GeV) = %0.3f \n\n", mu,
           PdfPtr->alphasQ2(mu * mu));

    for (std::size_t j = 0; j < NX + 1; ++j) {
      const double log10x = MINLOGX + j * stepX;
      const double x      = std::pow(10, log10x);

      printf("x = %0.5E \n", x);
      for (std::size_t k = 0; k < NQ2 + 1; ++k) {
        const double log10q2 = MINLOGQ2 + k * stepQ2;
        const double q2      = std::pow(10, log10q2);

        // Normal gluon pdf
        const double xf = xg_xQ2(x, q2);

        // Durham flux
        const double hxf = fg_xQ2Mu(x, q2, mu);

        printf(
            "(x = %0.3E, q2 = %0.2f, mu = %0.1f) : [gluon pdf: xg(x,q2), "
            "Durham flux: fg(x,q2,mu)] = (%0.2f,%0.2f) \n",
            x, q2, mu, xf, hxf);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

// Access QCD coupling alpha_s(Q^2) from LHAPDF
double MSudakov::AlphaS_Q2(double q2) const {
  if (q2 < pdf_q2_MIN * (1.0 - 1e-12)) {
    throw std::out_of_range("MSudakov::AlphaS_Q2: Q2 lies below the LHAPDF starting scale");
  }

  try {
    return PdfPtr->alphasQ2(q2);
  } catch (const std::exception &error) {
    throw std::invalid_argument("MSudakov::AlphaS_Q2: Problem with q2 = " +
                                std::to_string(q2) + ": " + error.what());
  }
}

// Calculate differentiation dxg/dQ2 via "Richardson's extrapolation"
// f'(x) = [4*D_0(h) - D_0(2h)] / 3 + O(h^4)
//
// Requires 4 evaluations of densities
//
// [REFERENCE: en.wikipedia.org/wiki/Richardson_extrapolation]
double MSudakov::diff_xg_xQ2_wrt_Q2(double x, double q2) const {
  const double h   = std::max(1E-5, 1E-4 * q2);
  const double hX2 = 2 * h;

  if (q2 - hX2 < pdf_q2_MIN) {
    const double f0 = xg_xQ2(x, q2);
    const double f1 = xg_xQ2(x, q2 + h);
    const double f2 = xg_xQ2(x, q2 + 2.0 * h);
    const double f3 = xg_xQ2(x, q2 + 3.0 * h);
    const double f4 = xg_xQ2(x, q2 + 4.0 * h);
    return (-25.0 * f0 + 48.0 * f1 - 36.0 * f2 + 16.0 * f3 - 3.0 * f4) /
           (12.0 * h);
  }

  const double D0_A = (xg_xQ2(x, q2 + h) - xg_xQ2(x, q2 - h)) / (2 * h);
  const double D0_B = (xg_xQ2(x, q2 + hX2) - xg_xQ2(x, q2 - hX2)) / (2 * hX2);

  return (4.0 * D0_A - D0_B) / 3.0;
}

// Durham flux (skewed gluon pdf)
// ~ Shuvaev transformed gluon pdf x Sudakov suppression
//
// f_g(x,x',qt^2,\mu)
// = \frac{\partial}{\partial \ln Q_t^2} [H_g(x/2,x/2,Q_t^2)\sqrt{T(Q_t^2,\mu)}]
//
double MSudakov::FluxDerivative(double q2, double Hg, double dHg, double Tg,
                                double dTg) const {
  if (Tg <= 1e-15) { return 0.0; }

  const double sqrt_tg = math::msqrt(Tg);
  return q2 * (dHg * sqrt_tg + Hg * dTg / (2.0 * sqrt_tg));
}

// Compute the perturbative integrated gluon and its logarithmic derivative
// G = H_g sqrt(T), f_g = dG/dln(Q_t^2)
MSudakov::MatchedFlux MSudakov::PerturbativeFlux(double x, double q2, double mu) const {
  if (!IsInsideGridDomain(q2, Numerics.q2_MIN, Numerics.q2_MAX, Numerics.SHUV_log_ON[0]) ||
      !IsInsideGridDomain(x, Numerics.x_MIN, Numerics.x_MAX, Numerics.SHUV_log_ON[1]) ||
      !IsInsideGridDomain(q2, Numerics.q2_MIN, Numerics.q2_MAX, Numerics.SUDA_log_ON[0]) || !SupportsMu(mu)) {
    return {};
  }

  const auto [Hg, dHg] = spdf.Interpolate2D(q2, x);
  const auto [Tg, dTg] = SudakovVeto(veto, q2, mu);
  if (!(Tg > 1e-15)) { return {}; }

  return {Hg * math::msqrt(Tg), FluxDerivative(q2, Hg, dHg, Tg, dTg)};
}

// Evaluate the signed matched Durham flux derivative
//
// PERTURBATIVE_ONLY treats q2 < Q0^2 as an excluded integration domain;
// the zero returned there is not a distributional continuation of G
double MSudakov::fg_xQ2Mu(double x, double q2, double mu) const {
  if (!(q2 > 0.0) || q2 > Numerics.q2_MAX ||
      !IsInsideGridDomain(x, Numerics.x_MIN, Numerics.x_MAX, Numerics.SHUV_log_ON[1]) || !SupportsMu(mu)) {
    return 0.0;
  }
  return (q2 >= Numerics.q2_MIN ? PerturbativeFlux(x, q2, mu)
                                : InfraredFlux(x, q2, mu)).flux;
}

// Compute the matched integrated amplitude-level gluon H_g sqrt(T)
double MSudakov::IntegratedFlux_xQ2Mu(double x, double q2, double mu) const {
  if (!(q2 > 0.0) || q2 > Numerics.q2_MAX ||
      !IsInsideGridDomain(x, Numerics.x_MIN, Numerics.x_MAX, Numerics.SHUV_log_ON[1]) || !SupportsMu(mu)) {
    return 0.0;
  }
  if (Model.Mode() == MSudakovIRMode::PERTURBATIVE_ONLY &&
      q2 < Numerics.q2_MIN) {
    throw std::domain_error(
        "MSudakov::IntegratedFlux_xQ2Mu: integrated flux is undefined below the strict perturbative domain");
  }
  return (q2 >= Numerics.q2_MIN ? PerturbativeFlux(x, q2, mu)
                                : InfraredFlux(x, q2, mu)).integrated;
}

// Prepare one event-local fixed-x and fixed-mu Durham flux evaluator
MSudakov::PreparedFlux MSudakov::PrepareFlux(double x, double mu) const {
  PreparedFlux prepared;
  prepared.owner = this;
  prepared.x = x;
  prepared.mu = mu;
  prepared.q2_min = Numerics.q2_MIN;
  prepared.q2_max = Numerics.q2_MAX;
  if (!IsInsideGridDomain(x, Numerics.x_MIN, Numerics.x_MAX, Numerics.SHUV_log_ON[1]) || !SupportsMu(mu)) {
    throw std::domain_error("MSudakov::PrepareFlux: x or mu is outside the gluon table domain");
  }

  prepared.shuvaev_x = spdf.PrepareSecondCoordinate(x);
  if (mu > Numerics.mu_MIN) { prepared.sudakov_mu = veto.PrepareSecondCoordinate(mu); }
  prepared.aligned_first_axis =
      spdf.islog[0] == veto.islog[0] && spdf.N[0] == veto.N[0] &&
      std::is_eq(spdf.MIN[0] <=> veto.MIN[0]) &&
      std::is_eq(spdf.MAX[0] <=> veto.MAX[0]) &&
      std::is_eq(spdf.STEP[0] <=> veto.STEP[0]);
  prepared.valid = true;
  return prepared;
}

// Populate one direct-mapped cell with fixed-x and fixed-mu interpolation data
MSudakov::PreparedFlux::PerturbativeCell &
MSudakov::PreparedFlux::PreparePerturbativeCell(int index) {
  PerturbativeCell &cell =
      perturbative_cache[static_cast<std::size_t>(index) % perturbative_cache.size()];
  if (cell.index == index) { return cell; }

  const auto prepare_field =
      [index](const IArray2D &grid,
              const IArray2D::PreparedSecondCoordinate &prepared_second) {
        InterpolationField field;
        const int j = prepared_second.index;
        const double ty = prepared_second.fraction;
        const double one_minus_ty = 1.0 - ty;
        field.value_1 =
            one_minus_ty * grid.F[index][j][2] +
            ty * grid.F[index][j + 1][2];
        field.value_2 =
            one_minus_ty * grid.F[index + 1][j][2] +
            ty * grid.F[index + 1][j + 1][2];
        field.slope_1 =
            (one_minus_ty * grid.F[index][j][3] +
             ty * grid.F[index][j + 1][3]) *
            (grid.islog[0] ? grid.first_axis_physical[index] : 1.0);
        field.slope_2 =
            (one_minus_ty * grid.F[index + 1][j][3] +
             ty * grid.F[index + 1][j + 1][3]) *
            (grid.islog[0] ? grid.first_axis_physical[index + 1] : 1.0);
        return field;
      };

  cell.index = index;
  cell.coordinate_1 = owner->spdf.F[index][shuvaev_x.index][0];
  cell.inverse_step =
      1.0 / (owner->spdf.F[index + 1][shuvaev_x.index][0] -
             cell.coordinate_1);
  cell.shuvaev = prepare_field(owner->spdf, shuvaev_x);
  if (mu > owner->Numerics.mu_MIN) { cell.sudakov = prepare_field(owner->veto, sudakov_mu); }
  return cell;
}

// Evaluate one cached cubic-Hermite field with the original interpolation algebra
std::pair<double, double>
MSudakov::PreparedFlux::EvaluatePerturbativeCell(
    const PerturbativeCell &cell, const InterpolationField &field,
    double input_q2, double coordinate_q2) const {
  const double tx = (coordinate_q2 - cell.coordinate_1) * cell.inverse_step;
  const double tx2 = tx * tx;
  const double tx3 = tx2 * tx;
  const double h00 = 2.0 * tx3 - 3.0 * tx2 + 1.0;
  const double h10 = tx3 - 2.0 * tx2 + tx;
  const double h01 = -2.0 * tx3 + 3.0 * tx2;
  const double h11 = tx3 - tx2;
  const double dh00 = 6.0 * tx2 - 6.0 * tx;
  const double dh10 = 3.0 * tx2 - 4.0 * tx + 1.0;
  const double dh01 = -6.0 * tx2 + 6.0 * tx;
  const double dh11 = 3.0 * tx2 - 2.0 * tx;
  const double step = 1.0 / cell.inverse_step;
  const double value =
      h00 * field.value_1 + h10 * step * field.slope_1 +
      h01 * field.value_2 + h11 * step * field.slope_2;
  const double derivative_coordinate =
      (dh00 * field.value_1 + dh10 * step * field.slope_1 +
       dh01 * field.value_2 + dh11 * step * field.slope_2) /
      step;
  const double derivative =
      derivative_coordinate / (owner->spdf.islog[0] ? input_q2 : 1.0);
  return {value, derivative};
}

// Compute f_g/Q2 from one event-local prepared Durham flux evaluator
double MSudakov::PreparedFlux::OverQ2(double q2) {
  if (!valid || owner == nullptr || !std::isfinite(q2) || q2 < 0.0 ||
      q2 > q2_max) {
    return 0.0;
  }

  if (q2 >= q2_min) {
    const double spdf_coordinate =
        owner->spdf.islog[0] ? std::log(q2) : q2;
    const double veto_coordinate =
        owner->veto.islog[0] == owner->spdf.islog[0]
            ? spdf_coordinate
            : (owner->veto.islog[0] ? std::log(q2) : q2);
    std::pair<double, double> shuvaev;
    std::pair<double, double> sudakov;
    if (aligned_first_axis) {
      int index =
          std::floor((spdf_coordinate - owner->spdf.MIN[0]) /
                     owner->spdf.STEP[0]);
      if (index < 0) { index = 0; }
      if (index >= static_cast<int>(owner->spdf.N[0])) {
        index = owner->spdf.N[0] - 1;
      }
      const PerturbativeCell &cell = PreparePerturbativeCell(index);
      shuvaev = EvaluatePerturbativeCell(
          cell, cell.shuvaev, q2, spdf_coordinate);
      if (q2 >= math::pow2(mu)) { return shuvaev.second; }
      sudakov = EvaluatePerturbativeCell(
          cell, cell.sudakov, q2, spdf_coordinate);
    } else {
      shuvaev = owner->spdf.InterpolatePreparedInDomain(
          q2, spdf_coordinate, shuvaev_x);
      if (q2 >= math::pow2(mu)) { return shuvaev.second; }
      sudakov = owner->veto.InterpolatePreparedInDomain(
          q2, veto_coordinate, sudakov_mu);
    }
    const auto [Hg, dHg] = shuvaev;
    const auto [Tg, dTg] = SudakovFromRadiator(q2, mu, sudakov);
    if (!(Tg > 1e-15)) { return 0.0; }
    const double sqrt_tg = math::msqrt(Tg);
    return dHg * sqrt_tg + Hg * dTg / (2.0 * sqrt_tg);
  }

  if (!infrared_ready) { PrepareInfrared(); }
  return InfraredOverQ2(q2);
}

// Compute the gluon kernel with alpha_s frozen at Q0 below the boundary
//
// The fixed coupling continuation is outside the coherent RKHS variation
double MSudakov::AlphaSFlux_xQ2Mu(double x, double q2, double mu,
                                  double alpha_scale2) const {
  if (!(alpha_scale2 > 0.0) || !std::isfinite(alpha_scale2)) { return 0.0; }
  const double matched_scale2 = std::max(alpha_scale2, Numerics.q2_MIN);
  const double flux = fg_xQ2Mu(x, q2, mu);
  return AlphaS_Q2(matched_scale2) * flux;
}

// Compute the perturbative Shuvaev-transformed gluon distribution
double MSudakov::PerturbativeShuvaevGluon_xQ2(double x, double q2) const {
  return spdf.Interpolate2D(q2, x).first;
}

// Compute the perturbative q2 derivative of the Shuvaev gluon
double MSudakov::PerturbativeShuvaevGluonDerivative_xQ2(double x, double q2) const {
  return spdf.Interpolate2D(q2, x).second;
}

// Compute the perturbative square-root Sudakov factor
double MSudakov::PerturbativeSqrtSudakov_Q2Mu(double q2, double mu) const {
  if (!IsInsideGridDomain(q2, Numerics.q2_MIN, Numerics.q2_MAX, Numerics.SUDA_log_ON[0]) || !SupportsMu(mu)) {
    return 0.0;
  }
  const double Tg = SudakovVeto(veto, q2, mu).first;
  return (Tg > 0.0) ? math::msqrt(Tg) : 0.0;
}

// Calculate Shuvaev integral transform
// ----------------------------------------------------------------------
// Identity:
//
// H_g(x, \xi -> 0) = xg(x)
//
// ----------------------------------------------------------------------
//
// Hg  = numerical integral transformed from standard pdf
// dHg = dHg/dq^2 (differentiated numerically)
//
// [REFERENCE: Harland-Lang, arxiv.org/abs/1306.6661]
//
std::pair<double, double> MSudakov::Shuvaev_H(double q2, double x) {
  return ShuvaevTransform(q2, x);
}

// Compute the full Shuvaev transform and its q2 derivative
std::pair<double, double> MSudakov::ShuvaevTransform(double q2, double x) {

  const double y_MIN  = x / 4.0;
  const double y_MAX  = 1.0;
  const double y_STEP = (y_MAX - y_MIN) / Numerics.ShuvaevIntegralN;

  double Hg  = 0.0;
  double dHg = 0.0;

  // Check that we are within valid domain (take into account floating points)
  const double EPS = 1e-5;
  if (x >= Numerics.x_MIN * (1 - EPS) && x <= Numerics.x_MAX * (1 + EPS) &&
      q2 >= Numerics.q2_MIN * (1 - EPS) && q2 <= Numerics.q2_MAX * (1 + EPS)) {
    // N+1!
    std::vector<double> fA(Numerics.ShuvaevIntegralN + 1, 0.0);
    std::vector<double> fB(Numerics.ShuvaevIntegralN + 1, 0.0);

    for (const auto &i : indices(fA)) {
      const double y        = y_MIN + i * y_STEP;
      const double argument = x / (4.0 * y);

      // H_g(x/2,x/2,Q^2) = 4x/\pi \int_{x/4}^1 dy
      // y^{1/2}(1-y)^{1/2}g(x/4y),q^2)
      //
      // => take into account that LHAPDF provides xg(), not g(), gives:
      const double factor = math::msqrt(math::pow3(y) * (1 - y));
      fA[i]               = factor * xg_xQ2(argument, q2);
      fB[i]               = factor * diff_xg_xQ2_wrt_Q2(argument, q2);
    }
    const double norm = 16.0 / math::PI;
    Hg                = norm * math::CS13Integral(fA, y_STEP);
    dHg               = norm * math::CS13Integral(fB, y_STEP);
  } else {
    // Fatal error
    throw std::invalid_argument("MSudakov::Shuvaev_H(q2,x): Input arguments out of domain: q2 = " +
                                std::to_string(q2) + ", x = " + std::to_string(x));
  }
  return {Hg, dHg};
}

// Compute the Durham soft emission cutoff Delta = kt / mu
// [REFERENCE: Coughlin, Forshaw, arxiv.org/abs/0912.3280v2, Eq. (74)]
double MSudakov::SudakovDelta(double kt2, double mu) const {
  if (!std::isfinite(kt2) || kt2 < 0.0 || !std::isfinite(mu) || mu <= 0.0) {
    throw std::invalid_argument("MSudakov::SudakovDelta: kt2 must be non-negative and mu positive");
  }
  return math::msqrt(kt2) / mu;
}

// Tabulate the smooth radiator K = -log(T)/L^2, L = log(mu^2/Q^2)
// Factoring L^2 imposes T = 1 and dT/dQ^2 = 0 exactly at Q^2 = mu^2
// [REFERENCE: Coughlin, Forshaw, arxiv.org/abs/0912.3280v2, Eq. (74)]
std::pair<double, double> MSudakov::SudakovRadiator(double qt2, double mu) {
  const double L = std::log(pow2(mu) / qt2);
  if (L <= std::sqrt(std::numeric_limits<double>::epsilon())) {
    return {AlphaS_Q2(pow2(mu)) * (6.0 + 0.5 * NumFlavor(pow2(mu))) / (8.0 * math::PI), 0.0};
  }
  const auto integrand = [&](double kt2, double flavour_scale) {
    const double delta = std::min(1.0, std::sqrt(kt2) / mu);
    return AlphaS_Q2(kt2) * (AP_gg(delta) + AP_qg(delta, flavour_scale)) / (2.0 * math::PI);
  };
  const double K = SudakovExponent(qt2, pow2(mu), pow2(charm_mass), pow2(bottom_mass), Numerics.SudakovIntegralN,
                                    integrand) / (L * L);
  return {K, (2.0 * K / L - integrand(qt2, qt2) / (L * L)) / qt2};
}

// Compute the direct Sudakov veto without interpolation
std::pair<double, double> MSudakov::Sudakov_T(double qt2, double mu) const {
  if (!(qt2 > 0.0) || !SupportsMu(mu)) { throw std::domain_error("MSudakov::Sudakov_T: invalid scales"); }
  if (qt2 >= pow2(mu)) { return {1.0, 0.0}; }
  const auto integrand = [&](double kt2, double flavour_scale) {
    const double delta = std::min(1.0, SudakovDelta(kt2, mu));
    return AlphaS_Q2(kt2) * (AP_gg(delta) + AP_qg(delta, flavour_scale)) / (2.0 * math::PI);
  };
  const double T = std::exp(-SudakovExponent(qt2, pow2(mu), pow2(charm_mass), pow2(bottom_mass),
                                             Numerics.SudakovIntegralN, integrand));
  return {T, T * integrand(qt2, qt2) / qt2};
}

// Differentiate the same interpolants used by the gluon flux at the matching scale
// d f/dln(Q^2) = f + Q^4 d^2[H_g sqrt(T)]/d(Q^2)^2
double MSudakov::PerturbativeFluxLogDerivative(double x, double q2, double mu) const {
  double ddH         = 0.0;
  double ddT         = 0.0;
  const auto [H, dH] = spdf.Interpolate2D(q2, x, &ddH);
  const auto [T, dT] = SudakovVeto(veto, q2, mu, &ddT);
  if (!(T > 1e-15)) { return 0.0; }
  const double root = std::sqrt(T);
  return FluxDerivative(q2, H, dH, T, dT) +
         q2 * q2 * root * (ddH + dH * dT / T + H * (0.5 * ddT / T - 0.25 * pow2(dT / T)));
}

// Altarelli-Parisi splitting function definite integral over z
//
// see standard literature, e.g
// [REFERENCE: B.R. Webber, CERN lectures 08, www.hep.phy.cam.ac.uk/theory/webber/QCDlect3.pdf]
// [REFERENCE: R.K. Ellis, W.J. Stirling, B.R. Webber, QCD and Collider Physics]
// [REFERENCE: http://pdg.lbl.gov/2018/download/db2018.pdf]
//
// \int z P_gg(z) dz
// \int z *{ 2 * CA * int [z/(1-z)_{+} + (1-z)/z + z*(1-z)] +
//           1/6*(11*CA - 4*Nf*TR)*deltafunc(1-z) } dz, z = 0 ... (1-delta) =
//           (below)
//
// First part: ...
// Second part: \int 1/6*(11*C - 4*N*T)*DiracDelta[1-x] dx, x = 0 ... 1 - delta,
// this vanishes with positive delta
//
// where + is the "plus-description"
//
double MSudakov::AP_gg(double delta) const {
  const double CA = 3.0;  // Structure constant
  return 2.0 * CA *
         (std::log(1.0 / delta) -
          math::pow2(1.0 - delta) * (3.0 * math::pow2(delta) - 2.0 * delta + 11.0) / 12.0);
}

// Altarelli-Parisi splitting function definite integral over z,
//
// \int \sum_q P_qz(z) dz
//
// with sum over active quark flavors with mass threshold qt2
// \sum_q TR * \int (z^2 + (1-z)^2) dz, z = 0 ... (1-delta)
//
double MSudakov::AP_qg(double delta, double qt2) const {
  const double TR = 0.5;  // Structure constant
  return TR * (-2.0 * math::pow3(delta) / 3.0 + math::pow2(delta) - delta + 2.0 / 3.0) *
         NumFlavor(qt2);
}

// Compute the number of quark flavors at scale q^2
//
double MSudakov::NumFlavor(double q2) const {
  if (q2 < math::pow2(charm_mass)) {
    return 3.0;
  } else if (q2 < math::pow2(bottom_mass)) {
    return 4.0;
  } else {
    return 5.0;
  }
}

// Constructs interpolation array values
void MSudakov::CalculateArray(IArray2D &arr,
                              std::pair<double, double> (MSudakov::*f)(double, double)) {
  MTimer timer;
  for (const auto &i : indices(arr.F)) {
    const double a = arr.MIN[0] + i * arr.STEP[0];

    // Transform input to linear if log stepping, for the function
    const double var1 = (arr.islog[0]) ? std::exp(a) : a;
    gra::aux::PrintProgress(i / static_cast<double>(arr.N[0] + 1));

    for (const auto &j : indices(arr.F[i])) {
      const double b = arr.MIN[1] + j * arr.STEP[1];

      // Transform input to linear if log stepping, for the function
      const double var2 = (arr.islog[1]) ? std::exp(b) : b;

      // Call function being pointed to
      const std::pair<double, double> output = (this->*f)(var1, var2);

      arr.F[i][j][0] = a;
      arr.F[i][j][1] = b;
      arr.F[i][j][2] = std::abs(output.first) < 1e-64 ? 0 : output.first;    // Underflow protection
      arr.F[i][j][3] = std::abs(output.second) < 1e-64 ? 0 : output.second;  //
    }
  }
  // Progressbar clearing
  gra::aux::ClearProgress();
  printf("MSudakov::CalculateArray: Time elapsed %0.1f sec \n", timer.ElapsedSec());
}

// Read a complete cache or rebuild it once while holding its interprocess lock
void MSudakov::LoadOrBuildArray(
    IArray2D &arr, const std::string &filename,
    std::pair<double, double> (MSudakov::*f)(double, double)) {
  aux::CreateDirectory(std::filesystem::path(filename).parent_path().string());
  if (arr.ReadArray(filename)) { return; }

  SudakovCacheFileLock lock(filename);
  if (arr.ReadArray(filename)) { return; }

  CalculateArray(arr, f);
  arr.WriteArray(filename, true);
  if (!arr.ReadArray(filename)) {
    throw std::runtime_error("MSudakov::LoadOrBuildArray: failed to validate cache " + filename);
  }
}

// Write the array through a same-directory temporary file and atomic rename
bool IArray2D::WriteArray(const std::string &filename, bool overwrite) const {
  if (!HasValidArrayShape(*this)) {
    throw std::invalid_argument("IArray2D::WriteArray: invalid array storage shape");
  }

  // Do not write if file exists already
  if (gra::aux::FileExist(filename) && !overwrite) { return true; }

  aux::CreateDirectory(std::filesystem::path(filename).parent_path().string());
  const std::string temporary = SudakovTemporaryFilename(filename);
  SudakovTemporaryFileGuard temporary_guard(temporary);
  std::ofstream             file(temporary, std::ios::out | std::ios::trunc);
  if (!file.is_open()) {
    std::string str = "IArray2D::WriteArray: Fatal IO-error with: " + temporary;
    throw std::invalid_argument(str);
  }

  std::cout << "IArray2D::WriteArray: " << std::flush;
  std::vector<double> values;
  values.reserve(F.size() * F[0].size() * 4);
  for (const auto &row : F) {
    for (const auto &cell : row) { values.insert(values.end(), cell.begin(), cell.end()); }
  }
  auto cache = ArrayMetadata(*this);
  cache["data"] = MJsonZip::CompressVector(values);
  file << cache.dump() << '\n';

  file.flush();
  if (!file.good()) {
    throw std::runtime_error("IArray2D::WriteArray: failed to flush " + temporary);
  }
  file.close();
  if (file.fail()) {
    throw std::runtime_error("IArray2D::WriteArray: failed to close " + temporary);
  }

  if (::rename(temporary.c_str(), filename.c_str()) != 0) {
    // Keep error text local because different PDF cache writers may fail concurrently
    const int error = errno;
    throw std::system_error(error, std::generic_category(), "IArray2D::WriteArray: failed to publish " + filename);
  }
  temporary_guard.Release();

  std::cout << rang::fg::green << "[DONE]" << rang::fg::reset << std::endl;
  return true;
}

// Read and validate exact compressed cache dimensions, coordinates and finite values
bool IArray2D::ReadArray(const std::string &filename) {
  if (!HasValidArrayShape(*this)) {
    throw std::invalid_argument("IArray2D::ReadArray: invalid array storage shape");
  }

  std::ifstream file(filename);
  if (!file.is_open()) {
    // Missing cache files are rebuilt by the caller
    return false;
  }

  std::cout << "IArray2D::ReadArray: ";
  try {
    nlohmann::json cache;
    file >> cache;
    std::vector<double> values;
    if (!MJsonZip::DecompressVector(cache.at("data"), values, F.size() * F[0].size() * 4)) {
      throw std::runtime_error("invalid compressed array");
    }
    cache.erase("data");
    if (cache != ArrayMetadata(*this)) { throw std::runtime_error("array metadata mismatch"); }
    for (const auto &i : indices(F)) {
      for (const auto &j : indices(F[i])) {
        const auto first = values.begin() + 4 * (i * F[i].size() + j);
        const double expected_a = MIN[0] + i * STEP[0];
        const double expected_b = MIN[1] + j * STEP[1];
        if (std::abs(first[0] - expected_a) > CacheCoordinateTolerance(expected_a) ||
            std::abs(first[1] - expected_b) > CacheCoordinateTolerance(expected_b)) {
          throw std::runtime_error("grid coordinate mismatch");
        }
        std::copy_n(first, 4, F[i][j].begin());
      }
    }
  } catch (const std::exception &e) {
    std::cout << rang::fg::red << "[CORRUPTED]" << rang::fg::reset << " " << filename
              << ": " << e.what() << std::endl;
    return false;
  }

  file.close();
  std::cout << rang::fg::green << "[DONE]" << rang::fg::reset << std::endl;
  return true;
}

// Prepare one repeatedly used physical coordinate on the second axis
IArray2D::PreparedSecondCoordinate IArray2D::PrepareSecondCoordinate(
    double input_b) const {
  if (!std::isfinite(input_b)) {
    throw std::invalid_argument(
        "IArray2D::PrepareSecondCoordinate: non-finite input");
  }
  if (islog[1] && !(input_b > 0.0)) {
    throw std::out_of_range(
        "IArray2D::PrepareSecondCoordinate: non-positive logarithmic input");
  }
  double b = islog[1] ? std::log(input_b) : input_b;
  const double tolerance_b = std::max(CacheCoordinateTolerance(MIN[1]),
                                      CacheCoordinateTolerance(MAX[1]));
  if (b < MIN[1] - tolerance_b || b > MAX[1] + tolerance_b) {
    const double min_b = islog[1] ? std::exp(MIN[1]) : MIN[1];
    const double max_b = islog[1] ? std::exp(MAX[1]) : MAX[1];
    std::ostringstream error;
    error << "IArray2D::PrepareSecondCoordinate: input out of grid domain: "
          << name[1] << "=" << input_b << " [" << min_b << ", " << max_b
          << "]";
    throw std::out_of_range(error.str());
  }
  b = std::clamp(b, MIN[1], MAX[1]);
  int j = std::floor((b - MIN[1]) / STEP[1]);
  if (j < 0) { j = 0; }
  if (j >= (int)N[1]) { j = N[1] - 1; }
  const double y1 = F[0][j][1];
  const double y2 = F[0][j + 1][1];
  return {this, j, (b - y1) / (y2 - y1)};
}

// Interpolate one field with prepared coordinates and cubic Hermite consistency
std::pair<double, double> IArray2D::InterpolatePrepared(double input_a, double a,
                                                        const PreparedSecondCoordinate& prepared_b,
                                                        double*                         second) const {
  if (prepared_b.owner != this || !std::isfinite(input_a) ||
      !std::isfinite(a)) {
    throw std::invalid_argument(
        "IArray2D::InterpolatePrepared: invalid prepared coordinates");
  }
  if (islog[0] && !(input_a > 0.0)) {
    throw std::out_of_range(
        "IArray2D::InterpolatePrepared: non-positive logarithmic input");
  }
  const double tolerance_a = std::max(CacheCoordinateTolerance(MIN[0]),
                                      CacheCoordinateTolerance(MAX[0]));
  if (a < MIN[0] - tolerance_a || a > MAX[0] + tolerance_a) {
    const double min_a = islog[0] ? std::exp(MIN[0]) : MIN[0];
    const double max_a = islog[0] ? std::exp(MAX[0]) : MAX[0];
    std::ostringstream error;
    error << "IArray2D::InterpolatePrepared: input out of grid domain: "
          << name[0] << "=" << input_a << " [" << min_a << ", " << max_a
          << "]";
    throw std::out_of_range(error.str());
  }
  return InterpolatePreparedInDomain(input_a, std::clamp(a, MIN[0], MAX[0]), prepared_b, second);
}

// Interpolate prepared coordinates with an already validated first-axis domain
std::pair<double, double> IArray2D::InterpolatePreparedInDomain(double input_a, double a,
                                                                const PreparedSecondCoordinate& prepared_b,
                                                                double*                         second) const {
  int i = std::floor((a - MIN[0]) / STEP[0]);
  if (i < 0) { i = 0; }
  if (i >= (int)N[0]) { i = N[0] - 1; }
  const int j = prepared_b.index;
  const double x1 = F[i][j][0];
  const double x2 = F[i + 1][j][0];
  const double tx = (a - x1) / (x2 - x1);
  const double ty = prepared_b.fraction;

  const double h00 = 2.0 * tx * tx * tx - 3.0 * tx * tx + 1.0;
  const double h10 = tx * tx * tx - 2.0 * tx * tx + tx;
  const double h01 = -2.0 * tx * tx * tx + 3.0 * tx * tx;
  const double h11 = tx * tx * tx - tx * tx;
  const double dh00 = 6.0 * tx * tx - 6.0 * tx;
  const double dh10 = 3.0 * tx * tx - 4.0 * tx + 1.0;
  const double dh01 = -6.0 * tx * tx + 6.0 * tx;
  const double dh11 = 3.0 * tx * tx - 2.0 * tx;
  const double xstep = x2 - x1;

  const double one_minus_ty = 1.0 - ty;
  const double value_1 =
      one_minus_ty * F[i][j][2] + ty * F[i][j + 1][2];
  const double value_2 =
      one_minus_ty * F[i + 1][j][2] + ty * F[i + 1][j + 1][2];
  const double derivative_scale_1 =
      islog[0] ? first_axis_physical[i] : 1.0;
  const double derivative_scale_2 =
      islog[0] ? first_axis_physical[i + 1] : 1.0;
  const double slope_1 =
      (one_minus_ty * F[i][j][3] + ty * F[i][j + 1][3]) *
      derivative_scale_1;
  const double slope_2 =
      (one_minus_ty * F[i + 1][j][3] + ty * F[i + 1][j + 1][3]) *
      derivative_scale_2;
  const double value =
      h00 * value_1 + h10 * xstep * slope_1 +
      h01 * value_2 + h11 * xstep * slope_2;
  const double derivative_coordinate =
      (dh00 * value_1 + dh10 * xstep * slope_1 +
       dh01 * value_2 + dh11 * xstep * slope_2) /
      xstep;
  const double derivative = derivative_coordinate / (islog[0] ? input_a : 1.0);
  if (second != nullptr) {
    const double curvature = ((12.0 * tx - 6.0) * value_1 + (6.0 * tx - 4.0) * xstep * slope_1 +
                              (6.0 - 12.0 * tx) * value_2 + (6.0 * tx - 2.0) * xstep * slope_2) /
                             (xstep * xstep);
    *second                = islog[0] ? (curvature - derivative_coordinate) / (input_a * input_a) : curvature;
  }
  return {value, derivative};
}

// Interpolate one field with cubic Hermite consistency along the derivative axis
std::pair<double, double> IArray2D::Interpolate2D(double a, double b, double* second) const {
  if (!std::isfinite(a) || !std::isfinite(b)) {
    throw std::invalid_argument("IArray2D::Interpolate2D: non-finite input");
  }
  if ((islog[0] && !(a > 0.0)) || (islog[1] && !(b > 0.0))) {
    throw std::out_of_range(
        "IArray2D::Interpolate2D: non-positive input for logarithmic grid");
  }
  const PreparedSecondCoordinate prepared_b = PrepareSecondCoordinate(b);
  const double coordinate_a = islog[0] ? std::log(a) : a;
  return InterpolatePrepared(a, coordinate_a, prepared_b, second);
}

}  // namespace gra
