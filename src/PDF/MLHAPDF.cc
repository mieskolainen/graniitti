// Shared LHAPDF access
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/PDF/MLHAPDF.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

namespace gra {
namespace {

// Initialize LHAPDF lazy caches before publishing a shared read-only member
void PrepareSharedPDF(LHAPDF::PDF &pdf) {
  const std::vector<int> &flavors = pdf.flavors();
  static_cast<void>(pdf.forcePositive());

  const LHAPDF::PDFInfo &info = pdf.info();
  const double x_min = info.get_entry_as<double>("XMin", 0.0);
  const double x_max = info.get_entry_as<double>("XMax", 1.0);
  const double q_min = info.get_entry_as<double>("QMin", 1.0);
  const double x_probe = 0.5 * (x_min + x_max);
  if (!flavors.empty() && std::isfinite(x_probe) && x_probe > 0.0 &&
      x_probe < 1.0 && std::isfinite(q_min) && q_min > 0.0) {
    static_cast<void>(pdf.xfxQ2(flavors.front(), x_probe, q_min * q_min));
  }

  if (!pdf.hasAlphaS()) {
    return;
  }

  double alpha_q = info.get_entry_as<double>("MZ", q_min);
  if (info.has_key("AlphaS_Qs")) {
    const std::vector<double> alpha_q_nodes =
        info.get_entry_as<std::vector<double>>("AlphaS_Qs");
    if (!alpha_q_nodes.empty()) {
      alpha_q = alpha_q_nodes[alpha_q_nodes.size() / 2];
    }
  }
  if (std::isfinite(alpha_q) && alpha_q > 0.0) {
    static_cast<void>(pdf.alphasQ2(alpha_q * alpha_q));
  }
}

} // namespace

// Compute one shared read-only LHAPDF member
MLHAPDFStore::PDFPtr MLHAPDFStore::GetPDF(const std::string &setname,
                                          int member) {
  aux::ValidateLHAPDFName(setname);
  if (member < 0) {
    throw std::invalid_argument("MLHAPDFStore::GetPDF: negative PDF member " +
                                std::to_string(member));
  }

  const Key key = {setname, member};
  return store.GetOrLoad(key, [this, &setname, member] {
    return LoadPDF(setname, member);
  });
}

// Load an LHAPDF member and download only when set data are missing
MLHAPDFStore::PDFPtr MLHAPDFStore::LoadPDF(const std::string &setname,
                                           int member) const {
  try {
    bool downloaded = false;
    if (LHAPDF::findpdfsetinfopath(setname).empty()) {
      aux::AutoDownloadLHAPDF(setname);
      downloaded = true;
    }

    // Check the member against set metadata before attempting a data download
    const LHAPDF::PDFSet set(setname);
    if (static_cast<std::size_t>(member) >= set.size()) {
      throw std::invalid_argument("member index outside set range");
    }
    if (!downloaded && LHAPDF::findpdfmempath(setname, member).empty()) {
      aux::AutoDownloadLHAPDF(setname);
    }
    std::shared_ptr<LHAPDF::PDF> pdf(
        LHAPDF::mkPDF(setname, static_cast<std::size_t>(member)));
    PrepareSharedPDF(*pdf);
    return pdf;
  } catch (const std::exception &e) {
    throw std::invalid_argument("MLHAPDFStore::LoadPDF: failed reading LHAPDF '" + setname +
                                "' member " + std::to_string(member) + ": " + e.what());
  }
}

} // namespace gra
