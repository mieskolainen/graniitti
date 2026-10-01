// Shared LHAPDF access
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MLHAPDF_H
#define MLHAPDF_H

// C++
#include <memory>
#include <string>
#include <utility>

// Own
#include "Graniitti/Tech/MFixedStore.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

namespace gra {

// Process-wide cache for read-only LHAPDF members
class MLHAPDFStore {
 public:
  using PDFPtr = std::shared_ptr<const LHAPDF::PDF>;

  // Compute one shared read-only LHAPDF member
  PDFPtr GetPDF(const std::string &setname, int member = 0);

 private:
  using Key = std::pair<std::string, int>;

  // Load an LHAPDF member and download only when set data are missing
  PDFPtr LoadPDF(const std::string &setname, int member) const;

  MFixedStore<Key, LHAPDF::PDF> store;
};

}  // namespace gra

#endif
