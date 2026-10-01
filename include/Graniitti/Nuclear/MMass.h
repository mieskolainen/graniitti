// Nuclear masses from evaluated tables and explicit theoretical completion
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARMASS_H
#define MNUCLEARMASS_H

#include <array>
#include <string>
#include <vector>

#include "Graniitti/Math/MMatrix.h"

namespace gra::nuclear {

enum class MassType { AME, FRDM };

struct MassTable {
  MassType    type = MassType::AME;
  std::string file;
};

struct MassParam {
  std::vector<MassTable> tables;
  std::array<double, 5>  bwm{};  // Volume, surface, Coulomb, asymmetry and pairing coefficients in MeV
  double                asymmetry = 0.0, pairing = 0.0;  // Mass-number scales
};

// Complete all daughter masses at initialization, preserving the ordered table priority
class MMass {
 public:
  // Load original mass tables and complete their missing isotopes with the configured BWM formula
  explicit MMass(const MassParam& param);

  // Compute a bare nuclear mass in GeV from the immutable completed table
  double Mass(unsigned int a, unsigned int z) const;

 private:
  MMatrix<double> mass_;
};

}  // namespace gra::nuclear

#endif
