// Container class for different type of histograms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MUSERHISTOGRAMS_H
#define MUSERHISTOGRAMS_H

// C++
#include <complex>
#include <map>
#include <string>
#include <vector>

// Own
#include "Graniitti/Tech/MH1.h"
#include "Graniitti/Tech/MH2.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra {
class MUserHistograms {
 public:
  // Constructor, destructor
  MUserHistograms() = default;
  ~MUserHistograms() = default;

  void InitHistograms();
  void FillHistograms(double totalweight, const gra::LORENTZSCALAR &scalar);
  void PrintHistograms();
  void SaveHistograms(const std::string filename);

  void SetHistograms(unsigned int in);
  void FillCosThetaPhi(double totalweight, const gra::LORENTZSCALAR &scalar);

  // Histograms indexed by name std::string
  std::map<std::string, MH1<double>> h1;
  std::map<std::string, MH2>         h2;

  unsigned int HIST = 0;  // Histogramming level
};

}  // namespace gra

#endif
