// Nuclear identity and normalized density model
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEUS_H
#define MNUCLEUS_H

#include "Graniitti/Nuclear/MDensity.h"
#include "Graniitti/Nuclear/MTypes.h"

namespace gra::nuclear {

// Compute whether a PDG code has a valid 10LZZZAAAI nuclear encoding
bool IsNuclearPDG(int pdg);

// Decode and validate one 10LZZZAAAI nuclear PDG code
NucleusID DecodeNuclearPDG(int pdg);

// Encode one validated 10LZZZAAAI nuclear PDG identity
int EncodeNuclearPDG(unsigned int a, unsigned int z, unsigned int lambda = 0, unsigned int isomer = 0,
                     bool anti = false);

// Compute tuned spherical density parameters for one encoded nucleus
NucleusParam DefaultNucleusParam(int pdg, double mass, const GeometryParam &geometry);

// Own one immutable nuclear identity and its charge and matter densities
class MNucleus {
 public:
  // Construct one nucleus from explicit physical inputs
  explicit MNucleus(const NucleusParam &param);

  // Compute the decoded nuclear identity
  const NucleusID &ID() const { return id_; }

  // Compute the validated physical input parameters
  const NucleusParam &Param() const { return param_; }

  // Compute the mass number
  unsigned int A() const { return id_.a; }

  // Compute the proton number
  unsigned int Z() const { return id_.z; }

  // Compute the neutron number for an ordinary nucleus
  unsigned int N() const;

  // Compute the signed electric charge in positron units
  int Charge() const;

  // Compute the bare nuclear mass in GeV
  double Mass() const { return param_.mass; }

  // Compute the normalized nuclear charge density
  const MDensity &ChargeDensity() const { return charge_; }

  // Compute the normalized nuclear matter density
  const MDensity &MatterDensity() const { return matter_; }

 private:
  NucleusParam param_;
  NucleusID    id_;
  MDensity     charge_;
  MDensity     matter_;
};

}  // namespace gra::nuclear

#endif
