#ifndef GRANIITTI_AMPLITUDE_MG5_PP_Z_PROCESSBASE_H
#define GRANIITTI_AMPLITUDE_MG5_PP_Z_PROCESSBASE_H

#include <string>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Model.h"

namespace MG5_PP_Z {

// Common base class for one generated MG5 subprocess
class ProcessBase {
 public:
  // Construct one generated subprocess matrix element
  ProcessBase() = default;
  ProcessBase(const ProcessBase &) = delete;
  ProcessBase &operator=(const ProcessBase &) = delete;
  ProcessBase(ProcessBase &&) = delete;
  ProcessBase &operator=(ProcessBase &&) = delete;

  // Destroy one generated subprocess instance
  virtual ~ProcessBase() {}

  // Initialize this subprocess from the shared parameter card
  virtual void initProc(std::string param_card_name) = 0;

  // Evaluate the kinematic matrix element for the current momenta
  virtual void sigmaKin() = 0;

  // Compute the flavour-filtered matrix element
  virtual double sigmaHat() = 0;

  // Compute spin-color averaged complex helicity amplitudes
  virtual std::vector<gra::mg5helas::HelicityComponent> helicityAmplitudes() = 0;

  // Initialize masses and couplings together from one model card
  virtual void InitParameters(SLHAReader slha) = 0;

  // Compute evaluated model particle parameters
  virtual gra::mg5::ParticleMap Particles() const = 0;

  // Compute the evaluated UFO electromagnetic coupling
  virtual double AlphaQED() const = 0;

  // Compute the maximum generated amplitude orders
  virtual int AlphaSPower() const = 0;
  virtual int AlphaQEDPower() const = 0;

  // Compute the current external mass table
  virtual const std::vector<double> &getMasses() const = 0;

  // Set the external momenta in MadGraph convention
  virtual void setMomenta(std::vector<double *> &momenta) = 0;

  // Set the current incoming PDG flavours
  virtual void setInitial(int inid1, int inid2) = 0;

  // Set alpha_s for event-dependent couplings
  virtual void setAlphaS(double in) = 0;

  // Select the fixed Thomson-limit QED coupling
  void setAlphaQEDZero(bool value) noexcept { alpha_qed_zero_ = value; }

 protected:
  // Compute whether the Thomson-limit QED coupling is selected
  bool alphaQEDZero() const noexcept { return alpha_qed_zero_; }

 private:
  bool alpha_qed_zero_ = true;
};

}  // namespace MG5_PP_Z

#endif
