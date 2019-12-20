// (Sub)-Processes and Amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSUBPROC_H
#define MSUBPROC_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MDurham.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MGamma.h"
#include "Graniitti/MPDG.h"
#include "Graniitti/MRegge.h"
#include "Graniitti/MTensorPomeron.h"


namespace gra {

using math::abs2;
using math::msqrt;


// Abstract process base class
//
// BASIC DESIGN:
// Use "late construction", that is, objects are pointers and initialized only
// just before the event generation.
//
//
class MProc {

public:
  MProc(const std::string& i, const std::string& j, const std::string& k) {
    ISTATE      = i;
    CHANNEL     = j;
    DESCRIPTION = k;
  };
  std::string ISTATE;
  std::string CHANNEL;
  std::string DESCRIPTION;

  virtual ~MProc() {} // Needs to be virtual

  virtual double Amp2(gra::LORENTZSCALAR& lts) = 0;


  // Apply non-collinear EPA fluxes at cross section level
  double ApplyktEPAfluxes(double amp2, gra::LORENTZSCALAR& lts) {

    const double gammaflux1 = lts.excite1
                                  ? gra::form::IncohFlux(lts.x1, lts.t1, lts.qt1, lts.pfinal[1].M2())
                                  : form::CohFlux(lts.x1, lts.t1, lts.qt1);
    const double gammaflux2 = lts.excite2
                                  ? gra::form::IncohFlux(lts.x2, lts.t2, lts.qt2, lts.pfinal[2].M2())
                                  : form::CohFlux(lts.x2, lts.t2, lts.qt2);

    // const double phasespace = lts.s / lts.s_hat;    // This gives problems at low masses
    const double phasespace = 1.0 / (lts.x1 * lts.x2);  // Consistent with kt-factorization

    // Combine all
    const double fluxes = gammaflux1 * gammaflux2 * phasespace;
    const double tot = fluxes * amp2;

    // --------------------------------------------------------------------
    // Apply fluxes to helicity amplitudes
    const double sqrt_fluxes = msqrt(fluxes);  // to "amplitude level"
    for (const auto &i : aux::indices(lts.hamp)) { lts.hamp[i] *= sqrt_fluxes; }
    // --------------------------------------------------------------------

    return tot;
  }

  // Gamma-Gamma collinear Drees-Zeppenfeld (coherent flux) at cross section level
  double ApplyDZfluxes(double amp2, gra::LORENTZSCALAR& lts) {

    // Evaluate gamma pdfs
    const double f1         = form::DZFlux(lts.x1);
    const double f2         = form::DZFlux(lts.x2);
    const double phasespace = 1.0 / (lts.x1 * lts.x2);
    const double tot        = f1 * f2 * amp2 * phasespace;

    return tot;
  }

  // Apply Gamma-Gamma LUX-pdf (use at \mu > 10 GeV) at cross section level
  double ApplyLUXfluxes(double amp2, gra::LORENTZSCALAR& lts) {

    // @@ MULTITHREADING LOCK NEEDED FOR INITIALIZATION @@
    gra::g_mutex.lock();

    // Not created yet
    if (lts.GlobalPdfPtr == nullptr) {
      std::string pdfname = lts.LHAPDFSET;

    retry:
      try {
        lts.GlobalPdfPtr = LHAPDF::mkPDF(pdfname, 0);
        lts.pdf_trials   = 0;  // fine
      } catch (...) {
        ++lts.pdf_trials;
        std::string str = "ApplyLUXfluxes: Problem with reading LHAPDF '" + pdfname + "'";
        aux::AutoDownloadLHAPDF(pdfname);  // Try autodownload
        gra::g_mutex.unlock();             // Remember before throw, otherwise deadlock
        if (lts.pdf_trials >= 2) {         // too many failures
          throw std::invalid_argument(str);
        } else {
          goto retry;
        }
      }
    }
    gra::g_mutex.unlock();
    // @@ MULTITHREADING UNLOCK @@

    // pdf factorization scale
    const double Q2 = lts.s_hat / 4.0;

    // Evaluate gamma pdfs
    double f1 = 0.0;
    double f2 = 0.0;
    try {
      // Divide x out
      f1 = lts.GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x1, Q2) / lts.x1;
      f2 = lts.GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x2, Q2) / lts.x2;
    } catch (...) { throw std::invalid_argument("ApplyLUXfluxes: Failed evaluating LHAPDF"); }

    const double phasespace = 1.0 / (lts.x1 * lts.x2);
    const double tot        = f1 * f2 * amp2 * phasespace;

    return tot;
  }

  // Processes usable with different fluxes
  // 2y ->
  double GammaGammaCON(gra::LORENTZSCALAR& lts) {

    double amp2 = 0.0;
    if (!AssertN(2, lts.decaytree.size())) {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2-body final state");
    }

    if        (AssertN({24,-24}, PDGlist(lts))) {
      amp2 = Gamma->AmpMG5_yy_ww.CalcAmp2(lts, 0.0);
    } else if (AssertLeptonOrQuarkPair(PDGlist(lts))) {
      amp2 = Gamma->yyffbar(lts);
    } else {
      ThrowUnknownFinalState();
    }

    return amp2;
  }

  // Process class containers
  std::unique_ptr<MGamma>  Gamma  = nullptr;
  std::unique_ptr<MDurham> Durham = nullptr;
  std::unique_ptr<MRegge>  Regge  = nullptr;
  std::unique_ptr<MTensorPomeron> Tensor = nullptr;

  // Delayed constructors
  void CallGamma(gra::LORENTZSCALAR& lts)  { if (Gamma  == nullptr) { Gamma  = make_unique<MGamma>(lts);      }}
  void CallDurham(gra::LORENTZSCALAR& lts) { if (Durham == nullptr) { Durham = make_unique<MDurham>();        }}
  void CallTensor(gra::LORENTZSCALAR& lts) { if (Tensor == nullptr) { Tensor = make_unique<MTensorPomeron>(); }}
  void CallRegge(gra::LORENTZSCALAR& lts)  { if (Regge  == nullptr) { Regge  = make_unique<MRegge>(lts);      }}


  // Assert the final state list
  bool AssertN(const std::vector<int>& reference, const std::vector<int>& input) const {
    return std::is_permutation(reference.begin(), reference.end(), input.begin());
  }

  // Assert the number of final states
  bool AssertN(int reference, int input) const {
    return (reference == input) ? true : false;
  }

  // Assert lepton or quark pair
  bool AssertLeptonOrQuarkPair(const std::vector<int>& input) const {

    // PDG identifiers
    static const MMatrix<int> X = {{11, -11}, {13, -13}, {15, -15}, {1, -1}, {2, -2}, {3, -3}, {4, -4}, {5, -5}, {6, -6}};

    for (std::size_t i = 0; i < X.size_row(); ++i) {
      if (AssertN({X[i][0], X[i][1]}, input)) { return true; }
    }
    return false;
  }

  // Get PDG list
  std::vector<int> PDGlist(const LORENTZSCALAR& lts) const {
    std::vector<int> L(lts.decaytree.size(), 0);
    for (const auto& i : aux::indices(lts.decaytree)) {
      L[i] = lts.decaytree[i].p.pdg;
    }
    return L;
  }

  void ThrowUnknownFinalState() const {
    throw std::invalid_argument(ISTATE + "[" + CHANNEL + "]" + " with unsupported final state");
  }

};


// -----------------------------------------------------------------------
// Gamma-Gamma processes
// -----------------------------------------------------------------------

class PROC_0 : public MProc {
public:
  PROC_0() : MProc("yy", "RES", "2xGamma to parametric resonance             [kt-EPA]") {}
  ~PROC_0() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallGamma(lts);
    const double amp2 = Gamma->yyX(lts, lts.RESONANCES.begin()->second);
    return ApplyktEPAfluxes(amp2, lts);
  }
};
class PROC_1 : public MProc {
public:
  PROC_1() : MProc("yy", "Higgs", "2xGamma to SM Higgs                         [kt-EPA]") {}
  ~PROC_1() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallGamma(lts);
    const double amp2 = Gamma->yyHiggs(lts);
    return ApplyktEPAfluxes(amp2, lts);
  }
};
class PROC_2 : public MProc {
public:
  PROC_2() : MProc("yy", "monopolium(0)", "2xGamma to Monopolium (J=0)                 [kt-EPA]") {}
  ~PROC_2() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallGamma(lts);
    const double amp2 = Gamma->yyMP(lts);
    return ApplyktEPAfluxes(amp2, lts);
  }
};
class PROC_3 : public MProc {
public:
  PROC_3() : MProc("yy", "CON", "2xGamma to l+l-, qqbar, W+W-, monopolepair  [kt-EPA]") {}
  ~PROC_3() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallGamma(lts);
    double amp2 = GammaGammaCON(lts);
    return ApplyktEPAfluxes(amp2, lts);
  }
};
class PROC_4 : public MProc {
public:
  PROC_4() : MProc("yy", "QED", "2xGamma to l+l-, qqbar                      [FULL QED]") {}
  ~PROC_4() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallTensor(lts);

    if (!AssertN(2, lts.decaytree.size())) {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2-body final state");
    }
    return Tensor->ME4(lts);
  }
};


// -----------------------------------------------------------------------
// Inclusive processes
// -----------------------------------------------------------------------

class PROC_5 : public MProc {
public:
  PROC_5() : MProc("X", "EL", "Elastic                                     [Eikonal Pomeron]     (Use with screening loop on)") {}
  ~PROC_5() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    return abs2( Regge->ME2(lts, 1) );
  }
};
class PROC_6 : public MProc {
public:
  PROC_6() : MProc("X", "SD", "Single Diffractive                          [Triple Pomeron]      (With TOY fragmentation)") {}
  ~PROC_6() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    return abs2( Regge->ME2(lts, 2) );
  }
};
class PROC_7 : public MProc {
public:
  PROC_7() : MProc("X", "DD", "Double Diffractive                          [Triple Pomeron]      (With TOY fragmentation)") {}
  ~PROC_7() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    return abs2( Regge->ME2(lts, 3) );
  }
};
class PROC_8 : public MProc {
public:
  PROC_8() : MProc("X", "ND", "Non-Diffractive                             [N-cut soft Pomerons] (With TOY fragmentation)") {}
  ~PROC_8() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    return 1.0;
  }
};

// -----------------------------------------------------------------------
// Pomeron-Pomeron processes
// -----------------------------------------------------------------------

class PROC_9 : public MProc {
public:
  PROC_9() : MProc("PP", "RES", "Regge parametric resonance                  [Pomeron]") {}
  ~PROC_9() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    std::complex<double> A = 0.0;

    // Coherent sum of Resonances (loop over)
    for (auto &x : lts.RESONANCES) {
      const int J = static_cast<int>(x.second.p.spinX2 / 2.0);

      // Gamma-Pomeron for vectors (could be Pomeron-Odderon too)
      if (J == 1 && x.second.p.P == -1) {
        A += Regge->PhotoME3(lts, x.second);
      }
      // Pomeron-Pomeron, J = 0,1,2,... all ok
      else {
        A += Regge->ME3(lts, x.second);
      }
    }

    // ------------------------------------------------------------------
    // Set for the screening loop
    lts.hamp = {A};
    // ------------------------------------------------------------------

    return abs2(A);
  }
};

class PROC_10 : public MProc {
public:
  PROC_10() : MProc("PP", "RESHEL", "Regge sliding helicity amplitudes           [Pomeron]             <-- DEVELOPER ONLY PROCESS!") {}
  ~PROC_10() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    std::complex<double> A = 0.0;
    A = Regge->ME3HEL(lts, lts.RESONANCES.begin()->second);
    return abs2(A);
  }
};

class PROC_11 : public MProc {
public:
  PROC_11() : MProc("PP", "RESTENSOR", "Regge resonance                             [Tensor Pomeron]") {}
  ~PROC_11() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallTensor(lts);
    return Tensor->ME3(lts);
  }
};

class PROC_12 : public MProc {
public:
  PROC_12() : MProc("PP", "CONTENSOR", "Regge continuum 2-body                      [Tensor Pomeron]") {}
  ~PROC_12() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallTensor(lts);

    if (!AssertN(2, lts.decaytree.size())) {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2-body final state");
    }
    return Tensor->ME4(lts);
  }
};

class PROC_13 : public MProc {
public:
  PROC_13() : MProc("PP", "CONTENSOR24", "Regge continuum 2-body > 4-body             [Tensor Pomeron]      <-- DEVELOPER ONLY PROCESS!") {}
  ~PROC_13() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallTensor(lts);

    if (!AssertN(2, lts.decaytree.size()) && !AssertN(4, lts.decaytree.size())) {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2 > {2x} or 4-body final state");
    }
    return Tensor->ME6(lts);
  }
};

class PROC_14 : public MProc {
public:
  PROC_14() : MProc("PP", "RES+CONTENSOR", "Regge resonances + continuum 2-body         [Tensor Pomeron / yP]") {}
  ~PROC_14() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallTensor(lts);

    if (!AssertN(2, lts.decaytree.size())) {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2-body final state");
    }

    // 1. Evaluate continuum matrix element -> helicity amplitudes to lts.hamp
    Tensor->ME4(lts);

    // Add to temp vector (init with zero!)
    std::vector<std::complex<double>> tempsum(lts.hamp.size(), 0.0);
    for (const auto &i : aux::indices(lts.hamp)) { tempsum[i] += lts.hamp[i]; }

    // 2. Evaluate resonance matrix elements -> helicity amplitudes to lts.hamp
    if (lts.RESONANCES.size() != 0) {
      // We loop over resonances inside ME3
      Tensor->ME3(lts);
      // Add to temp vector
      for (const auto &i : aux::indices(lts.hamp)) { tempsum[i] += lts.hamp[i]; }
    }
    // ------------------------------------------------------------------
    // Set helicity amplitudes for the screening loop
    lts.hamp = tempsum;
    // ------------------------------------------------------------------

    // Get total amplitude squared 1/4 \sum_h |A_h|^2
    double amp2 = 0.0;
    for (const auto &i : aux::indices(lts.hamp)) { amp2 += gra::math::abs2(lts.hamp[i]); }
    amp2 /= 4;  // Initial state helicity average

    return amp2;
  }
};


class PROC_15 : public MProc {
public:
  PROC_15() : MProc("PP", "CON", "Regge continuum 2/4/6-body                  [Pomeron]") {}
  ~PROC_15() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);

    std::complex<double> A = 0.0;
    if        (AssertN(2, lts.decaytree.size())) {
      A = Regge->ME4(lts, 1);
    } else if (AssertN(4, lts.decaytree.size())) {
      A = Regge->ME6(lts);
    } else if (AssertN(6, lts.decaytree.size())) {
      A = Regge->ME8(lts);
    } else {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2, 4 or 6-body final state");
    }
    return abs2(A);
  }
};


class PROC_16 : public MProc {
public:
  PROC_16() : MProc("PP", "CON-", "Regge continuum 2-body with [t-u] amplitude [Pomeron]") {}
  ~PROC_16() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);

    if (!AssertN(2, lts.decaytree.size())) {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2-body final state");
    }
    std::complex<double> A = Regge->ME4(lts, -1);  
    return abs2(A);
  }
};


class PROC_17 : public MProc {
public:
  PROC_17() : MProc("PP", "RES+CON", "Regge resonances + continuum 2-body         [Pomeron / yP]") {}
  ~PROC_17() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    std::complex<double> A = 0.0;

    if (!AssertN(2, lts.decaytree.size())) {
      throw std::invalid_argument(ISTATE + "[" + CHANNEL + "] requires 2-body final state");
    }

    // 1. Continuum matrix element
    A = Regge->ME4(lts, 1);

    // 2. Coherent sum of Resonances (loop over)
    for (auto &x : lts.RESONANCES) {
      const int J = static_cast<int>(x.second.p.spinX2 / 2.0);

      // Gamma-Pomeron for vectors
      if (J == 1 && x.second.p.P == -1) {
        A += Regge->PhotoME3(lts, x.second);
      }
      // Pomeron-Pomeron, J = 0,1,2,... all ok
      else {
        A += Regge->ME3(lts, x.second);
      }
    }

    // ------------------------------------------------------------------
    // Set helicity amplitudes for the screening loop
    lts.hamp = {A};
    // ------------------------------------------------------------------

    return abs2(A);
  }
};

// -----------------------------------------------------------------------
// Gamma(Odderon)-Pomeron processes
// -----------------------------------------------------------------------

class PROC_18 : public MProc {
public:
  PROC_18() : MProc("yP", "RES", "Photoproduced parametric resonance          [kt-EPA] x [Pomeron]") {}
  ~PROC_18() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    std::complex<double> A = 0.0;

    // Coherent sum of Resonances (loop over)
    for (auto &x : lts.RESONANCES) {
      const int J = static_cast<int>(x.second.p.spinX2 / 2.0);
      // Vectors only
      if (J == 1 && x.second.p.P == -1) {
        A += Regge->PhotoME3(lts, lts.RESONANCES.begin()->second);
      }
    }
    return abs2(A);  // amplitude squared
  }
};


class PROC_19 : public MProc {
public:
  PROC_19() : MProc("yP", "RESTENSOR", "Photoproduced resonance                     [QED] x [Tensor Pomeron]") {}
  ~PROC_19() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    return Tensor->ME3(lts);
  }
};


class PROC_20 : public MProc {
public:
  PROC_20() : MProc("OP", "RES", "Regge parametric vector resonance           [Odderon] x [Pomeron]") {}
  ~PROC_20() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallRegge(lts);
    std::complex<double> A = 0.0;
    // Coherent sum of Resonances (loop over)
    for (auto &x : lts.RESONANCES) {
      const int J = static_cast<int>(x.second.p.spinX2 / 2.0);
      // Vectors only
      if (J == 1 && x.second.p.P == -1) {
        A += Regge->ME3ODD(lts, x.second);
      }
    }
    return abs2(A);  // amplitude squared
  }
};


// -----------------------------------------------------------------------
// Durham QCD processes
// -----------------------------------------------------------------------

class PROC_21 : public MProc {
public:
  PROC_21() : MProc("gg", "chic(0)", "QCD resonance chic(0)                       [Durham QCD]") {}
  ~PROC_21() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallDurham(lts);
    double amp2 = 0.0;
    amp2 = Durham->DurhamQCD(lts, CHANNEL);
    return amp2;
  }
};

class PROC_22 : public MProc {
public:
  PROC_22() : MProc("gg", "CON", "QCD continuum to gg, 2 x pseudoscalar       [Durham QCD]          <-- UNDER VALIDATION!") {}
  ~PROC_22() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallDurham(lts);

    double amp2 = 0.0;
    if        (AssertN({21,21}, PDGlist(lts))) {
      amp2 = Durham->DurhamQCD(lts, "gg");
    } else if (AssertN(2, lts.decaytree.size())) {
      amp2 = Durham->DurhamQCD(lts, "MMbar");
    } else {
      ThrowUnknownFinalState();
    }
    return amp2;
  }
};


// -----------------------------------------------------------------------
// Gamma-Gamma processes
// -----------------------------------------------------------------------

class PROC_23 : public MProc {
public:
  PROC_23() : MProc("yy_LUX", "CON", "Collinear yy to l+l, qqbar, W+W- or monopolepair   [LUX-PDF]") {}
  ~PROC_23() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallGamma(lts);
    double amp2 = GammaGammaCON(lts);
    return ApplyLUXfluxes(amp2, lts);
  }
};

class PROC_24 : public MProc {
public:
  PROC_24() : MProc("yy_DZ", "CON", "Collinear yy to l+l, qqbar, W+W- or monopolepair   [Drees-Zeppenfeld EPA]") {}
  ~PROC_24() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallGamma(lts);
    double amp2 = GammaGammaCON(lts);
    return ApplyLUXfluxes(amp2, lts);
  }
};

class PROC_25 : public MProc {
public:
  PROC_25() : MProc("yy", "FLUX", "2xGamma with constant matrix element        [kt-EPA]") {}
  ~PROC_25() {}
  virtual double Amp2(gra::LORENTZSCALAR& lts) {
    CallGamma(lts);
    double amp2 = 1.0;
    return ApplyktEPAfluxes(amp2, lts);
  }
};


// Umbrella class
class MSubProc {

 public:
  MSubProc(const std::string &_ISTATE, const std::string &_CHANNEL);
  MSubProc(const std::vector<std::string> &first);
  void ConstructDescriptions(const std::string &first);

  MSubProc()  {}
  ~MSubProc();

  double GetBareAmplitude2(gra::LORENTZSCALAR &lts);

  std::string  ISTATE;      // "PP","yy","gg" etc.
  std::string  CHANNEL;     // "CON", "RES" etc.
  unsigned int LIPSDIM = 0; // Lorentz Invariant Phase Space Dimension
  
  // Available channels and their descriptions
  std::map<std::string, std::map<std::string, std::string>> descriptions;

 private:

  void CreateProcesses();
  void DeleteProcesses();
  void ActivateProcess();

  std::vector<MProc*> pr;

};

}  // namespace gra

#endif
