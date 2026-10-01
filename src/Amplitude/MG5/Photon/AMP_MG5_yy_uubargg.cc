// Generated MadGraph photon amplitude for yy_uubargg
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// @@@@ MadGraph to GRANIITTI conversion done @@@@
//==========================================================================

#include <cmath>
#include <span>

#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_uubargg.h"
#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/HelAmps_sm_lepton_masses.h"

using namespace MG5_sm_lepton_masses;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: a a > u u~ g g QCD<=2 QED<=2 @1

//--------------------------------------------------------------------------
// Initialize process.

void AMP_MG5_yy_uubargg::initProc(std::string param_card_name) {
  InitParameters(SLHAReader(param_card_name));
}

// Initialize model parameters before constructing external wavefunctions
void AMP_MG5_yy_uubargg::InitParameters(SLHAReader slha)
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm_lepton_masses();  // GRANIITTI
  pars.setIndependentParameters(slha);
  gra::mg5::ValidateModel(slha, pars.Particles());
  pars.setIndependentCouplings();
  // pars.printIndependentParameters();  // GRANIITTI
  // pars.printIndependentCouplings();  // GRANIITTI
  // Reinitialization must replace the external mass table
  mME.clear();
  // Set external particle masses for this matrix element
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  // Reinitialization must clear the fixed color-flow buffer
  jamp2[0].fill(0.0);
}


// Prepare one physical on-shell HELAS phase-space point
bool AMP_MG5_yy_uubargg::setup_kinematics(
    gra::LORENTZSCALAR &lts, gra::mg5helas::EPAHardFrame *frame) {
  // *** MADGRAPH CONVENTION IS [E,px,py,pz] ! ***

  // Reuse event buffers because this function runs at every screening node
  final_buffer.clear();

  final_buffer.reserve(lts.decaytree.size());
  for (std::size_t i = 0; i < lts.decaytree.size(); ++i) {
    final_buffer.push_back(lts.decaytree[i].p4);
  }
  if (!gra::mg5::OnShellFinal(final_buffer, mME)) { return false; }

  gra::M4Vec p1_;
  gra::M4Vec p2_;
  if (frame != nullptr) {
    if (!gra::mg5helas::PrepareEPAHardFrame(lts, final_buffer, *frame)) { return false; }
    p1_ = frame->incoming[0];
    p2_ = frame->incoming[1];
  } else if (!gra::mg5helas::PrepareOnShellKinematics(
                 lts, final_buffer, p1_, p2_)) {
    return false;
  }

  p.clear();
  momentum_buffer.clear();
  momentum_buffer.reserve(ninitial + final_buffer.size());

  momentum_buffer.push_back({p1_.E(), p1_.Px(), p1_.Py(), p1_.Pz()});
  momentum_buffer.push_back({p2_.E(), p2_.Px(), p2_.Py(), p2_.Pz()});

  for (const auto &p4 : final_buffer) {
    momentum_buffer.push_back({p4.E(), p4.Px(), p4.Py(), p4.Pz()});
  }

  for (auto &mom : momentum_buffer) { p.push_back(mom.data()); }
  return true;
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

gra::mg5helas::MatrixElementEvaluation AMP_MG5_yy_uubargg::Evaluate(gra::LORENTZSCALAR &lts, double alphas, bool coherent_epa)
{
  lts.hamp.clear();
  if (!std::isfinite(alphas) || alphas < 0.0) { return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0}; }
  // Set the parameters which change event by event
  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }
  // Reset color flows
  for(int i = 0; i < 2; i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 64;
  if (!setup_kinematics(lts)) {
    return {gra::mg5helas::EvaluationStatus::KinematicsFailure, 0.0};
  }
  bool goodhel[ncomb] = {};
  int ntry = 0, sum_hel = 0, ngood = 0;
  int igood[ncomb + 1] = {};
  int jhel = 0;
  std::complex<double> * * wfs;
  double t[nprocesses];
  // Helicities for the process
  const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {8};

  ntry = 1;  // GRANIITTI

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.;
  }
  // Define permutation
  int perm[nexternal];
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i;
  }

  goto SKIPLABEL;  // GRANIITTI: Skip this block
  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]);
        t[0] = matrix_1_aa_uuxgg();

        double tsum = 0;
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc];
          tsum += t[iproc];
        }
        // Store which helicities give non-zero result
        if (std::fpclassify(tsum) != FP_ZERO && !goodhel[ihel])
        {
          goodhel[ihel] = true;
          ngood++;
          igood[ngood] = ihel;
        }
      }
    }
    jhel = 0;
    sum_hel = std::min(sum_hel, ngood);
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++;
      if (jhel >= ngood)
        jhel = 0;
      double hwgt = double(ngood)/double(sum_hel);
      int ihel = igood[jhel];
      calculate_wavefunctions(perm, helicities[ihel]);
      t[0] = matrix_1_aa_uuxgg();

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt;
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i];

SKIPLABEL:

  // @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // Define permutation
  for (int i = 0; i < nexternal; ++i) { perm[i] = i; }

  // Loop over helicity combinations in the fixed GRANIITTI basis
  std::vector<gra::mg5helas::HelicityComponent> components;
  components.reserve(ncomb * ncolor);
  const auto color_denominators = ColorDenominators();
  const auto color_metric = ColorMetric();
  for (int ihel = 0; ihel < ncomb; ++ihel) {
    calculate_wavefunctions(perm, helicities[ihel]);
    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);
    // Contract the JAMP color basis with the full finite-Nc metric
    // [REFERENCE: Lifson and Mattelaer, Eur. Phys. J. C 82, 1144 (2022), arXiv:2210.07267]
    const std::vector<std::complex<double>> color_amplitudes =
        gra::mg5helas::ColorMetricAmplitudes(
            std::vector<std::complex<double>>(jamp, jamp + ncolor),
            color_denominators, color_metric);
    if (color_amplitudes.size() != ncolor) {
      lts.hamp.clear();
      return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
    }
    std::vector<int> outgoing;
    for (int leg = ninitial; leg < nexternal; ++leg) { outgoing.push_back(helicities[ihel][leg]); }
    for (int color = 0; color < ncolor; ++color) {
      components.push_back({{helicities[ihel][0], helicities[ihel][1]},
                            outgoing, static_cast<std::size_t>(color),
                            color_amplitudes[color]});
    }
  }

  lts.hamp.clear();
  if (coherent_epa) {
    const gra::M4Vec k1(p[0][1], p[0][2], p[0][3], p[0][0]);
    const gra::M4Vec k2(p[1][1], p[1][2], p[1][3], p[1][0]);
    lts.hamp = gra::mg5helas::ContractEPAPhotonSources(lts, components, k1, k2);
  } else {
    for (const auto &component : components) { lts.hamp.push_back(component.value); }
  }

  // Total amplitude squared over all helicity combinations individually
  const double final_symmetry_factor =
      gra::mg5helas::AppliedFinalStateSymmetryFactor(lts);
  const double normalization = std::sqrt(
      final_symmetry_factor / static_cast<double>(denominators[0]));
  gra::Scale(lts.hamp, normalization);
  const double amp2 = gra::SquaredNorm(lts.hamp);
  // Screening applies the incoming photon spin average after coherent summation
  gra::Scale(lts.hamp, 2.0);

  return {gra::mg5helas::EvaluationStatus::Success, amp2};
                // @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}


// Evaluate every generated color flow in the fixed helicity basis
void AMP_MG5_yy_uubargg::CalcColorFlowHelicity(gra::LORENTZSCALAR &lts, double alphas,
                                         ColorFlowHelicityMatrix &jamp_matrix) {
  jamp_matrix.clear();
  if (!std::isfinite(alphas) || alphas < 0.0) { return; }
  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }


  const int ncomb = 64;
  jamp_matrix.assign(ncomb, ColorFlowVector(ncolor, 0.0));
  if (!setup_kinematics(lts)) { return; }

  const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};

  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) { perm[i] = i; }

  for (int ihel = 0; ihel < ncomb; ++ihel) {
    calculate_wavefunctions(perm, helicities[ihel]);

    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);

    for (int icolor = 0; icolor < ncolor; ++icolor) {
      jamp_matrix[ihel][icolor] = jamp[icolor];
    }
  }
}

// Contract generated color flows directly into arbitrary hard-process projectors
gra::mg5helas::EvaluationStatus AMP_MG5_yy_uubargg::CalcColorProjectedHelicity(
    gra::LORENTZSCALAR &lts, double alphas,
    const std::complex<double> *color_projectors, int projector_count,
    std::complex<double> *projected, gra::M4Vec *hard_k1,
    gra::M4Vec *hard_k2) {
  const std::size_t output_size =
      projector_count > 0
          ? static_cast<std::size_t>(projector_count) * nhelicity
          : 0;
  if (projected != nullptr && output_size > 0) {
    std::fill(projected, projected + output_size,
              std::complex<double>(0.0));
  }
  if (hard_k1 != nullptr) { *hard_k1 = gra::M4Vec(); }
  if (hard_k2 != nullptr) { *hard_k2 = gra::M4Vec(); }
  if (color_projectors == nullptr || projected == nullptr ||
      projector_count <= 0) {
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }
  if (!std::isfinite(alphas) || alphas < 0.0) {
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }

  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }

  if (!setup_kinematics(lts)) {
    return gra::mg5helas::EvaluationStatus::KinematicsFailure;
  }

  const int ncomb = nhelicity;
  const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};

  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) { perm[i] = i; }

  for (int ihel = 0; ihel < nhelicity; ++ihel) {
    calculate_wavefunctions(perm, helicities[ihel]);

    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);
    const std::span<const std::complex<double>> jamp_view(jamp, ncolor);

    for (int projector = 0; projector < projector_count; ++projector) {
      const std::span<const std::complex<double>> projector_view(
          color_projectors + projector * ncolor, ncolor);
      const std::complex<double> value =
          gra::BilinearProduct(projector_view, jamp_view);
      projected[static_cast<std::size_t>(projector) * nhelicity + ihel] = value;
    }
  }
  if (hard_k1 != nullptr) {
    *hard_k1 = gra::M4Vec(p[0][1], p[0][2], p[0][3], p[0][0]);
  }
  if (hard_k2 != nullptr) {
    *hard_k2 = gra::M4Vec(p[1][1], p[1][2], p[1][3], p[1][0]);
  }
  return gra::mg5helas::EvaluationStatus::Success;
}

// Prepare the transfer-independent EPA hard tensor and its reusable frame
gra::mg5helas::EvaluationStatus
AMP_MG5_yy_uubargg::CalcColorProjectedEPAHardHelicity(
    gra::LORENTZSCALAR &lts, double alphas,
    const std::complex<double> *color_projectors, int projector_count,
    std::complex<double> *projected, gra::mg5helas::EPAHardFrame *frame) {
  const std::size_t output_size =
      projector_count > 0
          ? static_cast<std::size_t>(projector_count) * nhelicity
          : 0;
  if (projected != nullptr && output_size > 0) {
    std::fill(projected, projected + output_size,
              std::complex<double>(0.0));
  }
  if (frame != nullptr) { *frame = {}; }
  if (frame == nullptr ||
      color_projectors == nullptr || projected == nullptr ||
      projector_count <= 0) {
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }
  if (!std::isfinite(alphas) || alphas < 0.0) {
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }

  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }
  if (!setup_kinematics(lts, frame)) {
    return gra::mg5helas::EvaluationStatus::KinematicsFailure;
  }
  return CalcColorProjectedPrepared(
      color_projectors, projector_count, projected);
}

// Contract prepared HELAS amplitudes into arbitrary color projectors
gra::mg5helas::EvaluationStatus AMP_MG5_yy_uubargg::CalcColorProjectedPrepared(
    const std::complex<double> *color_projectors, int projector_count,
    std::complex<double> *projected) {
  const int ncomb = nhelicity;
  const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};

  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) { perm[i] = i; }

  for (int ihel = 0; ihel < nhelicity; ++ihel) {
    calculate_wavefunctions(perm, helicities[ihel]);

    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);
    const std::span<const std::complex<double>> jamp_view(jamp, ncolor);

    for (int projector = 0; projector < projector_count; ++projector) {
      const std::span<const std::complex<double>> projector_view(
          color_projectors + projector * ncolor, ncolor);
      const std::complex<double> value = gra::BilinearProduct(projector_view, jamp_view);
      projected[static_cast<std::size_t>(projector) * nhelicity + ihel] = value;
    }
  }
  return gra::mg5helas::EvaluationStatus::Success;
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double AMP_MG5_yy_uubargg::sigmaHat()
{
  // Select between the different processes
  if(id1 == 22 && id2 == 22)
  {
    // Add matrix elements for processes with beams (22, 22)
    return matrix_element[0];
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.;
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void AMP_MG5_yy_uubargg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j;

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]);
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]);
  FFV1_1(w[2], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[6]);
  FFV1_2(w[3], w[1], pars.GC_2, pars.ZERO, pars.ZERO, w[7]);
  FFV1_1(w[6], w[4], pars.GC_11, pars.ZERO, pars.ZERO, w[8]);
  FFV1_1(w[6], w[5], pars.GC_11, pars.ZERO, pars.ZERO, w[9]);
  VVV1P0_1(w[4], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[10]);
  FFV1_2(w[3], w[4], pars.GC_11, pars.ZERO, pars.ZERO, w[11]);
  FFV1_1(w[6], w[1], pars.GC_2, pars.ZERO, pars.ZERO, w[12]);
  FFV1_2(w[3], w[5], pars.GC_11, pars.ZERO, pars.ZERO, w[13]);
  FFV1_2(w[3], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[14]);
  FFV1_1(w[2], w[1], pars.GC_2, pars.ZERO, pars.ZERO, w[15]);
  FFV1_2(w[14], w[4], pars.GC_11, pars.ZERO, pars.ZERO, w[16]);
  FFV1_2(w[14], w[5], pars.GC_11, pars.ZERO, pars.ZERO, w[17]);
  FFV1_1(w[2], w[4], pars.GC_11, pars.ZERO, pars.ZERO, w[18]);
  FFV1_2(w[14], w[1], pars.GC_2, pars.ZERO, pars.ZERO, w[19]);
  FFV1_1(w[2], w[5], pars.GC_11, pars.ZERO, pars.ZERO, w[20]);
  FFV1_1(w[15], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[21]);
  FFV1_2(w[11], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[22]);
  FFV1_2(w[13], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[23]);
  FFV1_2(w[7], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[24]);
  FFV1_1(w[18], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[25]);
  FFV1_1(w[20], w[0], pars.GC_2, pars.ZERO, pars.ZERO, w[26]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[7], w[8], w[5], pars.GC_11, amp[0]);
  FFV1_0(w[7], w[9], w[4], pars.GC_11, amp[1]);
  FFV1_0(w[7], w[6], w[10], pars.GC_11, amp[2]);
  FFV1_0(w[11], w[12], w[5], pars.GC_11, amp[3]);
  FFV1_0(w[11], w[9], w[1], pars.GC_2, amp[4]);
  FFV1_0(w[13], w[12], w[4], pars.GC_11, amp[5]);
  FFV1_0(w[13], w[8], w[1], pars.GC_2, amp[6]);
  FFV1_0(w[3], w[12], w[10], pars.GC_11, amp[7]);
  FFV1_0(w[16], w[15], w[5], pars.GC_11, amp[8]);
  FFV1_0(w[17], w[15], w[4], pars.GC_11, amp[9]);
  FFV1_0(w[14], w[15], w[10], pars.GC_11, amp[10]);
  FFV1_0(w[19], w[18], w[5], pars.GC_11, amp[11]);
  FFV1_0(w[17], w[18], w[1], pars.GC_2, amp[12]);
  FFV1_0(w[19], w[20], w[4], pars.GC_11, amp[13]);
  FFV1_0(w[16], w[20], w[1], pars.GC_2, amp[14]);
  FFV1_0(w[19], w[2], w[10], pars.GC_11, amp[15]);
  FFV1_0(w[11], w[21], w[5], pars.GC_11, amp[16]);
  FFV1_0(w[22], w[15], w[5], pars.GC_11, amp[17]);
  FFV1_0(w[13], w[21], w[4], pars.GC_11, amp[18]);
  FFV1_0(w[23], w[15], w[4], pars.GC_11, amp[19]);
  FFV1_0(w[3], w[21], w[10], pars.GC_11, amp[20]);
  FFV1_0(w[24], w[18], w[5], pars.GC_11, amp[21]);
  FFV1_0(w[7], w[25], w[5], pars.GC_11, amp[22]);
  FFV1_0(w[24], w[20], w[4], pars.GC_11, amp[23]);
  FFV1_0(w[7], w[26], w[4], pars.GC_11, amp[24]);
  FFV1_0(w[24], w[2], w[10], pars.GC_11, amp[25]);
  FFV1_0(w[13], w[25], w[1], pars.GC_2, amp[26]);
  FFV1_0(w[23], w[18], w[1], pars.GC_2, amp[27]);
  FFV1_0(w[11], w[26], w[1], pars.GC_2, amp[28]);
  FFV1_0(w[22], w[20], w[1], pars.GC_2, amp[29]);

}

double AMP_MG5_yy_uubargg::matrix_1_aa_uuxgg() {
  int i, j;
  // Local variables
  const int            ngraphs = 30;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  const auto           denom   = ColorDenominators();
  const auto           cf      = ColorMetric();

  calculate_color_flows(jamp);

  // Sum and square the color flows to get the matrix element
  double matrix = 0;
  for (i = 0; i < ncolor; i++) {
    ztemp = 0.;
    for (j = 0; j < ncolor; j++) ztemp = ztemp + cf[i][j] * jamp[j];
    matrix = matrix + std::real(ztemp * std::conj(jamp[i])) / denom[i];
  }

  // Store the leading color flows for choice of color
  for (i = 0; i < ncolor; i++) jamp2[0][i] += std::real(jamp[i] * std::conj(jamp[i]));

  return matrix;
}

void AMP_MG5_yy_uubargg::calculate_color_flows(std::complex<double> jamp[ncolor]) const {
    jamp[0] = -amp[0] + std::complex<double> (0, 1) * amp[2] - amp[5] - amp[6] +
        std::complex<double> (0, 1) * amp[7] - amp[9] + std::complex<double> (0,
        1) * amp[10] - amp[11] - amp[12] + std::complex<double> (0, 1) * amp[15]
        - amp[18] - amp[19] + std::complex<double> (0, 1) * amp[20] - amp[21] -
        amp[22] + std::complex<double> (0, 1) * amp[25] - amp[26] - amp[27];
    jamp[1] = -amp[1] - std::complex<double> (0, 1) * amp[2] - amp[3] - amp[4] -
        std::complex<double> (0, 1) * amp[7] - amp[8] - std::complex<double> (0,
        1) * amp[10] - amp[13] - amp[14] - std::complex<double> (0, 1) * amp[15]
        - amp[16] - amp[17] - std::complex<double> (0, 1) * amp[20] - amp[23] -
        amp[24] - std::complex<double> (0, 1) * amp[25] - amp[28] - amp[29];
}

std::vector<double> AMP_MG5_yy_uubargg::ColorDenominators() const { return {3, 3}; }

std::vector<std::vector<double>> AMP_MG5_yy_uubargg::ColorMetric() const {
  return {{16, -2}, {-2, 16}};
}
