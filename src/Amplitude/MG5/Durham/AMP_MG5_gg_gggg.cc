//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// @@@@ MadGraph to GRANIITTI conversion done @@@@
//==========================================================================

#include <cmath>
#include <span>

#include "Graniitti/Amplitude/MG5/Durham/AMP_MG5_gg_gggg.h"
#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/HelAmps_sm_lepton_masses.h"

using namespace MG5_sm_lepton_masses;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > g g g g QED=0 @1

//--------------------------------------------------------------------------
// Initialize process.

void AMP_MG5_gg_gggg::initProc(std::string param_card_name) {
  InitParameters(SLHAReader(param_card_name));
}

// Initialize model parameters before constructing external wavefunctions
void AMP_MG5_gg_gggg::InitParameters(SLHAReader slha)
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
bool AMP_MG5_gg_gggg::setup_kinematics(gra::LORENTZSCALAR &lts) {
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
  if (!gra::mg5helas::PrepareOnShellKinematics(
          lts, final_buffer, p1_, p2_)) { return false; }

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

gra::mg5helas::MatrixElementEvaluation AMP_MG5_gg_gggg::Evaluate(gra::LORENTZSCALAR &lts, double alphas)
{
  lts.hamp.clear();
  if (!std::isfinite(alphas) || alphas < 0.0) { return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0}; }
  // Set the parameters which change event by event
  pars.setDependentParameters(alphas);
  pars.setDependentCouplings();
  // Reset color flows
  for(int i = 0; i < 120; i++ )
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
  const int denominators[nprocesses] = {6144};

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
        t[0] = matrix_1_gg_gggg();

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
      t[0] = matrix_1_gg_gggg();

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
  lts.hamp.clear();
  lts.hamp.reserve(ncomb * ncolor);
  const auto color_denominators = ColorDenominators();
  const auto color_metric = ColorMetric();
  for (int ihel = 0; ihel < ncomb; ++ihel) {
    calculate_wavefunctions(perm, helicities[ihel]);
    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);
    // Contract the JAMP color basis with the full finite-Nc metric
    // [REFERENCE: Lifson and Mattelaer, Eur. Phys. J. C 82, 1144 (2022), arXiv:2210.07267]
    const auto color_amplitudes = gra::mg5helas::ColorMetricAmplitudes(
        std::vector<std::complex<double>>(jamp, jamp + ncolor),
        color_denominators, color_metric);
    if (color_amplitudes.size() != ncolor) {
      lts.hamp.clear();
      return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
    }
    for (const auto amplitude : color_amplitudes) {
      lts.hamp.push_back(amplitude);
    }
  }

  // Total amplitude squared over all helicity combinations individually
  const double amp2 = gra::SquaredNorm(lts.hamp) / denominators[0];  // spin average matrix element squared

  return {gra::mg5helas::EvaluationStatus::Success, amp2};
                // @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}


// Evaluate every generated color flow in the fixed helicity basis
void AMP_MG5_gg_gggg::CalcColorFlowHelicity(gra::LORENTZSCALAR &lts, double alphas,
                                         ColorFlowHelicityMatrix &jamp_matrix) {
  jamp_matrix.clear();
  if (!std::isfinite(alphas) || alphas < 0.0) { return; }
  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();


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
gra::mg5helas::EvaluationStatus AMP_MG5_gg_gggg::CalcColorProjectedHelicity(
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

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double AMP_MG5_gg_gggg::sigmaHat()
{
  // Select between the different processes
  if(id1 == 21 && id2 == 21)
  {
    // Add matrix elements for processes with beams (21, 21)
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

void AMP_MG5_gg_gggg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j;

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  vxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  vxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]);
  VVV1P0_1(w[0], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[6]);
  VVV1P0_1(w[2], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[7]);
  VVV1P0_1(w[6], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[8]);
  VVV1P0_1(w[6], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[9]);
  VVV1P0_1(w[4], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[10]);
  VVV1P0_1(w[2], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[11]);
  VVV1P0_1(w[6], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[12]);
  VVV1P0_1(w[3], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[13]);
  VVV1P0_1(w[2], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[14]);
  VVV1P0_1(w[3], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[15]);
  VVV1P0_1(w[6], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[16]);
  VVVV1P0_1(w[2], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[17]);
  VVVV3P0_1(w[2], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[18]);
  VVVV4P0_1(w[2], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[19]);
  VVVV1P0_1(w[2], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[20]);
  VVVV3P0_1(w[2], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[21]);
  VVVV4P0_1(w[2], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[22]);
  VVVV1P0_1(w[2], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[23]);
  VVVV3P0_1(w[2], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[24]);
  VVVV4P0_1(w[2], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[25]);
  VVVV1P0_1(w[3], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[26]);
  VVVV3P0_1(w[3], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[27]);
  VVVV4P0_1(w[3], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[28]);
  VVV1P0_1(w[0], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[29]);
  VVV1P0_1(w[1], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[30]);
  VVV1P0_1(w[29], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[31]);
  VVV1P0_1(w[29], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[32]);
  VVV1P0_1(w[1], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[33]);
  VVV1P0_1(w[29], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[34]);
  VVV1P0_1(w[1], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[35]);
  VVV1P0_1(w[29], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[36]);
  VVVV1P0_1(w[1], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[37]);
  VVVV3P0_1(w[1], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[38]);
  VVVV4P0_1(w[1], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[39]);
  VVVV1P0_1(w[1], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[40]);
  VVVV3P0_1(w[1], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[41]);
  VVVV4P0_1(w[1], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[42]);
  VVVV1P0_1(w[1], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[43]);
  VVVV3P0_1(w[1], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[44]);
  VVVV4P0_1(w[1], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[45]);
  VVV1P0_1(w[0], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[46]);
  VVV1P0_1(w[1], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[47]);
  VVV1P0_1(w[46], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[48]);
  VVV1P0_1(w[46], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[49]);
  VVV1P0_1(w[46], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[50]);
  VVV1P0_1(w[46], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[51]);
  VVVV1P0_1(w[1], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[52]);
  VVVV3P0_1(w[1], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[53]);
  VVVV4P0_1(w[1], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[54]);
  VVVV1P0_1(w[1], w[2], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[55]);
  VVVV3P0_1(w[1], w[2], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[56]);
  VVVV4P0_1(w[1], w[2], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[57]);
  VVV1P0_1(w[0], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[58]);
  VVV1P0_1(w[58], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[59]);
  VVV1P0_1(w[58], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[60]);
  VVV1P0_1(w[58], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[61]);
  VVV1P0_1(w[58], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[62]);
  VVVV1P0_1(w[1], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[63]);
  VVVV3P0_1(w[1], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[64]);
  VVVV4P0_1(w[1], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[65]);
  VVV1P0_1(w[0], w[5], pars.GC_10, pars.ZERO, pars.ZERO, w[66]);
  VVV1P0_1(w[66], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[67]);
  VVV1P0_1(w[66], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[68]);
  VVV1P0_1(w[66], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[69]);
  VVV1P0_1(w[66], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[70]);
  VVV1P0_1(w[0], w[47], pars.GC_10, pars.ZERO, pars.ZERO, w[71]);
  VVV1P0_1(w[0], w[15], pars.GC_10, pars.ZERO, pars.ZERO, w[72]);
  VVV1P0_1(w[0], w[13], pars.GC_10, pars.ZERO, pars.ZERO, w[73]);
  VVV1P0_1(w[0], w[10], pars.GC_10, pars.ZERO, pars.ZERO, w[74]);
  VVV1P0_1(w[0], w[30], pars.GC_10, pars.ZERO, pars.ZERO, w[75]);
  VVV1P0_1(w[0], w[11], pars.GC_10, pars.ZERO, pars.ZERO, w[76]);
  VVV1P0_1(w[0], w[14], pars.GC_10, pars.ZERO, pars.ZERO, w[77]);
  VVV1P0_1(w[0], w[33], pars.GC_10, pars.ZERO, pars.ZERO, w[78]);
  VVV1P0_1(w[0], w[7], pars.GC_10, pars.ZERO, pars.ZERO, w[79]);
  VVV1P0_1(w[0], w[35], pars.GC_10, pars.ZERO, pars.ZERO, w[80]);
  VVVV1P0_1(w[0], w[1], w[2], pars.GC_12, pars.ZERO, pars.ZERO, w[81]);
  VVVV3P0_1(w[0], w[1], w[2], pars.GC_12, pars.ZERO, pars.ZERO, w[82]);
  VVVV4P0_1(w[0], w[1], w[2], pars.GC_12, pars.ZERO, pars.ZERO, w[83]);
  VVVV1P0_1(w[0], w[1], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[84]);
  VVVV3P0_1(w[0], w[1], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[85]);
  VVVV4P0_1(w[0], w[1], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[86]);
  VVVV1P0_1(w[0], w[1], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[87]);
  VVVV3P0_1(w[0], w[1], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[88]);
  VVVV4P0_1(w[0], w[1], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[89]);
  VVVV1P0_1(w[0], w[1], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[90]);
  VVVV3P0_1(w[0], w[1], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[91]);
  VVVV4P0_1(w[0], w[1], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[92]);
  VVVV1P0_1(w[0], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[93]);
  VVVV3P0_1(w[0], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[94]);
  VVVV4P0_1(w[0], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[95]);
  VVVV1P0_1(w[0], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[96]);
  VVVV3P0_1(w[0], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[97]);
  VVVV4P0_1(w[0], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[98]);
  VVVV1P0_1(w[0], w[2], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[99]);
  VVVV3P0_1(w[0], w[2], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[100]);
  VVVV4P0_1(w[0], w[2], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[101]);
  VVVV1P0_1(w[0], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[102]);
  VVVV3P0_1(w[0], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[103]);
  VVVV4P0_1(w[0], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[104]);
  VVVV1P0_1(w[0], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[105]);
  VVVV3P0_1(w[0], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[106]);
  VVVV4P0_1(w[0], w[3], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[107]);
  VVVV1P0_1(w[0], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[108]);
  VVVV3P0_1(w[0], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[109]);
  VVVV4P0_1(w[0], w[4], w[5], pars.GC_12, pars.ZERO, pars.ZERO, w[110]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVVV1_0(w[6], w[7], w[4], w[5], pars.GC_12, amp[0]);
  VVVV3_0(w[6], w[7], w[4], w[5], pars.GC_12, amp[1]);
  VVVV4_0(w[6], w[7], w[4], w[5], pars.GC_12, amp[2]);
  VVV1_0(w[7], w[5], w[8], pars.GC_10, amp[3]);
  VVV1_0(w[7], w[4], w[9], pars.GC_10, amp[4]);
  VVV1_0(w[6], w[7], w[10], pars.GC_10, amp[5]);
  VVVV1_0(w[6], w[11], w[3], w[5], pars.GC_12, amp[6]);
  VVVV3_0(w[6], w[11], w[3], w[5], pars.GC_12, amp[7]);
  VVVV4_0(w[6], w[11], w[3], w[5], pars.GC_12, amp[8]);
  VVV1_0(w[11], w[5], w[12], pars.GC_10, amp[9]);
  VVV1_0(w[11], w[3], w[9], pars.GC_10, amp[10]);
  VVV1_0(w[6], w[11], w[13], pars.GC_10, amp[11]);
  VVVV1_0(w[6], w[14], w[3], w[4], pars.GC_12, amp[12]);
  VVVV3_0(w[6], w[14], w[3], w[4], pars.GC_12, amp[13]);
  VVVV4_0(w[6], w[14], w[3], w[4], pars.GC_12, amp[14]);
  VVV1_0(w[14], w[4], w[12], pars.GC_10, amp[15]);
  VVV1_0(w[14], w[3], w[8], pars.GC_10, amp[16]);
  VVV1_0(w[6], w[14], w[15], pars.GC_10, amp[17]);
  VVVV1_0(w[6], w[2], w[15], w[5], pars.GC_12, amp[18]);
  VVVV3_0(w[6], w[2], w[15], w[5], pars.GC_12, amp[19]);
  VVVV4_0(w[6], w[2], w[15], w[5], pars.GC_12, amp[20]);
  VVV1_0(w[15], w[5], w[16], pars.GC_10, amp[21]);
  VVV1_0(w[2], w[15], w[9], pars.GC_10, amp[22]);
  VVVV1_0(w[6], w[2], w[13], w[4], pars.GC_12, amp[23]);
  VVVV3_0(w[6], w[2], w[13], w[4], pars.GC_12, amp[24]);
  VVVV4_0(w[6], w[2], w[13], w[4], pars.GC_12, amp[25]);
  VVV1_0(w[13], w[4], w[16], pars.GC_10, amp[26]);
  VVV1_0(w[2], w[13], w[8], pars.GC_10, amp[27]);
  VVVV1_0(w[6], w[2], w[3], w[10], pars.GC_12, amp[28]);
  VVVV3_0(w[6], w[2], w[3], w[10], pars.GC_12, amp[29]);
  VVVV4_0(w[6], w[2], w[3], w[10], pars.GC_12, amp[30]);
  VVV1_0(w[3], w[10], w[16], pars.GC_10, amp[31]);
  VVV1_0(w[2], w[10], w[12], pars.GC_10, amp[32]);
  VVV1_0(w[6], w[17], w[5], pars.GC_10, amp[33]);
  VVV1_0(w[6], w[18], w[5], pars.GC_10, amp[34]);
  VVV1_0(w[6], w[19], w[5], pars.GC_10, amp[35]);
  VVV1_0(w[6], w[20], w[4], pars.GC_10, amp[36]);
  VVV1_0(w[6], w[21], w[4], pars.GC_10, amp[37]);
  VVV1_0(w[6], w[22], w[4], pars.GC_10, amp[38]);
  VVV1_0(w[6], w[23], w[3], pars.GC_10, amp[39]);
  VVV1_0(w[6], w[24], w[3], pars.GC_10, amp[40]);
  VVV1_0(w[6], w[25], w[3], pars.GC_10, amp[41]);
  VVV1_0(w[6], w[2], w[26], pars.GC_10, amp[42]);
  VVV1_0(w[6], w[2], w[27], pars.GC_10, amp[43]);
  VVV1_0(w[6], w[2], w[28], pars.GC_10, amp[44]);
  VVVV1_0(w[29], w[30], w[4], w[5], pars.GC_12, amp[45]);
  VVVV3_0(w[29], w[30], w[4], w[5], pars.GC_12, amp[46]);
  VVVV4_0(w[29], w[30], w[4], w[5], pars.GC_12, amp[47]);
  VVV1_0(w[30], w[5], w[31], pars.GC_10, amp[48]);
  VVV1_0(w[30], w[4], w[32], pars.GC_10, amp[49]);
  VVV1_0(w[29], w[30], w[10], pars.GC_10, amp[50]);
  VVVV1_0(w[29], w[33], w[3], w[5], pars.GC_12, amp[51]);
  VVVV3_0(w[29], w[33], w[3], w[5], pars.GC_12, amp[52]);
  VVVV4_0(w[29], w[33], w[3], w[5], pars.GC_12, amp[53]);
  VVV1_0(w[33], w[5], w[34], pars.GC_10, amp[54]);
  VVV1_0(w[33], w[3], w[32], pars.GC_10, amp[55]);
  VVV1_0(w[29], w[33], w[13], pars.GC_10, amp[56]);
  VVVV1_0(w[29], w[35], w[3], w[4], pars.GC_12, amp[57]);
  VVVV3_0(w[29], w[35], w[3], w[4], pars.GC_12, amp[58]);
  VVVV4_0(w[29], w[35], w[3], w[4], pars.GC_12, amp[59]);
  VVV1_0(w[35], w[4], w[34], pars.GC_10, amp[60]);
  VVV1_0(w[35], w[3], w[31], pars.GC_10, amp[61]);
  VVV1_0(w[29], w[35], w[15], pars.GC_10, amp[62]);
  VVVV1_0(w[29], w[1], w[15], w[5], pars.GC_12, amp[63]);
  VVVV3_0(w[29], w[1], w[15], w[5], pars.GC_12, amp[64]);
  VVVV4_0(w[29], w[1], w[15], w[5], pars.GC_12, amp[65]);
  VVV1_0(w[15], w[5], w[36], pars.GC_10, amp[66]);
  VVV1_0(w[1], w[15], w[32], pars.GC_10, amp[67]);
  VVVV1_0(w[29], w[1], w[13], w[4], pars.GC_12, amp[68]);
  VVVV3_0(w[29], w[1], w[13], w[4], pars.GC_12, amp[69]);
  VVVV4_0(w[29], w[1], w[13], w[4], pars.GC_12, amp[70]);
  VVV1_0(w[13], w[4], w[36], pars.GC_10, amp[71]);
  VVV1_0(w[1], w[13], w[31], pars.GC_10, amp[72]);
  VVVV1_0(w[29], w[1], w[3], w[10], pars.GC_12, amp[73]);
  VVVV3_0(w[29], w[1], w[3], w[10], pars.GC_12, amp[74]);
  VVVV4_0(w[29], w[1], w[3], w[10], pars.GC_12, amp[75]);
  VVV1_0(w[3], w[10], w[36], pars.GC_10, amp[76]);
  VVV1_0(w[1], w[10], w[34], pars.GC_10, amp[77]);
  VVV1_0(w[29], w[37], w[5], pars.GC_10, amp[78]);
  VVV1_0(w[29], w[38], w[5], pars.GC_10, amp[79]);
  VVV1_0(w[29], w[39], w[5], pars.GC_10, amp[80]);
  VVV1_0(w[29], w[40], w[4], pars.GC_10, amp[81]);
  VVV1_0(w[29], w[41], w[4], pars.GC_10, amp[82]);
  VVV1_0(w[29], w[42], w[4], pars.GC_10, amp[83]);
  VVV1_0(w[29], w[43], w[3], pars.GC_10, amp[84]);
  VVV1_0(w[29], w[44], w[3], pars.GC_10, amp[85]);
  VVV1_0(w[29], w[45], w[3], pars.GC_10, amp[86]);
  VVV1_0(w[29], w[1], w[26], pars.GC_10, amp[87]);
  VVV1_0(w[29], w[1], w[27], pars.GC_10, amp[88]);
  VVV1_0(w[29], w[1], w[28], pars.GC_10, amp[89]);
  VVVV1_0(w[46], w[47], w[4], w[5], pars.GC_12, amp[90]);
  VVVV3_0(w[46], w[47], w[4], w[5], pars.GC_12, amp[91]);
  VVVV4_0(w[46], w[47], w[4], w[5], pars.GC_12, amp[92]);
  VVV1_0(w[47], w[5], w[48], pars.GC_10, amp[93]);
  VVV1_0(w[47], w[4], w[49], pars.GC_10, amp[94]);
  VVV1_0(w[46], w[47], w[10], pars.GC_10, amp[95]);
  VVVV1_0(w[46], w[33], w[2], w[5], pars.GC_12, amp[96]);
  VVVV3_0(w[46], w[33], w[2], w[5], pars.GC_12, amp[97]);
  VVVV4_0(w[46], w[33], w[2], w[5], pars.GC_12, amp[98]);
  VVV1_0(w[33], w[5], w[50], pars.GC_10, amp[99]);
  VVV1_0(w[33], w[2], w[49], pars.GC_10, amp[100]);
  VVV1_0(w[46], w[33], w[14], pars.GC_10, amp[101]);
  VVVV1_0(w[46], w[35], w[2], w[4], pars.GC_12, amp[102]);
  VVVV3_0(w[46], w[35], w[2], w[4], pars.GC_12, amp[103]);
  VVVV4_0(w[46], w[35], w[2], w[4], pars.GC_12, amp[104]);
  VVV1_0(w[35], w[4], w[50], pars.GC_10, amp[105]);
  VVV1_0(w[35], w[2], w[48], pars.GC_10, amp[106]);
  VVV1_0(w[46], w[35], w[11], pars.GC_10, amp[107]);
  VVVV1_0(w[46], w[1], w[11], w[5], pars.GC_12, amp[108]);
  VVVV3_0(w[46], w[1], w[11], w[5], pars.GC_12, amp[109]);
  VVVV4_0(w[46], w[1], w[11], w[5], pars.GC_12, amp[110]);
  VVV1_0(w[11], w[5], w[51], pars.GC_10, amp[111]);
  VVV1_0(w[1], w[11], w[49], pars.GC_10, amp[112]);
  VVVV1_0(w[46], w[1], w[14], w[4], pars.GC_12, amp[113]);
  VVVV3_0(w[46], w[1], w[14], w[4], pars.GC_12, amp[114]);
  VVVV4_0(w[46], w[1], w[14], w[4], pars.GC_12, amp[115]);
  VVV1_0(w[14], w[4], w[51], pars.GC_10, amp[116]);
  VVV1_0(w[1], w[14], w[48], pars.GC_10, amp[117]);
  VVVV1_0(w[46], w[1], w[2], w[10], pars.GC_12, amp[118]);
  VVVV3_0(w[46], w[1], w[2], w[10], pars.GC_12, amp[119]);
  VVVV4_0(w[46], w[1], w[2], w[10], pars.GC_12, amp[120]);
  VVV1_0(w[2], w[10], w[51], pars.GC_10, amp[121]);
  VVV1_0(w[1], w[10], w[50], pars.GC_10, amp[122]);
  VVV1_0(w[46], w[52], w[5], pars.GC_10, amp[123]);
  VVV1_0(w[46], w[53], w[5], pars.GC_10, amp[124]);
  VVV1_0(w[46], w[54], w[5], pars.GC_10, amp[125]);
  VVV1_0(w[46], w[55], w[4], pars.GC_10, amp[126]);
  VVV1_0(w[46], w[56], w[4], pars.GC_10, amp[127]);
  VVV1_0(w[46], w[57], w[4], pars.GC_10, amp[128]);
  VVV1_0(w[46], w[43], w[2], pars.GC_10, amp[129]);
  VVV1_0(w[46], w[44], w[2], pars.GC_10, amp[130]);
  VVV1_0(w[46], w[45], w[2], pars.GC_10, amp[131]);
  VVV1_0(w[46], w[1], w[23], pars.GC_10, amp[132]);
  VVV1_0(w[46], w[1], w[24], pars.GC_10, amp[133]);
  VVV1_0(w[46], w[1], w[25], pars.GC_10, amp[134]);
  VVVV1_0(w[58], w[47], w[3], w[5], pars.GC_12, amp[135]);
  VVVV3_0(w[58], w[47], w[3], w[5], pars.GC_12, amp[136]);
  VVVV4_0(w[58], w[47], w[3], w[5], pars.GC_12, amp[137]);
  VVV1_0(w[47], w[5], w[59], pars.GC_10, amp[138]);
  VVV1_0(w[47], w[3], w[60], pars.GC_10, amp[139]);
  VVV1_0(w[58], w[47], w[13], pars.GC_10, amp[140]);
  VVVV1_0(w[58], w[30], w[2], w[5], pars.GC_12, amp[141]);
  VVVV3_0(w[58], w[30], w[2], w[5], pars.GC_12, amp[142]);
  VVVV4_0(w[58], w[30], w[2], w[5], pars.GC_12, amp[143]);
  VVV1_0(w[30], w[5], w[61], pars.GC_10, amp[144]);
  VVV1_0(w[30], w[2], w[60], pars.GC_10, amp[145]);
  VVV1_0(w[58], w[30], w[14], pars.GC_10, amp[146]);
  VVVV1_0(w[58], w[35], w[2], w[3], pars.GC_12, amp[147]);
  VVVV3_0(w[58], w[35], w[2], w[3], pars.GC_12, amp[148]);
  VVVV4_0(w[58], w[35], w[2], w[3], pars.GC_12, amp[149]);
  VVV1_0(w[35], w[3], w[61], pars.GC_10, amp[150]);
  VVV1_0(w[35], w[2], w[59], pars.GC_10, amp[151]);
  VVV1_0(w[58], w[35], w[7], pars.GC_10, amp[152]);
  VVVV1_0(w[58], w[1], w[7], w[5], pars.GC_12, amp[153]);
  VVVV3_0(w[58], w[1], w[7], w[5], pars.GC_12, amp[154]);
  VVVV4_0(w[58], w[1], w[7], w[5], pars.GC_12, amp[155]);
  VVV1_0(w[7], w[5], w[62], pars.GC_10, amp[156]);
  VVV1_0(w[1], w[7], w[60], pars.GC_10, amp[157]);
  VVVV1_0(w[58], w[1], w[14], w[3], pars.GC_12, amp[158]);
  VVVV3_0(w[58], w[1], w[14], w[3], pars.GC_12, amp[159]);
  VVVV4_0(w[58], w[1], w[14], w[3], pars.GC_12, amp[160]);
  VVV1_0(w[14], w[3], w[62], pars.GC_10, amp[161]);
  VVV1_0(w[1], w[14], w[59], pars.GC_10, amp[162]);
  VVVV1_0(w[58], w[1], w[2], w[13], pars.GC_12, amp[163]);
  VVVV3_0(w[58], w[1], w[2], w[13], pars.GC_12, amp[164]);
  VVVV4_0(w[58], w[1], w[2], w[13], pars.GC_12, amp[165]);
  VVV1_0(w[2], w[13], w[62], pars.GC_10, amp[166]);
  VVV1_0(w[1], w[13], w[61], pars.GC_10, amp[167]);
  VVV1_0(w[58], w[63], w[5], pars.GC_10, amp[168]);
  VVV1_0(w[58], w[64], w[5], pars.GC_10, amp[169]);
  VVV1_0(w[58], w[65], w[5], pars.GC_10, amp[170]);
  VVV1_0(w[58], w[55], w[3], pars.GC_10, amp[171]);
  VVV1_0(w[58], w[56], w[3], pars.GC_10, amp[172]);
  VVV1_0(w[58], w[57], w[3], pars.GC_10, amp[173]);
  VVV1_0(w[58], w[40], w[2], pars.GC_10, amp[174]);
  VVV1_0(w[58], w[41], w[2], pars.GC_10, amp[175]);
  VVV1_0(w[58], w[42], w[2], pars.GC_10, amp[176]);
  VVV1_0(w[58], w[1], w[20], pars.GC_10, amp[177]);
  VVV1_0(w[58], w[1], w[21], pars.GC_10, amp[178]);
  VVV1_0(w[58], w[1], w[22], pars.GC_10, amp[179]);
  VVVV1_0(w[66], w[47], w[3], w[4], pars.GC_12, amp[180]);
  VVVV3_0(w[66], w[47], w[3], w[4], pars.GC_12, amp[181]);
  VVVV4_0(w[66], w[47], w[3], w[4], pars.GC_12, amp[182]);
  VVV1_0(w[47], w[4], w[67], pars.GC_10, amp[183]);
  VVV1_0(w[47], w[3], w[68], pars.GC_10, amp[184]);
  VVV1_0(w[66], w[47], w[15], pars.GC_10, amp[185]);
  VVVV1_0(w[66], w[30], w[2], w[4], pars.GC_12, amp[186]);
  VVVV3_0(w[66], w[30], w[2], w[4], pars.GC_12, amp[187]);
  VVVV4_0(w[66], w[30], w[2], w[4], pars.GC_12, amp[188]);
  VVV1_0(w[30], w[4], w[69], pars.GC_10, amp[189]);
  VVV1_0(w[30], w[2], w[68], pars.GC_10, amp[190]);
  VVV1_0(w[66], w[30], w[11], pars.GC_10, amp[191]);
  VVVV1_0(w[66], w[33], w[2], w[3], pars.GC_12, amp[192]);
  VVVV3_0(w[66], w[33], w[2], w[3], pars.GC_12, amp[193]);
  VVVV4_0(w[66], w[33], w[2], w[3], pars.GC_12, amp[194]);
  VVV1_0(w[33], w[3], w[69], pars.GC_10, amp[195]);
  VVV1_0(w[33], w[2], w[67], pars.GC_10, amp[196]);
  VVV1_0(w[66], w[33], w[7], pars.GC_10, amp[197]);
  VVVV1_0(w[66], w[1], w[7], w[4], pars.GC_12, amp[198]);
  VVVV3_0(w[66], w[1], w[7], w[4], pars.GC_12, amp[199]);
  VVVV4_0(w[66], w[1], w[7], w[4], pars.GC_12, amp[200]);
  VVV1_0(w[7], w[4], w[70], pars.GC_10, amp[201]);
  VVV1_0(w[1], w[7], w[68], pars.GC_10, amp[202]);
  VVVV1_0(w[66], w[1], w[11], w[3], pars.GC_12, amp[203]);
  VVVV3_0(w[66], w[1], w[11], w[3], pars.GC_12, amp[204]);
  VVVV4_0(w[66], w[1], w[11], w[3], pars.GC_12, amp[205]);
  VVV1_0(w[11], w[3], w[70], pars.GC_10, amp[206]);
  VVV1_0(w[1], w[11], w[67], pars.GC_10, amp[207]);
  VVVV1_0(w[66], w[1], w[2], w[15], pars.GC_12, amp[208]);
  VVVV3_0(w[66], w[1], w[2], w[15], pars.GC_12, amp[209]);
  VVVV4_0(w[66], w[1], w[2], w[15], pars.GC_12, amp[210]);
  VVV1_0(w[2], w[15], w[70], pars.GC_10, amp[211]);
  VVV1_0(w[1], w[15], w[69], pars.GC_10, amp[212]);
  VVV1_0(w[66], w[63], w[4], pars.GC_10, amp[213]);
  VVV1_0(w[66], w[64], w[4], pars.GC_10, amp[214]);
  VVV1_0(w[66], w[65], w[4], pars.GC_10, amp[215]);
  VVV1_0(w[66], w[52], w[3], pars.GC_10, amp[216]);
  VVV1_0(w[66], w[53], w[3], pars.GC_10, amp[217]);
  VVV1_0(w[66], w[54], w[3], pars.GC_10, amp[218]);
  VVV1_0(w[66], w[37], w[2], pars.GC_10, amp[219]);
  VVV1_0(w[66], w[38], w[2], pars.GC_10, amp[220]);
  VVV1_0(w[66], w[39], w[2], pars.GC_10, amp[221]);
  VVV1_0(w[66], w[1], w[17], pars.GC_10, amp[222]);
  VVV1_0(w[66], w[1], w[18], pars.GC_10, amp[223]);
  VVV1_0(w[66], w[1], w[19], pars.GC_10, amp[224]);
  VVVV1_0(w[0], w[47], w[15], w[5], pars.GC_12, amp[225]);
  VVVV3_0(w[0], w[47], w[15], w[5], pars.GC_12, amp[226]);
  VVVV4_0(w[0], w[47], w[15], w[5], pars.GC_12, amp[227]);
  VVV1_0(w[15], w[5], w[71], pars.GC_10, amp[228]);
  VVV1_0(w[47], w[5], w[72], pars.GC_10, amp[229]);
  VVVV1_0(w[0], w[47], w[13], w[4], pars.GC_12, amp[230]);
  VVVV3_0(w[0], w[47], w[13], w[4], pars.GC_12, amp[231]);
  VVVV4_0(w[0], w[47], w[13], w[4], pars.GC_12, amp[232]);
  VVV1_0(w[13], w[4], w[71], pars.GC_10, amp[233]);
  VVV1_0(w[47], w[4], w[73], pars.GC_10, amp[234]);
  VVVV1_0(w[0], w[47], w[3], w[10], pars.GC_12, amp[235]);
  VVVV3_0(w[0], w[47], w[3], w[10], pars.GC_12, amp[236]);
  VVVV4_0(w[0], w[47], w[3], w[10], pars.GC_12, amp[237]);
  VVV1_0(w[3], w[10], w[71], pars.GC_10, amp[238]);
  VVV1_0(w[47], w[3], w[74], pars.GC_10, amp[239]);
  VVV1_0(w[0], w[47], w[26], pars.GC_10, amp[240]);
  VVV1_0(w[0], w[47], w[27], pars.GC_10, amp[241]);
  VVV1_0(w[0], w[47], w[28], pars.GC_10, amp[242]);
  VVVV1_0(w[0], w[30], w[11], w[5], pars.GC_12, amp[243]);
  VVVV3_0(w[0], w[30], w[11], w[5], pars.GC_12, amp[244]);
  VVVV4_0(w[0], w[30], w[11], w[5], pars.GC_12, amp[245]);
  VVV1_0(w[11], w[5], w[75], pars.GC_10, amp[246]);
  VVV1_0(w[30], w[5], w[76], pars.GC_10, amp[247]);
  VVVV1_0(w[0], w[30], w[14], w[4], pars.GC_12, amp[248]);
  VVVV3_0(w[0], w[30], w[14], w[4], pars.GC_12, amp[249]);
  VVVV4_0(w[0], w[30], w[14], w[4], pars.GC_12, amp[250]);
  VVV1_0(w[14], w[4], w[75], pars.GC_10, amp[251]);
  VVV1_0(w[30], w[4], w[77], pars.GC_10, amp[252]);
  VVVV1_0(w[0], w[30], w[2], w[10], pars.GC_12, amp[253]);
  VVVV3_0(w[0], w[30], w[2], w[10], pars.GC_12, amp[254]);
  VVVV4_0(w[0], w[30], w[2], w[10], pars.GC_12, amp[255]);
  VVV1_0(w[2], w[10], w[75], pars.GC_10, amp[256]);
  VVV1_0(w[30], w[2], w[74], pars.GC_10, amp[257]);
  VVV1_0(w[0], w[30], w[23], pars.GC_10, amp[258]);
  VVV1_0(w[0], w[30], w[24], pars.GC_10, amp[259]);
  VVV1_0(w[0], w[30], w[25], pars.GC_10, amp[260]);
  VVVV1_0(w[0], w[33], w[7], w[5], pars.GC_12, amp[261]);
  VVVV3_0(w[0], w[33], w[7], w[5], pars.GC_12, amp[262]);
  VVVV4_0(w[0], w[33], w[7], w[5], pars.GC_12, amp[263]);
  VVV1_0(w[7], w[5], w[78], pars.GC_10, amp[264]);
  VVV1_0(w[33], w[5], w[79], pars.GC_10, amp[265]);
  VVVV1_0(w[0], w[33], w[14], w[3], pars.GC_12, amp[266]);
  VVVV3_0(w[0], w[33], w[14], w[3], pars.GC_12, amp[267]);
  VVVV4_0(w[0], w[33], w[14], w[3], pars.GC_12, amp[268]);
  VVV1_0(w[14], w[3], w[78], pars.GC_10, amp[269]);
  VVV1_0(w[33], w[3], w[77], pars.GC_10, amp[270]);
  VVVV1_0(w[0], w[33], w[2], w[13], pars.GC_12, amp[271]);
  VVVV3_0(w[0], w[33], w[2], w[13], pars.GC_12, amp[272]);
  VVVV4_0(w[0], w[33], w[2], w[13], pars.GC_12, amp[273]);
  VVV1_0(w[2], w[13], w[78], pars.GC_10, amp[274]);
  VVV1_0(w[33], w[2], w[73], pars.GC_10, amp[275]);
  VVV1_0(w[0], w[33], w[20], pars.GC_10, amp[276]);
  VVV1_0(w[0], w[33], w[21], pars.GC_10, amp[277]);
  VVV1_0(w[0], w[33], w[22], pars.GC_10, amp[278]);
  VVVV1_0(w[0], w[35], w[7], w[4], pars.GC_12, amp[279]);
  VVVV3_0(w[0], w[35], w[7], w[4], pars.GC_12, amp[280]);
  VVVV4_0(w[0], w[35], w[7], w[4], pars.GC_12, amp[281]);
  VVV1_0(w[7], w[4], w[80], pars.GC_10, amp[282]);
  VVV1_0(w[35], w[4], w[79], pars.GC_10, amp[283]);
  VVVV1_0(w[0], w[35], w[11], w[3], pars.GC_12, amp[284]);
  VVVV3_0(w[0], w[35], w[11], w[3], pars.GC_12, amp[285]);
  VVVV4_0(w[0], w[35], w[11], w[3], pars.GC_12, amp[286]);
  VVV1_0(w[11], w[3], w[80], pars.GC_10, amp[287]);
  VVV1_0(w[35], w[3], w[76], pars.GC_10, amp[288]);
  VVVV1_0(w[0], w[35], w[2], w[15], pars.GC_12, amp[289]);
  VVVV3_0(w[0], w[35], w[2], w[15], pars.GC_12, amp[290]);
  VVVV4_0(w[0], w[35], w[2], w[15], pars.GC_12, amp[291]);
  VVV1_0(w[2], w[15], w[80], pars.GC_10, amp[292]);
  VVV1_0(w[35], w[2], w[72], pars.GC_10, amp[293]);
  VVV1_0(w[0], w[35], w[17], pars.GC_10, amp[294]);
  VVV1_0(w[0], w[35], w[18], pars.GC_10, amp[295]);
  VVV1_0(w[0], w[35], w[19], pars.GC_10, amp[296]);
  VVVV1_0(w[0], w[1], w[7], w[10], pars.GC_12, amp[297]);
  VVVV3_0(w[0], w[1], w[7], w[10], pars.GC_12, amp[298]);
  VVVV4_0(w[0], w[1], w[7], w[10], pars.GC_12, amp[299]);
  VVV1_0(w[1], w[10], w[79], pars.GC_10, amp[300]);
  VVV1_0(w[1], w[7], w[74], pars.GC_10, amp[301]);
  VVVV1_0(w[0], w[1], w[11], w[13], pars.GC_12, amp[302]);
  VVVV3_0(w[0], w[1], w[11], w[13], pars.GC_12, amp[303]);
  VVVV4_0(w[0], w[1], w[11], w[13], pars.GC_12, amp[304]);
  VVV1_0(w[1], w[13], w[76], pars.GC_10, amp[305]);
  VVV1_0(w[1], w[11], w[73], pars.GC_10, amp[306]);
  VVVV1_0(w[0], w[1], w[14], w[15], pars.GC_12, amp[307]);
  VVVV3_0(w[0], w[1], w[14], w[15], pars.GC_12, amp[308]);
  VVVV4_0(w[0], w[1], w[14], w[15], pars.GC_12, amp[309]);
  VVV1_0(w[1], w[15], w[77], pars.GC_10, amp[310]);
  VVV1_0(w[1], w[14], w[72], pars.GC_10, amp[311]);
  VVV1_0(w[81], w[15], w[5], pars.GC_10, amp[312]);
  VVV1_0(w[82], w[15], w[5], pars.GC_10, amp[313]);
  VVV1_0(w[83], w[15], w[5], pars.GC_10, amp[314]);
  VVV1_0(w[81], w[13], w[4], pars.GC_10, amp[315]);
  VVV1_0(w[82], w[13], w[4], pars.GC_10, amp[316]);
  VVV1_0(w[83], w[13], w[4], pars.GC_10, amp[317]);
  VVV1_0(w[81], w[3], w[10], pars.GC_10, amp[318]);
  VVV1_0(w[82], w[3], w[10], pars.GC_10, amp[319]);
  VVV1_0(w[83], w[3], w[10], pars.GC_10, amp[320]);
  VVVV1_0(w[3], w[4], w[5], w[81], pars.GC_12, amp[321]);
  VVVV3_0(w[3], w[4], w[5], w[81], pars.GC_12, amp[322]);
  VVVV4_0(w[3], w[4], w[5], w[81], pars.GC_12, amp[323]);
  VVVV1_0(w[3], w[4], w[5], w[82], pars.GC_12, amp[324]);
  VVVV3_0(w[3], w[4], w[5], w[82], pars.GC_12, amp[325]);
  VVVV4_0(w[3], w[4], w[5], w[82], pars.GC_12, amp[326]);
  VVVV1_0(w[3], w[4], w[5], w[83], pars.GC_12, amp[327]);
  VVVV3_0(w[3], w[4], w[5], w[83], pars.GC_12, amp[328]);
  VVVV4_0(w[3], w[4], w[5], w[83], pars.GC_12, amp[329]);
  VVV1_0(w[84], w[11], w[5], pars.GC_10, amp[330]);
  VVV1_0(w[85], w[11], w[5], pars.GC_10, amp[331]);
  VVV1_0(w[86], w[11], w[5], pars.GC_10, amp[332]);
  VVV1_0(w[84], w[14], w[4], pars.GC_10, amp[333]);
  VVV1_0(w[85], w[14], w[4], pars.GC_10, amp[334]);
  VVV1_0(w[86], w[14], w[4], pars.GC_10, amp[335]);
  VVV1_0(w[84], w[2], w[10], pars.GC_10, amp[336]);
  VVV1_0(w[85], w[2], w[10], pars.GC_10, amp[337]);
  VVV1_0(w[86], w[2], w[10], pars.GC_10, amp[338]);
  VVVV1_0(w[2], w[4], w[5], w[84], pars.GC_12, amp[339]);
  VVVV3_0(w[2], w[4], w[5], w[84], pars.GC_12, amp[340]);
  VVVV4_0(w[2], w[4], w[5], w[84], pars.GC_12, amp[341]);
  VVVV1_0(w[2], w[4], w[5], w[85], pars.GC_12, amp[342]);
  VVVV3_0(w[2], w[4], w[5], w[85], pars.GC_12, amp[343]);
  VVVV4_0(w[2], w[4], w[5], w[85], pars.GC_12, amp[344]);
  VVVV1_0(w[2], w[4], w[5], w[86], pars.GC_12, amp[345]);
  VVVV3_0(w[2], w[4], w[5], w[86], pars.GC_12, amp[346]);
  VVVV4_0(w[2], w[4], w[5], w[86], pars.GC_12, amp[347]);
  VVV1_0(w[87], w[7], w[5], pars.GC_10, amp[348]);
  VVV1_0(w[88], w[7], w[5], pars.GC_10, amp[349]);
  VVV1_0(w[89], w[7], w[5], pars.GC_10, amp[350]);
  VVV1_0(w[87], w[14], w[3], pars.GC_10, amp[351]);
  VVV1_0(w[88], w[14], w[3], pars.GC_10, amp[352]);
  VVV1_0(w[89], w[14], w[3], pars.GC_10, amp[353]);
  VVV1_0(w[87], w[2], w[13], pars.GC_10, amp[354]);
  VVV1_0(w[88], w[2], w[13], pars.GC_10, amp[355]);
  VVV1_0(w[89], w[2], w[13], pars.GC_10, amp[356]);
  VVVV1_0(w[2], w[3], w[5], w[87], pars.GC_12, amp[357]);
  VVVV3_0(w[2], w[3], w[5], w[87], pars.GC_12, amp[358]);
  VVVV4_0(w[2], w[3], w[5], w[87], pars.GC_12, amp[359]);
  VVVV1_0(w[2], w[3], w[5], w[88], pars.GC_12, amp[360]);
  VVVV3_0(w[2], w[3], w[5], w[88], pars.GC_12, amp[361]);
  VVVV4_0(w[2], w[3], w[5], w[88], pars.GC_12, amp[362]);
  VVVV1_0(w[2], w[3], w[5], w[89], pars.GC_12, amp[363]);
  VVVV3_0(w[2], w[3], w[5], w[89], pars.GC_12, amp[364]);
  VVVV4_0(w[2], w[3], w[5], w[89], pars.GC_12, amp[365]);
  VVV1_0(w[90], w[7], w[4], pars.GC_10, amp[366]);
  VVV1_0(w[91], w[7], w[4], pars.GC_10, amp[367]);
  VVV1_0(w[92], w[7], w[4], pars.GC_10, amp[368]);
  VVV1_0(w[90], w[11], w[3], pars.GC_10, amp[369]);
  VVV1_0(w[91], w[11], w[3], pars.GC_10, amp[370]);
  VVV1_0(w[92], w[11], w[3], pars.GC_10, amp[371]);
  VVV1_0(w[90], w[2], w[15], pars.GC_10, amp[372]);
  VVV1_0(w[91], w[2], w[15], pars.GC_10, amp[373]);
  VVV1_0(w[92], w[2], w[15], pars.GC_10, amp[374]);
  VVVV1_0(w[2], w[3], w[4], w[90], pars.GC_12, amp[375]);
  VVVV3_0(w[2], w[3], w[4], w[90], pars.GC_12, amp[376]);
  VVVV4_0(w[2], w[3], w[4], w[90], pars.GC_12, amp[377]);
  VVVV1_0(w[2], w[3], w[4], w[91], pars.GC_12, amp[378]);
  VVVV3_0(w[2], w[3], w[4], w[91], pars.GC_12, amp[379]);
  VVVV4_0(w[2], w[3], w[4], w[91], pars.GC_12, amp[380]);
  VVVV1_0(w[2], w[3], w[4], w[92], pars.GC_12, amp[381]);
  VVVV3_0(w[2], w[3], w[4], w[92], pars.GC_12, amp[382]);
  VVVV4_0(w[2], w[3], w[4], w[92], pars.GC_12, amp[383]);
  VVV1_0(w[93], w[33], w[5], pars.GC_10, amp[384]);
  VVV1_0(w[94], w[33], w[5], pars.GC_10, amp[385]);
  VVV1_0(w[95], w[33], w[5], pars.GC_10, amp[386]);
  VVV1_0(w[93], w[35], w[4], pars.GC_10, amp[387]);
  VVV1_0(w[94], w[35], w[4], pars.GC_10, amp[388]);
  VVV1_0(w[95], w[35], w[4], pars.GC_10, amp[389]);
  VVV1_0(w[93], w[1], w[10], pars.GC_10, amp[390]);
  VVV1_0(w[94], w[1], w[10], pars.GC_10, amp[391]);
  VVV1_0(w[95], w[1], w[10], pars.GC_10, amp[392]);
  VVVV1_0(w[1], w[4], w[5], w[93], pars.GC_12, amp[393]);
  VVVV3_0(w[1], w[4], w[5], w[93], pars.GC_12, amp[394]);
  VVVV4_0(w[1], w[4], w[5], w[93], pars.GC_12, amp[395]);
  VVVV1_0(w[1], w[4], w[5], w[94], pars.GC_12, amp[396]);
  VVVV3_0(w[1], w[4], w[5], w[94], pars.GC_12, amp[397]);
  VVVV4_0(w[1], w[4], w[5], w[94], pars.GC_12, amp[398]);
  VVVV1_0(w[1], w[4], w[5], w[95], pars.GC_12, amp[399]);
  VVVV3_0(w[1], w[4], w[5], w[95], pars.GC_12, amp[400]);
  VVVV4_0(w[1], w[4], w[5], w[95], pars.GC_12, amp[401]);
  VVV1_0(w[96], w[30], w[5], pars.GC_10, amp[402]);
  VVV1_0(w[97], w[30], w[5], pars.GC_10, amp[403]);
  VVV1_0(w[98], w[30], w[5], pars.GC_10, amp[404]);
  VVV1_0(w[96], w[35], w[3], pars.GC_10, amp[405]);
  VVV1_0(w[97], w[35], w[3], pars.GC_10, amp[406]);
  VVV1_0(w[98], w[35], w[3], pars.GC_10, amp[407]);
  VVV1_0(w[96], w[1], w[13], pars.GC_10, amp[408]);
  VVV1_0(w[97], w[1], w[13], pars.GC_10, amp[409]);
  VVV1_0(w[98], w[1], w[13], pars.GC_10, amp[410]);
  VVVV1_0(w[1], w[3], w[5], w[96], pars.GC_12, amp[411]);
  VVVV3_0(w[1], w[3], w[5], w[96], pars.GC_12, amp[412]);
  VVVV4_0(w[1], w[3], w[5], w[96], pars.GC_12, amp[413]);
  VVVV1_0(w[1], w[3], w[5], w[97], pars.GC_12, amp[414]);
  VVVV3_0(w[1], w[3], w[5], w[97], pars.GC_12, amp[415]);
  VVVV4_0(w[1], w[3], w[5], w[97], pars.GC_12, amp[416]);
  VVVV1_0(w[1], w[3], w[5], w[98], pars.GC_12, amp[417]);
  VVVV3_0(w[1], w[3], w[5], w[98], pars.GC_12, amp[418]);
  VVVV4_0(w[1], w[3], w[5], w[98], pars.GC_12, amp[419]);
  VVV1_0(w[99], w[30], w[4], pars.GC_10, amp[420]);
  VVV1_0(w[100], w[30], w[4], pars.GC_10, amp[421]);
  VVV1_0(w[101], w[30], w[4], pars.GC_10, amp[422]);
  VVV1_0(w[99], w[33], w[3], pars.GC_10, amp[423]);
  VVV1_0(w[100], w[33], w[3], pars.GC_10, amp[424]);
  VVV1_0(w[101], w[33], w[3], pars.GC_10, amp[425]);
  VVV1_0(w[99], w[1], w[15], pars.GC_10, amp[426]);
  VVV1_0(w[100], w[1], w[15], pars.GC_10, amp[427]);
  VVV1_0(w[101], w[1], w[15], pars.GC_10, amp[428]);
  VVVV1_0(w[1], w[3], w[4], w[99], pars.GC_12, amp[429]);
  VVVV3_0(w[1], w[3], w[4], w[99], pars.GC_12, amp[430]);
  VVVV4_0(w[1], w[3], w[4], w[99], pars.GC_12, amp[431]);
  VVVV1_0(w[1], w[3], w[4], w[100], pars.GC_12, amp[432]);
  VVVV3_0(w[1], w[3], w[4], w[100], pars.GC_12, amp[433]);
  VVVV4_0(w[1], w[3], w[4], w[100], pars.GC_12, amp[434]);
  VVVV1_0(w[1], w[3], w[4], w[101], pars.GC_12, amp[435]);
  VVVV3_0(w[1], w[3], w[4], w[101], pars.GC_12, amp[436]);
  VVVV4_0(w[1], w[3], w[4], w[101], pars.GC_12, amp[437]);
  VVV1_0(w[102], w[47], w[5], pars.GC_10, amp[438]);
  VVV1_0(w[103], w[47], w[5], pars.GC_10, amp[439]);
  VVV1_0(w[104], w[47], w[5], pars.GC_10, amp[440]);
  VVV1_0(w[102], w[35], w[2], pars.GC_10, amp[441]);
  VVV1_0(w[103], w[35], w[2], pars.GC_10, amp[442]);
  VVV1_0(w[104], w[35], w[2], pars.GC_10, amp[443]);
  VVV1_0(w[102], w[1], w[14], pars.GC_10, amp[444]);
  VVV1_0(w[103], w[1], w[14], pars.GC_10, amp[445]);
  VVV1_0(w[104], w[1], w[14], pars.GC_10, amp[446]);
  VVVV1_0(w[1], w[2], w[5], w[102], pars.GC_12, amp[447]);
  VVVV3_0(w[1], w[2], w[5], w[102], pars.GC_12, amp[448]);
  VVVV4_0(w[1], w[2], w[5], w[102], pars.GC_12, amp[449]);
  VVVV1_0(w[1], w[2], w[5], w[103], pars.GC_12, amp[450]);
  VVVV3_0(w[1], w[2], w[5], w[103], pars.GC_12, amp[451]);
  VVVV4_0(w[1], w[2], w[5], w[103], pars.GC_12, amp[452]);
  VVVV1_0(w[1], w[2], w[5], w[104], pars.GC_12, amp[453]);
  VVVV3_0(w[1], w[2], w[5], w[104], pars.GC_12, amp[454]);
  VVVV4_0(w[1], w[2], w[5], w[104], pars.GC_12, amp[455]);
  VVV1_0(w[105], w[47], w[4], pars.GC_10, amp[456]);
  VVV1_0(w[106], w[47], w[4], pars.GC_10, amp[457]);
  VVV1_0(w[107], w[47], w[4], pars.GC_10, amp[458]);
  VVV1_0(w[105], w[33], w[2], pars.GC_10, amp[459]);
  VVV1_0(w[106], w[33], w[2], pars.GC_10, amp[460]);
  VVV1_0(w[107], w[33], w[2], pars.GC_10, amp[461]);
  VVV1_0(w[105], w[1], w[11], pars.GC_10, amp[462]);
  VVV1_0(w[106], w[1], w[11], pars.GC_10, amp[463]);
  VVV1_0(w[107], w[1], w[11], pars.GC_10, amp[464]);
  VVVV1_0(w[1], w[2], w[4], w[105], pars.GC_12, amp[465]);
  VVVV3_0(w[1], w[2], w[4], w[105], pars.GC_12, amp[466]);
  VVVV4_0(w[1], w[2], w[4], w[105], pars.GC_12, amp[467]);
  VVVV1_0(w[1], w[2], w[4], w[106], pars.GC_12, amp[468]);
  VVVV3_0(w[1], w[2], w[4], w[106], pars.GC_12, amp[469]);
  VVVV4_0(w[1], w[2], w[4], w[106], pars.GC_12, amp[470]);
  VVVV1_0(w[1], w[2], w[4], w[107], pars.GC_12, amp[471]);
  VVVV3_0(w[1], w[2], w[4], w[107], pars.GC_12, amp[472]);
  VVVV4_0(w[1], w[2], w[4], w[107], pars.GC_12, amp[473]);
  VVV1_0(w[108], w[47], w[3], pars.GC_10, amp[474]);
  VVV1_0(w[109], w[47], w[3], pars.GC_10, amp[475]);
  VVV1_0(w[110], w[47], w[3], pars.GC_10, amp[476]);
  VVV1_0(w[108], w[30], w[2], pars.GC_10, amp[477]);
  VVV1_0(w[109], w[30], w[2], pars.GC_10, amp[478]);
  VVV1_0(w[110], w[30], w[2], pars.GC_10, amp[479]);
  VVV1_0(w[108], w[1], w[7], pars.GC_10, amp[480]);
  VVV1_0(w[109], w[1], w[7], pars.GC_10, amp[481]);
  VVV1_0(w[110], w[1], w[7], pars.GC_10, amp[482]);
  VVVV1_0(w[1], w[2], w[3], w[108], pars.GC_12, amp[483]);
  VVVV3_0(w[1], w[2], w[3], w[108], pars.GC_12, amp[484]);
  VVVV4_0(w[1], w[2], w[3], w[108], pars.GC_12, amp[485]);
  VVVV1_0(w[1], w[2], w[3], w[109], pars.GC_12, amp[486]);
  VVVV3_0(w[1], w[2], w[3], w[109], pars.GC_12, amp[487]);
  VVVV4_0(w[1], w[2], w[3], w[109], pars.GC_12, amp[488]);
  VVVV1_0(w[1], w[2], w[3], w[110], pars.GC_12, amp[489]);
  VVVV3_0(w[1], w[2], w[3], w[110], pars.GC_12, amp[490]);
  VVVV4_0(w[1], w[2], w[3], w[110], pars.GC_12, amp[491]);
  VVV1_0(w[0], w[63], w[10], pars.GC_10, amp[492]);
  VVV1_0(w[0], w[64], w[10], pars.GC_10, amp[493]);
  VVV1_0(w[0], w[65], w[10], pars.GC_10, amp[494]);
  VVV1_0(w[0], w[52], w[13], pars.GC_10, amp[495]);
  VVV1_0(w[0], w[53], w[13], pars.GC_10, amp[496]);
  VVV1_0(w[0], w[54], w[13], pars.GC_10, amp[497]);
  VVV1_0(w[0], w[55], w[15], pars.GC_10, amp[498]);
  VVV1_0(w[0], w[56], w[15], pars.GC_10, amp[499]);
  VVV1_0(w[0], w[57], w[15], pars.GC_10, amp[500]);
  VVV1_0(w[0], w[37], w[14], pars.GC_10, amp[501]);
  VVV1_0(w[0], w[38], w[14], pars.GC_10, amp[502]);
  VVV1_0(w[0], w[39], w[14], pars.GC_10, amp[503]);
  VVV1_0(w[0], w[40], w[11], pars.GC_10, amp[504]);
  VVV1_0(w[0], w[41], w[11], pars.GC_10, amp[505]);
  VVV1_0(w[0], w[42], w[11], pars.GC_10, amp[506]);
  VVV1_0(w[0], w[43], w[7], pars.GC_10, amp[507]);
  VVV1_0(w[0], w[44], w[7], pars.GC_10, amp[508]);
  VVV1_0(w[0], w[45], w[7], pars.GC_10, amp[509]);

}

double AMP_MG5_gg_gggg::matrix_1_gg_gggg() {
  int i, j;
  // Local variables
  const int            ngraphs = 510;
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

void AMP_MG5_gg_gggg::calculate_color_flows(std::complex<double> jamp[ncolor]) const {
    jamp[0] = +2. * (+amp[0] - amp[2] - amp[4] + amp[5] + amp[18] - amp[20] +
        amp[21] - amp[22] + amp[28] - amp[30] + amp[31] + amp[33] - amp[35] +
        amp[42] - amp[44] - amp[180] + amp[182] + amp[184] - amp[185] - amp[198]
        + amp[200] - amp[201] + amp[202] - amp[208] + amp[210] - amp[211] -
        amp[213] + amp[215] - amp[222] + amp[224] + amp[225] - amp[227] +
        amp[228] + amp[235] - amp[237] + amp[238] - amp[239] + amp[240] -
        amp[242] + amp[297] - amp[299] - amp[301] + amp[312] - amp[314] +
        amp[318] - amp[320] + amp[321] - amp[327] - amp[323] + amp[329] +
        amp[367] + amp[368] + amp[373] + amp[374] + amp[378] + amp[381] -
        amp[380] - amp[383] - amp[474] - amp[475] - amp[480] - amp[481] -
        amp[483] - amp[486] + amp[485] + amp[488] + amp[492] - amp[494]);
    jamp[1] = +2. * (-amp[0] - amp[1] - amp[3] - amp[5] + amp[23] - amp[25] +
        amp[26] - amp[27] - amp[28] + amp[30] - amp[31] + amp[36] - amp[38] -
        amp[42] - amp[43] - amp[135] + amp[137] + amp[139] - amp[140] - amp[153]
        + amp[155] - amp[156] + amp[157] - amp[163] + amp[165] - amp[166] -
        amp[168] + amp[170] - amp[177] + amp[179] + amp[230] - amp[232] +
        amp[233] - amp[235] + amp[237] - amp[238] + amp[239] - amp[240] -
        amp[241] - amp[297] + amp[299] + amp[301] + amp[315] - amp[317] -
        amp[318] + amp[320] + amp[322] - amp[328] + amp[323] - amp[329] +
        amp[349] + amp[350] + amp[355] + amp[356] + amp[360] + amp[363] -
        amp[362] - amp[365] + amp[474] - amp[476] + amp[480] - amp[482] +
        amp[483] - amp[489] - amp[485] + amp[491] - amp[492] + amp[494]);
    jamp[2] = +2. * (+amp[6] - amp[8] - amp[10] + amp[11] - amp[18] + amp[20] -
        amp[21] + amp[22] - amp[23] - amp[24] - amp[26] - amp[33] - amp[34] +
        amp[43] + amp[44] + amp[180] + amp[181] + amp[183] + amp[185] - amp[203]
        + amp[205] - amp[206] + amp[207] + amp[208] - amp[210] + amp[211] -
        amp[216] + amp[218] + amp[222] + amp[223] - amp[225] + amp[227] -
        amp[228] - amp[230] - amp[231] - amp[233] - amp[234] + amp[241] +
        amp[242] + amp[302] - amp[304] - amp[306] - amp[312] + amp[314] -
        amp[315] + amp[317] - amp[321] + amp[327] - amp[322] + amp[328] +
        amp[370] + amp[371] - amp[373] - amp[374] + amp[379] + amp[382] +
        amp[380] + amp[383] - amp[456] - amp[457] - amp[462] - amp[463] -
        amp[465] - amp[468] + amp[467] + amp[470] + amp[495] - amp[497]);
    jamp[3] = +2. * (-amp[6] - amp[7] - amp[9] - amp[11] + amp[23] + amp[24] +
        amp[26] - amp[28] - amp[29] - amp[31] - amp[32] + amp[39] - amp[41] -
        amp[42] - amp[43] - amp[90] + amp[92] + amp[94] - amp[95] - amp[108] +
        amp[110] - amp[111] + amp[112] - amp[118] + amp[120] - amp[121] -
        amp[123] + amp[125] - amp[132] + amp[134] + amp[230] + amp[231] +
        amp[233] + amp[234] - amp[235] - amp[236] - amp[238] - amp[240] -
        amp[241] - amp[302] + amp[304] + amp[306] + amp[315] - amp[317] -
        amp[318] + amp[320] + amp[322] - amp[328] + amp[323] - amp[329] +
        amp[331] + amp[332] + amp[337] + amp[338] + amp[342] + amp[345] -
        amp[344] - amp[347] + amp[456] - amp[458] + amp[462] - amp[464] +
        amp[465] - amp[471] - amp[467] + amp[473] - amp[495] + amp[497]);
    jamp[4] = +2. * (+amp[12] - amp[14] - amp[16] + amp[17] - amp[18] - amp[19] -
        amp[21] - amp[23] + amp[25] - amp[26] + amp[27] - amp[36] - amp[37] +
        amp[43] + amp[44] + amp[135] + amp[136] + amp[138] + amp[140] - amp[158]
        + amp[160] - amp[161] + amp[162] + amp[163] - amp[165] + amp[166] -
        amp[171] + amp[173] + amp[177] + amp[178] - amp[225] - amp[226] -
        amp[228] - amp[229] - amp[230] + amp[232] - amp[233] + amp[241] +
        amp[242] + amp[307] - amp[309] - amp[311] - amp[312] + amp[314] -
        amp[315] + amp[317] - amp[321] + amp[327] - amp[322] + amp[328] +
        amp[352] + amp[353] - amp[355] - amp[356] + amp[361] + amp[364] +
        amp[362] + amp[365] - amp[438] - amp[439] - amp[444] - amp[445] -
        amp[447] - amp[450] + amp[449] + amp[452] + amp[498] - amp[500]);
    jamp[5] = +2. * (-amp[12] - amp[13] - amp[15] - amp[17] + amp[18] + amp[19] +
        amp[21] + amp[28] + amp[29] + amp[31] + amp[32] - amp[39] - amp[40] +
        amp[42] - amp[44] + amp[90] + amp[91] + amp[93] + amp[95] - amp[113] +
        amp[115] - amp[116] + amp[117] + amp[118] - amp[120] + amp[121] -
        amp[126] + amp[128] + amp[132] + amp[133] + amp[225] + amp[226] +
        amp[228] + amp[229] + amp[235] + amp[236] + amp[238] + amp[240] -
        amp[242] - amp[307] + amp[309] + amp[311] + amp[312] - amp[314] +
        amp[318] - amp[320] + amp[321] - amp[327] - amp[323] + amp[329] +
        amp[334] + amp[335] - amp[337] - amp[338] + amp[343] + amp[346] +
        amp[344] + amp[347] + amp[438] - amp[440] + amp[444] - amp[446] +
        amp[447] - amp[453] - amp[449] + amp[455] - amp[498] + amp[500]);
    jamp[6] = +2. * (-amp[0] + amp[2] + amp[4] - amp[5] + amp[7] + amp[8] +
        amp[9] + amp[10] + amp[29] + amp[30] + amp[32] + amp[34] + amp[35] -
        amp[39] + amp[41] - amp[186] + amp[188] + amp[190] - amp[191] + amp[198]
        - amp[200] + amp[201] - amp[202] + amp[203] + amp[204] + amp[206] +
        amp[213] + amp[214] - amp[223] - amp[224] + amp[243] - amp[245] +
        amp[246] + amp[253] - amp[255] + amp[256] - amp[257] + amp[258] -
        amp[260] - amp[297] + amp[299] + amp[301] + amp[330] - amp[332] +
        amp[336] - amp[338] + amp[339] - amp[345] - amp[341] + amp[347] -
        amp[367] - amp[368] - amp[370] - amp[371] - amp[378] - amp[381] -
        amp[379] - amp[382] - amp[477] - amp[478] + amp[480] + amp[481] -
        amp[484] - amp[487] - amp[485] - amp[488] - amp[492] - amp[493]);
    jamp[7] = +2. * (+amp[0] + amp[1] + amp[3] + amp[5] + amp[13] + amp[14] +
        amp[15] + amp[16] - amp[29] - amp[30] - amp[32] + amp[37] + amp[38] +
        amp[39] + amp[40] - amp[141] + amp[143] + amp[145] - amp[146] + amp[153]
        - amp[155] + amp[156] - amp[157] + amp[158] + amp[159] + amp[161] +
        amp[168] + amp[169] - amp[178] - amp[179] + amp[248] - amp[250] +
        amp[251] - amp[253] + amp[255] - amp[256] + amp[257] - amp[258] -
        amp[259] + amp[297] - amp[299] - amp[301] + amp[333] - amp[335] -
        amp[336] + amp[338] + amp[340] - amp[346] + amp[341] - amp[347] -
        amp[349] - amp[350] - amp[352] - amp[353] - amp[360] - amp[363] -
        amp[361] - amp[364] + amp[477] - amp[479] - amp[480] + amp[482] +
        amp[484] - amp[490] + amp[485] - amp[491] + amp[492] + amp[493]);
    jamp[8] = +2. * (-amp[7] - amp[8] - amp[9] - amp[10] - amp[12] - amp[13] -
        amp[15] - amp[17] + amp[19] + amp[20] + amp[22] - amp[33] - amp[34] -
        amp[40] - amp[41] + amp[186] + amp[187] + amp[189] + amp[191] - amp[203]
        - amp[204] - amp[206] + amp[208] + amp[209] + amp[211] + amp[212] -
        amp[219] + amp[221] + amp[222] + amp[223] - amp[243] + amp[245] -
        amp[246] - amp[248] - amp[249] - amp[251] - amp[252] + amp[259] +
        amp[260] - amp[307] - amp[308] - amp[310] - amp[330] + amp[332] -
        amp[333] + amp[335] - amp[339] + amp[345] - amp[340] + amp[346] +
        amp[370] + amp[371] - amp[373] - amp[374] + amp[379] + amp[382] +
        amp[380] + amp[383] - amp[420] - amp[421] - amp[426] - amp[427] -
        amp[429] - amp[432] + amp[431] + amp[434] + amp[501] - amp[503]);
    jamp[9] = +2. * (+amp[12] + amp[13] + amp[15] + amp[17] - amp[18] - amp[19] -
        amp[21] - amp[28] - amp[29] - amp[31] - amp[32] + amp[39] + amp[40] -
        amp[42] + amp[44] - amp[45] + amp[47] + amp[49] - amp[50] - amp[63] +
        amp[65] - amp[66] + amp[67] - amp[73] + amp[75] - amp[76] - amp[78] +
        amp[80] - amp[87] + amp[89] + amp[248] + amp[249] + amp[251] + amp[252] -
        amp[253] - amp[254] - amp[256] - amp[258] - amp[259] + amp[307] +
        amp[308] + amp[310] + amp[313] + amp[314] + amp[319] + amp[320] +
        amp[324] + amp[327] - amp[326] - amp[329] + amp[333] - amp[335] -
        amp[336] + amp[338] + amp[340] - amp[346] + amp[341] - amp[347] +
        amp[420] - amp[422] + amp[426] - amp[428] + amp[429] - amp[435] -
        amp[431] + amp[437] - amp[501] + amp[503]);
    jamp[10] = +2. * (-amp[6] - amp[7] - amp[9] - amp[11] - amp[13] - amp[14] -
        amp[15] - amp[16] + amp[24] + amp[25] + amp[27] - amp[36] - amp[37] -
        amp[40] - amp[41] + amp[141] + amp[142] + amp[144] + amp[146] - amp[158]
        - amp[159] - amp[161] + amp[163] + amp[164] + amp[166] + amp[167] -
        amp[174] + amp[176] + amp[177] + amp[178] - amp[243] - amp[244] -
        amp[246] - amp[247] - amp[248] + amp[250] - amp[251] + amp[259] +
        amp[260] - amp[302] - amp[303] - amp[305] - amp[330] + amp[332] -
        amp[333] + amp[335] - amp[339] + amp[345] - amp[340] + amp[346] +
        amp[352] + amp[353] - amp[355] - amp[356] + amp[361] + amp[364] +
        amp[362] + amp[365] - amp[402] - amp[403] - amp[408] - amp[409] -
        amp[411] - amp[414] + amp[413] + amp[416] + amp[504] - amp[506]);
    jamp[11] = +2. * (+amp[6] + amp[7] + amp[9] + amp[11] - amp[23] - amp[24] -
        amp[26] + amp[28] + amp[29] + amp[31] + amp[32] - amp[39] + amp[41] +
        amp[42] + amp[43] + amp[45] + amp[46] + amp[48] + amp[50] - amp[68] +
        amp[70] - amp[71] + amp[72] + amp[73] - amp[75] + amp[76] - amp[81] +
        amp[83] + amp[87] + amp[88] + amp[243] + amp[244] + amp[246] + amp[247] +
        amp[253] + amp[254] + amp[256] + amp[258] - amp[260] + amp[302] +
        amp[303] + amp[305] + amp[316] + amp[317] - amp[319] - amp[320] +
        amp[325] + amp[328] + amp[326] + amp[329] + amp[330] - amp[332] +
        amp[336] - amp[338] + amp[339] - amp[345] - amp[341] + amp[347] +
        amp[402] - amp[404] + amp[408] - amp[410] + amp[411] - amp[417] -
        amp[413] + amp[419] - amp[504] + amp[506]);
    jamp[12] = +2. * (+amp[1] + amp[2] + amp[3] + amp[4] - amp[6] + amp[8] +
        amp[10] - amp[11] + amp[24] + amp[25] + amp[27] + amp[34] + amp[35] -
        amp[36] + amp[38] - amp[192] + amp[194] + amp[196] - amp[197] + amp[198]
        + amp[199] + amp[201] + amp[203] - amp[205] + amp[206] - amp[207] +
        amp[216] + amp[217] - amp[223] - amp[224] + amp[261] - amp[263] +
        amp[264] + amp[271] - amp[273] + amp[274] - amp[275] + amp[276] -
        amp[278] - amp[302] + amp[304] + amp[306] + amp[348] - amp[350] +
        amp[354] - amp[356] + amp[357] - amp[363] - amp[359] + amp[365] -
        amp[367] - amp[368] - amp[370] - amp[371] - amp[378] - amp[381] -
        amp[379] - amp[382] - amp[459] - amp[460] + amp[462] + amp[463] -
        amp[466] - amp[469] - amp[467] - amp[470] - amp[495] - amp[496]);
    jamp[13] = +2. * (+amp[6] + amp[7] + amp[9] + amp[11] + amp[13] + amp[14] +
        amp[15] + amp[16] - amp[24] - amp[25] - amp[27] + amp[36] + amp[37] +
        amp[40] + amp[41] - amp[96] + amp[98] + amp[100] - amp[101] + amp[108] -
        amp[110] + amp[111] - amp[112] + amp[113] + amp[114] + amp[116] +
        amp[123] + amp[124] - amp[133] - amp[134] + amp[266] - amp[268] +
        amp[269] - amp[271] + amp[273] - amp[274] + amp[275] - amp[276] -
        amp[277] + amp[302] - amp[304] - amp[306] - amp[331] - amp[332] -
        amp[334] - amp[335] - amp[342] - amp[345] - amp[343] - amp[346] +
        amp[351] - amp[353] - amp[354] + amp[356] + amp[358] - amp[364] +
        amp[359] - amp[365] + amp[459] - amp[461] - amp[462] + amp[464] +
        amp[466] - amp[472] + amp[467] - amp[473] + amp[495] + amp[496]);
    jamp[14] = +2. * (-amp[1] - amp[2] - amp[3] - amp[4] + amp[12] - amp[14] -
        amp[16] + amp[17] - amp[19] - amp[20] - amp[22] + amp[33] - amp[35] -
        amp[37] - amp[38] + amp[192] + amp[193] + amp[195] + amp[197] - amp[198]
        - amp[199] - amp[201] - amp[208] - amp[209] - amp[211] - amp[212] +
        amp[219] + amp[220] - amp[222] + amp[224] - amp[261] + amp[263] -
        amp[264] - amp[266] - amp[267] - amp[269] - amp[270] + amp[277] +
        amp[278] + amp[307] + amp[308] + amp[310] - amp[348] + amp[350] -
        amp[351] + amp[353] - amp[357] + amp[363] - amp[358] + amp[364] +
        amp[367] + amp[368] + amp[373] + amp[374] + amp[378] + amp[381] -
        amp[380] - amp[383] - amp[423] - amp[424] + amp[426] + amp[427] -
        amp[430] - amp[433] - amp[431] - amp[434] - amp[501] - amp[502]);
    jamp[15] = +2. * (-amp[12] + amp[14] + amp[16] - amp[17] + amp[18] + amp[19]
        + amp[21] + amp[23] - amp[25] + amp[26] - amp[27] + amp[36] + amp[37] -
        amp[43] - amp[44] - amp[51] + amp[53] + amp[55] - amp[56] + amp[63] -
        amp[65] + amp[66] - amp[67] + amp[68] + amp[69] + amp[71] + amp[78] +
        amp[79] - amp[88] - amp[89] + amp[266] + amp[267] + amp[269] + amp[270] -
        amp[271] - amp[272] - amp[274] - amp[276] - amp[277] - amp[307] -
        amp[308] - amp[310] - amp[313] - amp[314] - amp[316] - amp[317] -
        amp[324] - amp[327] - amp[325] - amp[328] + amp[351] - amp[353] -
        amp[354] + amp[356] + amp[358] - amp[364] + amp[359] - amp[365] +
        amp[423] - amp[425] - amp[426] + amp[428] + amp[430] - amp[436] +
        amp[431] - amp[437] + amp[501] + amp[502]);
    jamp[16] = +2. * (-amp[0] - amp[1] - amp[3] - amp[5] - amp[13] - amp[14] -
        amp[15] - amp[16] + amp[29] + amp[30] + amp[32] - amp[37] - amp[38] -
        amp[39] - amp[40] + amp[96] + amp[97] + amp[99] + amp[101] - amp[113] -
        amp[114] - amp[116] + amp[118] + amp[119] + amp[121] + amp[122] -
        amp[129] + amp[131] + amp[132] + amp[133] - amp[261] - amp[262] -
        amp[264] - amp[265] - amp[266] + amp[268] - amp[269] + amp[277] +
        amp[278] - amp[297] - amp[298] - amp[300] + amp[334] + amp[335] -
        amp[337] - amp[338] + amp[343] + amp[346] + amp[344] + amp[347] -
        amp[348] + amp[350] - amp[351] + amp[353] - amp[357] + amp[363] -
        amp[358] + amp[364] - amp[384] - amp[385] - amp[390] - amp[391] -
        amp[393] - amp[396] + amp[395] + amp[398] + amp[507] - amp[509]);
    jamp[17] = +2. * (+amp[0] + amp[1] + amp[3] + amp[5] - amp[23] + amp[25] -
        amp[26] + amp[27] + amp[28] - amp[30] + amp[31] - amp[36] + amp[38] +
        amp[42] + amp[43] + amp[51] + amp[52] + amp[54] + amp[56] - amp[68] -
        amp[69] - amp[71] + amp[73] + amp[74] + amp[76] + amp[77] - amp[84] +
        amp[86] + amp[87] + amp[88] + amp[261] + amp[262] + amp[264] + amp[265] +
        amp[271] + amp[272] + amp[274] + amp[276] - amp[278] + amp[297] +
        amp[298] + amp[300] + amp[316] + amp[317] - amp[319] - amp[320] +
        amp[325] + amp[328] + amp[326] + amp[329] + amp[348] - amp[350] +
        amp[354] - amp[356] + amp[357] - amp[363] - amp[359] + amp[365] +
        amp[384] - amp[386] + amp[390] - amp[392] + amp[393] - amp[399] -
        amp[395] + amp[401] - amp[507] + amp[509]);
    jamp[18] = +2. * (+amp[1] + amp[2] + amp[3] + amp[4] - amp[12] + amp[14] +
        amp[16] - amp[17] + amp[19] + amp[20] + amp[22] - amp[33] + amp[35] +
        amp[37] + amp[38] - amp[147] + amp[149] + amp[151] - amp[152] + amp[153]
        + amp[154] + amp[156] + amp[158] - amp[160] + amp[161] - amp[162] +
        amp[171] + amp[172] - amp[178] - amp[179] + amp[279] - amp[281] +
        amp[282] + amp[289] - amp[291] + amp[292] - amp[293] + amp[294] -
        amp[296] - amp[307] + amp[309] + amp[311] - amp[349] - amp[350] -
        amp[352] - amp[353] - amp[360] - amp[363] - amp[361] - amp[364] +
        amp[366] - amp[368] + amp[372] - amp[374] + amp[375] - amp[381] -
        amp[377] + amp[383] - amp[441] - amp[442] + amp[444] + amp[445] -
        amp[448] - amp[451] - amp[449] - amp[452] - amp[498] - amp[499]);
    jamp[19] = +2. * (+amp[7] + amp[8] + amp[9] + amp[10] + amp[12] + amp[13] +
        amp[15] + amp[17] - amp[19] - amp[20] - amp[22] + amp[33] + amp[34] +
        amp[40] + amp[41] - amp[102] + amp[104] + amp[106] - amp[107] + amp[108]
        + amp[109] + amp[111] + amp[113] - amp[115] + amp[116] - amp[117] +
        amp[126] + amp[127] - amp[133] - amp[134] + amp[284] - amp[286] +
        amp[287] - amp[289] + amp[291] - amp[292] + amp[293] - amp[294] -
        amp[295] + amp[307] - amp[309] - amp[311] - amp[331] - amp[332] -
        amp[334] - amp[335] - amp[342] - amp[345] - amp[343] - amp[346] +
        amp[369] - amp[371] - amp[372] + amp[374] + amp[376] - amp[382] +
        amp[377] - amp[383] + amp[441] - amp[443] - amp[444] + amp[446] +
        amp[448] - amp[454] + amp[449] - amp[455] + amp[498] + amp[499]);
    jamp[20] = +2. * (-amp[1] - amp[2] - amp[3] - amp[4] + amp[6] - amp[8] -
        amp[10] + amp[11] - amp[24] - amp[25] - amp[27] - amp[34] - amp[35] +
        amp[36] - amp[38] + amp[147] + amp[148] + amp[150] + amp[152] - amp[153]
        - amp[154] - amp[156] - amp[163] - amp[164] - amp[166] - amp[167] +
        amp[174] + amp[175] - amp[177] + amp[179] - amp[279] + amp[281] -
        amp[282] - amp[284] - amp[285] - amp[287] - amp[288] + amp[295] +
        amp[296] + amp[302] + amp[303] + amp[305] + amp[349] + amp[350] +
        amp[355] + amp[356] + amp[360] + amp[363] - amp[362] - amp[365] -
        amp[366] + amp[368] - amp[369] + amp[371] - amp[375] + amp[381] -
        amp[376] + amp[382] - amp[405] - amp[406] + amp[408] + amp[409] -
        amp[412] - amp[415] - amp[413] - amp[416] - amp[504] - amp[505]);
    jamp[21] = +2. * (-amp[6] + amp[8] + amp[10] - amp[11] + amp[18] - amp[20] +
        amp[21] - amp[22] + amp[23] + amp[24] + amp[26] + amp[33] + amp[34] -
        amp[43] - amp[44] - amp[57] + amp[59] + amp[61] - amp[62] + amp[63] +
        amp[64] + amp[66] + amp[68] - amp[70] + amp[71] - amp[72] + amp[81] +
        amp[82] - amp[88] - amp[89] + amp[284] + amp[285] + amp[287] + amp[288] -
        amp[289] - amp[290] - amp[292] - amp[294] - amp[295] - amp[302] -
        amp[303] - amp[305] - amp[313] - amp[314] - amp[316] - amp[317] -
        amp[324] - amp[327] - amp[325] - amp[328] + amp[369] - amp[371] -
        amp[372] + amp[374] + amp[376] - amp[382] + amp[377] - amp[383] +
        amp[405] - amp[407] - amp[408] + amp[410] + amp[412] - amp[418] +
        amp[413] - amp[419] + amp[504] + amp[505]);
    jamp[22] = +2. * (+amp[0] - amp[2] - amp[4] + amp[5] - amp[7] - amp[8] -
        amp[9] - amp[10] - amp[29] - amp[30] - amp[32] - amp[34] - amp[35] +
        amp[39] - amp[41] + amp[102] + amp[103] + amp[105] + amp[107] - amp[108]
        - amp[109] - amp[111] - amp[118] - amp[119] - amp[121] - amp[122] +
        amp[129] + amp[130] - amp[132] + amp[134] - amp[279] - amp[280] -
        amp[282] - amp[283] - amp[284] + amp[286] - amp[287] + amp[295] +
        amp[296] + amp[297] + amp[298] + amp[300] + amp[331] + amp[332] +
        amp[337] + amp[338] + amp[342] + amp[345] - amp[344] - amp[347] -
        amp[366] + amp[368] - amp[369] + amp[371] - amp[375] + amp[381] -
        amp[376] + amp[382] - amp[387] - amp[388] + amp[390] + amp[391] -
        amp[394] - amp[397] - amp[395] - amp[398] - amp[507] - amp[508]);
    jamp[23] = +2. * (-amp[0] + amp[2] + amp[4] - amp[5] - amp[18] + amp[20] -
        amp[21] + amp[22] - amp[28] + amp[30] - amp[31] - amp[33] + amp[35] -
        amp[42] + amp[44] + amp[57] + amp[58] + amp[60] + amp[62] - amp[63] -
        amp[64] - amp[66] - amp[73] - amp[74] - amp[76] - amp[77] + amp[84] +
        amp[85] - amp[87] + amp[89] + amp[279] + amp[280] + amp[282] + amp[283] +
        amp[289] + amp[290] + amp[292] + amp[294] - amp[296] - amp[297] -
        amp[298] - amp[300] + amp[313] + amp[314] + amp[319] + amp[320] +
        amp[324] + amp[327] - amp[326] - amp[329] + amp[366] - amp[368] +
        amp[372] - amp[374] + amp[375] - amp[381] - amp[377] + amp[383] +
        amp[387] - amp[389] - amp[390] + amp[392] + amp[394] - amp[400] +
        amp[395] - amp[401] + amp[507] + amp[508]);
    jamp[24] = +2. * (+amp[45] - amp[47] - amp[49] + amp[50] + amp[63] - amp[65]
        + amp[66] - amp[67] + amp[73] - amp[75] + amp[76] + amp[78] - amp[80] +
        amp[87] - amp[89] + amp[180] - amp[182] - amp[184] + amp[185] - amp[187]
        - amp[188] - amp[189] - amp[190] - amp[209] - amp[210] - amp[212] -
        amp[214] - amp[215] + amp[219] - amp[221] - amp[225] + amp[227] -
        amp[228] - amp[235] + amp[237] - amp[238] + amp[239] - amp[240] +
        amp[242] + amp[254] + amp[255] + amp[257] - amp[312] - amp[313] -
        amp[318] - amp[319] - amp[321] - amp[324] + amp[323] + amp[326] +
        amp[421] + amp[422] + amp[427] + amp[428] + amp[432] + amp[435] -
        amp[434] - amp[437] + amp[474] + amp[475] + amp[477] + amp[478] +
        amp[483] + amp[486] + amp[484] + amp[487] + amp[493] + amp[494]);
    jamp[25] = +2. * (-amp[45] - amp[46] - amp[48] - amp[50] + amp[68] - amp[70]
        + amp[71] - amp[72] - amp[73] + amp[75] - amp[76] + amp[81] - amp[83] -
        amp[87] - amp[88] + amp[135] - amp[137] - amp[139] + amp[140] - amp[142]
        - amp[143] - amp[144] - amp[145] - amp[164] - amp[165] - amp[167] -
        amp[169] - amp[170] + amp[174] - amp[176] - amp[230] + amp[232] -
        amp[233] + amp[235] - amp[237] + amp[238] - amp[239] + amp[240] +
        amp[241] - amp[254] - amp[255] - amp[257] - amp[315] - amp[316] +
        amp[318] + amp[319] - amp[322] - amp[325] - amp[323] - amp[326] +
        amp[403] + amp[404] + amp[409] + amp[410] + amp[414] + amp[417] -
        amp[416] - amp[419] - amp[474] + amp[476] - amp[477] + amp[479] -
        amp[483] + amp[489] - amp[484] + amp[490] - amp[493] - amp[494]);
    jamp[26] = +2. * (+amp[51] - amp[53] - amp[55] + amp[56] - amp[63] + amp[65]
        - amp[66] + amp[67] - amp[68] - amp[69] - amp[71] - amp[78] - amp[79] +
        amp[88] + amp[89] - amp[180] - amp[181] - amp[183] - amp[185] - amp[193]
        - amp[194] - amp[195] - amp[196] + amp[209] + amp[210] + amp[212] -
        amp[217] - amp[218] - amp[219] - amp[220] + amp[225] - amp[227] +
        amp[228] + amp[230] + amp[231] + amp[233] + amp[234] - amp[241] -
        amp[242] + amp[272] + amp[273] + amp[275] + amp[312] + amp[313] +
        amp[315] + amp[316] + amp[321] + amp[324] + amp[322] + amp[325] +
        amp[424] + amp[425] - amp[427] - amp[428] + amp[433] + amp[436] +
        amp[434] + amp[437] + amp[456] + amp[457] + amp[459] + amp[460] +
        amp[465] + amp[468] + amp[466] + amp[469] + amp[496] + amp[497]);
    jamp[27] = +2. * (-amp[51] - amp[52] - amp[54] - amp[56] + amp[68] + amp[69]
        + amp[71] - amp[73] - amp[74] - amp[76] - amp[77] + amp[84] - amp[86] -
        amp[87] - amp[88] + amp[90] - amp[92] - amp[94] + amp[95] - amp[97] -
        amp[98] - amp[99] - amp[100] - amp[119] - amp[120] - amp[122] - amp[124]
        - amp[125] + amp[129] - amp[131] - amp[230] - amp[231] - amp[233] -
        amp[234] + amp[235] + amp[236] + amp[238] + amp[240] + amp[241] -
        amp[272] - amp[273] - amp[275] - amp[315] - amp[316] + amp[318] +
        amp[319] - amp[322] - amp[325] - amp[323] - amp[326] + amp[385] +
        amp[386] + amp[391] + amp[392] + amp[396] + amp[399] - amp[398] -
        amp[401] - amp[456] + amp[458] - amp[459] + amp[461] - amp[465] +
        amp[471] - amp[466] + amp[472] - amp[496] - amp[497]);
    jamp[28] = +2. * (+amp[57] - amp[59] - amp[61] + amp[62] - amp[63] - amp[64]
        - amp[66] - amp[68] + amp[70] - amp[71] + amp[72] - amp[81] - amp[82] +
        amp[88] + amp[89] - amp[135] - amp[136] - amp[138] - amp[140] - amp[148]
        - amp[149] - amp[150] - amp[151] + amp[164] + amp[165] + amp[167] -
        amp[172] - amp[173] - amp[174] - amp[175] + amp[225] + amp[226] +
        amp[228] + amp[229] + amp[230] - amp[232] + amp[233] - amp[241] -
        amp[242] + amp[290] + amp[291] + amp[293] + amp[312] + amp[313] +
        amp[315] + amp[316] + amp[321] + amp[324] + amp[322] + amp[325] +
        amp[406] + amp[407] - amp[409] - amp[410] + amp[415] + amp[418] +
        amp[416] + amp[419] + amp[438] + amp[439] + amp[441] + amp[442] +
        amp[447] + amp[450] + amp[448] + amp[451] + amp[499] + amp[500]);
    jamp[29] = +2. * (-amp[57] - amp[58] - amp[60] - amp[62] + amp[63] + amp[64]
        + amp[66] + amp[73] + amp[74] + amp[76] + amp[77] - amp[84] - amp[85] +
        amp[87] - amp[89] - amp[90] - amp[91] - amp[93] - amp[95] - amp[103] -
        amp[104] - amp[105] - amp[106] + amp[119] + amp[120] + amp[122] -
        amp[127] - amp[128] - amp[129] - amp[130] - amp[225] - amp[226] -
        amp[228] - amp[229] - amp[235] - amp[236] - amp[238] - amp[240] +
        amp[242] - amp[290] - amp[291] - amp[293] - amp[312] - amp[313] -
        amp[318] - amp[319] - amp[321] - amp[324] + amp[323] + amp[326] +
        amp[388] + amp[389] - amp[391] - amp[392] + amp[397] + amp[400] +
        amp[398] + amp[401] - amp[438] + amp[440] - amp[441] + amp[443] -
        amp[447] + amp[453] - amp[448] + amp[454] - amp[499] - amp[500]);
    jamp[30] = +2. * (-amp[45] + amp[47] + amp[49] - amp[50] + amp[52] + amp[53]
        + amp[54] + amp[55] + amp[74] + amp[75] + amp[77] + amp[79] + amp[80] -
        amp[84] + amp[86] + amp[187] + amp[188] + amp[189] + amp[190] + amp[192]
        + amp[193] + amp[195] + amp[197] - amp[199] - amp[200] - amp[202] +
        amp[213] + amp[214] + amp[220] + amp[221] - amp[254] - amp[255] -
        amp[257] + amp[262] + amp[263] + amp[265] + amp[298] + amp[299] +
        amp[300] + amp[301] + amp[384] - amp[386] + amp[390] - amp[392] +
        amp[393] - amp[399] - amp[395] + amp[401] - amp[421] - amp[422] -
        amp[424] - amp[425] - amp[432] - amp[435] - amp[433] - amp[436] -
        amp[477] - amp[478] + amp[480] + amp[481] - amp[484] - amp[487] -
        amp[485] - amp[488] - amp[492] - amp[493] - amp[507] + amp[509]);
    jamp[31] = +2. * (+amp[45] + amp[46] + amp[48] + amp[50] + amp[58] + amp[59]
        + amp[60] + amp[61] - amp[74] - amp[75] - amp[77] + amp[82] + amp[83] +
        amp[84] + amp[85] + amp[142] + amp[143] + amp[144] + amp[145] + amp[147]
        + amp[148] + amp[150] + amp[152] - amp[154] - amp[155] - amp[157] +
        amp[168] + amp[169] + amp[175] + amp[176] + amp[254] + amp[255] +
        amp[257] + amp[280] + amp[281] + amp[283] - amp[298] - amp[299] -
        amp[300] - amp[301] + amp[387] - amp[389] - amp[390] + amp[392] +
        amp[394] - amp[400] + amp[395] - amp[401] - amp[403] - amp[404] -
        amp[406] - amp[407] - amp[414] - amp[417] - amp[415] - amp[418] +
        amp[477] - amp[479] - amp[480] + amp[482] + amp[484] - amp[490] +
        amp[485] - amp[491] + amp[492] + amp[493] + amp[507] + amp[508]);
    jamp[32] = +2. * (-amp[52] - amp[53] - amp[54] - amp[55] - amp[57] - amp[58]
        - amp[60] - amp[62] + amp[64] + amp[65] + amp[67] - amp[78] - amp[79] -
        amp[85] - amp[86] - amp[192] - amp[193] - amp[195] - amp[197] + amp[198]
        + amp[199] + amp[201] + amp[208] + amp[209] + amp[211] + amp[212] -
        amp[219] - amp[220] + amp[222] - amp[224] - amp[262] - amp[263] -
        amp[265] - amp[279] - amp[280] - amp[282] - amp[283] - amp[289] -
        amp[290] - amp[292] - amp[294] + amp[296] - amp[366] - amp[367] -
        amp[372] - amp[373] - amp[375] - amp[378] + amp[377] + amp[380] -
        amp[384] + amp[386] - amp[387] + amp[389] - amp[393] + amp[399] -
        amp[394] + amp[400] + amp[424] + amp[425] - amp[427] - amp[428] +
        amp[433] + amp[436] + amp[434] + amp[437] - amp[508] - amp[509]);
    jamp[33] = +2. * (-amp[0] + amp[2] + amp[4] - amp[5] - amp[18] + amp[20] -
        amp[21] + amp[22] - amp[28] + amp[30] - amp[31] - amp[33] + amp[35] -
        amp[42] + amp[44] + amp[57] + amp[58] + amp[60] + amp[62] - amp[63] -
        amp[64] - amp[66] - amp[73] - amp[74] - amp[76] - amp[77] + amp[84] +
        amp[85] - amp[87] + amp[89] + amp[279] + amp[280] + amp[282] + amp[283] +
        amp[289] + amp[290] + amp[292] + amp[294] - amp[296] - amp[297] -
        amp[298] - amp[300] + amp[313] + amp[314] + amp[319] + amp[320] +
        amp[324] + amp[327] - amp[326] - amp[329] + amp[366] - amp[368] +
        amp[372] - amp[374] + amp[375] - amp[381] - amp[377] + amp[383] +
        amp[387] - amp[389] - amp[390] + amp[392] + amp[394] - amp[400] +
        amp[395] - amp[401] + amp[507] + amp[508]);
    jamp[34] = +2. * (-amp[51] - amp[52] - amp[54] - amp[56] - amp[58] - amp[59]
        - amp[60] - amp[61] + amp[69] + amp[70] + amp[72] - amp[81] - amp[82] -
        amp[85] - amp[86] - amp[147] - amp[148] - amp[150] - amp[152] + amp[153]
        + amp[154] + amp[156] + amp[163] + amp[164] + amp[166] + amp[167] -
        amp[174] - amp[175] + amp[177] - amp[179] - amp[261] - amp[262] -
        amp[264] - amp[265] - amp[271] - amp[272] - amp[274] - amp[276] +
        amp[278] - amp[280] - amp[281] - amp[283] - amp[348] - amp[349] -
        amp[354] - amp[355] - amp[357] - amp[360] + amp[359] + amp[362] -
        amp[384] + amp[386] - amp[387] + amp[389] - amp[393] + amp[399] -
        amp[394] + amp[400] + amp[406] + amp[407] - amp[409] - amp[410] +
        amp[415] + amp[418] + amp[416] + amp[419] - amp[508] - amp[509]);
    jamp[35] = +2. * (+amp[0] + amp[1] + amp[3] + amp[5] - amp[23] + amp[25] -
        amp[26] + amp[27] + amp[28] - amp[30] + amp[31] - amp[36] + amp[38] +
        amp[42] + amp[43] + amp[51] + amp[52] + amp[54] + amp[56] - amp[68] -
        amp[69] - amp[71] + amp[73] + amp[74] + amp[76] + amp[77] - amp[84] +
        amp[86] + amp[87] + amp[88] + amp[261] + amp[262] + amp[264] + amp[265] +
        amp[271] + amp[272] + amp[274] + amp[276] - amp[278] + amp[297] +
        amp[298] + amp[300] + amp[316] + amp[317] - amp[319] - amp[320] +
        amp[325] + amp[328] + amp[326] + amp[329] + amp[348] - amp[350] +
        amp[354] - amp[356] + amp[357] - amp[363] - amp[359] + amp[365] +
        amp[384] - amp[386] + amp[390] - amp[392] + amp[393] - amp[399] -
        amp[395] + amp[401] - amp[507] + amp[509]);
    jamp[36] = +2. * (+amp[46] + amp[47] + amp[48] + amp[49] - amp[51] + amp[53]
        + amp[55] - amp[56] + amp[69] + amp[70] + amp[72] + amp[79] + amp[80] -
        amp[81] + amp[83] + amp[186] + amp[187] + amp[189] + amp[191] + amp[193]
        + amp[194] + amp[195] + amp[196] - amp[204] - amp[205] - amp[207] +
        amp[216] + amp[217] + amp[220] + amp[221] + amp[244] + amp[245] +
        amp[247] - amp[272] - amp[273] - amp[275] + amp[303] + amp[304] +
        amp[305] + amp[306] + amp[402] - amp[404] + amp[408] - amp[410] +
        amp[411] - amp[417] - amp[413] + amp[419] - amp[421] - amp[422] -
        amp[424] - amp[425] - amp[432] - amp[435] - amp[433] - amp[436] -
        amp[459] - amp[460] + amp[462] + amp[463] - amp[466] - amp[469] -
        amp[467] - amp[470] - amp[495] - amp[496] - amp[504] + amp[506]);
    jamp[37] = +2. * (+amp[51] + amp[52] + amp[54] + amp[56] + amp[58] + amp[59]
        + amp[60] + amp[61] - amp[69] - amp[70] - amp[72] + amp[81] + amp[82] +
        amp[85] + amp[86] + amp[97] + amp[98] + amp[99] + amp[100] + amp[102] +
        amp[103] + amp[105] + amp[107] - amp[109] - amp[110] - amp[112] +
        amp[123] + amp[124] + amp[130] + amp[131] + amp[272] + amp[273] +
        amp[275] + amp[285] + amp[286] + amp[288] - amp[303] - amp[304] -
        amp[305] - amp[306] - amp[385] - amp[386] - amp[388] - amp[389] -
        amp[396] - amp[399] - amp[397] - amp[400] + amp[405] - amp[407] -
        amp[408] + amp[410] + amp[412] - amp[418] + amp[413] - amp[419] +
        amp[459] - amp[461] - amp[462] + amp[464] + amp[466] - amp[472] +
        amp[467] - amp[473] + amp[495] + amp[496] + amp[504] + amp[505]);
    jamp[38] = +2. * (-amp[46] - amp[47] - amp[48] - amp[49] + amp[57] - amp[59]
        - amp[61] + amp[62] - amp[64] - amp[65] - amp[67] + amp[78] - amp[80] -
        amp[82] - amp[83] - amp[186] - amp[187] - amp[189] - amp[191] + amp[203]
        + amp[204] + amp[206] - amp[208] - amp[209] - amp[211] - amp[212] +
        amp[219] - amp[221] - amp[222] - amp[223] - amp[244] - amp[245] -
        amp[247] - amp[284] - amp[285] - amp[287] - amp[288] + amp[289] +
        amp[290] + amp[292] + amp[294] + amp[295] - amp[369] - amp[370] +
        amp[372] + amp[373] - amp[376] - amp[379] - amp[377] - amp[380] -
        amp[402] + amp[404] - amp[405] + amp[407] - amp[411] + amp[417] -
        amp[412] + amp[418] + amp[421] + amp[422] + amp[427] + amp[428] +
        amp[432] + amp[435] - amp[434] - amp[437] - amp[505] - amp[506]);
    jamp[39] = +2. * (-amp[6] + amp[8] + amp[10] - amp[11] + amp[18] - amp[20] +
        amp[21] - amp[22] + amp[23] + amp[24] + amp[26] + amp[33] + amp[34] -
        amp[43] - amp[44] - amp[57] + amp[59] + amp[61] - amp[62] + amp[63] +
        amp[64] + amp[66] + amp[68] - amp[70] + amp[71] - amp[72] + amp[81] +
        amp[82] - amp[88] - amp[89] + amp[284] + amp[285] + amp[287] + amp[288] -
        amp[289] - amp[290] - amp[292] - amp[294] - amp[295] - amp[302] -
        amp[303] - amp[305] - amp[313] - amp[314] - amp[316] - amp[317] -
        amp[324] - amp[327] - amp[325] - amp[328] + amp[369] - amp[371] -
        amp[372] + amp[374] + amp[376] - amp[382] + amp[377] - amp[383] +
        amp[405] - amp[407] - amp[408] + amp[410] + amp[412] - amp[418] +
        amp[413] - amp[419] + amp[504] + amp[505]);
    jamp[40] = +2. * (-amp[45] - amp[46] - amp[48] - amp[50] - amp[58] - amp[59]
        - amp[60] - amp[61] + amp[74] + amp[75] + amp[77] - amp[82] - amp[83] -
        amp[84] - amp[85] - amp[102] - amp[103] - amp[105] - amp[107] + amp[108]
        + amp[109] + amp[111] + amp[118] + amp[119] + amp[121] + amp[122] -
        amp[129] - amp[130] + amp[132] - amp[134] - amp[243] - amp[244] -
        amp[246] - amp[247] - amp[253] - amp[254] - amp[256] - amp[258] +
        amp[260] - amp[285] - amp[286] - amp[288] - amp[330] - amp[331] -
        amp[336] - amp[337] - amp[339] - amp[342] + amp[341] + amp[344] +
        amp[388] + amp[389] - amp[391] - amp[392] + amp[397] + amp[400] +
        amp[398] + amp[401] - amp[402] + amp[404] - amp[405] + amp[407] -
        amp[411] + amp[417] - amp[412] + amp[418] - amp[505] - amp[506]);
    jamp[41] = +2. * (+amp[6] + amp[7] + amp[9] + amp[11] - amp[23] - amp[24] -
        amp[26] + amp[28] + amp[29] + amp[31] + amp[32] - amp[39] + amp[41] +
        amp[42] + amp[43] + amp[45] + amp[46] + amp[48] + amp[50] - amp[68] +
        amp[70] - amp[71] + amp[72] + amp[73] - amp[75] + amp[76] - amp[81] +
        amp[83] + amp[87] + amp[88] + amp[243] + amp[244] + amp[246] + amp[247] +
        amp[253] + amp[254] + amp[256] + amp[258] - amp[260] + amp[302] +
        amp[303] + amp[305] + amp[316] + amp[317] - amp[319] - amp[320] +
        amp[325] + amp[328] + amp[326] + amp[329] + amp[330] - amp[332] +
        amp[336] - amp[338] + amp[339] - amp[345] - amp[341] + amp[347] +
        amp[402] - amp[404] + amp[408] - amp[410] + amp[411] - amp[417] -
        amp[413] + amp[419] - amp[504] + amp[506]);
    jamp[42] = +2. * (+amp[46] + amp[47] + amp[48] + amp[49] - amp[57] + amp[59]
        + amp[61] - amp[62] + amp[64] + amp[65] + amp[67] - amp[78] + amp[80] +
        amp[82] + amp[83] + amp[141] + amp[142] + amp[144] + amp[146] + amp[148]
        + amp[149] + amp[150] + amp[151] - amp[159] - amp[160] - amp[162] +
        amp[171] + amp[172] + amp[175] + amp[176] + amp[249] + amp[250] +
        amp[252] - amp[290] - amp[291] - amp[293] + amp[308] + amp[309] +
        amp[310] + amp[311] - amp[403] - amp[404] - amp[406] - amp[407] -
        amp[414] - amp[417] - amp[415] - amp[418] + amp[420] - amp[422] +
        amp[426] - amp[428] + amp[429] - amp[435] - amp[431] + amp[437] -
        amp[441] - amp[442] + amp[444] + amp[445] - amp[448] - amp[451] -
        amp[449] - amp[452] - amp[498] - amp[499] - amp[501] + amp[503]);
    jamp[43] = +2. * (+amp[52] + amp[53] + amp[54] + amp[55] + amp[57] + amp[58]
        + amp[60] + amp[62] - amp[64] - amp[65] - amp[67] + amp[78] + amp[79] +
        amp[85] + amp[86] + amp[96] + amp[97] + amp[99] + amp[101] + amp[103] +
        amp[104] + amp[105] + amp[106] - amp[114] - amp[115] - amp[117] +
        amp[126] + amp[127] + amp[130] + amp[131] + amp[267] + amp[268] +
        amp[270] + amp[290] + amp[291] + amp[293] - amp[308] - amp[309] -
        amp[310] - amp[311] - amp[385] - amp[386] - amp[388] - amp[389] -
        amp[396] - amp[399] - amp[397] - amp[400] + amp[423] - amp[425] -
        amp[426] + amp[428] + amp[430] - amp[436] + amp[431] - amp[437] +
        amp[441] - amp[443] - amp[444] + amp[446] + amp[448] - amp[454] +
        amp[449] - amp[455] + amp[498] + amp[499] + amp[501] + amp[502]);
    jamp[44] = +2. * (-amp[46] - amp[47] - amp[48] - amp[49] + amp[51] - amp[53]
        - amp[55] + amp[56] - amp[69] - amp[70] - amp[72] - amp[79] - amp[80] +
        amp[81] - amp[83] - amp[141] - amp[142] - amp[144] - amp[146] + amp[158]
        + amp[159] + amp[161] - amp[163] - amp[164] - amp[166] - amp[167] +
        amp[174] - amp[176] - amp[177] - amp[178] - amp[249] - amp[250] -
        amp[252] - amp[266] - amp[267] - amp[269] - amp[270] + amp[271] +
        amp[272] + amp[274] + amp[276] + amp[277] - amp[351] - amp[352] +
        amp[354] + amp[355] - amp[358] - amp[361] - amp[359] - amp[362] +
        amp[403] + amp[404] + amp[409] + amp[410] + amp[414] + amp[417] -
        amp[416] - amp[419] - amp[420] + amp[422] - amp[423] + amp[425] -
        amp[429] + amp[435] - amp[430] + amp[436] - amp[502] - amp[503]);
    jamp[45] = +2. * (-amp[12] + amp[14] + amp[16] - amp[17] + amp[18] + amp[19]
        + amp[21] + amp[23] - amp[25] + amp[26] - amp[27] + amp[36] + amp[37] -
        amp[43] - amp[44] - amp[51] + amp[53] + amp[55] - amp[56] + amp[63] -
        amp[65] + amp[66] - amp[67] + amp[68] + amp[69] + amp[71] + amp[78] +
        amp[79] - amp[88] - amp[89] + amp[266] + amp[267] + amp[269] + amp[270] -
        amp[271] - amp[272] - amp[274] - amp[276] - amp[277] - amp[307] -
        amp[308] - amp[310] - amp[313] - amp[314] - amp[316] - amp[317] -
        amp[324] - amp[327] - amp[325] - amp[328] + amp[351] - amp[353] -
        amp[354] + amp[356] + amp[358] - amp[364] + amp[359] - amp[365] +
        amp[423] - amp[425] - amp[426] + amp[428] + amp[430] - amp[436] +
        amp[431] - amp[437] + amp[501] + amp[502]);
    jamp[46] = +2. * (+amp[45] - amp[47] - amp[49] + amp[50] - amp[52] - amp[53]
        - amp[54] - amp[55] - amp[74] - amp[75] - amp[77] - amp[79] - amp[80] +
        amp[84] - amp[86] - amp[96] - amp[97] - amp[99] - amp[101] + amp[113] +
        amp[114] + amp[116] - amp[118] - amp[119] - amp[121] - amp[122] +
        amp[129] - amp[131] - amp[132] - amp[133] - amp[248] - amp[249] -
        amp[251] - amp[252] + amp[253] + amp[254] + amp[256] + amp[258] +
        amp[259] - amp[267] - amp[268] - amp[270] - amp[333] - amp[334] +
        amp[336] + amp[337] - amp[340] - amp[343] - amp[341] - amp[344] +
        amp[385] + amp[386] + amp[391] + amp[392] + amp[396] + amp[399] -
        amp[398] - amp[401] - amp[420] + amp[422] - amp[423] + amp[425] -
        amp[429] + amp[435] - amp[430] + amp[436] - amp[502] - amp[503]);
    jamp[47] = +2. * (+amp[12] + amp[13] + amp[15] + amp[17] - amp[18] - amp[19]
        - amp[21] - amp[28] - amp[29] - amp[31] - amp[32] + amp[39] + amp[40] -
        amp[42] + amp[44] - amp[45] + amp[47] + amp[49] - amp[50] - amp[63] +
        amp[65] - amp[66] + amp[67] - amp[73] + amp[75] - amp[76] - amp[78] +
        amp[80] - amp[87] + amp[89] + amp[248] + amp[249] + amp[251] + amp[252] -
        amp[253] - amp[254] - amp[256] - amp[258] - amp[259] + amp[307] +
        amp[308] + amp[310] + amp[313] + amp[314] + amp[319] + amp[320] +
        amp[324] + amp[327] - amp[326] - amp[329] + amp[333] - amp[335] -
        amp[336] + amp[338] + amp[340] - amp[346] + amp[341] - amp[347] +
        amp[420] - amp[422] + amp[426] - amp[428] + amp[429] - amp[435] -
        amp[431] + amp[437] - amp[501] + amp[503]);
    jamp[48] = +2. * (+amp[90] - amp[92] - amp[94] + amp[95] + amp[108] -
        amp[110] + amp[111] - amp[112] + amp[118] - amp[120] + amp[121] +
        amp[123] - amp[125] + amp[132] - amp[134] - amp[181] - amp[182] -
        amp[183] - amp[184] + amp[186] - amp[188] - amp[190] + amp[191] -
        amp[204] - amp[205] - amp[207] - amp[214] - amp[215] + amp[216] -
        amp[218] + amp[236] + amp[237] + amp[239] - amp[243] + amp[245] -
        amp[246] - amp[253] + amp[255] - amp[256] + amp[257] - amp[258] +
        amp[260] - amp[330] - amp[331] - amp[336] - amp[337] - amp[339] -
        amp[342] + amp[341] + amp[344] + amp[457] + amp[458] + amp[463] +
        amp[464] + amp[468] + amp[471] - amp[470] - amp[473] + amp[474] +
        amp[475] + amp[477] + amp[478] + amp[483] + amp[486] + amp[484] +
        amp[487] + amp[493] + amp[494]);
    jamp[49] = +2. * (-amp[90] - amp[91] - amp[93] - amp[95] + amp[113] -
        amp[115] + amp[116] - amp[117] - amp[118] + amp[120] - amp[121] +
        amp[126] - amp[128] - amp[132] - amp[133] - amp[136] - amp[137] -
        amp[138] - amp[139] + amp[141] - amp[143] - amp[145] + amp[146] -
        amp[159] - amp[160] - amp[162] - amp[169] - amp[170] + amp[171] -
        amp[173] - amp[236] - amp[237] - amp[239] - amp[248] + amp[250] -
        amp[251] + amp[253] - amp[255] + amp[256] - amp[257] + amp[258] +
        amp[259] - amp[333] - amp[334] + amp[336] + amp[337] - amp[340] -
        amp[343] - amp[341] - amp[344] + amp[439] + amp[440] + amp[445] +
        amp[446] + amp[450] + amp[453] - amp[452] - amp[455] - amp[474] +
        amp[476] - amp[477] + amp[479] - amp[483] + amp[489] - amp[484] +
        amp[490] - amp[493] - amp[494]);
    jamp[50] = +2. * (+amp[96] - amp[98] - amp[100] + amp[101] - amp[108] +
        amp[110] - amp[111] + amp[112] - amp[113] - amp[114] - amp[116] -
        amp[123] - amp[124] + amp[133] + amp[134] - amp[186] - amp[187] -
        amp[189] - amp[191] - amp[193] - amp[194] - amp[195] - amp[196] +
        amp[204] + amp[205] + amp[207] - amp[216] - amp[217] - amp[220] -
        amp[221] + amp[243] - amp[245] + amp[246] + amp[248] + amp[249] +
        amp[251] + amp[252] - amp[259] - amp[260] + amp[267] + amp[268] +
        amp[270] + amp[330] + amp[331] + amp[333] + amp[334] + amp[339] +
        amp[342] + amp[340] + amp[343] + amp[420] + amp[421] + amp[423] +
        amp[424] + amp[429] + amp[432] + amp[430] + amp[433] + amp[460] +
        amp[461] - amp[463] - amp[464] + amp[469] + amp[472] + amp[470] +
        amp[473] + amp[502] + amp[503]);
    jamp[51] = +2. * (+amp[45] - amp[47] - amp[49] + amp[50] - amp[52] - amp[53]
        - amp[54] - amp[55] - amp[74] - amp[75] - amp[77] - amp[79] - amp[80] +
        amp[84] - amp[86] - amp[96] - amp[97] - amp[99] - amp[101] + amp[113] +
        amp[114] + amp[116] - amp[118] - amp[119] - amp[121] - amp[122] +
        amp[129] - amp[131] - amp[132] - amp[133] - amp[248] - amp[249] -
        amp[251] - amp[252] + amp[253] + amp[254] + amp[256] + amp[258] +
        amp[259] - amp[267] - amp[268] - amp[270] - amp[333] - amp[334] +
        amp[336] + amp[337] - amp[340] - amp[343] - amp[341] - amp[344] +
        amp[385] + amp[386] + amp[391] + amp[392] + amp[396] + amp[399] -
        amp[398] - amp[401] - amp[420] + amp[422] - amp[423] + amp[425] -
        amp[429] + amp[435] - amp[430] + amp[436] - amp[502] - amp[503]);
    jamp[52] = +2. * (+amp[102] - amp[104] - amp[106] + amp[107] - amp[108] -
        amp[109] - amp[111] - amp[113] + amp[115] - amp[116] + amp[117] -
        amp[126] - amp[127] + amp[133] + amp[134] - amp[141] - amp[142] -
        amp[144] - amp[146] - amp[148] - amp[149] - amp[150] - amp[151] +
        amp[159] + amp[160] + amp[162] - amp[171] - amp[172] - amp[175] -
        amp[176] + amp[243] + amp[244] + amp[246] + amp[247] + amp[248] -
        amp[250] + amp[251] - amp[259] - amp[260] + amp[285] + amp[286] +
        amp[288] + amp[330] + amp[331] + amp[333] + amp[334] + amp[339] +
        amp[342] + amp[340] + amp[343] + amp[402] + amp[403] + amp[405] +
        amp[406] + amp[411] + amp[414] + amp[412] + amp[415] + amp[442] +
        amp[443] - amp[445] - amp[446] + amp[451] + amp[454] + amp[452] +
        amp[455] + amp[505] + amp[506]);
    jamp[53] = +2. * (-amp[45] - amp[46] - amp[48] - amp[50] - amp[58] - amp[59]
        - amp[60] - amp[61] + amp[74] + amp[75] + amp[77] - amp[82] - amp[83] -
        amp[84] - amp[85] - amp[102] - amp[103] - amp[105] - amp[107] + amp[108]
        + amp[109] + amp[111] + amp[118] + amp[119] + amp[121] + amp[122] -
        amp[129] - amp[130] + amp[132] - amp[134] - amp[243] - amp[244] -
        amp[246] - amp[247] - amp[253] - amp[254] - amp[256] - amp[258] +
        amp[260] - amp[285] - amp[286] - amp[288] - amp[330] - amp[331] -
        amp[336] - amp[337] - amp[339] - amp[342] + amp[341] + amp[344] +
        amp[388] + amp[389] - amp[391] - amp[392] + amp[397] + amp[400] +
        amp[398] + amp[401] - amp[402] + amp[404] - amp[405] + amp[407] -
        amp[411] + amp[417] - amp[412] + amp[418] - amp[505] - amp[506]);
    jamp[54] = +2. * (-amp[90] + amp[92] + amp[94] - amp[95] + amp[97] + amp[98]
        + amp[99] + amp[100] + amp[119] + amp[120] + amp[122] + amp[124] +
        amp[125] - amp[129] + amp[131] + amp[181] + amp[182] + amp[183] +
        amp[184] - amp[192] + amp[194] + amp[196] - amp[197] + amp[199] +
        amp[200] + amp[202] - amp[213] + amp[215] + amp[217] + amp[218] -
        amp[236] - amp[237] - amp[239] - amp[262] - amp[263] - amp[265] -
        amp[298] - amp[299] - amp[300] - amp[301] - amp[384] - amp[385] -
        amp[390] - amp[391] - amp[393] - amp[396] + amp[395] + amp[398] -
        amp[457] - amp[458] - amp[460] - amp[461] - amp[468] - amp[471] -
        amp[469] - amp[472] - amp[474] - amp[475] - amp[480] - amp[481] -
        amp[483] - amp[486] + amp[485] + amp[488] + amp[492] - amp[494] +
        amp[507] - amp[509]);
    jamp[55] = +2. * (+amp[90] + amp[91] + amp[93] + amp[95] + amp[103] +
        amp[104] + amp[105] + amp[106] - amp[119] - amp[120] - amp[122] +
        amp[127] + amp[128] + amp[129] + amp[130] + amp[136] + amp[137] +
        amp[138] + amp[139] - amp[147] + amp[149] + amp[151] - amp[152] +
        amp[154] + amp[155] + amp[157] - amp[168] + amp[170] + amp[172] +
        amp[173] + amp[236] + amp[237] + amp[239] - amp[280] - amp[281] -
        amp[283] + amp[298] + amp[299] + amp[300] + amp[301] - amp[387] -
        amp[388] + amp[390] + amp[391] - amp[394] - amp[397] - amp[395] -
        amp[398] - amp[439] - amp[440] - amp[442] - amp[443] - amp[450] -
        amp[453] - amp[451] - amp[454] + amp[474] - amp[476] + amp[480] -
        amp[482] + amp[483] - amp[489] - amp[485] + amp[491] - amp[492] +
        amp[494] - amp[507] - amp[508]);
    jamp[56] = +2. * (-amp[97] - amp[98] - amp[99] - amp[100] - amp[102] -
        amp[103] - amp[105] - amp[107] + amp[109] + amp[110] + amp[112] -
        amp[123] - amp[124] - amp[130] - amp[131] + amp[192] - amp[194] -
        amp[196] + amp[197] - amp[198] - amp[199] - amp[201] - amp[203] +
        amp[205] - amp[206] + amp[207] - amp[216] - amp[217] + amp[223] +
        amp[224] + amp[262] + amp[263] + amp[265] + amp[279] + amp[280] +
        amp[282] + amp[283] + amp[284] - amp[286] + amp[287] - amp[295] -
        amp[296] + amp[366] + amp[367] + amp[369] + amp[370] + amp[375] +
        amp[378] + amp[376] + amp[379] + amp[384] + amp[385] + amp[387] +
        amp[388] + amp[393] + amp[396] + amp[394] + amp[397] + amp[460] +
        amp[461] - amp[463] - amp[464] + amp[469] + amp[472] + amp[470] +
        amp[473] + amp[508] + amp[509]);
    jamp[57] = +2. * (+amp[0] - amp[2] - amp[4] + amp[5] - amp[7] - amp[8] -
        amp[9] - amp[10] - amp[29] - amp[30] - amp[32] - amp[34] - amp[35] +
        amp[39] - amp[41] + amp[102] + amp[103] + amp[105] + amp[107] - amp[108]
        - amp[109] - amp[111] - amp[118] - amp[119] - amp[121] - amp[122] +
        amp[129] + amp[130] - amp[132] + amp[134] - amp[279] - amp[280] -
        amp[282] - amp[283] - amp[284] + amp[286] - amp[287] + amp[295] +
        amp[296] + amp[297] + amp[298] + amp[300] + amp[331] + amp[332] +
        amp[337] + amp[338] + amp[342] + amp[345] - amp[344] - amp[347] -
        amp[366] + amp[368] - amp[369] + amp[371] - amp[375] + amp[381] -
        amp[376] + amp[382] - amp[387] - amp[388] + amp[390] + amp[391] -
        amp[394] - amp[397] - amp[395] - amp[398] - amp[507] - amp[508]);
    jamp[58] = +2. * (-amp[96] - amp[97] - amp[99] - amp[101] - amp[103] -
        amp[104] - amp[105] - amp[106] + amp[114] + amp[115] + amp[117] -
        amp[126] - amp[127] - amp[130] - amp[131] + amp[147] - amp[149] -
        amp[151] + amp[152] - amp[153] - amp[154] - amp[156] - amp[158] +
        amp[160] - amp[161] + amp[162] - amp[171] - amp[172] + amp[178] +
        amp[179] + amp[261] + amp[262] + amp[264] + amp[265] + amp[266] -
        amp[268] + amp[269] - amp[277] - amp[278] + amp[280] + amp[281] +
        amp[283] + amp[348] + amp[349] + amp[351] + amp[352] + amp[357] +
        amp[360] + amp[358] + amp[361] + amp[384] + amp[385] + amp[387] +
        amp[388] + amp[393] + amp[396] + amp[394] + amp[397] + amp[442] +
        amp[443] - amp[445] - amp[446] + amp[451] + amp[454] + amp[452] +
        amp[455] + amp[508] + amp[509]);
    jamp[59] = +2. * (-amp[0] - amp[1] - amp[3] - amp[5] - amp[13] - amp[14] -
        amp[15] - amp[16] + amp[29] + amp[30] + amp[32] - amp[37] - amp[38] -
        amp[39] - amp[40] + amp[96] + amp[97] + amp[99] + amp[101] - amp[113] -
        amp[114] - amp[116] + amp[118] + amp[119] + amp[121] + amp[122] -
        amp[129] + amp[131] + amp[132] + amp[133] - amp[261] - amp[262] -
        amp[264] - amp[265] - amp[266] + amp[268] - amp[269] + amp[277] +
        amp[278] - amp[297] - amp[298] - amp[300] + amp[334] + amp[335] -
        amp[337] - amp[338] + amp[343] + amp[346] + amp[344] + amp[347] -
        amp[348] + amp[350] - amp[351] + amp[353] - amp[357] + amp[363] -
        amp[358] + amp[364] - amp[384] - amp[385] - amp[390] - amp[391] -
        amp[393] - amp[396] + amp[395] + amp[398] + amp[507] - amp[509]);
    jamp[60] = +2. * (+amp[91] + amp[92] + amp[93] + amp[94] - amp[96] + amp[98]
        + amp[100] - amp[101] + amp[114] + amp[115] + amp[117] + amp[124] +
        amp[125] - amp[126] + amp[128] + amp[180] + amp[181] + amp[183] +
        amp[185] + amp[193] + amp[194] + amp[195] + amp[196] - amp[209] -
        amp[210] - amp[212] + amp[217] + amp[218] + amp[219] + amp[220] +
        amp[226] + amp[227] + amp[229] - amp[267] - amp[268] - amp[270] +
        amp[308] + amp[309] + amp[310] + amp[311] - amp[423] - amp[424] +
        amp[426] + amp[427] - amp[430] - amp[433] - amp[431] - amp[434] +
        amp[438] - amp[440] + amp[444] - amp[446] + amp[447] - amp[453] -
        amp[449] + amp[455] - amp[457] - amp[458] - amp[460] - amp[461] -
        amp[468] - amp[471] - amp[469] - amp[472] - amp[498] + amp[500] -
        amp[501] - amp[502]);
    jamp[61] = +2. * (+amp[52] + amp[53] + amp[54] + amp[55] + amp[57] + amp[58]
        + amp[60] + amp[62] - amp[64] - amp[65] - amp[67] + amp[78] + amp[79] +
        amp[85] + amp[86] + amp[96] + amp[97] + amp[99] + amp[101] + amp[103] +
        amp[104] + amp[105] + amp[106] - amp[114] - amp[115] - amp[117] +
        amp[126] + amp[127] + amp[130] + amp[131] + amp[267] + amp[268] +
        amp[270] + amp[290] + amp[291] + amp[293] - amp[308] - amp[309] -
        amp[310] - amp[311] - amp[385] - amp[386] - amp[388] - amp[389] -
        amp[396] - amp[399] - amp[397] - amp[400] + amp[423] - amp[425] -
        amp[426] + amp[428] + amp[430] - amp[436] + amp[431] - amp[437] +
        amp[441] - amp[443] - amp[444] + amp[446] + amp[448] - amp[454] +
        amp[449] - amp[455] + amp[498] + amp[499] + amp[501] + amp[502]);
    jamp[62] = +2. * (-amp[91] - amp[92] - amp[93] - amp[94] + amp[102] -
        amp[104] - amp[106] + amp[107] - amp[109] - amp[110] - amp[112] +
        amp[123] - amp[125] - amp[127] - amp[128] - amp[180] - amp[181] -
        amp[183] - amp[185] + amp[203] - amp[205] + amp[206] - amp[207] -
        amp[208] + amp[210] - amp[211] + amp[216] - amp[218] - amp[222] -
        amp[223] - amp[226] - amp[227] - amp[229] - amp[284] + amp[286] -
        amp[287] + amp[289] - amp[291] + amp[292] - amp[293] + amp[294] +
        amp[295] - amp[369] - amp[370] + amp[372] + amp[373] - amp[376] -
        amp[379] - amp[377] - amp[380] - amp[438] + amp[440] - amp[441] +
        amp[443] - amp[447] + amp[453] - amp[448] + amp[454] + amp[457] +
        amp[458] + amp[463] + amp[464] + amp[468] + amp[471] - amp[470] -
        amp[473] - amp[499] - amp[500]);
    jamp[63] = +2. * (+amp[7] + amp[8] + amp[9] + amp[10] + amp[12] + amp[13] +
        amp[15] + amp[17] - amp[19] - amp[20] - amp[22] + amp[33] + amp[34] +
        amp[40] + amp[41] - amp[102] + amp[104] + amp[106] - amp[107] + amp[108]
        + amp[109] + amp[111] + amp[113] - amp[115] + amp[116] - amp[117] +
        amp[126] + amp[127] - amp[133] - amp[134] + amp[284] - amp[286] +
        amp[287] - amp[289] + amp[291] - amp[292] + amp[293] - amp[294] -
        amp[295] + amp[307] - amp[309] - amp[311] - amp[331] - amp[332] -
        amp[334] - amp[335] - amp[342] - amp[345] - amp[343] - amp[346] +
        amp[369] - amp[371] - amp[372] + amp[374] + amp[376] - amp[382] +
        amp[377] - amp[383] + amp[441] - amp[443] - amp[444] + amp[446] +
        amp[448] - amp[454] + amp[449] - amp[455] + amp[498] + amp[499]);
    jamp[64] = +2. * (-amp[57] - amp[58] - amp[60] - amp[62] + amp[63] + amp[64]
        + amp[66] + amp[73] + amp[74] + amp[76] + amp[77] - amp[84] - amp[85] +
        amp[87] - amp[89] - amp[90] - amp[91] - amp[93] - amp[95] - amp[103] -
        amp[104] - amp[105] - amp[106] + amp[119] + amp[120] + amp[122] -
        amp[127] - amp[128] - amp[129] - amp[130] - amp[225] - amp[226] -
        amp[228] - amp[229] - amp[235] - amp[236] - amp[238] - amp[240] +
        amp[242] - amp[290] - amp[291] - amp[293] - amp[312] - amp[313] -
        amp[318] - amp[319] - amp[321] - amp[324] + amp[323] + amp[326] +
        amp[388] + amp[389] - amp[391] - amp[392] + amp[397] + amp[400] +
        amp[398] + amp[401] - amp[438] + amp[440] - amp[441] + amp[443] -
        amp[447] + amp[453] - amp[448] + amp[454] - amp[499] - amp[500]);
    jamp[65] = +2. * (-amp[12] - amp[13] - amp[15] - amp[17] + amp[18] + amp[19]
        + amp[21] + amp[28] + amp[29] + amp[31] + amp[32] - amp[39] - amp[40] +
        amp[42] - amp[44] + amp[90] + amp[91] + amp[93] + amp[95] - amp[113] +
        amp[115] - amp[116] + amp[117] + amp[118] - amp[120] + amp[121] -
        amp[126] + amp[128] + amp[132] + amp[133] + amp[225] + amp[226] +
        amp[228] + amp[229] + amp[235] + amp[236] + amp[238] + amp[240] -
        amp[242] - amp[307] + amp[309] + amp[311] + amp[312] - amp[314] +
        amp[318] - amp[320] + amp[321] - amp[327] - amp[323] + amp[329] +
        amp[334] + amp[335] - amp[337] - amp[338] + amp[343] + amp[346] +
        amp[344] + amp[347] + amp[438] - amp[440] + amp[444] - amp[446] +
        amp[447] - amp[453] - amp[449] + amp[455] - amp[498] + amp[500]);
    jamp[66] = +2. * (+amp[91] + amp[92] + amp[93] + amp[94] - amp[102] +
        amp[104] + amp[106] - amp[107] + amp[109] + amp[110] + amp[112] -
        amp[123] + amp[125] + amp[127] + amp[128] + amp[135] + amp[136] +
        amp[138] + amp[140] + amp[148] + amp[149] + amp[150] + amp[151] -
        amp[164] - amp[165] - amp[167] + amp[172] + amp[173] + amp[174] +
        amp[175] + amp[231] + amp[232] + amp[234] - amp[285] - amp[286] -
        amp[288] + amp[303] + amp[304] + amp[305] + amp[306] - amp[405] -
        amp[406] + amp[408] + amp[409] - amp[412] - amp[415] - amp[413] -
        amp[416] - amp[439] - amp[440] - amp[442] - amp[443] - amp[450] -
        amp[453] - amp[451] - amp[454] + amp[456] - amp[458] + amp[462] -
        amp[464] + amp[465] - amp[471] - amp[467] + amp[473] - amp[495] +
        amp[497] - amp[504] - amp[505]);
    jamp[67] = +2. * (+amp[51] + amp[52] + amp[54] + amp[56] + amp[58] + amp[59]
        + amp[60] + amp[61] - amp[69] - amp[70] - amp[72] + amp[81] + amp[82] +
        amp[85] + amp[86] + amp[97] + amp[98] + amp[99] + amp[100] + amp[102] +
        amp[103] + amp[105] + amp[107] - amp[109] - amp[110] - amp[112] +
        amp[123] + amp[124] + amp[130] + amp[131] + amp[272] + amp[273] +
        amp[275] + amp[285] + amp[286] + amp[288] - amp[303] - amp[304] -
        amp[305] - amp[306] - amp[385] - amp[386] - amp[388] - amp[389] -
        amp[396] - amp[399] - amp[397] - amp[400] + amp[405] - amp[407] -
        amp[408] + amp[410] + amp[412] - amp[418] + amp[413] - amp[419] +
        amp[459] - amp[461] - amp[462] + amp[464] + amp[466] - amp[472] +
        amp[467] - amp[473] + amp[495] + amp[496] + amp[504] + amp[505]);
    jamp[68] = +2. * (-amp[91] - amp[92] - amp[93] - amp[94] + amp[96] - amp[98]
        - amp[100] + amp[101] - amp[114] - amp[115] - amp[117] - amp[124] -
        amp[125] + amp[126] - amp[128] - amp[135] - amp[136] - amp[138] -
        amp[140] + amp[158] - amp[160] + amp[161] - amp[162] - amp[163] +
        amp[165] - amp[166] + amp[171] - amp[173] - amp[177] - amp[178] -
        amp[231] - amp[232] - amp[234] - amp[266] + amp[268] - amp[269] +
        amp[271] - amp[273] + amp[274] - amp[275] + amp[276] + amp[277] -
        amp[351] - amp[352] + amp[354] + amp[355] - amp[358] - amp[361] -
        amp[359] - amp[362] + amp[439] + amp[440] + amp[445] + amp[446] +
        amp[450] + amp[453] - amp[452] - amp[455] - amp[456] + amp[458] -
        amp[459] + amp[461] - amp[465] + amp[471] - amp[466] + amp[472] -
        amp[496] - amp[497]);
    jamp[69] = +2. * (+amp[6] + amp[7] + amp[9] + amp[11] + amp[13] + amp[14] +
        amp[15] + amp[16] - amp[24] - amp[25] - amp[27] + amp[36] + amp[37] +
        amp[40] + amp[41] - amp[96] + amp[98] + amp[100] - amp[101] + amp[108] -
        amp[110] + amp[111] - amp[112] + amp[113] + amp[114] + amp[116] +
        amp[123] + amp[124] - amp[133] - amp[134] + amp[266] - amp[268] +
        amp[269] - amp[271] + amp[273] - amp[274] + amp[275] - amp[276] -
        amp[277] + amp[302] - amp[304] - amp[306] - amp[331] - amp[332] -
        amp[334] - amp[335] - amp[342] - amp[345] - amp[343] - amp[346] +
        amp[351] - amp[353] - amp[354] + amp[356] + amp[358] - amp[364] +
        amp[359] - amp[365] + amp[459] - amp[461] - amp[462] + amp[464] +
        amp[466] - amp[472] + amp[467] - amp[473] + amp[495] + amp[496]);
    jamp[70] = +2. * (-amp[51] - amp[52] - amp[54] - amp[56] + amp[68] + amp[69]
        + amp[71] - amp[73] - amp[74] - amp[76] - amp[77] + amp[84] - amp[86] -
        amp[87] - amp[88] + amp[90] - amp[92] - amp[94] + amp[95] - amp[97] -
        amp[98] - amp[99] - amp[100] - amp[119] - amp[120] - amp[122] - amp[124]
        - amp[125] + amp[129] - amp[131] - amp[230] - amp[231] - amp[233] -
        amp[234] + amp[235] + amp[236] + amp[238] + amp[240] + amp[241] -
        amp[272] - amp[273] - amp[275] - amp[315] - amp[316] + amp[318] +
        amp[319] - amp[322] - amp[325] - amp[323] - amp[326] + amp[385] +
        amp[386] + amp[391] + amp[392] + amp[396] + amp[399] - amp[398] -
        amp[401] - amp[456] + amp[458] - amp[459] + amp[461] - amp[465] +
        amp[471] - amp[466] + amp[472] - amp[496] - amp[497]);
    jamp[71] = +2. * (-amp[6] - amp[7] - amp[9] - amp[11] + amp[23] + amp[24] +
        amp[26] - amp[28] - amp[29] - amp[31] - amp[32] + amp[39] - amp[41] -
        amp[42] - amp[43] - amp[90] + amp[92] + amp[94] - amp[95] - amp[108] +
        amp[110] - amp[111] + amp[112] - amp[118] + amp[120] - amp[121] -
        amp[123] + amp[125] - amp[132] + amp[134] + amp[230] + amp[231] +
        amp[233] + amp[234] - amp[235] - amp[236] - amp[238] - amp[240] -
        amp[241] - amp[302] + amp[304] + amp[306] + amp[315] - amp[317] -
        amp[318] + amp[320] + amp[322] - amp[328] + amp[323] - amp[329] +
        amp[331] + amp[332] + amp[337] + amp[338] + amp[342] + amp[345] -
        amp[344] - amp[347] + amp[456] - amp[458] + amp[462] - amp[464] +
        amp[465] - amp[471] - amp[467] + amp[473] - amp[495] + amp[497]);
    jamp[72] = +2. * (+amp[135] - amp[137] - amp[139] + amp[140] + amp[153] -
        amp[155] + amp[156] - amp[157] + amp[163] - amp[165] + amp[166] +
        amp[168] - amp[170] + amp[177] - amp[179] - amp[181] - amp[182] -
        amp[183] - amp[184] + amp[192] - amp[194] - amp[196] + amp[197] -
        amp[199] - amp[200] - amp[202] + amp[213] - amp[215] - amp[217] -
        amp[218] + amp[231] + amp[232] + amp[234] - amp[261] + amp[263] -
        amp[264] - amp[271] + amp[273] - amp[274] + amp[275] - amp[276] +
        amp[278] - amp[348] - amp[349] - amp[354] - amp[355] - amp[357] -
        amp[360] + amp[359] + amp[362] + amp[456] + amp[457] + amp[459] +
        amp[460] + amp[465] + amp[468] + amp[466] + amp[469] + amp[475] +
        amp[476] + amp[481] + amp[482] + amp[486] + amp[489] - amp[488] -
        amp[491] + amp[496] + amp[497]);
    jamp[73] = +2. * (-amp[91] - amp[92] - amp[93] - amp[94] + amp[96] - amp[98]
        - amp[100] + amp[101] - amp[114] - amp[115] - amp[117] - amp[124] -
        amp[125] + amp[126] - amp[128] - amp[135] - amp[136] - amp[138] -
        amp[140] + amp[158] - amp[160] + amp[161] - amp[162] - amp[163] +
        amp[165] - amp[166] + amp[171] - amp[173] - amp[177] - amp[178] -
        amp[231] - amp[232] - amp[234] - amp[266] + amp[268] - amp[269] +
        amp[271] - amp[273] + amp[274] - amp[275] + amp[276] + amp[277] -
        amp[351] - amp[352] + amp[354] + amp[355] - amp[358] - amp[361] -
        amp[359] - amp[362] + amp[439] + amp[440] + amp[445] + amp[446] +
        amp[450] + amp[453] - amp[452] - amp[455] - amp[456] + amp[458] -
        amp[459] + amp[461] - amp[465] + amp[471] - amp[466] + amp[472] -
        amp[496] - amp[497]);
    jamp[74] = +2. * (+amp[141] - amp[143] - amp[145] + amp[146] - amp[153] +
        amp[155] - amp[156] + amp[157] - amp[158] - amp[159] - amp[161] -
        amp[168] - amp[169] + amp[178] + amp[179] - amp[187] - amp[188] -
        amp[189] - amp[190] - amp[192] - amp[193] - amp[195] - amp[197] +
        amp[199] + amp[200] + amp[202] - amp[213] - amp[214] - amp[220] -
        amp[221] + amp[249] + amp[250] + amp[252] + amp[261] - amp[263] +
        amp[264] + amp[266] + amp[267] + amp[269] + amp[270] - amp[277] -
        amp[278] + amp[348] + amp[349] + amp[351] + amp[352] + amp[357] +
        amp[360] + amp[358] + amp[361] + amp[420] + amp[421] + amp[423] +
        amp[424] + amp[429] + amp[432] + amp[430] + amp[433] + amp[478] +
        amp[479] - amp[481] - amp[482] + amp[487] + amp[490] + amp[488] +
        amp[491] + amp[502] + amp[503]);
    jamp[75] = +2. * (-amp[46] - amp[47] - amp[48] - amp[49] + amp[51] - amp[53]
        - amp[55] + amp[56] - amp[69] - amp[70] - amp[72] - amp[79] - amp[80] +
        amp[81] - amp[83] - amp[141] - amp[142] - amp[144] - amp[146] + amp[158]
        + amp[159] + amp[161] - amp[163] - amp[164] - amp[166] - amp[167] +
        amp[174] - amp[176] - amp[177] - amp[178] - amp[249] - amp[250] -
        amp[252] - amp[266] - amp[267] - amp[269] - amp[270] + amp[271] +
        amp[272] + amp[274] + amp[276] + amp[277] - amp[351] - amp[352] +
        amp[354] + amp[355] - amp[358] - amp[361] - amp[359] - amp[362] +
        amp[403] + amp[404] + amp[409] + amp[410] + amp[414] + amp[417] -
        amp[416] - amp[419] - amp[420] + amp[422] - amp[423] + amp[425] -
        amp[429] + amp[435] - amp[430] + amp[436] - amp[502] - amp[503]);
    jamp[76] = +2. * (-amp[96] - amp[97] - amp[99] - amp[101] - amp[103] -
        amp[104] - amp[105] - amp[106] + amp[114] + amp[115] + amp[117] -
        amp[126] - amp[127] - amp[130] - amp[131] + amp[147] - amp[149] -
        amp[151] + amp[152] - amp[153] - amp[154] - amp[156] - amp[158] +
        amp[160] - amp[161] + amp[162] - amp[171] - amp[172] + amp[178] +
        amp[179] + amp[261] + amp[262] + amp[264] + amp[265] + amp[266] -
        amp[268] + amp[269] - amp[277] - amp[278] + amp[280] + amp[281] +
        amp[283] + amp[348] + amp[349] + amp[351] + amp[352] + amp[357] +
        amp[360] + amp[358] + amp[361] + amp[384] + amp[385] + amp[387] +
        amp[388] + amp[393] + amp[396] + amp[394] + amp[397] + amp[442] +
        amp[443] - amp[445] - amp[446] + amp[451] + amp[454] + amp[452] +
        amp[455] + amp[508] + amp[509]);
    jamp[77] = +2. * (-amp[51] - amp[52] - amp[54] - amp[56] - amp[58] - amp[59]
        - amp[60] - amp[61] + amp[69] + amp[70] + amp[72] - amp[81] - amp[82] -
        amp[85] - amp[86] - amp[147] - amp[148] - amp[150] - amp[152] + amp[153]
        + amp[154] + amp[156] + amp[163] + amp[164] + amp[166] + amp[167] -
        amp[174] - amp[175] + amp[177] - amp[179] - amp[261] - amp[262] -
        amp[264] - amp[265] - amp[271] - amp[272] - amp[274] - amp[276] +
        amp[278] - amp[280] - amp[281] - amp[283] - amp[348] - amp[349] -
        amp[354] - amp[355] - amp[357] - amp[360] + amp[359] + amp[362] -
        amp[384] + amp[386] - amp[387] + amp[389] - amp[393] + amp[399] -
        amp[394] + amp[400] + amp[406] + amp[407] - amp[409] - amp[410] +
        amp[415] + amp[418] + amp[416] + amp[419] - amp[508] - amp[509]);
    jamp[78] = +2. * (-amp[135] + amp[137] + amp[139] - amp[140] + amp[142] +
        amp[143] + amp[144] + amp[145] + amp[164] + amp[165] + amp[167] +
        amp[169] + amp[170] - amp[174] + amp[176] + amp[181] + amp[182] +
        amp[183] + amp[184] - amp[186] + amp[188] + amp[190] - amp[191] +
        amp[204] + amp[205] + amp[207] + amp[214] + amp[215] - amp[216] +
        amp[218] - amp[231] - amp[232] - amp[234] - amp[244] - amp[245] -
        amp[247] - amp[303] - amp[304] - amp[305] - amp[306] - amp[402] -
        amp[403] - amp[408] - amp[409] - amp[411] - amp[414] + amp[413] +
        amp[416] - amp[456] - amp[457] - amp[462] - amp[463] - amp[465] -
        amp[468] + amp[467] + amp[470] - amp[475] - amp[476] - amp[478] -
        amp[479] - amp[486] - amp[489] - amp[487] - amp[490] + amp[495] -
        amp[497] + amp[504] - amp[506]);
    jamp[79] = +2. * (+amp[91] + amp[92] + amp[93] + amp[94] - amp[102] +
        amp[104] + amp[106] - amp[107] + amp[109] + amp[110] + amp[112] -
        amp[123] + amp[125] + amp[127] + amp[128] + amp[135] + amp[136] +
        amp[138] + amp[140] + amp[148] + amp[149] + amp[150] + amp[151] -
        amp[164] - amp[165] - amp[167] + amp[172] + amp[173] + amp[174] +
        amp[175] + amp[231] + amp[232] + amp[234] - amp[285] - amp[286] -
        amp[288] + amp[303] + amp[304] + amp[305] + amp[306] - amp[405] -
        amp[406] + amp[408] + amp[409] - amp[412] - amp[415] - amp[413] -
        amp[416] - amp[439] - amp[440] - amp[442] - amp[443] - amp[450] -
        amp[453] - amp[451] - amp[454] + amp[456] - amp[458] + amp[462] -
        amp[464] + amp[465] - amp[471] - amp[467] + amp[473] - amp[495] +
        amp[497] - amp[504] - amp[505]);
    jamp[80] = +2. * (-amp[142] - amp[143] - amp[144] - amp[145] - amp[147] -
        amp[148] - amp[150] - amp[152] + amp[154] + amp[155] + amp[157] -
        amp[168] - amp[169] - amp[175] - amp[176] + amp[186] - amp[188] -
        amp[190] + amp[191] - amp[198] + amp[200] - amp[201] + amp[202] -
        amp[203] - amp[204] - amp[206] - amp[213] - amp[214] + amp[223] +
        amp[224] + amp[244] + amp[245] + amp[247] + amp[279] - amp[281] +
        amp[282] + amp[284] + amp[285] + amp[287] + amp[288] - amp[295] -
        amp[296] + amp[366] + amp[367] + amp[369] + amp[370] + amp[375] +
        amp[378] + amp[376] + amp[379] + amp[402] + amp[403] + amp[405] +
        amp[406] + amp[411] + amp[414] + amp[412] + amp[415] + amp[478] +
        amp[479] - amp[481] - amp[482] + amp[487] + amp[490] + amp[488] +
        amp[491] + amp[505] + amp[506]);
    jamp[81] = +2. * (-amp[1] - amp[2] - amp[3] - amp[4] + amp[6] - amp[8] -
        amp[10] + amp[11] - amp[24] - amp[25] - amp[27] - amp[34] - amp[35] +
        amp[36] - amp[38] + amp[147] + amp[148] + amp[150] + amp[152] - amp[153]
        - amp[154] - amp[156] - amp[163] - amp[164] - amp[166] - amp[167] +
        amp[174] + amp[175] - amp[177] + amp[179] - amp[279] + amp[281] -
        amp[282] - amp[284] - amp[285] - amp[287] - amp[288] + amp[295] +
        amp[296] + amp[302] + amp[303] + amp[305] + amp[349] + amp[350] +
        amp[355] + amp[356] + amp[360] + amp[363] - amp[362] - amp[365] -
        amp[366] + amp[368] - amp[369] + amp[371] - amp[375] + amp[381] -
        amp[376] + amp[382] - amp[405] - amp[406] + amp[408] + amp[409] -
        amp[412] - amp[415] - amp[413] - amp[416] - amp[504] - amp[505]);
    jamp[82] = +2. * (+amp[102] - amp[104] - amp[106] + amp[107] - amp[108] -
        amp[109] - amp[111] - amp[113] + amp[115] - amp[116] + amp[117] -
        amp[126] - amp[127] + amp[133] + amp[134] - amp[141] - amp[142] -
        amp[144] - amp[146] - amp[148] - amp[149] - amp[150] - amp[151] +
        amp[159] + amp[160] + amp[162] - amp[171] - amp[172] - amp[175] -
        amp[176] + amp[243] + amp[244] + amp[246] + amp[247] + amp[248] -
        amp[250] + amp[251] - amp[259] - amp[260] + amp[285] + amp[286] +
        amp[288] + amp[330] + amp[331] + amp[333] + amp[334] + amp[339] +
        amp[342] + amp[340] + amp[343] + amp[402] + amp[403] + amp[405] +
        amp[406] + amp[411] + amp[414] + amp[412] + amp[415] + amp[442] +
        amp[443] - amp[445] - amp[446] + amp[451] + amp[454] + amp[452] +
        amp[455] + amp[505] + amp[506]);
    jamp[83] = +2. * (-amp[6] - amp[7] - amp[9] - amp[11] - amp[13] - amp[14] -
        amp[15] - amp[16] + amp[24] + amp[25] + amp[27] - amp[36] - amp[37] -
        amp[40] - amp[41] + amp[141] + amp[142] + amp[144] + amp[146] - amp[158]
        - amp[159] - amp[161] + amp[163] + amp[164] + amp[166] + amp[167] -
        amp[174] + amp[176] + amp[177] + amp[178] - amp[243] - amp[244] -
        amp[246] - amp[247] - amp[248] + amp[250] - amp[251] + amp[259] +
        amp[260] - amp[302] - amp[303] - amp[305] - amp[330] + amp[332] -
        amp[333] + amp[335] - amp[339] + amp[345] - amp[340] + amp[346] +
        amp[352] + amp[353] - amp[355] - amp[356] + amp[361] + amp[364] +
        amp[362] + amp[365] - amp[402] - amp[403] - amp[408] - amp[409] -
        amp[411] - amp[414] + amp[413] + amp[416] + amp[504] - amp[506]);
    jamp[84] = +2. * (+amp[136] + amp[137] + amp[138] + amp[139] - amp[141] +
        amp[143] + amp[145] - amp[146] + amp[159] + amp[160] + amp[162] +
        amp[169] + amp[170] - amp[171] + amp[173] - amp[180] + amp[182] +
        amp[184] - amp[185] + amp[187] + amp[188] + amp[189] + amp[190] +
        amp[209] + amp[210] + amp[212] + amp[214] + amp[215] - amp[219] +
        amp[221] - amp[226] - amp[227] - amp[229] - amp[249] - amp[250] -
        amp[252] - amp[308] - amp[309] - amp[310] - amp[311] - amp[420] -
        amp[421] - amp[426] - amp[427] - amp[429] - amp[432] + amp[431] +
        amp[434] - amp[438] - amp[439] - amp[444] - amp[445] - amp[447] -
        amp[450] + amp[449] + amp[452] - amp[475] - amp[476] - amp[478] -
        amp[479] - amp[486] - amp[489] - amp[487] - amp[490] + amp[498] -
        amp[500] + amp[501] - amp[503]);
    jamp[85] = +2. * (+amp[46] + amp[47] + amp[48] + amp[49] - amp[57] + amp[59]
        + amp[61] - amp[62] + amp[64] + amp[65] + amp[67] - amp[78] + amp[80] +
        amp[82] + amp[83] + amp[141] + amp[142] + amp[144] + amp[146] + amp[148]
        + amp[149] + amp[150] + amp[151] - amp[159] - amp[160] - amp[162] +
        amp[171] + amp[172] + amp[175] + amp[176] + amp[249] + amp[250] +
        amp[252] - amp[290] - amp[291] - amp[293] + amp[308] + amp[309] +
        amp[310] + amp[311] - amp[403] - amp[404] - amp[406] - amp[407] -
        amp[414] - amp[417] - amp[415] - amp[418] + amp[420] - amp[422] +
        amp[426] - amp[428] + amp[429] - amp[435] - amp[431] + amp[437] -
        amp[441] - amp[442] + amp[444] + amp[445] - amp[448] - amp[451] -
        amp[449] - amp[452] - amp[498] - amp[499] - amp[501] + amp[503]);
    jamp[86] = +2. * (-amp[136] - amp[137] - amp[138] - amp[139] + amp[147] -
        amp[149] - amp[151] + amp[152] - amp[154] - amp[155] - amp[157] +
        amp[168] - amp[170] - amp[172] - amp[173] + amp[180] - amp[182] -
        amp[184] + amp[185] + amp[198] - amp[200] + amp[201] - amp[202] +
        amp[208] - amp[210] + amp[211] + amp[213] - amp[215] + amp[222] -
        amp[224] + amp[226] + amp[227] + amp[229] - amp[279] + amp[281] -
        amp[282] - amp[289] + amp[291] - amp[292] + amp[293] - amp[294] +
        amp[296] - amp[366] - amp[367] - amp[372] - amp[373] - amp[375] -
        amp[378] + amp[377] + amp[380] + amp[438] + amp[439] + amp[441] +
        amp[442] + amp[447] + amp[450] + amp[448] + amp[451] + amp[475] +
        amp[476] + amp[481] + amp[482] + amp[486] + amp[489] - amp[488] -
        amp[491] + amp[499] + amp[500]);
    jamp[87] = +2. * (+amp[1] + amp[2] + amp[3] + amp[4] - amp[12] + amp[14] +
        amp[16] - amp[17] + amp[19] + amp[20] + amp[22] - amp[33] + amp[35] +
        amp[37] + amp[38] - amp[147] + amp[149] + amp[151] - amp[152] + amp[153]
        + amp[154] + amp[156] + amp[158] - amp[160] + amp[161] - amp[162] +
        amp[171] + amp[172] - amp[178] - amp[179] + amp[279] - amp[281] +
        amp[282] + amp[289] - amp[291] + amp[292] - amp[293] + amp[294] -
        amp[296] - amp[307] + amp[309] + amp[311] - amp[349] - amp[350] -
        amp[352] - amp[353] - amp[360] - amp[363] - amp[361] - amp[364] +
        amp[366] - amp[368] + amp[372] - amp[374] + amp[375] - amp[381] -
        amp[377] + amp[383] - amp[441] - amp[442] + amp[444] + amp[445] -
        amp[448] - amp[451] - amp[449] - amp[452] - amp[498] - amp[499]);
    jamp[88] = +2. * (+amp[57] - amp[59] - amp[61] + amp[62] - amp[63] - amp[64]
        - amp[66] - amp[68] + amp[70] - amp[71] + amp[72] - amp[81] - amp[82] +
        amp[88] + amp[89] - amp[135] - amp[136] - amp[138] - amp[140] - amp[148]
        - amp[149] - amp[150] - amp[151] + amp[164] + amp[165] + amp[167] -
        amp[172] - amp[173] - amp[174] - amp[175] + amp[225] + amp[226] +
        amp[228] + amp[229] + amp[230] - amp[232] + amp[233] - amp[241] -
        amp[242] + amp[290] + amp[291] + amp[293] + amp[312] + amp[313] +
        amp[315] + amp[316] + amp[321] + amp[324] + amp[322] + amp[325] +
        amp[406] + amp[407] - amp[409] - amp[410] + amp[415] + amp[418] +
        amp[416] + amp[419] + amp[438] + amp[439] + amp[441] + amp[442] +
        amp[447] + amp[450] + amp[448] + amp[451] + amp[499] + amp[500]);
    jamp[89] = +2. * (+amp[12] - amp[14] - amp[16] + amp[17] - amp[18] - amp[19]
        - amp[21] - amp[23] + amp[25] - amp[26] + amp[27] - amp[36] - amp[37] +
        amp[43] + amp[44] + amp[135] + amp[136] + amp[138] + amp[140] - amp[158]
        + amp[160] - amp[161] + amp[162] + amp[163] - amp[165] + amp[166] -
        amp[171] + amp[173] + amp[177] + amp[178] - amp[225] - amp[226] -
        amp[228] - amp[229] - amp[230] + amp[232] - amp[233] + amp[241] +
        amp[242] + amp[307] - amp[309] - amp[311] - amp[312] + amp[314] -
        amp[315] + amp[317] - amp[321] + amp[327] - amp[322] + amp[328] +
        amp[352] + amp[353] - amp[355] - amp[356] + amp[361] + amp[364] +
        amp[362] + amp[365] - amp[438] - amp[439] - amp[444] - amp[445] -
        amp[447] - amp[450] + amp[449] + amp[452] + amp[498] - amp[500]);
    jamp[90] = +2. * (+amp[90] + amp[91] + amp[93] + amp[95] + amp[103] +
        amp[104] + amp[105] + amp[106] - amp[119] - amp[120] - amp[122] +
        amp[127] + amp[128] + amp[129] + amp[130] + amp[136] + amp[137] +
        amp[138] + amp[139] - amp[147] + amp[149] + amp[151] - amp[152] +
        amp[154] + amp[155] + amp[157] - amp[168] + amp[170] + amp[172] +
        amp[173] + amp[236] + amp[237] + amp[239] - amp[280] - amp[281] -
        amp[283] + amp[298] + amp[299] + amp[300] + amp[301] - amp[387] -
        amp[388] + amp[390] + amp[391] - amp[394] - amp[397] - amp[395] -
        amp[398] - amp[439] - amp[440] - amp[442] - amp[443] - amp[450] -
        amp[453] - amp[451] - amp[454] + amp[474] - amp[476] + amp[480] -
        amp[482] + amp[483] - amp[489] - amp[485] + amp[491] - amp[492] +
        amp[494] - amp[507] - amp[508]);
    jamp[91] = +2. * (+amp[45] + amp[46] + amp[48] + amp[50] + amp[58] + amp[59]
        + amp[60] + amp[61] - amp[74] - amp[75] - amp[77] + amp[82] + amp[83] +
        amp[84] + amp[85] + amp[142] + amp[143] + amp[144] + amp[145] + amp[147]
        + amp[148] + amp[150] + amp[152] - amp[154] - amp[155] - amp[157] +
        amp[168] + amp[169] + amp[175] + amp[176] + amp[254] + amp[255] +
        amp[257] + amp[280] + amp[281] + amp[283] - amp[298] - amp[299] -
        amp[300] - amp[301] + amp[387] - amp[389] - amp[390] + amp[392] +
        amp[394] - amp[400] + amp[395] - amp[401] - amp[403] - amp[404] -
        amp[406] - amp[407] - amp[414] - amp[417] - amp[415] - amp[418] +
        amp[477] - amp[479] - amp[480] + amp[482] + amp[484] - amp[490] +
        amp[485] - amp[491] + amp[492] + amp[493] + amp[507] + amp[508]);
    jamp[92] = +2. * (-amp[90] - amp[91] - amp[93] - amp[95] + amp[113] -
        amp[115] + amp[116] - amp[117] - amp[118] + amp[120] - amp[121] +
        amp[126] - amp[128] - amp[132] - amp[133] - amp[136] - amp[137] -
        amp[138] - amp[139] + amp[141] - amp[143] - amp[145] + amp[146] -
        amp[159] - amp[160] - amp[162] - amp[169] - amp[170] + amp[171] -
        amp[173] - amp[236] - amp[237] - amp[239] - amp[248] + amp[250] -
        amp[251] + amp[253] - amp[255] + amp[256] - amp[257] + amp[258] +
        amp[259] - amp[333] - amp[334] + amp[336] + amp[337] - amp[340] -
        amp[343] - amp[341] - amp[344] + amp[439] + amp[440] + amp[445] +
        amp[446] + amp[450] + amp[453] - amp[452] - amp[455] - amp[474] +
        amp[476] - amp[477] + amp[479] - amp[483] + amp[489] - amp[484] +
        amp[490] - amp[493] - amp[494]);
    jamp[93] = +2. * (+amp[0] + amp[1] + amp[3] + amp[5] + amp[13] + amp[14] +
        amp[15] + amp[16] - amp[29] - amp[30] - amp[32] + amp[37] + amp[38] +
        amp[39] + amp[40] - amp[141] + amp[143] + amp[145] - amp[146] + amp[153]
        - amp[155] + amp[156] - amp[157] + amp[158] + amp[159] + amp[161] +
        amp[168] + amp[169] - amp[178] - amp[179] + amp[248] - amp[250] +
        amp[251] - amp[253] + amp[255] - amp[256] + amp[257] - amp[258] -
        amp[259] + amp[297] - amp[299] - amp[301] + amp[333] - amp[335] -
        amp[336] + amp[338] + amp[340] - amp[346] + amp[341] - amp[347] -
        amp[349] - amp[350] - amp[352] - amp[353] - amp[360] - amp[363] -
        amp[361] - amp[364] + amp[477] - amp[479] - amp[480] + amp[482] +
        amp[484] - amp[490] + amp[485] - amp[491] + amp[492] + amp[493]);
    jamp[94] = +2. * (-amp[45] - amp[46] - amp[48] - amp[50] + amp[68] - amp[70]
        + amp[71] - amp[72] - amp[73] + amp[75] - amp[76] + amp[81] - amp[83] -
        amp[87] - amp[88] + amp[135] - amp[137] - amp[139] + amp[140] - amp[142]
        - amp[143] - amp[144] - amp[145] - amp[164] - amp[165] - amp[167] -
        amp[169] - amp[170] + amp[174] - amp[176] - amp[230] + amp[232] -
        amp[233] + amp[235] - amp[237] + amp[238] - amp[239] + amp[240] +
        amp[241] - amp[254] - amp[255] - amp[257] - amp[315] - amp[316] +
        amp[318] + amp[319] - amp[322] - amp[325] - amp[323] - amp[326] +
        amp[403] + amp[404] + amp[409] + amp[410] + amp[414] + amp[417] -
        amp[416] - amp[419] - amp[474] + amp[476] - amp[477] + amp[479] -
        amp[483] + amp[489] - amp[484] + amp[490] - amp[493] - amp[494]);
    jamp[95] = +2. * (-amp[0] - amp[1] - amp[3] - amp[5] + amp[23] - amp[25] +
        amp[26] - amp[27] - amp[28] + amp[30] - amp[31] + amp[36] - amp[38] -
        amp[42] - amp[43] - amp[135] + amp[137] + amp[139] - amp[140] - amp[153]
        + amp[155] - amp[156] + amp[157] - amp[163] + amp[165] - amp[166] -
        amp[168] + amp[170] - amp[177] + amp[179] + amp[230] - amp[232] +
        amp[233] - amp[235] + amp[237] - amp[238] + amp[239] - amp[240] -
        amp[241] - amp[297] + amp[299] + amp[301] + amp[315] - amp[317] -
        amp[318] + amp[320] + amp[322] - amp[328] + amp[323] - amp[329] +
        amp[349] + amp[350] + amp[355] + amp[356] + amp[360] + amp[363] -
        amp[362] - amp[365] + amp[474] - amp[476] + amp[480] - amp[482] +
        amp[483] - amp[489] - amp[485] + amp[491] - amp[492] + amp[494]);
    jamp[96] = +2. * (-amp[136] - amp[137] - amp[138] - amp[139] + amp[147] -
        amp[149] - amp[151] + amp[152] - amp[154] - amp[155] - amp[157] +
        amp[168] - amp[170] - amp[172] - amp[173] + amp[180] - amp[182] -
        amp[184] + amp[185] + amp[198] - amp[200] + amp[201] - amp[202] +
        amp[208] - amp[210] + amp[211] + amp[213] - amp[215] + amp[222] -
        amp[224] + amp[226] + amp[227] + amp[229] - amp[279] + amp[281] -
        amp[282] - amp[289] + amp[291] - amp[292] + amp[293] - amp[294] +
        amp[296] - amp[366] - amp[367] - amp[372] - amp[373] - amp[375] -
        amp[378] + amp[377] + amp[380] + amp[438] + amp[439] + amp[441] +
        amp[442] + amp[447] + amp[450] + amp[448] + amp[451] + amp[475] +
        amp[476] + amp[481] + amp[482] + amp[486] + amp[489] - amp[488] -
        amp[491] + amp[499] + amp[500]);
    jamp[97] = +2. * (-amp[91] - amp[92] - amp[93] - amp[94] + amp[102] -
        amp[104] - amp[106] + amp[107] - amp[109] - amp[110] - amp[112] +
        amp[123] - amp[125] - amp[127] - amp[128] - amp[180] - amp[181] -
        amp[183] - amp[185] + amp[203] - amp[205] + amp[206] - amp[207] -
        amp[208] + amp[210] - amp[211] + amp[216] - amp[218] - amp[222] -
        amp[223] - amp[226] - amp[227] - amp[229] - amp[284] + amp[286] -
        amp[287] + amp[289] - amp[291] + amp[292] - amp[293] + amp[294] +
        amp[295] - amp[369] - amp[370] + amp[372] + amp[373] - amp[376] -
        amp[379] - amp[377] - amp[380] - amp[438] + amp[440] - amp[441] +
        amp[443] - amp[447] + amp[453] - amp[448] + amp[454] + amp[457] +
        amp[458] + amp[463] + amp[464] + amp[468] + amp[471] - amp[470] -
        amp[473] - amp[499] - amp[500]);
    jamp[98] = +2. * (-amp[142] - amp[143] - amp[144] - amp[145] - amp[147] -
        amp[148] - amp[150] - amp[152] + amp[154] + amp[155] + amp[157] -
        amp[168] - amp[169] - amp[175] - amp[176] + amp[186] - amp[188] -
        amp[190] + amp[191] - amp[198] + amp[200] - amp[201] + amp[202] -
        amp[203] - amp[204] - amp[206] - amp[213] - amp[214] + amp[223] +
        amp[224] + amp[244] + amp[245] + amp[247] + amp[279] - amp[281] +
        amp[282] + amp[284] + amp[285] + amp[287] + amp[288] - amp[295] -
        amp[296] + amp[366] + amp[367] + amp[369] + amp[370] + amp[375] +
        amp[378] + amp[376] + amp[379] + amp[402] + amp[403] + amp[405] +
        amp[406] + amp[411] + amp[414] + amp[412] + amp[415] + amp[478] +
        amp[479] - amp[481] - amp[482] + amp[487] + amp[490] + amp[488] +
        amp[491] + amp[505] + amp[506]);
    jamp[99] = +2. * (-amp[46] - amp[47] - amp[48] - amp[49] + amp[57] - amp[59]
        - amp[61] + amp[62] - amp[64] - amp[65] - amp[67] + amp[78] - amp[80] -
        amp[82] - amp[83] - amp[186] - amp[187] - amp[189] - amp[191] + amp[203]
        + amp[204] + amp[206] - amp[208] - amp[209] - amp[211] - amp[212] +
        amp[219] - amp[221] - amp[222] - amp[223] - amp[244] - amp[245] -
        amp[247] - amp[284] - amp[285] - amp[287] - amp[288] + amp[289] +
        amp[290] + amp[292] + amp[294] + amp[295] - amp[369] - amp[370] +
        amp[372] + amp[373] - amp[376] - amp[379] - amp[377] - amp[380] -
        amp[402] + amp[404] - amp[405] + amp[407] - amp[411] + amp[417] -
        amp[412] + amp[418] + amp[421] + amp[422] + amp[427] + amp[428] +
        amp[432] + amp[435] - amp[434] - amp[437] - amp[505] - amp[506]);
    jamp[100] = +2. * (-amp[97] - amp[98] - amp[99] - amp[100] - amp[102] -
        amp[103] - amp[105] - amp[107] + amp[109] + amp[110] + amp[112] -
        amp[123] - amp[124] - amp[130] - amp[131] + amp[192] - amp[194] -
        amp[196] + amp[197] - amp[198] - amp[199] - amp[201] - amp[203] +
        amp[205] - amp[206] + amp[207] - amp[216] - amp[217] + amp[223] +
        amp[224] + amp[262] + amp[263] + amp[265] + amp[279] + amp[280] +
        amp[282] + amp[283] + amp[284] - amp[286] + amp[287] - amp[295] -
        amp[296] + amp[366] + amp[367] + amp[369] + amp[370] + amp[375] +
        amp[378] + amp[376] + amp[379] + amp[384] + amp[385] + amp[387] +
        amp[388] + amp[393] + amp[396] + amp[394] + amp[397] + amp[460] +
        amp[461] - amp[463] - amp[464] + amp[469] + amp[472] + amp[470] +
        amp[473] + amp[508] + amp[509]);
    jamp[101] = +2. * (-amp[52] - amp[53] - amp[54] - amp[55] - amp[57] - amp[58]
        - amp[60] - amp[62] + amp[64] + amp[65] + amp[67] - amp[78] - amp[79] -
        amp[85] - amp[86] - amp[192] - amp[193] - amp[195] - amp[197] + amp[198]
        + amp[199] + amp[201] + amp[208] + amp[209] + amp[211] + amp[212] -
        amp[219] - amp[220] + amp[222] - amp[224] - amp[262] - amp[263] -
        amp[265] - amp[279] - amp[280] - amp[282] - amp[283] - amp[289] -
        amp[290] - amp[292] - amp[294] + amp[296] - amp[366] - amp[367] -
        amp[372] - amp[373] - amp[375] - amp[378] + amp[377] + amp[380] -
        amp[384] + amp[386] - amp[387] + amp[389] - amp[393] + amp[399] -
        amp[394] + amp[400] + amp[424] + amp[425] - amp[427] - amp[428] +
        amp[433] + amp[436] + amp[434] + amp[437] - amp[508] - amp[509]);
    jamp[102] = +2. * (+amp[136] + amp[137] + amp[138] + amp[139] - amp[141] +
        amp[143] + amp[145] - amp[146] + amp[159] + amp[160] + amp[162] +
        amp[169] + amp[170] - amp[171] + amp[173] - amp[180] + amp[182] +
        amp[184] - amp[185] + amp[187] + amp[188] + amp[189] + amp[190] +
        amp[209] + amp[210] + amp[212] + amp[214] + amp[215] - amp[219] +
        amp[221] - amp[226] - amp[227] - amp[229] - amp[249] - amp[250] -
        amp[252] - amp[308] - amp[309] - amp[310] - amp[311] - amp[420] -
        amp[421] - amp[426] - amp[427] - amp[429] - amp[432] + amp[431] +
        amp[434] - amp[438] - amp[439] - amp[444] - amp[445] - amp[447] -
        amp[450] + amp[449] + amp[452] - amp[475] - amp[476] - amp[478] -
        amp[479] - amp[486] - amp[489] - amp[487] - amp[490] + amp[498] -
        amp[500] + amp[501] - amp[503]);
    jamp[103] = +2. * (+amp[91] + amp[92] + amp[93] + amp[94] - amp[96] + amp[98]
        + amp[100] - amp[101] + amp[114] + amp[115] + amp[117] + amp[124] +
        amp[125] - amp[126] + amp[128] + amp[180] + amp[181] + amp[183] +
        amp[185] + amp[193] + amp[194] + amp[195] + amp[196] - amp[209] -
        amp[210] - amp[212] + amp[217] + amp[218] + amp[219] + amp[220] +
        amp[226] + amp[227] + amp[229] - amp[267] - amp[268] - amp[270] +
        amp[308] + amp[309] + amp[310] + amp[311] - amp[423] - amp[424] +
        amp[426] + amp[427] - amp[430] - amp[433] - amp[431] - amp[434] +
        amp[438] - amp[440] + amp[444] - amp[446] + amp[447] - amp[453] -
        amp[449] + amp[455] - amp[457] - amp[458] - amp[460] - amp[461] -
        amp[468] - amp[471] - amp[469] - amp[472] - amp[498] + amp[500] -
        amp[501] - amp[502]);
    jamp[104] = +2. * (+amp[141] - amp[143] - amp[145] + amp[146] - amp[153] +
        amp[155] - amp[156] + amp[157] - amp[158] - amp[159] - amp[161] -
        amp[168] - amp[169] + amp[178] + amp[179] - amp[187] - amp[188] -
        amp[189] - amp[190] - amp[192] - amp[193] - amp[195] - amp[197] +
        amp[199] + amp[200] + amp[202] - amp[213] - amp[214] - amp[220] -
        amp[221] + amp[249] + amp[250] + amp[252] + amp[261] - amp[263] +
        amp[264] + amp[266] + amp[267] + amp[269] + amp[270] - amp[277] -
        amp[278] + amp[348] + amp[349] + amp[351] + amp[352] + amp[357] +
        amp[360] + amp[358] + amp[361] + amp[420] + amp[421] + amp[423] +
        amp[424] + amp[429] + amp[432] + amp[430] + amp[433] + amp[478] +
        amp[479] - amp[481] - amp[482] + amp[487] + amp[490] + amp[488] +
        amp[491] + amp[502] + amp[503]);
    jamp[105] = +2. * (-amp[1] - amp[2] - amp[3] - amp[4] + amp[12] - amp[14] -
        amp[16] + amp[17] - amp[19] - amp[20] - amp[22] + amp[33] - amp[35] -
        amp[37] - amp[38] + amp[192] + amp[193] + amp[195] + amp[197] - amp[198]
        - amp[199] - amp[201] - amp[208] - amp[209] - amp[211] - amp[212] +
        amp[219] + amp[220] - amp[222] + amp[224] - amp[261] + amp[263] -
        amp[264] - amp[266] - amp[267] - amp[269] - amp[270] + amp[277] +
        amp[278] + amp[307] + amp[308] + amp[310] - amp[348] + amp[350] -
        amp[351] + amp[353] - amp[357] + amp[363] - amp[358] + amp[364] +
        amp[367] + amp[368] + amp[373] + amp[374] + amp[378] + amp[381] -
        amp[380] - amp[383] - amp[423] - amp[424] + amp[426] + amp[427] -
        amp[430] - amp[433] - amp[431] - amp[434] - amp[501] - amp[502]);
    jamp[106] = +2. * (+amp[96] - amp[98] - amp[100] + amp[101] - amp[108] +
        amp[110] - amp[111] + amp[112] - amp[113] - amp[114] - amp[116] -
        amp[123] - amp[124] + amp[133] + amp[134] - amp[186] - amp[187] -
        amp[189] - amp[191] - amp[193] - amp[194] - amp[195] - amp[196] +
        amp[204] + amp[205] + amp[207] - amp[216] - amp[217] - amp[220] -
        amp[221] + amp[243] - amp[245] + amp[246] + amp[248] + amp[249] +
        amp[251] + amp[252] - amp[259] - amp[260] + amp[267] + amp[268] +
        amp[270] + amp[330] + amp[331] + amp[333] + amp[334] + amp[339] +
        amp[342] + amp[340] + amp[343] + amp[420] + amp[421] + amp[423] +
        amp[424] + amp[429] + amp[432] + amp[430] + amp[433] + amp[460] +
        amp[461] - amp[463] - amp[464] + amp[469] + amp[472] + amp[470] +
        amp[473] + amp[502] + amp[503]);
    jamp[107] = +2. * (-amp[7] - amp[8] - amp[9] - amp[10] - amp[12] - amp[13] -
        amp[15] - amp[17] + amp[19] + amp[20] + amp[22] - amp[33] - amp[34] -
        amp[40] - amp[41] + amp[186] + amp[187] + amp[189] + amp[191] - amp[203]
        - amp[204] - amp[206] + amp[208] + amp[209] + amp[211] + amp[212] -
        amp[219] + amp[221] + amp[222] + amp[223] - amp[243] + amp[245] -
        amp[246] - amp[248] - amp[249] - amp[251] - amp[252] + amp[259] +
        amp[260] - amp[307] - amp[308] - amp[310] - amp[330] + amp[332] -
        amp[333] + amp[335] - amp[339] + amp[345] - amp[340] + amp[346] +
        amp[370] + amp[371] - amp[373] - amp[374] + amp[379] + amp[382] +
        amp[380] + amp[383] - amp[420] - amp[421] - amp[426] - amp[427] -
        amp[429] - amp[432] + amp[431] + amp[434] + amp[501] - amp[503]);
    jamp[108] = +2. * (-amp[135] + amp[137] + amp[139] - amp[140] + amp[142] +
        amp[143] + amp[144] + amp[145] + amp[164] + amp[165] + amp[167] +
        amp[169] + amp[170] - amp[174] + amp[176] + amp[181] + amp[182] +
        amp[183] + amp[184] - amp[186] + amp[188] + amp[190] - amp[191] +
        amp[204] + amp[205] + amp[207] + amp[214] + amp[215] - amp[216] +
        amp[218] - amp[231] - amp[232] - amp[234] - amp[244] - amp[245] -
        amp[247] - amp[303] - amp[304] - amp[305] - amp[306] - amp[402] -
        amp[403] - amp[408] - amp[409] - amp[411] - amp[414] + amp[413] +
        amp[416] - amp[456] - amp[457] - amp[462] - amp[463] - amp[465] -
        amp[468] + amp[467] + amp[470] - amp[475] - amp[476] - amp[478] -
        amp[479] - amp[486] - amp[489] - amp[487] - amp[490] + amp[495] -
        amp[497] + amp[504] - amp[506]);
    jamp[109] = +2. * (+amp[46] + amp[47] + amp[48] + amp[49] - amp[51] + amp[53]
        + amp[55] - amp[56] + amp[69] + amp[70] + amp[72] + amp[79] + amp[80] -
        amp[81] + amp[83] + amp[186] + amp[187] + amp[189] + amp[191] + amp[193]
        + amp[194] + amp[195] + amp[196] - amp[204] - amp[205] - amp[207] +
        amp[216] + amp[217] + amp[220] + amp[221] + amp[244] + amp[245] +
        amp[247] - amp[272] - amp[273] - amp[275] + amp[303] + amp[304] +
        amp[305] + amp[306] + amp[402] - amp[404] + amp[408] - amp[410] +
        amp[411] - amp[417] - amp[413] + amp[419] - amp[421] - amp[422] -
        amp[424] - amp[425] - amp[432] - amp[435] - amp[433] - amp[436] -
        amp[459] - amp[460] + amp[462] + amp[463] - amp[466] - amp[469] -
        amp[467] - amp[470] - amp[495] - amp[496] - amp[504] + amp[506]);
    jamp[110] = +2. * (+amp[135] - amp[137] - amp[139] + amp[140] + amp[153] -
        amp[155] + amp[156] - amp[157] + amp[163] - amp[165] + amp[166] +
        amp[168] - amp[170] + amp[177] - amp[179] - amp[181] - amp[182] -
        amp[183] - amp[184] + amp[192] - amp[194] - amp[196] + amp[197] -
        amp[199] - amp[200] - amp[202] + amp[213] - amp[215] - amp[217] -
        amp[218] + amp[231] + amp[232] + amp[234] - amp[261] + amp[263] -
        amp[264] - amp[271] + amp[273] - amp[274] + amp[275] - amp[276] +
        amp[278] - amp[348] - amp[349] - amp[354] - amp[355] - amp[357] -
        amp[360] + amp[359] + amp[362] + amp[456] + amp[457] + amp[459] +
        amp[460] + amp[465] + amp[468] + amp[466] + amp[469] + amp[475] +
        amp[476] + amp[481] + amp[482] + amp[486] + amp[489] - amp[488] -
        amp[491] + amp[496] + amp[497]);
    jamp[111] = +2. * (+amp[1] + amp[2] + amp[3] + amp[4] - amp[6] + amp[8] +
        amp[10] - amp[11] + amp[24] + amp[25] + amp[27] + amp[34] + amp[35] -
        amp[36] + amp[38] - amp[192] + amp[194] + amp[196] - amp[197] + amp[198]
        + amp[199] + amp[201] + amp[203] - amp[205] + amp[206] - amp[207] +
        amp[216] + amp[217] - amp[223] - amp[224] + amp[261] - amp[263] +
        amp[264] + amp[271] - amp[273] + amp[274] - amp[275] + amp[276] -
        amp[278] - amp[302] + amp[304] + amp[306] + amp[348] - amp[350] +
        amp[354] - amp[356] + amp[357] - amp[363] - amp[359] + amp[365] -
        amp[367] - amp[368] - amp[370] - amp[371] - amp[378] - amp[381] -
        amp[379] - amp[382] - amp[459] - amp[460] + amp[462] + amp[463] -
        amp[466] - amp[469] - amp[467] - amp[470] - amp[495] - amp[496]);
    jamp[112] = +2. * (+amp[51] - amp[53] - amp[55] + amp[56] - amp[63] + amp[65]
        - amp[66] + amp[67] - amp[68] - amp[69] - amp[71] - amp[78] - amp[79] +
        amp[88] + amp[89] - amp[180] - amp[181] - amp[183] - amp[185] - amp[193]
        - amp[194] - amp[195] - amp[196] + amp[209] + amp[210] + amp[212] -
        amp[217] - amp[218] - amp[219] - amp[220] + amp[225] - amp[227] +
        amp[228] + amp[230] + amp[231] + amp[233] + amp[234] - amp[241] -
        amp[242] + amp[272] + amp[273] + amp[275] + amp[312] + amp[313] +
        amp[315] + amp[316] + amp[321] + amp[324] + amp[322] + amp[325] +
        amp[424] + amp[425] - amp[427] - amp[428] + amp[433] + amp[436] +
        amp[434] + amp[437] + amp[456] + amp[457] + amp[459] + amp[460] +
        amp[465] + amp[468] + amp[466] + amp[469] + amp[496] + amp[497]);
    jamp[113] = +2. * (+amp[6] - amp[8] - amp[10] + amp[11] - amp[18] + amp[20] -
        amp[21] + amp[22] - amp[23] - amp[24] - amp[26] - amp[33] - amp[34] +
        amp[43] + amp[44] + amp[180] + amp[181] + amp[183] + amp[185] - amp[203]
        + amp[205] - amp[206] + amp[207] + amp[208] - amp[210] + amp[211] -
        amp[216] + amp[218] + amp[222] + amp[223] - amp[225] + amp[227] -
        amp[228] - amp[230] - amp[231] - amp[233] - amp[234] + amp[241] +
        amp[242] + amp[302] - amp[304] - amp[306] - amp[312] + amp[314] -
        amp[315] + amp[317] - amp[321] + amp[327] - amp[322] + amp[328] +
        amp[370] + amp[371] - amp[373] - amp[374] + amp[379] + amp[382] +
        amp[380] + amp[383] - amp[456] - amp[457] - amp[462] - amp[463] -
        amp[465] - amp[468] + amp[467] + amp[470] + amp[495] - amp[497]);
    jamp[114] = +2. * (-amp[90] + amp[92] + amp[94] - amp[95] + amp[97] + amp[98]
        + amp[99] + amp[100] + amp[119] + amp[120] + amp[122] + amp[124] +
        amp[125] - amp[129] + amp[131] + amp[181] + amp[182] + amp[183] +
        amp[184] - amp[192] + amp[194] + amp[196] - amp[197] + amp[199] +
        amp[200] + amp[202] - amp[213] + amp[215] + amp[217] + amp[218] -
        amp[236] - amp[237] - amp[239] - amp[262] - amp[263] - amp[265] -
        amp[298] - amp[299] - amp[300] - amp[301] - amp[384] - amp[385] -
        amp[390] - amp[391] - amp[393] - amp[396] + amp[395] + amp[398] -
        amp[457] - amp[458] - amp[460] - amp[461] - amp[468] - amp[471] -
        amp[469] - amp[472] - amp[474] - amp[475] - amp[480] - amp[481] -
        amp[483] - amp[486] + amp[485] + amp[488] + amp[492] - amp[494] +
        amp[507] - amp[509]);
    jamp[115] = +2. * (-amp[45] + amp[47] + amp[49] - amp[50] + amp[52] + amp[53]
        + amp[54] + amp[55] + amp[74] + amp[75] + amp[77] + amp[79] + amp[80] -
        amp[84] + amp[86] + amp[187] + amp[188] + amp[189] + amp[190] + amp[192]
        + amp[193] + amp[195] + amp[197] - amp[199] - amp[200] - amp[202] +
        amp[213] + amp[214] + amp[220] + amp[221] - amp[254] - amp[255] -
        amp[257] + amp[262] + amp[263] + amp[265] + amp[298] + amp[299] +
        amp[300] + amp[301] + amp[384] - amp[386] + amp[390] - amp[392] +
        amp[393] - amp[399] - amp[395] + amp[401] - amp[421] - amp[422] -
        amp[424] - amp[425] - amp[432] - amp[435] - amp[433] - amp[436] -
        amp[477] - amp[478] + amp[480] + amp[481] - amp[484] - amp[487] -
        amp[485] - amp[488] - amp[492] - amp[493] - amp[507] + amp[509]);
    jamp[116] = +2. * (+amp[90] - amp[92] - amp[94] + amp[95] + amp[108] -
        amp[110] + amp[111] - amp[112] + amp[118] - amp[120] + amp[121] +
        amp[123] - amp[125] + amp[132] - amp[134] - amp[181] - amp[182] -
        amp[183] - amp[184] + amp[186] - amp[188] - amp[190] + amp[191] -
        amp[204] - amp[205] - amp[207] - amp[214] - amp[215] + amp[216] -
        amp[218] + amp[236] + amp[237] + amp[239] - amp[243] + amp[245] -
        amp[246] - amp[253] + amp[255] - amp[256] + amp[257] - amp[258] +
        amp[260] - amp[330] - amp[331] - amp[336] - amp[337] - amp[339] -
        amp[342] + amp[341] + amp[344] + amp[457] + amp[458] + amp[463] +
        amp[464] + amp[468] + amp[471] - amp[470] - amp[473] + amp[474] +
        amp[475] + amp[477] + amp[478] + amp[483] + amp[486] + amp[484] +
        amp[487] + amp[493] + amp[494]);
    jamp[117] = +2. * (-amp[0] + amp[2] + amp[4] - amp[5] + amp[7] + amp[8] +
        amp[9] + amp[10] + amp[29] + amp[30] + amp[32] + amp[34] + amp[35] -
        amp[39] + amp[41] - amp[186] + amp[188] + amp[190] - amp[191] + amp[198]
        - amp[200] + amp[201] - amp[202] + amp[203] + amp[204] + amp[206] +
        amp[213] + amp[214] - amp[223] - amp[224] + amp[243] - amp[245] +
        amp[246] + amp[253] - amp[255] + amp[256] - amp[257] + amp[258] -
        amp[260] - amp[297] + amp[299] + amp[301] + amp[330] - amp[332] +
        amp[336] - amp[338] + amp[339] - amp[345] - amp[341] + amp[347] -
        amp[367] - amp[368] - amp[370] - amp[371] - amp[378] - amp[381] -
        amp[379] - amp[382] - amp[477] - amp[478] + amp[480] + amp[481] -
        amp[484] - amp[487] - amp[485] - amp[488] - amp[492] - amp[493]);
    jamp[118] = +2. * (+amp[45] - amp[47] - amp[49] + amp[50] + amp[63] - amp[65]
        + amp[66] - amp[67] + amp[73] - amp[75] + amp[76] + amp[78] - amp[80] +
        amp[87] - amp[89] + amp[180] - amp[182] - amp[184] + amp[185] - amp[187]
        - amp[188] - amp[189] - amp[190] - amp[209] - amp[210] - amp[212] -
        amp[214] - amp[215] + amp[219] - amp[221] - amp[225] + amp[227] -
        amp[228] - amp[235] + amp[237] - amp[238] + amp[239] - amp[240] +
        amp[242] + amp[254] + amp[255] + amp[257] - amp[312] - amp[313] -
        amp[318] - amp[319] - amp[321] - amp[324] + amp[323] + amp[326] +
        amp[421] + amp[422] + amp[427] + amp[428] + amp[432] + amp[435] -
        amp[434] - amp[437] + amp[474] + amp[475] + amp[477] + amp[478] +
        amp[483] + amp[486] + amp[484] + amp[487] + amp[493] + amp[494]);
    jamp[119] = +2. * (+amp[0] - amp[2] - amp[4] + amp[5] + amp[18] - amp[20] +
        amp[21] - amp[22] + amp[28] - amp[30] + amp[31] + amp[33] - amp[35] +
        amp[42] - amp[44] - amp[180] + amp[182] + amp[184] - amp[185] - amp[198]
        + amp[200] - amp[201] + amp[202] - amp[208] + amp[210] - amp[211] -
        amp[213] + amp[215] - amp[222] + amp[224] + amp[225] - amp[227] +
        amp[228] + amp[235] - amp[237] + amp[238] - amp[239] + amp[240] -
        amp[242] + amp[297] - amp[299] - amp[301] + amp[312] - amp[314] +
        amp[318] - amp[320] + amp[321] - amp[327] - amp[323] + amp[329] +
        amp[367] + amp[368] + amp[373] + amp[374] + amp[378] + amp[381] -
        amp[380] - amp[383] - amp[474] - amp[475] - amp[480] - amp[481] -
        amp[483] - amp[486] + amp[485] + amp[488] + amp[492] - amp[494]);
}

std::vector<double> AMP_MG5_gg_gggg::ColorDenominators() const { return {648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648,
      648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648, 648}; }

std::vector<std::vector<double>> AMP_MG5_gg_gggg::ColorMetric() const {
  return {{3641, -454, -454, 50, 50, 572,
      -454, 68, 50, 50, -13, -76, 50, -13, 572, -76, 554, 428, 50, -76, -76,
      428, 428, -184, -454, 68, 68, -76, -76, -40, 50, -76, 50, -454, -13, 68,
      -13, 5, -76, 68, -13, -76, -13, -58, 5, -76, -58, -40, 50, -76, -13, -13,
      5, -58, 572, -40, -76, 68, -58, 104, 554, -13, 428, -76, 50, 572, -13, 5,
      -58, -58, -76, -40, 50, -13, -76, 5, -13, -58, -76, -58, 428, -76, 5,
      -58, 428, -58, -184, -40, -76, -40, 50, -76, -76, 428, 428, -184, -454,
      68, 68, -76, -76, -40, 68, 104, -76, 572, -58, -40, -76, -58, -40, -40,
      428, -184, 572, -40, -40, -184, -184, 176}, {-454, 3641, 50, 572, -454,
      50, 68, -454, -13, -76, 50, 50, 50, -76, -76, 428, 428, -184, 50, -13,
      572, -76, 554, 428, 68, -454, -76, -40, 68, -76, -76, 50, -13, 68, 50,
      -454, -13, -58, 5, -76, -58, -40, -13, 5, -76, 68, -13, -76, -76, 50, 5,
      -58, -13, -13, -40, 572, -58, 104, -76, 68, -13, 5, -58, -58, -76, -40,
      554, -13, 428, -76, 50, 572, -454, 68, 68, -76, -76, -40, 68, 104, -76,
      572, -58, -40, -76, -58, -40, -40, 428, -184, 572, -40, -40, -184, -184,
      176, 50, -13, -76, 5, -13, -58, -76, -58, 428, -76, 5, -58, 428, -58,
      -184, -40, -76, -40, 50, -76, -76, 428, 428, -184}, {-454, 50, 3641,
      -454, 572, 50, 50, -13, 572, -76, 554, 428, -454, 68, 50, 50, -13, -76,
      -76, 50, 428, -184, -76, 428, 68, -76, -454, 68, -40, -76, -13, 5, -76,
      68, -13, -76, 50, -76, 50, -454, -13, 68, -58, -13, -58, -40, 5, -76, 50,
      -13, -76, 5, -13, -58, -76, -58, 428, -76, 5, -58, 428, -58, -184, -40,
      -76, -40, 50, -76, -76, 428, 428, -184, 50, -76, -13, -13, 5, -58, 572,
      -40, -76, 68, -58, 104, 554, -13, 428, -76, 50, 572, -13, 5, -58, -58,
      -76, -40, 68, -454, -76, -40, 68, -76, 104, 68, -58, -40, -76, 572, 572,
      -40, -40, -184, -184, 176, -76, -58, -40, -40, 428, -184}, {50, 572,
      -454, 3641, 50, -454, 50, -76, -76, 428, 428, -184, 68, -454, -13, -76,
      50, 50, -13, 50, 554, 428, 572, -76, -76, -40, 68, -454, -76, 68, -13,
      -58, 5, -76, -58, -40, -76, 50, -13, 68, 50, -454, 5, -13, -13, -76, -76,
      68, -454, 68, 68, -76, -76, -40, 68, 104, -76, 572, -58, -40, -76, -58,
      -40, -40, 428, -184, 572, -40, -40, -184, -184, 176, -76, 50, 5, -58,
      -13, -13, -40, 572, -58, 104, -76, 68, -13, 5, -58, -58, -76, -40, 554,
      -13, 428, -76, 50, 572, -13, 50, -13, -58, -76, 5, -58, -76, 5, -58, 428,
      -76, 50, -76, -76, 428, 428, -184, 428, -58, -184, -40, -76, -40}, {50,
      -454, 572, 50, 3641, -454, -13, 50, 554, 428, 572, -76, -76, 50, 428,
      -184, -76, 428, -454, 68, 50, 50, -13, -76, -76, 68, -40, -76, -454, 68,
      5, -13, -13, -76, -76, 68, -58, -13, -58, -40, 5, -76, 50, -76, 50, -454,
      -13, 68, -13, 50, -13, -58, -76, 5, -58, -76, 5, -58, 428, -76, 50, -76,
      -76, 428, 428, -184, 428, -58, -184, -40, -76, -40, 68, -454, -76, -40,
      68, -76, 104, 68, -58, -40, -76, 572, 572, -40, -40, -184, -184, 176,
      -76, -58, -40, -40, 428, -184, 50, -76, -13, -13, 5, -58, 572, -40, -76,
      68, -58, 104, 554, -13, 428, -76, 50, 572, -13, 5, -58, -58, -76, -40},
      {572, 50, 50, -454, -454, 3641, -76, 50, 428, -184, -76, 428, -13, 50,
      554, 428, 572, -76, 68, -454, -13, -76, 50, 50, -40, -76, -76, 68, 68,
      -454, -58, -13, -58, -40, 5, -76, 5, -13, -13, -76, -76, 68, -76, 50,
      -13, 68, 50, -454, 68, -454, -76, -40, 68, -76, 104, 68, -58, -40, -76,
      572, 572, -40, -40, -184, -184, 176, -76, -58, -40, -40, 428, -184, -13,
      50, -13, -58, -76, 5, -58, -76, 5, -58, 428, -76, 50, -76, -76, 428, 428,
      -184, 428, -58, -184, -40, -76, -40, -76, 50, 5, -58, -13, -13, -40, 572,
      -58, 104, -76, 68, -13, 5, -58, -58, -76, -40, 554, -13, 428, -76, 50,
      572}, {-454, 68, 50, 50, -13, -76, 3641, -454, -454, 50, 50, 572, 572,
      -76, 50, -13, 428, 554, -76, 428, 50, -76, -184, 428, 50, -76, -13, -13,
      5, -58, 572, -40, -76, 68, -58, 104, 554, -13, 428, -76, 50, 572, -13, 5,
      -58, -58, -76, -40, -454, 68, 68, -76, -76, -40, 50, -76, 50, -454, -13,
      68, -13, 5, -76, 68, -13, -76, -13, -58, 5, -76, -58, -40, -76, 5, 50,
      -13, -58, -13, 428, -58, -184, -40, -76, -40, -76, -58, 428, -76, 5, -58,
      -76, 50, 428, -184, -76, 428, 68, -76, -454, 68, -40, -76, -76, -58, -40,
      -40, 428, -184, 68, 104, -76, 572, -58, -40, -40, 572, -184, 176, -40,
      -184}, {68, -454, -13, -76, 50, 50, -454, 3641, 50, 572, -454, 50, -76,
      428, 50, -76, -184, 428, 572, -76, 50, -13, 428, 554, -76, 50, 5, -58,
      -13, -13, -40, 572, -58, 104, -76, 68, -13, 5, -58, -58, -76, -40, 554,
      -13, 428, -76, 50, 572, 68, -454, -76, -40, 68, -76, -76, 50, -13, 68,
      50, -454, -13, -58, 5, -76, -58, -40, -13, 5, -76, 68, -13, -76, 68, -76,
      -454, 68, -40, -76, -76, -58, -40, -40, 428, -184, 68, 104, -76, 572,
      -58, -40, -40, 572, -184, 176, -40, -184, -76, 5, 50, -13, -58, -13, 428,
      -58, -184, -40, -76, -40, -76, -58, 428, -76, 5, -58, -76, 50, 428, -184,
      -76, 428}, {50, -13, 572, -76, 554, 428, -454, 50, 3641, -454, 572, 50,
      50, 50, -454, 68, -76, -13, 428, -184, -76, 50, 428, -76, 50, -13, -76,
      5, -13, -58, -76, -58, 428, -76, 5, -58, 428, -58, -184, -40, -76, -40,
      50, -76, -76, 428, 428, -184, 68, -76, -454, 68, -40, -76, -13, 5, -76,
      68, -13, -76, 50, -76, 50, -454, -13, 68, -58, -13, -58, -40, 5, -76,
      -13, -13, 50, -76, -58, 5, 554, -13, 428, -76, 50, 572, 572, -40, -76,
      68, -58, 104, 5, -13, -76, -40, -58, -58, -76, -40, 68, -454, -76, 68,
      572, -40, -40, -184, -184, 176, 104, 68, -58, -40, -76, 572, -58, -76,
      428, -184, -40, -40}, {50, -76, -76, 428, 428, -184, 50, 572, -454, 3641,
      50, -454, -13, -76, 68, -454, 50, 50, 554, 428, -13, 50, -76, 572, -454,
      68, 68, -76, -76, -40, 68, 104, -76, 572, -58, -40, -76, -58, -40, -40,
      428, -184, 572, -40, -40, -184, -184, 176, -76, -40, 68, -454, -76, 68,
      -13, -58, 5, -76, -58, -40, -76, 50, -13, 68, 50, -454, 5, -13, -13, -76,
      -76, 68, 5, -58, -76, 50, -13, -13, -13, 5, -58, -58, -76, -40, -40, 572,
      -58, 104, -76, 68, -13, 554, 50, 572, 428, -76, -13, -58, -13, 50, 5,
      -76, 50, -76, -76, 428, 428, -184, -58, -76, 5, -58, 428, -76, -58, 428,
      -76, -40, -184, -40}, {-13, 50, 554, 428, 572, -76, 50, -454, 572, 50,
      3641, -454, 428, -184, -76, 50, 428, -76, 50, 50, -454, 68, -76, -13,
      -13, 50, -13, -58, -76, 5, -58, -76, 5, -58, 428, -76, 50, -76, -76, 428,
      428, -184, 428, -58, -184, -40, -76, -40, -76, 68, -40, -76, -454, 68, 5,
      -13, -13, -76, -76, 68, -58, -13, -58, -40, 5, -76, 50, -76, 50, -454,
      -13, 68, -76, -40, 68, -454, -76, 68, 572, -40, -40, -184, -184, 176,
      104, 68, -58, -40, -76, 572, -58, -76, 428, -184, -40, -40, -13, -13, 50,
      -76, -58, 5, 554, -13, 428, -76, 50, 572, 572, -40, -76, 68, -58, 104, 5,
      -13, -76, -40, -58, -58}, {-76, 50, 428, -184, -76, 428, 572, 50, 50,
      -454, -454, 3641, 554, 428, -13, 50, -76, 572, -13, -76, 68, -454, 50,
      50, 68, -454, -76, -40, 68, -76, 104, 68, -58, -40, -76, 572, 572, -40,
      -40, -184, -184, 176, -76, -58, -40, -40, 428, -184, -40, -76, -76, 68,
      68, -454, -58, -13, -58, -40, 5, -76, 5, -13, -13, -76, -76, 68, -76, 50,
      -13, 68, 50, -454, -13, -58, -13, 50, 5, -76, 50, -76, -76, 428, 428,
      -184, -58, -76, 5, -58, 428, -76, -58, 428, -76, -40, -184, -40, 5, -58,
      -76, 50, -13, -13, -13, 5, -58, -58, -76, -40, -40, 572, -58, 104, -76,
      68, -13, 554, 50, 572, 428, -76}, {50, 50, -454, 68, -76, -13, 572, -76,
      50, -13, 428, 554, 3641, -454, -454, 50, 50, 572, 428, -76, -184, 428,
      50, -76, -13, -13, 50, -76, -58, 5, 554, -13, 428, -76, 50, 572, 572,
      -40, -76, 68, -58, 104, 5, -13, -76, -40, -58, -58, -76, 5, 50, -13, -58,
      -13, 428, -58, -184, -40, -76, -40, -76, -58, 428, -76, 5, -58, -76, 50,
      428, -184, -76, 428, -454, 68, 68, -76, -76, -40, 50, -76, 50, -454, -13,
      68, -13, 5, -76, 68, -13, -76, -13, -58, 5, -76, -58, -40, -76, 68, -40,
      -76, -454, 68, -58, -76, 428, -184, -40, -40, -40, 572, -184, 176, -40,
      -184, 68, 104, -76, 572, -58, -40}, {-13, -76, 68, -454, 50, 50, -76,
      428, 50, -76, -184, 428, -454, 3641, 50, 572, -454, 50, -76, 572, 428,
      554, 50, -13, 5, -58, -76, 50, -13, -13, -13, 5, -58, -58, -76, -40, -40,
      572, -58, 104, -76, 68, -13, 554, 50, 572, 428, -76, 68, -76, -454, 68,
      -40, -76, -76, -58, -40, -40, 428, -184, 68, 104, -76, 572, -58, -40,
      -40, 572, -184, 176, -40, -184, 68, -454, -76, -40, 68, -76, -76, 50,
      -13, 68, 50, -454, -13, -58, 5, -76, -58, -40, -13, 5, -76, 68, -13, -76,
      5, -76, -58, -13, 50, -13, -58, 428, -76, -40, -184, -40, -76, 50, 428,
      -184, -76, 428, -76, -58, 428, -76, 5, -58}, {572, -76, 50, -13, 428,
      554, 50, 50, -454, 68, -76, -13, -454, 50, 3641, -454, 572, 50, -184,
      428, 428, -76, -76, 50, -76, 5, 50, -13, -58, -13, 428, -58, -184, -40,
      -76, -40, -76, -58, 428, -76, 5, -58, -76, 50, 428, -184, -76, 428, -13,
      -13, 50, -76, -58, 5, 554, -13, 428, -76, 50, 572, 572, -40, -76, 68,
      -58, 104, 5, -13, -76, -40, -58, -58, 68, -76, -454, 68, -40, -76, -13,
      5, -76, 68, -13, -76, 50, -76, 50, -454, -13, 68, -58, -13, -58, -40, 5,
      -76, -40, -76, -76, 68, 68, -454, -40, 572, -184, 176, -40, -184, -58,
      -76, 428, -184, -40, -40, 104, 68, -58, -40, -76, 572}, {-76, 428, 50,
      -76, -184, 428, -13, -76, 68, -454, 50, 50, 50, 572, -454, 3641, 50,
      -454, 428, 554, -76, 572, -13, 50, 68, -76, -454, 68, -40, -76, -76, -58,
      -40, -40, 428, -184, 68, 104, -76, 572, -58, -40, -40, 572, -184, 176,
      -40, -184, 5, -58, -76, 50, -13, -13, -13, 5, -58, -58, -76, -40, -40,
      572, -58, 104, -76, 68, -13, 554, 50, 572, 428, -76, -76, -40, 68, -454,
      -76, 68, -13, -58, 5, -76, -58, -40, -76, 50, -13, 68, 50, -454, 5, -13,
      -13, -76, -76, 68, -58, -13, 5, -76, -13, 50, -76, 50, 428, -184, -76,
      428, -58, 428, -76, -40, -184, -40, -58, -76, 5, -58, 428, -76}, {554,
      428, -13, 50, -76, 572, 428, -184, -76, 50, 428, -76, 50, -454, 572, 50,
      3641, -454, 50, 50, -76, -13, -454, 68, -13, -58, -13, 50, 5, -76, 50,
      -76, -76, 428, 428, -184, -58, -76, 5, -58, 428, -76, -58, 428, -76, -40,
      -184, -40, -76, -40, 68, -454, -76, 68, 572, -40, -40, -184, -184, 176,
      104, 68, -58, -40, -76, 572, -58, -76, 428, -184, -40, -40, -76, 68, -40,
      -76, -454, 68, 5, -13, -13, -76, -76, 68, -58, -13, -58, -40, 5, -76, 50,
      -76, 50, -454, -13, 68, -13, -13, -58, 5, 50, -76, -13, 554, 50, 572,
      428, -76, 5, -13, -76, -40, -58, -58, 572, -40, -76, 68, -58, 104}, {428,
      -184, -76, 50, 428, -76, 554, 428, -13, 50, -76, 572, 572, 50, 50, -454,
      -454, 3641, -76, -13, 50, 50, 68, -454, -76, -40, 68, -454, -76, 68, 572,
      -40, -40, -184, -184, 176, 104, 68, -58, -40, -76, 572, -58, -76, 428,
      -184, -40, -40, -13, -58, -13, 50, 5, -76, 50, -76, -76, 428, 428, -184,
      -58, -76, 5, -58, 428, -76, -58, 428, -76, -40, -184, -40, -40, -76, -76,
      68, 68, -454, -58, -13, -58, -40, 5, -76, 5, -13, -13, -76, -76, 68, -76,
      50, -13, 68, 50, -454, -58, 5, -13, -13, -76, 50, 5, -13, -76, -40, -58,
      -58, -13, 554, 50, 572, 428, -76, -40, 572, -58, 104, -76, 68}, {50, 50,
      -76, -13, -454, 68, -76, 572, 428, 554, 50, -13, 428, -76, -184, 428, 50,
      -76, 3641, -454, -454, 50, 50, 572, -13, -13, -58, 5, 50, -76, -13, 554,
      50, 572, 428, -76, 5, -13, -76, -40, -58, -58, 572, -40, -76, 68, -58,
      104, 5, -76, -58, -13, 50, -13, -58, 428, -76, -40, -184, -40, -76, 50,
      428, -184, -76, 428, -76, -58, 428, -76, 5, -58, -76, 68, -40, -76, -454,
      68, -58, -76, 428, -184, -40, -40, -40, 572, -184, 176, -40, -184, 68,
      104, -76, 572, -58, -40, -454, 68, 68, -76, -76, -40, 50, -76, 50, -454,
      -13, 68, -13, 5, -76, 68, -13, -76, -13, -58, 5, -76, -58, -40}, {-76,
      -13, 50, 50, 68, -454, 428, -76, -184, 428, 50, -76, -76, 572, 428, 554,
      50, -13, -454, 3641, 50, 572, -454, 50, -58, 5, -13, -13, -76, 50, 5,
      -13, -76, -40, -58, -58, -13, 554, 50, 572, 428, -76, -40, 572, -58, 104,
      -76, 68, -76, 68, -40, -76, -454, 68, -58, -76, 428, -184, -40, -40, -40,
      572, -184, 176, -40, -184, 68, 104, -76, 572, -58, -40, 5, -76, -58, -13,
      50, -13, -58, 428, -76, -40, -184, -40, -76, 50, 428, -184, -76, 428,
      -76, -58, 428, -76, 5, -58, 68, -454, -76, -40, 68, -76, -76, 50, -13,
      68, 50, -454, -13, -58, 5, -76, -58, -40, -13, 5, -76, 68, -13, -76},
      {-76, 572, 428, 554, 50, -13, 50, 50, -76, -13, -454, 68, -184, 428, 428,
      -76, -76, 50, -454, 50, 3641, -454, 572, 50, 5, -76, -58, -13, 50, -13,
      -58, 428, -76, -40, -184, -40, -76, 50, 428, -184, -76, 428, -76, -58,
      428, -76, 5, -58, -13, -13, -58, 5, 50, -76, -13, 554, 50, 572, 428, -76,
      5, -13, -76, -40, -58, -58, 572, -40, -76, 68, -58, 104, -40, -76, -76,
      68, 68, -454, -40, 572, -184, 176, -40, -184, -58, -76, 428, -184, -40,
      -40, 104, 68, -58, -40, -76, 572, 68, -76, -454, 68, -40, -76, -13, 5,
      -76, 68, -13, -76, 50, -76, 50, -454, -13, 68, -58, -13, -58, -40, 5,
      -76}, {428, -76, -184, 428, 50, -76, -76, -13, 50, 50, 68, -454, 428,
      554, -76, 572, -13, 50, 50, 572, -454, 3641, 50, -454, -76, 68, -40, -76,
      -454, 68, -58, -76, 428, -184, -40, -40, -40, 572, -184, 176, -40, -184,
      68, 104, -76, 572, -58, -40, -58, 5, -13, -13, -76, 50, 5, -13, -76, -40,
      -58, -58, -13, 554, 50, 572, 428, -76, -40, 572, -58, 104, -76, 68, -58,
      -13, 5, -76, -13, 50, -76, 50, 428, -184, -76, 428, -58, 428, -76, -40,
      -184, -40, -58, -76, 5, -58, 428, -76, -76, -40, 68, -454, -76, 68, -13,
      -58, 5, -76, -58, -40, -76, 50, -13, 68, 50, -454, 5, -13, -13, -76, -76,
      68}, {428, 554, -76, 572, -13, 50, -184, 428, 428, -76, -76, 50, 50, 50,
      -76, -13, -454, 68, 50, -454, 572, 50, 3641, -454, -58, -13, 5, -76, -13,
      50, -76, 50, 428, -184, -76, 428, -58, 428, -76, -40, -184, -40, -58,
      -76, 5, -58, 428, -76, -40, -76, -76, 68, 68, -454, -40, 572, -184, 176,
      -40, -184, -58, -76, 428, -184, -40, -40, 104, 68, -58, -40, -76, 572,
      -13, -13, -58, 5, 50, -76, -13, 554, 50, 572, 428, -76, 5, -13, -76, -40,
      -58, -58, 572, -40, -76, 68, -58, 104, -76, 68, -40, -76, -454, 68, 5,
      -13, -13, -76, -76, 68, -58, -13, -58, -40, 5, -76, 50, -76, 50, -454,
      -13, 68}, {-184, 428, 428, -76, -76, 50, 428, 554, -76, 572, -13, 50,
      -76, -13, 50, 50, 68, -454, 572, 50, 50, -454, -454, 3641, -40, -76, -76,
      68, 68, -454, -40, 572, -184, 176, -40, -184, -58, -76, 428, -184, -40,
      -40, 104, 68, -58, -40, -76, 572, -58, -13, 5, -76, -13, 50, -76, 50,
      428, -184, -76, 428, -58, 428, -76, -40, -184, -40, -58, -76, 5, -58,
      428, -76, -58, 5, -13, -13, -76, 50, 5, -13, -76, -40, -58, -58, -13,
      554, 50, 572, 428, -76, -40, 572, -58, 104, -76, 68, -40, -76, -76, 68,
      68, -454, -58, -13, -58, -40, 5, -76, 5, -13, -13, -76, -76, 68, -76, 50,
      -13, 68, 50, -454}, {-454, 68, 68, -76, -76, -40, 50, -76, 50, -454, -13,
      68, -13, 5, -76, 68, -13, -76, -13, -58, 5, -76, -58, -40, 3641, -454,
      -454, 50, 50, 572, -454, 68, 50, 50, -13, -76, 50, -13, 572, -76, 554,
      428, 50, -76, -76, 428, 428, -184, 572, -40, -76, 68, -58, 104, 50, -76,
      -13, -13, 5, -58, 428, -76, 554, -13, 572, 50, -58, -58, -13, 5, -40,
      -76, -76, -58, 428, -76, 5, -58, 50, -13, -76, 5, -13, -58, -184, -40,
      428, -58, -40, -76, -76, 428, 50, -76, -184, 428, 68, 104, -76, 572, -58,
      -40, -454, 68, 68, -76, -76, -40, -40, -40, -76, -58, -184, 428, -40,
      -184, 572, -40, 176, -184}, {68, -454, -76, -40, 68, -76, -76, 50, -13,
      68, 50, -454, -13, -58, 5, -76, -58, -40, -13, 5, -76, 68, -13, -76,
      -454, 3641, 50, 572, -454, 50, 68, -454, -13, -76, 50, 50, 50, -76, -76,
      428, 428, -184, 50, -13, 572, -76, 554, 428, -40, 572, -58, 104, -76, 68,
      -76, 50, 5, -58, -13, -13, -58, -58, -13, 5, -40, -76, 428, -76, 554,
      -13, 572, 50, 68, 104, -76, 572, -58, -40, -454, 68, 68, -76, -76, -40,
      -40, -40, -76, -58, -184, 428, -40, -184, 572, -40, 176, -184, -76, -58,
      428, -76, 5, -58, 50, -13, -76, 5, -13, -58, -184, -40, 428, -58, -40,
      -76, -76, 428, 50, -76, -184, 428}, {68, -76, -454, 68, -40, -76, -13, 5,
      -76, 68, -13, -76, 50, -76, 50, -454, -13, 68, -58, -13, -58, -40, 5,
      -76, -454, 50, 3641, -454, 572, 50, 50, -13, 572, -76, 554, 428, -454,
      68, 50, 50, -13, -76, -76, 50, 428, -184, -76, 428, -76, -58, 428, -76,
      5, -58, 50, -13, -76, 5, -13, -58, -184, -40, 428, -58, -40, -76, -76,
      428, 50, -76, -184, 428, 572, -40, -76, 68, -58, 104, 50, -76, -13, -13,
      5, -58, 428, -76, 554, -13, 572, 50, -58, -58, -13, 5, -40, -76, 104, 68,
      -58, -40, -76, 572, 68, -454, -76, -40, 68, -76, -40, -184, 572, -40,
      176, -184, -40, -40, -76, -58, -184, 428}, {-76, -40, 68, -454, -76, 68,
      -13, -58, 5, -76, -58, -40, -76, 50, -13, 68, 50, -454, 5, -13, -13, -76,
      -76, 68, 50, 572, -454, 3641, 50, -454, 50, -76, -76, 428, 428, -184, 68,
      -454, -13, -76, 50, 50, -13, 50, 554, 428, 572, -76, 68, 104, -76, 572,
      -58, -40, -454, 68, 68, -76, -76, -40, -40, -40, -76, -58, -184, 428,
      -40, -184, 572, -40, 176, -184, -40, 572, -58, 104, -76, 68, -76, 50, 5,
      -58, -13, -13, -58, -58, -13, 5, -40, -76, 428, -76, 554, -13, 572, 50,
      -58, -76, 5, -58, 428, -76, -13, 50, -13, -58, -76, 5, -76, 428, 50, -76,
      -184, 428, -184, -40, 428, -58, -40, -76}, {-76, 68, -40, -76, -454, 68,
      5, -13, -13, -76, -76, 68, -58, -13, -58, -40, 5, -76, 50, -76, 50, -454,
      -13, 68, 50, -454, 572, 50, 3641, -454, -13, 50, 554, 428, 572, -76, -76,
      50, 428, -184, -76, 428, -454, 68, 50, 50, -13, -76, -58, -76, 5, -58,
      428, -76, -13, 50, -13, -58, -76, 5, -76, 428, 50, -76, -184, 428, -184,
      -40, 428, -58, -40, -76, 104, 68, -58, -40, -76, 572, 68, -454, -76, -40,
      68, -76, -40, -184, 572, -40, 176, -184, -40, -40, -76, -58, -184, 428,
      572, -40, -76, 68, -58, 104, 50, -76, -13, -13, 5, -58, 428, -76, 554,
      -13, 572, 50, -58, -58, -13, 5, -40, -76}, {-40, -76, -76, 68, 68, -454,
      -58, -13, -58, -40, 5, -76, 5, -13, -13, -76, -76, 68, -76, 50, -13, 68,
      50, -454, 572, 50, 50, -454, -454, 3641, -76, 50, 428, -184, -76, 428,
      -13, 50, 554, 428, 572, -76, 68, -454, -13, -76, 50, 50, 104, 68, -58,
      -40, -76, 572, 68, -454, -76, -40, 68, -76, -40, -184, 572, -40, 176,
      -184, -40, -40, -76, -58, -184, 428, -58, -76, 5, -58, 428, -76, -13, 50,
      -13, -58, -76, 5, -76, 428, 50, -76, -184, 428, -184, -40, 428, -58, -40,
      -76, -40, 572, -58, 104, -76, 68, -76, 50, 5, -58, -13, -13, -58, -58,
      -13, 5, -40, -76, 428, -76, 554, -13, 572, 50}, {50, -76, -13, -13, 5,
      -58, 572, -40, -76, 68, -58, 104, 554, -13, 428, -76, 50, 572, -13, 5,
      -58, -58, -76, -40, -454, 68, 50, 50, -13, -76, 3641, -454, -454, 50, 50,
      572, 572, -76, 50, -13, 428, 554, -76, 428, 50, -76, -184, 428, 50, -76,
      50, -454, -13, 68, -454, 68, 68, -76, -76, -40, -76, 68, -13, 5, -76,
      -13, 5, -76, -13, -58, -40, -58, 428, -58, -184, -40, -76, -40, -76, 5,
      50, -13, -58, -13, 428, -76, -76, -58, -58, 5, 428, -184, -76, 50, 428,
      -76, -76, -58, -40, -40, 428, -184, 68, -76, -454, 68, -40, -76, -76,
      572, 68, 104, -40, -58, -184, 176, -40, 572, -184, -40}, {-76, 50, 5,
      -58, -13, -13, -40, 572, -58, 104, -76, 68, -13, 5, -58, -58, -76, -40,
      554, -13, 428, -76, 50, 572, 68, -454, -13, -76, 50, 50, -454, 3641, 50,
      572, -454, 50, -76, 428, 50, -76, -184, 428, 572, -76, 50, -13, 428, 554,
      -76, 50, -13, 68, 50, -454, 68, -454, -76, -40, 68, -76, 5, -76, -13,
      -58, -40, -58, -76, 68, -13, 5, -76, -13, -76, -58, -40, -40, 428, -184,
      68, -76, -454, 68, -40, -76, -76, 572, 68, 104, -40, -58, -184, 176, -40,
      572, -184, -40, 428, -58, -184, -40, -76, -40, -76, 5, 50, -13, -58, -13,
      428, -76, -76, -58, -58, 5, 428, -184, -76, 50, 428, -76}, {50, -13, -76,
      5, -13, -58, -76, -58, 428, -76, 5, -58, 428, -58, -184, -40, -76, -40,
      50, -76, -76, 428, 428, -184, 50, -13, 572, -76, 554, 428, -454, 50,
      3641, -454, 572, 50, 50, 50, -454, 68, -76, -13, 428, -184, -76, 50, 428,
      -76, -13, 5, -76, 68, -13, -76, 68, -76, -454, 68, -40, -76, 50, -454,
      50, -76, 68, -13, -58, -40, -58, -13, -76, 5, 554, -13, 428, -76, 50,
      572, -13, -13, 50, -76, -58, 5, -76, 68, 572, -40, 104, -58, -76, -40, 5,
      -13, -58, -58, 572, -40, -40, -184, -184, 176, -76, -40, 68, -454, -76,
      68, -58, -40, 104, 68, 572, -76, 428, -184, -58, -76, -40, -40}, {-454,
      68, 68, -76, -76, -40, 68, 104, -76, 572, -58, -40, -76, -58, -40, -40,
      428, -184, 572, -40, -40, -184, -184, 176, 50, -76, -76, 428, 428, -184,
      50, 572, -454, 3641, 50, -454, -13, -76, 68, -454, 50, 50, 554, 428, -13,
      50, -76, 572, -13, -58, 5, -76, -58, -40, -76, -40, 68, -454, -76, 68,
      -13, 68, -76, 50, -454, 50, -13, -76, 5, -13, 68, -76, -13, 5, -58, -58,
      -76, -40, 5, -58, -76, 50, -13, -13, -58, 104, -40, 572, 68, -76, 50,
      572, -13, 554, -76, 428, 50, -76, -76, 428, 428, -184, -13, -58, -13, 50,
      5, -76, 5, -58, -58, -76, -76, 428, -76, -40, -58, 428, -40, -184}, {-13,
      50, -13, -58, -76, 5, -58, -76, 5, -58, 428, -76, 50, -76, -76, 428, 428,
      -184, 428, -58, -184, -40, -76, -40, -13, 50, 554, 428, 572, -76, 50,
      -454, 572, 50, 3641, -454, 428, -184, -76, 50, 428, -76, 50, 50, -454,
      68, -76, -13, 5, -13, -13, -76, -76, 68, -76, 68, -40, -76, -454, 68,
      -58, -40, -58, -13, -76, 5, 50, -454, 50, -76, 68, -13, 572, -40, -40,
      -184, -184, 176, -76, -40, 68, -454, -76, 68, -58, -40, 104, 68, 572,
      -76, 428, -184, -58, -76, -40, -40, 554, -13, 428, -76, 50, 572, -13,
      -13, 50, -76, -58, 5, -76, 68, 572, -40, 104, -58, -76, -40, 5, -13, -58,
      -58}, {68, -454, -76, -40, 68, -76, 104, 68, -58, -40, -76, 572, 572,
      -40, -40, -184, -184, 176, -76, -58, -40, -40, 428, -184, -76, 50, 428,
      -184, -76, 428, 572, 50, 50, -454, -454, 3641, 554, 428, -13, 50, -76,
      572, -13, -76, 68, -454, 50, 50, -58, -13, -58, -40, 5, -76, -40, -76,
      -76, 68, 68, -454, -13, -76, 5, -13, 68, -76, -13, 68, -76, 50, -454, 50,
      50, -76, -76, 428, 428, -184, -13, -58, -13, 50, 5, -76, 5, -58, -58,
      -76, -76, 428, -76, -40, -58, 428, -40, -184, -13, 5, -58, -58, -76, -40,
      5, -58, -76, 50, -13, -13, -58, 104, -40, 572, 68, -76, 50, 572, -13,
      554, -76, 428}, {-13, -13, 50, -76, -58, 5, 554, -13, 428, -76, 50, 572,
      572, -40, -76, 68, -58, 104, 5, -13, -76, -40, -58, -58, 50, 50, -454,
      68, -76, -13, 572, -76, 50, -13, 428, 554, 3641, -454, -454, 50, 50, 572,
      428, -76, -184, 428, 50, -76, 428, -58, -184, -40, -76, -40, -76, 5, 50,
      -13, -58, -13, 428, -76, -76, -58, -58, 5, 428, -184, -76, 50, 428, -76,
      50, -76, 50, -454, -13, 68, -454, 68, 68, -76, -76, -40, -76, 68, -13, 5,
      -76, -13, 5, -76, -13, -58, -40, -58, -58, -76, 428, -184, -40, -40, -76,
      68, -40, -76, -454, 68, -184, 176, -40, 572, -184, -40, -76, 572, 68,
      104, -40, -58}, {5, -58, -76, 50, -13, -13, -13, 5, -58, -58, -76, -40,
      -40, 572, -58, 104, -76, 68, -13, 554, 50, 572, 428, -76, -13, -76, 68,
      -454, 50, 50, -76, 428, 50, -76, -184, 428, -454, 3641, 50, 572, -454,
      50, -76, 572, 428, 554, 50, -13, -76, -58, -40, -40, 428, -184, 68, -76,
      -454, 68, -40, -76, -76, 572, 68, 104, -40, -58, -184, 176, -40, 572,
      -184, -40, -76, 50, -13, 68, 50, -454, 68, -454, -76, -40, 68, -76, 5,
      -76, -13, -58, -40, -58, -76, 68, -13, 5, -76, -13, -58, 428, -76, -40,
      -184, -40, 5, -76, -58, -13, 50, -13, 428, -184, -76, 50, 428, -76, 428,
      -76, -76, -58, -58, 5}, {-76, 5, 50, -13, -58, -13, 428, -58, -184, -40,
      -76, -40, -76, -58, 428, -76, 5, -58, -76, 50, 428, -184, -76, 428, 572,
      -76, 50, -13, 428, 554, 50, 50, -454, 68, -76, -13, -454, 50, 3641, -454,
      572, 50, -184, 428, 428, -76, -76, 50, 554, -13, 428, -76, 50, 572, -13,
      -13, 50, -76, -58, 5, -76, 68, 572, -40, 104, -58, -76, -40, 5, -13, -58,
      -58, -13, 5, -76, 68, -13, -76, 68, -76, -454, 68, -40, -76, 50, -454,
      50, -76, 68, -13, -58, -40, -58, -13, -76, 5, -40, 572, -184, 176, -40,
      -184, -40, -76, -76, 68, 68, -454, 428, -184, -58, -76, -40, -40, -58,
      -40, 104, 68, 572, -76}, {68, -76, -454, 68, -40, -76, -76, -58, -40,
      -40, 428, -184, 68, 104, -76, 572, -58, -40, -40, 572, -184, 176, -40,
      -184, -76, 428, 50, -76, -184, 428, -13, -76, 68, -454, 50, 50, 50, 572,
      -454, 3641, 50, -454, 428, 554, -76, 572, -13, 50, -13, 5, -58, -58, -76,
      -40, 5, -58, -76, 50, -13, -13, -58, 104, -40, 572, 68, -76, 50, 572,
      -13, 554, -76, 428, -13, -58, 5, -76, -58, -40, -76, -40, 68, -454, -76,
      68, -13, 68, -76, 50, -454, 50, -13, -76, 5, -13, 68, -76, -76, 50, 428,
      -184, -76, 428, -58, -13, 5, -76, -13, 50, -76, -40, -58, 428, -40, -184,
      5, -58, -58, -76, -76, 428}, {-13, -58, -13, 50, 5, -76, 50, -76, -76,
      428, 428, -184, -58, -76, 5, -58, 428, -76, -58, 428, -76, -40, -184,
      -40, 554, 428, -13, 50, -76, 572, 428, -184, -76, 50, 428, -76, 50, -454,
      572, 50, 3641, -454, 50, 50, -76, -13, -454, 68, 572, -40, -40, -184,
      -184, 176, -76, -40, 68, -454, -76, 68, -58, -40, 104, 68, 572, -76, 428,
      -184, -58, -76, -40, -40, 5, -13, -13, -76, -76, 68, -76, 68, -40, -76,
      -454, 68, -58, -40, -58, -13, -76, 5, 50, -454, 50, -76, 68, -13, -13,
      554, 50, 572, 428, -76, -13, -13, -58, 5, 50, -76, -76, -40, 5, -13, -58,
      -58, -76, 68, 572, -40, 104, -58}, {-76, -40, 68, -454, -76, 68, 572,
      -40, -40, -184, -184, 176, 104, 68, -58, -40, -76, 572, -58, -76, 428,
      -184, -40, -40, 428, -184, -76, 50, 428, -76, 554, 428, -13, 50, -76,
      572, 572, 50, 50, -454, -454, 3641, -76, -13, 50, 50, 68, -454, 50, -76,
      -76, 428, 428, -184, -13, -58, -13, 50, 5, -76, 5, -58, -58, -76, -76,
      428, -76, -40, -58, 428, -40, -184, -58, -13, -58, -40, 5, -76, -40, -76,
      -76, 68, 68, -454, -13, -76, 5, -13, 68, -76, -13, 68, -76, 50, -454, 50,
      5, -13, -76, -40, -58, -58, -58, 5, -13, -13, -76, 50, 50, 572, -13, 554,
      -76, 428, -58, 104, -40, 572, 68, -76}, {-13, -13, -58, 5, 50, -76, -13,
      554, 50, 572, 428, -76, 5, -13, -76, -40, -58, -58, 572, -40, -76, 68,
      -58, 104, 50, 50, -76, -13, -454, 68, -76, 572, 428, 554, 50, -13, 428,
      -76, -184, 428, 50, -76, 3641, -454, -454, 50, 50, 572, -58, 428, -76,
      -40, -184, -40, 5, -76, -58, -13, 50, -13, 428, -184, -76, 50, 428, -76,
      428, -76, -76, -58, -58, 5, -58, -76, 428, -184, -40, -40, -76, 68, -40,
      -76, -454, 68, -184, 176, -40, 572, -184, -40, -76, 572, 68, 104, -40,
      -58, 50, -76, 50, -454, -13, 68, -454, 68, 68, -76, -76, -40, -76, 68,
      -13, 5, -76, -13, 5, -76, -13, -58, -40, -58}, {-58, 5, -13, -13, -76,
      50, 5, -13, -76, -40, -58, -58, -13, 554, 50, 572, 428, -76, -40, 572,
      -58, 104, -76, 68, -76, -13, 50, 50, 68, -454, 428, -76, -184, 428, 50,
      -76, -76, 572, 428, 554, 50, -13, -454, 3641, 50, 572, -454, 50, -58,
      -76, 428, -184, -40, -40, -76, 68, -40, -76, -454, 68, -184, 176, -40,
      572, -184, -40, -76, 572, 68, 104, -40, -58, -58, 428, -76, -40, -184,
      -40, 5, -76, -58, -13, 50, -13, 428, -184, -76, 50, 428, -76, 428, -76,
      -76, -58, -58, 5, -76, 50, -13, 68, 50, -454, 68, -454, -76, -40, 68,
      -76, 5, -76, -13, -58, -40, -58, -76, 68, -13, 5, -76, -13}, {5, -76,
      -58, -13, 50, -13, -58, 428, -76, -40, -184, -40, -76, 50, 428, -184,
      -76, 428, -76, -58, 428, -76, 5, -58, -76, 572, 428, 554, 50, -13, 50,
      50, -76, -13, -454, 68, -184, 428, 428, -76, -76, 50, -454, 50, 3641,
      -454, 572, 50, -13, 554, 50, 572, 428, -76, -13, -13, -58, 5, 50, -76,
      -76, -40, 5, -13, -58, -58, -76, 68, 572, -40, 104, -58, -40, 572, -184,
      176, -40, -184, -40, -76, -76, 68, 68, -454, 428, -184, -58, -76, -40,
      -40, -58, -40, 104, 68, 572, -76, -13, 5, -76, 68, -13, -76, 68, -76,
      -454, 68, -40, -76, 50, -454, 50, -76, 68, -13, -58, -40, -58, -13, -76,
      5}, {-76, 68, -40, -76, -454, 68, -58, -76, 428, -184, -40, -40, -40,
      572, -184, 176, -40, -184, 68, 104, -76, 572, -58, -40, 428, -76, -184,
      428, 50, -76, -76, -13, 50, 50, 68, -454, 428, 554, -76, 572, -13, 50,
      50, 572, -454, 3641, 50, -454, 5, -13, -76, -40, -58, -58, -58, 5, -13,
      -13, -76, 50, 50, 572, -13, 554, -76, 428, -58, 104, -40, 572, 68, -76,
      -76, 50, 428, -184, -76, 428, -58, -13, 5, -76, -13, 50, -76, -40, -58,
      428, -40, -184, 5, -58, -58, -76, -76, 428, -13, -58, 5, -76, -58, -40,
      -76, -40, 68, -454, -76, 68, -13, 68, -76, 50, -454, 50, -13, -76, 5,
      -13, 68, -76}, {-58, -13, 5, -76, -13, 50, -76, 50, 428, -184, -76, 428,
      -58, 428, -76, -40, -184, -40, -58, -76, 5, -58, 428, -76, 428, 554, -76,
      572, -13, 50, -184, 428, 428, -76, -76, 50, 50, 50, -76, -13, -454, 68,
      50, -454, 572, 50, 3641, -454, -40, 572, -184, 176, -40, -184, -40, -76,
      -76, 68, 68, -454, 428, -184, -58, -76, -40, -40, -58, -40, 104, 68, 572,
      -76, -13, 554, 50, 572, 428, -76, -13, -13, -58, 5, 50, -76, -76, -40, 5,
      -13, -58, -58, -76, 68, 572, -40, 104, -58, 5, -13, -13, -76, -76, 68,
      -76, 68, -40, -76, -454, 68, -58, -40, -58, -13, -76, 5, 50, -454, 50,
      -76, 68, -13}, {-40, -76, -76, 68, 68, -454, -40, 572, -184, 176, -40,
      -184, -58, -76, 428, -184, -40, -40, 104, 68, -58, -40, -76, 572, -184,
      428, 428, -76, -76, 50, 428, 554, -76, 572, -13, 50, -76, -13, 50, 50,
      68, -454, 572, 50, 50, -454, -454, 3641, -76, 50, 428, -184, -76, 428,
      -58, -13, 5, -76, -13, 50, -76, -40, -58, 428, -40, -184, 5, -58, -58,
      -76, -76, 428, 5, -13, -76, -40, -58, -58, -58, 5, -13, -13, -76, 50, 50,
      572, -13, 554, -76, 428, -58, 104, -40, 572, 68, -76, -58, -13, -58, -40,
      5, -76, -40, -76, -76, 68, 68, -454, -13, -76, 5, -13, 68, -76, -13, 68,
      -76, 50, -454, 50}, {50, -76, 50, -454, -13, 68, -454, 68, 68, -76, -76,
      -40, -76, 68, -13, 5, -76, -13, 5, -76, -13, -58, -40, -58, 572, -40,
      -76, 68, -58, 104, 50, -76, -13, -13, 5, -58, 428, -76, 554, -13, 572,
      50, -58, -58, -13, 5, -40, -76, 3641, -454, -454, 50, 50, 572, -454, 68,
      50, 50, -13, -76, 50, -13, 572, -76, 554, 428, 50, -76, -76, 428, 428,
      -184, 428, -76, -76, -58, -58, 5, -184, -40, 428, -58, -40, -76, 50, -13,
      -76, 5, -13, -58, 428, -76, -184, 428, 50, -76, -76, 572, 68, 104, -40,
      -58, -40, -40, -76, -58, -184, 428, -454, 68, 68, -76, -76, -40, -184,
      -40, 176, -184, 572, -40}, {-76, 50, -13, 68, 50, -454, 68, -454, -76,
      -40, 68, -76, 5, -76, -13, -58, -40, -58, -76, 68, -13, 5, -76, -13, -40,
      572, -58, 104, -76, 68, -76, 50, 5, -58, -13, -13, -58, -58, -13, 5, -40,
      -76, 428, -76, 554, -13, 572, 50, -454, 3641, 50, 572, -454, 50, 68,
      -454, -13, -76, 50, 50, 50, -76, -76, 428, 428, -184, 50, -13, 572, -76,
      554, 428, -76, 572, 68, 104, -40, -58, -40, -40, -76, -58, -184, 428,
      -454, 68, 68, -76, -76, -40, -184, -40, 176, -184, 572, -40, 428, -76,
      -76, -58, -58, 5, -184, -40, 428, -58, -40, -76, 50, -13, -76, 5, -13,
      -58, 428, -76, -184, 428, 50, -76}, {-13, 5, -76, 68, -13, -76, 68, -76,
      -454, 68, -40, -76, 50, -454, 50, -76, 68, -13, -58, -40, -58, -13, -76,
      5, -76, -58, 428, -76, 5, -58, 50, -13, -76, 5, -13, -58, -184, -40, 428,
      -58, -40, -76, -76, 428, 50, -76, -184, 428, -454, 50, 3641, -454, 572,
      50, 50, -13, 572, -76, 554, 428, -454, 68, 50, 50, -13, -76, -76, 50,
      428, -184, -76, 428, -76, 68, 572, -40, 104, -58, 428, -76, 554, -13,
      572, 50, 50, -76, -13, -13, 5, -58, -58, -58, -40, -76, -13, 5, -58, -40,
      104, 68, 572, -76, -40, -184, 572, -40, 176, -184, 68, -454, -76, -40,
      68, -76, -40, -40, -184, 428, -76, -58}, {-13, -58, 5, -76, -58, -40,
      -76, -40, 68, -454, -76, 68, -13, 68, -76, 50, -454, 50, -13, -76, 5,
      -13, 68, -76, 68, 104, -76, 572, -58, -40, -454, 68, 68, -76, -76, -40,
      -40, -40, -76, -58, -184, 428, -40, -184, 572, -40, 176, -184, 50, 572,
      -454, 3641, 50, -454, 50, -76, -76, 428, 428, -184, 68, -454, -13, -76,
      50, 50, -13, 50, 554, 428, 572, -76, -58, 104, -40, 572, 68, -76, -58,
      -58, -13, 5, -40, -76, -76, 50, 5, -58, -13, -13, -76, 428, 572, 50, 554,
      -13, 5, -58, -58, -76, -76, 428, -76, 428, 50, -76, -184, 428, -13, 50,
      -13, -58, -76, 5, -40, -184, -40, -76, 428, -58}, {5, -13, -13, -76, -76,
      68, -76, 68, -40, -76, -454, 68, -58, -40, -58, -13, -76, 5, 50, -454,
      50, -76, 68, -13, -58, -76, 5, -58, 428, -76, -13, 50, -13, -58, -76, 5,
      -76, 428, 50, -76, -184, 428, -184, -40, 428, -58, -40, -76, 50, -454,
      572, 50, 3641, -454, -13, 50, 554, 428, 572, -76, -76, 50, 428, -184,
      -76, 428, -454, 68, 50, 50, -13, -76, -58, -40, 104, 68, 572, -76, -40,
      -184, 572, -40, 176, -184, 68, -454, -76, -40, 68, -76, -40, -40, -184,
      428, -76, -58, -76, 68, 572, -40, 104, -58, 428, -76, 554, -13, 572, 50,
      50, -76, -13, -13, 5, -58, -58, -58, -40, -76, -13, 5}, {-58, -13, -58,
      -40, 5, -76, -40, -76, -76, 68, 68, -454, -13, -76, 5, -13, 68, -76, -13,
      68, -76, 50, -454, 50, 104, 68, -58, -40, -76, 572, 68, -454, -76, -40,
      68, -76, -40, -184, 572, -40, 176, -184, -40, -40, -76, -58, -184, 428,
      572, 50, 50, -454, -454, 3641, -76, 50, 428, -184, -76, 428, -13, 50,
      554, 428, 572, -76, 68, -454, -13, -76, 50, 50, 5, -58, -58, -76, -76,
      428, -76, 428, 50, -76, -184, 428, -13, 50, -13, -58, -76, 5, -40, -184,
      -40, -76, 428, -58, -58, 104, -40, 572, 68, -76, -58, -58, -13, 5, -40,
      -76, -76, 50, 5, -58, -13, -13, -76, 428, 572, 50, 554, -13}, {572, -40,
      -76, 68, -58, 104, 50, -76, -13, -13, 5, -58, 428, -76, 554, -13, 572,
      50, -58, -58, -13, 5, -40, -76, 50, -76, 50, -454, -13, 68, -454, 68, 68,
      -76, -76, -40, -76, 68, -13, 5, -76, -13, 5, -76, -13, -58, -40, -58,
      -454, 68, 50, 50, -13, -76, 3641, -454, -454, 50, 50, 572, 572, -76, 50,
      -13, 428, 554, -76, 428, 50, -76, -184, 428, -184, -40, 428, -58, -40,
      -76, 428, -76, -76, -58, -58, 5, -76, 5, 50, -13, -58, -13, -184, 428,
      428, -76, -76, 50, -40, -40, -76, -58, -184, 428, -76, 572, 68, 104, -40,
      -58, 68, -76, -454, 68, -40, -76, 176, -184, -184, -40, -40, 572}, {-40,
      572, -58, 104, -76, 68, -76, 50, 5, -58, -13, -13, -58, -58, -13, 5, -40,
      -76, 428, -76, 554, -13, 572, 50, -76, 50, -13, 68, 50, -454, 68, -454,
      -76, -40, 68, -76, 5, -76, -13, -58, -40, -58, -76, 68, -13, 5, -76, -13,
      68, -454, -13, -76, 50, 50, -454, 3641, 50, 572, -454, 50, -76, 428, 50,
      -76, -184, 428, 572, -76, 50, -13, 428, 554, -40, -40, -76, -58, -184,
      428, -76, 572, 68, 104, -40, -58, 68, -76, -454, 68, -40, -76, 176, -184,
      -184, -40, -40, 572, -184, -40, 428, -58, -40, -76, 428, -76, -76, -58,
      -58, 5, -76, 5, 50, -13, -58, -13, -184, 428, 428, -76, -76, 50}, {-76,
      -58, 428, -76, 5, -58, 50, -13, -76, 5, -13, -58, -184, -40, 428, -58,
      -40, -76, -76, 428, 50, -76, -184, 428, -13, 5, -76, 68, -13, -76, 68,
      -76, -454, 68, -40, -76, 50, -454, 50, -76, 68, -13, -58, -40, -58, -13,
      -76, 5, 50, -13, 572, -76, 554, 428, -454, 50, 3641, -454, 572, 50, 50,
      50, -454, 68, -76, -13, 428, -184, -76, 50, 428, -76, 428, -76, 554, -13,
      572, 50, -76, 68, 572, -40, 104, -58, -13, -13, 50, -76, -58, 5, -40,
      -76, -58, -58, 5, -13, -40, -184, 572, -40, 176, -184, -58, -40, 104, 68,
      572, -76, -76, -40, 68, -454, -76, 68, -184, 428, -40, -40, -58, -76},
      {68, 104, -76, 572, -58, -40, -454, 68, 68, -76, -76, -40, -40, -40, -76,
      -58, -184, 428, -40, -184, 572, -40, 176, -184, -13, -58, 5, -76, -58,
      -40, -76, -40, 68, -454, -76, 68, -13, 68, -76, 50, -454, 50, -13, -76,
      5, -13, 68, -76, 50, -76, -76, 428, 428, -184, 50, 572, -454, 3641, 50,
      -454, -13, -76, 68, -454, 50, 50, 554, 428, -13, 50, -76, 572, -58, -58,
      -13, 5, -40, -76, -58, 104, -40, 572, 68, -76, 5, -58, -76, 50, -13, -13,
      572, 50, -76, 428, -13, 554, -76, 428, 50, -76, -184, 428, 5, -58, -58,
      -76, -76, 428, -13, -58, -13, 50, 5, -76, -40, -76, -40, -184, -58, 428},
      {-58, -76, 5, -58, 428, -76, -13, 50, -13, -58, -76, 5, -76, 428, 50,
      -76, -184, 428, -184, -40, 428, -58, -40, -76, 5, -13, -13, -76, -76, 68,
      -76, 68, -40, -76, -454, 68, -58, -40, -58, -13, -76, 5, 50, -454, 50,
      -76, 68, -13, -13, 50, 554, 428, 572, -76, 50, -454, 572, 50, 3641, -454,
      428, -184, -76, 50, 428, -76, 50, 50, -454, 68, -76, -13, -40, -184, 572,
      -40, 176, -184, -58, -40, 104, 68, 572, -76, -76, -40, 68, -454, -76, 68,
      -184, 428, -40, -40, -58, -76, 428, -76, 554, -13, 572, 50, -76, 68, 572,
      -40, 104, -58, -13, -13, 50, -76, -58, 5, -40, -76, -58, -58, 5, -13},
      {104, 68, -58, -40, -76, 572, 68, -454, -76, -40, 68, -76, -40, -184,
      572, -40, 176, -184, -40, -40, -76, -58, -184, 428, -58, -13, -58, -40,
      5, -76, -40, -76, -76, 68, 68, -454, -13, -76, 5, -13, 68, -76, -13, 68,
      -76, 50, -454, 50, -76, 50, 428, -184, -76, 428, 572, 50, 50, -454, -454,
      3641, 554, 428, -13, 50, -76, 572, -13, -76, 68, -454, 50, 50, -76, 428,
      50, -76, -184, 428, 5, -58, -58, -76, -76, 428, -13, -58, -13, 50, 5,
      -76, -40, -76, -40, -184, -58, 428, -58, -58, -13, 5, -40, -76, -58, 104,
      -40, 572, 68, -76, 5, -58, -76, 50, -13, -13, 572, 50, -76, 428, -13,
      554}, {554, -13, 428, -76, 50, 572, -13, -13, 50, -76, -58, 5, -76, 68,
      572, -40, 104, -58, -76, -40, 5, -13, -58, -58, 428, -58, -184, -40, -76,
      -40, -76, 5, 50, -13, -58, -13, 428, -76, -76, -58, -58, 5, 428, -184,
      -76, 50, 428, -76, 50, 50, -454, 68, -76, -13, 572, -76, 50, -13, 428,
      554, 3641, -454, -454, 50, 50, 572, 428, -76, -184, 428, 50, -76, 50,
      -454, 50, -76, 68, -13, -76, 68, -13, 5, -76, -13, -454, 68, 68, -76,
      -76, -40, -76, 5, -40, -58, -13, -58, 428, -184, -58, -76, -40, -40,
      -184, 176, -40, 572, -184, -40, -76, 68, -40, -76, -454, 68, 572, -76,
      -40, -58, 68, 104}, {-13, 5, -58, -58, -76, -40, 5, -58, -76, 50, -13,
      -13, -58, 104, -40, 572, 68, -76, 50, 572, -13, 554, -76, 428, -76, -58,
      -40, -40, 428, -184, 68, -76, -454, 68, -40, -76, -76, 572, 68, 104, -40,
      -58, -184, 176, -40, 572, -184, -40, -13, -76, 68, -454, 50, 50, -76,
      428, 50, -76, -184, 428, -454, 3641, 50, 572, -454, 50, -76, 572, 428,
      554, 50, -13, -13, 68, -76, 50, -454, 50, 5, -76, -13, -58, -40, -58, 68,
      -454, -76, -40, 68, -76, 68, -76, -76, -13, -13, 5, -76, -40, -58, 428,
      -40, -184, 428, -184, -76, 50, 428, -76, 5, -76, -58, -13, 50, -13, -76,
      428, -58, 5, -76, -58}, {428, -58, -184, -40, -76, -40, -76, 5, 50, -13,
      -58, -13, 428, -76, -76, -58, -58, 5, 428, -184, -76, 50, 428, -76, 554,
      -13, 428, -76, 50, 572, -13, -13, 50, -76, -58, 5, -76, 68, 572, -40,
      104, -58, -76, -40, 5, -13, -58, -58, 572, -76, 50, -13, 428, 554, 50,
      50, -454, 68, -76, -13, -454, 50, 3641, -454, 572, 50, -184, 428, 428,
      -76, -76, 50, -76, 68, -13, 5, -76, -13, 50, -454, 50, -76, 68, -13, 68,
      -76, -454, 68, -40, -76, -40, -58, -76, 5, -58, -13, -184, 176, -40, 572,
      -184, -40, 428, -184, -58, -76, -40, -40, -40, -76, -76, 68, 68, -454,
      -40, -58, 572, -76, 104, 68}, {-76, -58, -40, -40, 428, -184, 68, -76,
      -454, 68, -40, -76, -76, 572, 68, 104, -40, -58, -184, 176, -40, 572,
      -184, -40, -13, 5, -58, -58, -76, -40, 5, -58, -76, 50, -13, -13, -58,
      104, -40, 572, 68, -76, 50, 572, -13, 554, -76, 428, -76, 428, 50, -76,
      -184, 428, -13, -76, 68, -454, 50, 50, 50, 572, -454, 3641, 50, -454,
      428, 554, -76, 572, -13, 50, 5, -76, -13, -58, -40, -58, -13, 68, -76,
      50, -454, 50, -76, -40, 68, -454, -76, 68, -76, -13, 68, -76, 5, -13,
      428, -184, -76, 50, 428, -76, -76, -40, -58, 428, -40, -184, -58, -13, 5,
      -76, -13, 50, -58, 5, -76, 428, -58, -76}, {50, -76, -76, 428, 428, -184,
      -13, -58, -13, 50, 5, -76, 5, -58, -58, -76, -76, 428, -76, -40, -58,
      428, -40, -184, 572, -40, -40, -184, -184, 176, -76, -40, 68, -454, -76,
      68, -58, -40, 104, 68, 572, -76, 428, -184, -58, -76, -40, -40, 554, 428,
      -13, 50, -76, 572, 428, -184, -76, 50, 428, -76, 50, -454, 572, 50, 3641,
      -454, 50, 50, -76, -13, -454, 68, -13, -76, 5, -13, 68, -76, -58, -40,
      -58, -13, -76, 5, -76, 68, -40, -76, -454, 68, -454, 50, 68, -13, 50,
      -76, 50, 572, -13, 554, -76, 428, -76, -40, 5, -13, -58, -58, -13, -13,
      -58, 5, 50, -76, 68, -76, 104, -58, 572, -40}, {572, -40, -40, -184,
      -184, 176, -76, -40, 68, -454, -76, 68, -58, -40, 104, 68, 572, -76, 428,
      -184, -58, -76, -40, -40, 50, -76, -76, 428, 428, -184, -13, -58, -13,
      50, 5, -76, 5, -58, -58, -76, -76, 428, -76, -40, -58, 428, -40, -184,
      428, -184, -76, 50, 428, -76, 554, 428, -13, 50, -76, 572, 572, 50, 50,
      -454, -454, 3641, -76, -13, 50, 50, 68, -454, -58, -40, -58, -13, -76, 5,
      -13, -76, 5, -13, 68, -76, -40, -76, -76, 68, 68, -454, 68, -13, -454,
      50, -76, 50, -76, -40, 5, -13, -58, -58, 50, 572, -13, 554, -76, 428,
      -58, 5, -13, -13, -76, 50, 104, -58, 68, -76, -40, 572}, {-13, 554, 50,
      572, 428, -76, -13, -13, -58, 5, 50, -76, -76, -40, 5, -13, -58, -58,
      -76, 68, 572, -40, 104, -58, -58, 428, -76, -40, -184, -40, 5, -76, -58,
      -13, 50, -13, 428, -184, -76, 50, 428, -76, 428, -76, -76, -58, -58, 5,
      50, 50, -76, -13, -454, 68, -76, 572, 428, 554, 50, -13, 428, -76, -184,
      428, 50, -76, 3641, -454, -454, 50, 50, 572, 428, -184, -58, -76, -40,
      -40, -184, 176, -40, 572, -184, -40, -76, 68, -40, -76, -454, 68, 572,
      -76, -40, -58, 68, 104, 50, -454, 50, -76, 68, -13, -76, 68, -13, 5, -76,
      -13, -454, 68, 68, -76, -76, -40, -76, 5, -40, -58, -13, -58}, {5, -13,
      -76, -40, -58, -58, -58, 5, -13, -13, -76, 50, 50, 572, -13, 554, -76,
      428, -58, 104, -40, 572, 68, -76, -58, -76, 428, -184, -40, -40, -76, 68,
      -40, -76, -454, 68, -184, 176, -40, 572, -184, -40, -76, 572, 68, 104,
      -40, -58, -76, -13, 50, 50, 68, -454, 428, -76, -184, 428, 50, -76, -76,
      572, 428, 554, 50, -13, -454, 3641, 50, 572, -454, 50, -76, -40, -58,
      428, -40, -184, 428, -184, -76, 50, 428, -76, 5, -76, -58, -13, 50, -13,
      -76, 428, -58, 5, -76, -58, -13, 68, -76, 50, -454, 50, 5, -76, -13, -58,
      -40, -58, 68, -454, -76, -40, 68, -76, 68, -76, -76, -13, -13, 5}, {-58,
      428, -76, -40, -184, -40, 5, -76, -58, -13, 50, -13, 428, -184, -76, 50,
      428, -76, 428, -76, -76, -58, -58, 5, -13, 554, 50, 572, 428, -76, -13,
      -13, -58, 5, 50, -76, -76, -40, 5, -13, -58, -58, -76, 68, 572, -40, 104,
      -58, -76, 572, 428, 554, 50, -13, 50, 50, -76, -13, -454, 68, -184, 428,
      428, -76, -76, 50, -454, 50, 3641, -454, 572, 50, -184, 176, -40, 572,
      -184, -40, 428, -184, -58, -76, -40, -40, -40, -76, -76, 68, 68, -454,
      -40, -58, 572, -76, 104, 68, -76, 68, -13, 5, -76, -13, 50, -454, 50,
      -76, 68, -13, 68, -76, -454, 68, -40, -76, -40, -58, -76, 5, -58, -13},
      {-58, -76, 428, -184, -40, -40, -76, 68, -40, -76, -454, 68, -184, 176,
      -40, 572, -184, -40, -76, 572, 68, 104, -40, -58, 5, -13, -76, -40, -58,
      -58, -58, 5, -13, -13, -76, 50, 50, 572, -13, 554, -76, 428, -58, 104,
      -40, 572, 68, -76, 428, -76, -184, 428, 50, -76, -76, -13, 50, 50, 68,
      -454, 428, 554, -76, 572, -13, 50, 50, 572, -454, 3641, 50, -454, 428,
      -184, -76, 50, 428, -76, -76, -40, -58, 428, -40, -184, -58, -13, 5, -76,
      -13, 50, -58, 5, -76, 428, -58, -76, 5, -76, -13, -58, -40, -58, -13, 68,
      -76, 50, -454, 50, -76, -40, 68, -454, -76, 68, -76, -13, 68, -76, 5,
      -13}, {-76, 50, 428, -184, -76, 428, -58, -13, 5, -76, -13, 50, -76, -40,
      -58, 428, -40, -184, 5, -58, -58, -76, -76, 428, -40, 572, -184, 176,
      -40, -184, -40, -76, -76, 68, 68, -454, 428, -184, -58, -76, -40, -40,
      -58, -40, 104, 68, 572, -76, 428, 554, -76, 572, -13, 50, -184, 428, 428,
      -76, -76, 50, 50, 50, -76, -13, -454, 68, 50, -454, 572, 50, 3641, -454,
      50, 572, -13, 554, -76, 428, -76, -40, 5, -13, -58, -58, -13, -13, -58,
      5, 50, -76, 68, -76, 104, -58, 572, -40, -13, -76, 5, -13, 68, -76, -58,
      -40, -58, -13, -76, 5, -76, 68, -40, -76, -454, 68, -454, 50, 68, -13,
      50, -76}, {-40, 572, -184, 176, -40, -184, -40, -76, -76, 68, 68, -454,
      428, -184, -58, -76, -40, -40, -58, -40, 104, 68, 572, -76, -76, 50, 428,
      -184, -76, 428, -58, -13, 5, -76, -13, 50, -76, -40, -58, 428, -40, -184,
      5, -58, -58, -76, -76, 428, -184, 428, 428, -76, -76, 50, 428, 554, -76,
      572, -13, 50, -76, -13, 50, 50, 68, -454, 572, 50, 50, -454, -454, 3641,
      -76, -40, 5, -13, -58, -58, 50, 572, -13, 554, -76, 428, -58, 5, -13,
      -13, -76, 50, 104, -58, 68, -76, -40, 572, -58, -40, -58, -13, -76, 5,
      -13, -76, 5, -13, 68, -76, -40, -76, -76, 68, 68, -454, 68, -13, -454,
      50, -76, 50}, {50, -454, 50, -76, 68, -13, -76, 68, -13, 5, -76, -13,
      -454, 68, 68, -76, -76, -40, -76, 5, -40, -58, -13, -58, -76, 68, 572,
      -40, 104, -58, 428, -76, 554, -13, 572, 50, 50, -76, -13, -13, 5, -58,
      -58, -58, -40, -76, -13, 5, 428, -76, -76, -58, -58, 5, -184, -40, 428,
      -58, -40, -76, 50, -13, -76, 5, -13, -58, 428, -76, -184, 428, 50, -76,
      3641, -454, -454, 50, 50, 572, -454, 68, 50, 50, -13, -76, 50, -13, 572,
      -76, 554, 428, 50, -76, -76, 428, 428, -184, 572, -76, -40, -58, 68, 104,
      -40, -40, -184, 428, -76, -58, -184, -40, 176, -184, 572, -40, -454, 68,
      68, -76, -76, -40}, {-13, 68, -76, 50, -454, 50, 5, -76, -13, -58, -40,
      -58, 68, -454, -76, -40, 68, -76, 68, -76, -76, -13, -13, 5, -58, 104,
      -40, 572, 68, -76, -58, -58, -13, 5, -40, -76, -76, 50, 5, -58, -13, -13,
      -76, 428, 572, 50, 554, -13, -76, 572, 68, 104, -40, -58, -40, -40, -76,
      -58, -184, 428, -454, 68, 68, -76, -76, -40, -184, -40, 176, -184, 572,
      -40, -454, 3641, 50, 572, -454, 50, 68, -454, -13, -76, 50, 50, 50, -76,
      -76, 428, 428, -184, 50, -13, 572, -76, 554, 428, -76, 428, -58, 5, -76,
      -58, -40, -184, -40, -76, 428, -58, 428, -76, -184, 428, 50, -76, 50,
      -13, -76, 5, -13, -58}, {-76, 68, -13, 5, -76, -13, 50, -454, 50, -76,
      68, -13, 68, -76, -454, 68, -40, -76, -40, -58, -76, 5, -58, -13, 428,
      -76, -76, -58, -58, 5, -184, -40, 428, -58, -40, -76, 50, -13, -76, 5,
      -13, -58, 428, -76, -184, 428, 50, -76, -76, 68, 572, -40, 104, -58, 428,
      -76, 554, -13, 572, 50, 50, -76, -13, -13, 5, -58, -58, -58, -40, -76,
      -13, 5, -454, 50, 3641, -454, 572, 50, 50, -13, 572, -76, 554, 428, -454,
      68, 50, 50, -13, -76, -76, 50, 428, -184, -76, 428, -40, -58, 572, -76,
      104, 68, -184, -40, 176, -184, 572, -40, -40, -40, -184, 428, -76, -58,
      68, -454, -76, -40, 68, -76}, {5, -76, -13, -58, -40, -58, -13, 68, -76,
      50, -454, 50, -76, -40, 68, -454, -76, 68, -76, -13, 68, -76, 5, -13,
      -76, 572, 68, 104, -40, -58, -40, -40, -76, -58, -184, 428, -454, 68, 68,
      -76, -76, -40, -184, -40, 176, -184, 572, -40, -58, 104, -40, 572, 68,
      -76, -58, -58, -13, 5, -40, -76, -76, 50, 5, -58, -13, -13, -76, 428,
      572, 50, 554, -13, 50, 572, -454, 3641, 50, -454, 50, -76, -76, 428, 428,
      -184, 68, -454, -13, -76, 50, 50, -13, 50, 554, 428, 572, -76, -58, 5,
      -76, 428, -58, -76, 428, -76, -184, 428, 50, -76, -40, -184, -40, -76,
      428, -58, -13, 50, -13, -58, -76, 5}, {-13, -76, 5, -13, 68, -76, -58,
      -40, -58, -13, -76, 5, -76, 68, -40, -76, -454, 68, -454, 50, 68, -13,
      50, -76, 5, -58, -58, -76, -76, 428, -76, 428, 50, -76, -184, 428, -13,
      50, -13, -58, -76, 5, -40, -184, -40, -76, 428, -58, -58, -40, 104, 68,
      572, -76, -40, -184, 572, -40, 176, -184, 68, -454, -76, -40, 68, -76,
      -40, -40, -184, 428, -76, -58, 50, -454, 572, 50, 3641, -454, -13, 50,
      554, 428, 572, -76, -76, 50, 428, -184, -76, 428, -454, 68, 50, 50, -13,
      -76, 68, -76, 104, -58, 572, -40, -76, 428, 572, 50, 554, -13, -58, -58,
      -40, -76, -13, 5, 50, -76, -13, -13, 5, -58}, {-58, -40, -58, -13, -76,
      5, -13, -76, 5, -13, 68, -76, -40, -76, -76, 68, 68, -454, 68, -13, -454,
      50, -76, 50, -58, -40, 104, 68, 572, -76, -40, -184, 572, -40, 176, -184,
      68, -454, -76, -40, 68, -76, -40, -40, -184, 428, -76, -58, 5, -58, -58,
      -76, -76, 428, -76, 428, 50, -76, -184, 428, -13, 50, -13, -58, -76, 5,
      -40, -184, -40, -76, 428, -58, 572, 50, 50, -454, -454, 3641, -76, 50,
      428, -184, -76, 428, -13, 50, 554, 428, 572, -76, 68, -454, -13, -76, 50,
      50, 104, -58, 68, -76, -40, 572, -58, -58, -40, -76, -13, 5, -76, 428,
      572, 50, 554, -13, -76, 50, 5, -58, -13, -13}, {-76, 68, 572, -40, 104,
      -58, 428, -76, 554, -13, 572, 50, 50, -76, -13, -13, 5, -58, -58, -58,
      -40, -76, -13, 5, 50, -454, 50, -76, 68, -13, -76, 68, -13, 5, -76, -13,
      -454, 68, 68, -76, -76, -40, -76, 5, -40, -58, -13, -58, -184, -40, 428,
      -58, -40, -76, 428, -76, -76, -58, -58, 5, -76, 5, 50, -13, -58, -13,
      -184, 428, 428, -76, -76, 50, -454, 68, 50, 50, -13, -76, 3641, -454,
      -454, 50, 50, 572, 572, -76, 50, -13, 428, 554, -76, 428, 50, -76, -184,
      428, -40, -40, -184, 428, -76, -58, 572, -76, -40, -58, 68, 104, 176,
      -184, -184, -40, -40, 572, 68, -76, -454, 68, -40, -76}, {-58, 104, -40,
      572, 68, -76, -58, -58, -13, 5, -40, -76, -76, 50, 5, -58, -13, -13, -76,
      428, 572, 50, 554, -13, -13, 68, -76, 50, -454, 50, 5, -76, -13, -58,
      -40, -58, 68, -454, -76, -40, 68, -76, 68, -76, -76, -13, -13, 5, -40,
      -40, -76, -58, -184, 428, -76, 572, 68, 104, -40, -58, 68, -76, -454, 68,
      -40, -76, 176, -184, -184, -40, -40, 572, 68, -454, -13, -76, 50, 50,
      -454, 3641, 50, 572, -454, 50, -76, 428, 50, -76, -184, 428, 572, -76,
      50, -13, 428, 554, -40, -184, -40, -76, 428, -58, -76, 428, -58, 5, -76,
      -58, -184, 428, 428, -76, -76, 50, -76, 5, 50, -13, -58, -13}, {428, -76,
      -76, -58, -58, 5, -184, -40, 428, -58, -40, -76, 50, -13, -76, 5, -13,
      -58, 428, -76, -184, 428, 50, -76, -76, 68, -13, 5, -76, -13, 50, -454,
      50, -76, 68, -13, 68, -76, -454, 68, -40, -76, -40, -58, -76, 5, -58,
      -13, 428, -76, 554, -13, 572, 50, -76, 68, 572, -40, 104, -58, -13, -13,
      50, -76, -58, 5, -40, -76, -58, -58, 5, -13, 50, -13, 572, -76, 554, 428,
      -454, 50, 3641, -454, 572, 50, 50, 50, -454, 68, -76, -13, 428, -184,
      -76, 50, 428, -76, -184, -40, 176, -184, 572, -40, -40, -58, 572, -76,
      104, 68, -184, 428, -40, -40, -58, -76, -76, -40, 68, -454, -76, 68},
      {-76, 572, 68, 104, -40, -58, -40, -40, -76, -58, -184, 428, -454, 68,
      68, -76, -76, -40, -184, -40, 176, -184, 572, -40, 5, -76, -13, -58, -40,
      -58, -13, 68, -76, 50, -454, 50, -76, -40, 68, -454, -76, 68, -76, -13,
      68, -76, 5, -13, -58, -58, -13, 5, -40, -76, -58, 104, -40, 572, 68, -76,
      5, -58, -76, 50, -13, -13, 572, 50, -76, 428, -13, 554, 50, -76, -76,
      428, 428, -184, 50, 572, -454, 3641, 50, -454, -13, -76, 68, -454, 50,
      50, 554, 428, -13, 50, -76, 572, 428, -76, -184, 428, 50, -76, -58, 5,
      -76, 428, -58, -76, -40, -76, -40, -184, -58, 428, -13, -58, -13, 50, 5,
      -76}, {5, -58, -58, -76, -76, 428, -76, 428, 50, -76, -184, 428, -13, 50,
      -13, -58, -76, 5, -40, -184, -40, -76, 428, -58, -13, -76, 5, -13, 68,
      -76, -58, -40, -58, -13, -76, 5, -76, 68, -40, -76, -454, 68, -454, 50,
      68, -13, 50, -76, -40, -184, 572, -40, 176, -184, -58, -40, 104, 68, 572,
      -76, -76, -40, 68, -454, -76, 68, -184, 428, -40, -40, -58, -76, -13, 50,
      554, 428, 572, -76, 50, -454, 572, 50, 3641, -454, 428, -184, -76, 50,
      428, -76, 50, 50, -454, 68, -76, -13, -76, 428, 572, 50, 554, -13, 68,
      -76, 104, -58, 572, -40, -40, -76, -58, -58, 5, -13, -13, -13, 50, -76,
      -58, 5}, {-58, -40, 104, 68, 572, -76, -40, -184, 572, -40, 176, -184,
      68, -454, -76, -40, 68, -76, -40, -40, -184, 428, -76, -58, -58, -40,
      -58, -13, -76, 5, -13, -76, 5, -13, 68, -76, -40, -76, -76, 68, 68, -454,
      68, -13, -454, 50, -76, 50, -76, 428, 50, -76, -184, 428, 5, -58, -58,
      -76, -76, 428, -13, -58, -13, 50, 5, -76, -40, -76, -40, -184, -58, 428,
      -76, 50, 428, -184, -76, 428, 572, 50, 50, -454, -454, 3641, 554, 428,
      -13, 50, -76, 572, -13, -76, 68, -454, 50, 50, -58, -58, -40, -76, -13,
      5, 104, -58, 68, -76, -40, 572, 572, 50, -76, 428, -13, 554, 5, -58, -76,
      50, -13, -13}, {428, -76, 554, -13, 572, 50, -76, 68, 572, -40, 104, -58,
      -13, -13, 50, -76, -58, 5, -40, -76, -58, -58, 5, -13, -184, -40, 428,
      -58, -40, -76, 428, -76, -76, -58, -58, 5, -76, 5, 50, -13, -58, -13,
      -184, 428, 428, -76, -76, 50, 50, -454, 50, -76, 68, -13, -76, 68, -13,
      5, -76, -13, -454, 68, 68, -76, -76, -40, -76, 5, -40, -58, -13, -58, 50,
      50, -454, 68, -76, -13, 572, -76, 50, -13, 428, 554, 3641, -454, -454,
      50, 50, 572, 428, -76, -184, 428, 50, -76, -184, 428, -40, -40, -58, -76,
      176, -184, -184, -40, -40, 572, 572, -76, -40, -58, 68, 104, -76, 68,
      -40, -76, -454, 68}, {-58, -58, -13, 5, -40, -76, -58, 104, -40, 572, 68,
      -76, 5, -58, -76, 50, -13, -13, 572, 50, -76, 428, -13, 554, -40, -40,
      -76, -58, -184, 428, -76, 572, 68, 104, -40, -58, 68, -76, -454, 68, -40,
      -76, 176, -184, -184, -40, -40, 572, -13, 68, -76, 50, -454, 50, 5, -76,
      -13, -58, -40, -58, 68, -454, -76, -40, 68, -76, 68, -76, -76, -13, -13,
      5, -13, -76, 68, -454, 50, 50, -76, 428, 50, -76, -184, 428, -454, 3641,
      50, 572, -454, 50, -76, 572, 428, 554, 50, -13, -40, -76, -40, -184, -58,
      428, -184, 428, 428, -76, -76, 50, -76, 428, -58, 5, -76, -58, 5, -76,
      -58, -13, 50, -13}, {-184, -40, 428, -58, -40, -76, 428, -76, -76, -58,
      -58, 5, -76, 5, 50, -13, -58, -13, -184, 428, 428, -76, -76, 50, 428,
      -76, 554, -13, 572, 50, -76, 68, 572, -40, 104, -58, -13, -13, 50, -76,
      -58, 5, -40, -76, -58, -58, 5, -13, -76, 68, -13, 5, -76, -13, 50, -454,
      50, -76, 68, -13, 68, -76, -454, 68, -40, -76, -40, -58, -76, 5, -58,
      -13, 572, -76, 50, -13, 428, 554, 50, 50, -454, 68, -76, -13, -454, 50,
      3641, -454, 572, 50, -184, 428, 428, -76, -76, 50, 176, -184, -184, -40,
      -40, 572, -184, 428, -40, -40, -58, -76, -40, -58, 572, -76, 104, 68,
      -40, -76, -76, 68, 68, -454}, {-40, -40, -76, -58, -184, 428, -76, 572,
      68, 104, -40, -58, 68, -76, -454, 68, -40, -76, 176, -184, -184, -40,
      -40, 572, -58, -58, -13, 5, -40, -76, -58, 104, -40, 572, 68, -76, 5,
      -58, -76, 50, -13, -13, 572, 50, -76, 428, -13, 554, 5, -76, -13, -58,
      -40, -58, -13, 68, -76, 50, -454, 50, -76, -40, 68, -454, -76, 68, -76,
      -13, 68, -76, 5, -13, -76, 428, 50, -76, -184, 428, -13, -76, 68, -454,
      50, 50, 50, 572, -454, 3641, 50, -454, 428, 554, -76, 572, -13, 50, -184,
      428, 428, -76, -76, 50, -40, -76, -40, -184, -58, 428, -58, 5, -76, 428,
      -58, -76, -58, -13, 5, -76, -13, 50}, {-76, 428, 50, -76, -184, 428, 5,
      -58, -58, -76, -76, 428, -13, -58, -13, 50, 5, -76, -40, -76, -40, -184,
      -58, 428, -40, -184, 572, -40, 176, -184, -58, -40, 104, 68, 572, -76,
      -76, -40, 68, -454, -76, 68, -184, 428, -40, -40, -58, -76, -13, -76, 5,
      -13, 68, -76, -58, -40, -58, -13, -76, 5, -76, 68, -40, -76, -454, 68,
      -454, 50, 68, -13, 50, -76, 554, 428, -13, 50, -76, 572, 428, -184, -76,
      50, 428, -76, 50, -454, 572, 50, 3641, -454, 50, 50, -76, -13, -454, 68,
      572, 50, -76, 428, -13, 554, -40, -76, -58, -58, 5, -13, 68, -76, 104,
      -58, 572, -40, -13, -13, -58, 5, 50, -76}, {-40, -184, 572, -40, 176,
      -184, -58, -40, 104, 68, 572, -76, -76, -40, 68, -454, -76, 68, -184,
      428, -40, -40, -58, -76, -76, 428, 50, -76, -184, 428, 5, -58, -58, -76,
      -76, 428, -13, -58, -13, 50, 5, -76, -40, -76, -40, -184, -58, 428, -58,
      -40, -58, -13, -76, 5, -13, -76, 5, -13, 68, -76, -40, -76, -76, 68, 68,
      -454, 68, -13, -454, 50, -76, 50, 428, -184, -76, 50, 428, -76, 554, 428,
      -13, 50, -76, 572, 572, 50, 50, -454, -454, 3641, -76, -13, 50, 50, 68,
      -454, -40, -76, -58, -58, 5, -13, 572, 50, -76, 428, -13, 554, 104, -58,
      68, -76, -40, 572, -58, 5, -13, -13, -76, 50}, {50, 572, -13, 554, -76,
      428, -76, -40, 5, -13, -58, -58, -13, -13, -58, 5, 50, -76, 68, -76, 104,
      -58, 572, -40, -76, -40, -58, 428, -40, -184, 428, -184, -76, 50, 428,
      -76, 5, -76, -58, -13, 50, -13, -76, 428, -58, 5, -76, -58, 428, -184,
      -58, -76, -40, -40, -184, 176, -40, 572, -184, -40, -76, 68, -40, -76,
      -454, 68, 572, -76, -40, -58, 68, 104, 50, 50, -76, -13, -454, 68, -76,
      572, 428, 554, 50, -13, 428, -76, -184, 428, 50, -76, 3641, -454, -454,
      50, 50, 572, -454, 50, 68, -13, 50, -76, 68, -76, -76, -13, -13, 5, -76,
      5, -40, -58, -13, -58, -454, 68, 68, -76, -76, -40}, {-76, -40, 5, -13,
      -58, -58, 50, 572, -13, 554, -76, 428, -58, 5, -13, -13, -76, 50, 104,
      -58, 68, -76, -40, 572, 428, -184, -58, -76, -40, -40, -184, 176, -40,
      572, -184, -40, -76, 68, -40, -76, -454, 68, 572, -76, -40, -58, 68, 104,
      -76, -40, -58, 428, -40, -184, 428, -184, -76, 50, 428, -76, 5, -76, -58,
      -13, 50, -13, -76, 428, -58, 5, -76, -58, -76, -13, 50, 50, 68, -454,
      428, -76, -184, 428, 50, -76, -76, 572, 428, 554, 50, -13, -454, 3641,
      50, 572, -454, 50, 68, -13, -454, 50, -76, 50, -76, 5, -40, -58, -13,
      -58, 68, -76, -76, -13, -13, 5, 68, -454, -76, -40, 68, -76}, {-76, -40,
      -58, 428, -40, -184, 428, -184, -76, 50, 428, -76, 5, -76, -58, -13, 50,
      -13, -76, 428, -58, 5, -76, -58, 50, 572, -13, 554, -76, 428, -76, -40,
      5, -13, -58, -58, -13, -13, -58, 5, 50, -76, 68, -76, 104, -58, 572, -40,
      -184, 176, -40, 572, -184, -40, 428, -184, -58, -76, -40, -40, -40, -76,
      -76, 68, 68, -454, -40, -58, 572, -76, 104, 68, -76, 572, 428, 554, 50,
      -13, 50, 50, -76, -13, -454, 68, -184, 428, 428, -76, -76, 50, -454, 50,
      3641, -454, 572, 50, 68, -76, -76, -13, -13, 5, -454, 50, 68, -13, 50,
      -76, -40, -58, -76, 5, -58, -13, 68, -76, -454, 68, -40, -76}, {428,
      -184, -58, -76, -40, -40, -184, 176, -40, 572, -184, -40, -76, 68, -40,
      -76, -454, 68, 572, -76, -40, -58, 68, 104, -76, -40, 5, -13, -58, -58,
      50, 572, -13, 554, -76, 428, -58, 5, -13, -13, -76, 50, 104, -58, 68,
      -76, -40, 572, 428, -184, -76, 50, 428, -76, -76, -40, -58, 428, -40,
      -184, -58, -13, 5, -76, -13, 50, -58, 5, -76, 428, -58, -76, 428, -76,
      -184, 428, 50, -76, -76, -13, 50, 50, 68, -454, 428, 554, -76, 572, -13,
      50, 50, 572, -454, 3641, 50, -454, -76, 5, -40, -58, -13, -58, 68, -13,
      -454, 50, -76, 50, -76, -13, 68, -76, 5, -13, -76, -40, 68, -454, -76,
      68}, {428, -184, -76, 50, 428, -76, -76, -40, -58, 428, -40, -184, -58,
      -13, 5, -76, -13, 50, -58, 5, -76, 428, -58, -76, -184, 176, -40, 572,
      -184, -40, 428, -184, -58, -76, -40, -40, -40, -76, -76, 68, 68, -454,
      -40, -58, 572, -76, 104, 68, 50, 572, -13, 554, -76, 428, -76, -40, 5,
      -13, -58, -58, -13, -13, -58, 5, 50, -76, 68, -76, 104, -58, 572, -40,
      428, 554, -76, 572, -13, 50, -184, 428, 428, -76, -76, 50, 50, 50, -76,
      -13, -454, 68, 50, -454, 572, 50, 3641, -454, -76, -13, 68, -76, 5, -13,
      -40, -58, -76, 5, -58, -13, -454, 50, 68, -13, 50, -76, -76, 68, -40,
      -76, -454, 68}, {-184, 176, -40, 572, -184, -40, 428, -184, -58, -76,
      -40, -40, -40, -76, -76, 68, 68, -454, -40, -58, 572, -76, 104, 68, 428,
      -184, -76, 50, 428, -76, -76, -40, -58, 428, -40, -184, -58, -13, 5, -76,
      -13, 50, -58, 5, -76, 428, -58, -76, -76, -40, 5, -13, -58, -58, 50, 572,
      -13, 554, -76, 428, -58, 5, -13, -13, -76, 50, 104, -58, 68, -76, -40,
      572, -184, 428, 428, -76, -76, 50, 428, 554, -76, 572, -13, 50, -76, -13,
      50, 50, 68, -454, 572, 50, 50, -454, -454, 3641, -40, -58, -76, 5, -58,
      -13, -76, -13, 68, -76, 5, -13, 68, -13, -454, 50, -76, 50, -40, -76,
      -76, 68, 68, -454}, {-454, 50, 68, -13, 50, -76, 68, -76, -76, -13, -13,
      5, -76, 5, -40, -58, -13, -58, -454, 68, 68, -76, -76, -40, 68, -76, 104,
      -58, 572, -40, -76, 428, 572, 50, 554, -13, -58, -58, -40, -76, -13, 5,
      50, -76, -13, -13, 5, -58, -76, 428, -58, 5, -76, -58, -40, -184, -40,
      -76, 428, -58, 428, -76, -184, 428, 50, -76, 50, -13, -76, 5, -13, -58,
      572, -76, -40, -58, 68, 104, -40, -40, -184, 428, -76, -58, -184, -40,
      176, -184, 572, -40, -454, 68, 68, -76, -76, -40, 3641, -454, -454, 50,
      50, 572, -454, 68, 50, 50, -13, -76, 50, -13, 572, -76, 554, 428, 50,
      -76, -76, 428, 428, -184}, {68, -13, -454, 50, -76, 50, -76, 5, -40, -58,
      -13, -58, 68, -76, -76, -13, -13, 5, 68, -454, -76, -40, 68, -76, 104,
      -58, 68, -76, -40, 572, -58, -58, -40, -76, -13, 5, -76, 428, 572, 50,
      554, -13, -76, 50, 5, -58, -13, -13, 572, -76, -40, -58, 68, 104, -40,
      -40, -184, 428, -76, -58, -184, -40, 176, -184, 572, -40, -454, 68, 68,
      -76, -76, -40, -76, 428, -58, 5, -76, -58, -40, -184, -40, -76, 428, -58,
      428, -76, -184, 428, 50, -76, 50, -13, -76, 5, -13, -58, -454, 3641, 50,
      572, -454, 50, 68, -454, -13, -76, 50, 50, 50, -76, -76, 428, 428, -184,
      50, -13, 572, -76, 554, 428}, {68, -76, -76, -13, -13, 5, -454, 50, 68,
      -13, 50, -76, -40, -58, -76, 5, -58, -13, 68, -76, -454, 68, -40, -76,
      -76, 428, -58, 5, -76, -58, -40, -184, -40, -76, 428, -58, 428, -76,
      -184, 428, 50, -76, 50, -13, -76, 5, -13, -58, 68, -76, 104, -58, 572,
      -40, -76, 428, 572, 50, 554, -13, -58, -58, -40, -76, -13, 5, 50, -76,
      -13, -13, 5, -58, -40, -58, 572, -76, 104, 68, -184, -40, 176, -184, 572,
      -40, -40, -40, -184, 428, -76, -58, 68, -454, -76, -40, 68, -76, -454,
      50, 3641, -454, 572, 50, 50, -13, 572, -76, 554, 428, -454, 68, 50, 50,
      -13, -76, -76, 50, 428, -184, -76, 428}, {-76, 5, -40, -58, -13, -58, 68,
      -13, -454, 50, -76, 50, -76, -13, 68, -76, 5, -13, -76, -40, 68, -454,
      -76, 68, 572, -76, -40, -58, 68, 104, -40, -40, -184, 428, -76, -58,
      -184, -40, 176, -184, 572, -40, -454, 68, 68, -76, -76, -40, 104, -58,
      68, -76, -40, 572, -58, -58, -40, -76, -13, 5, -76, 428, 572, 50, 554,
      -13, -76, 50, 5, -58, -13, -13, -58, 5, -76, 428, -58, -76, 428, -76,
      -184, 428, 50, -76, -40, -184, -40, -76, 428, -58, -13, 50, -13, -58,
      -76, 5, 50, 572, -454, 3641, 50, -454, 50, -76, -76, 428, 428, -184, 68,
      -454, -13, -76, 50, 50, -13, 50, 554, 428, 572, -76}, {-76, -13, 68, -76,
      5, -13, -40, -58, -76, 5, -58, -13, -454, 50, 68, -13, 50, -76, -76, 68,
      -40, -76, -454, 68, -58, 5, -76, 428, -58, -76, 428, -76, -184, 428, 50,
      -76, -40, -184, -40, -76, 428, -58, -13, 50, -13, -58, -76, 5, -40, -58,
      572, -76, 104, 68, -184, -40, 176, -184, 572, -40, -40, -40, -184, 428,
      -76, -58, 68, -454, -76, -40, 68, -76, 68, -76, 104, -58, 572, -40, -76,
      428, 572, 50, 554, -13, -58, -58, -40, -76, -13, 5, 50, -76, -13, -13, 5,
      -58, 50, -454, 572, 50, 3641, -454, -13, 50, 554, 428, 572, -76, -76, 50,
      428, -184, -76, 428, -454, 68, 50, 50, -13, -76}, {-40, -58, -76, 5, -58,
      -13, -76, -13, 68, -76, 5, -13, 68, -13, -454, 50, -76, 50, -40, -76,
      -76, 68, 68, -454, -40, -58, 572, -76, 104, 68, -184, -40, 176, -184,
      572, -40, -40, -40, -184, 428, -76, -58, 68, -454, -76, -40, 68, -76,
      -58, 5, -76, 428, -58, -76, 428, -76, -184, 428, 50, -76, -40, -184, -40,
      -76, 428, -58, -13, 50, -13, -58, -76, 5, 104, -58, 68, -76, -40, 572,
      -58, -58, -40, -76, -13, 5, -76, 428, 572, 50, 554, -13, -76, 50, 5, -58,
      -13, -13, 572, 50, 50, -454, -454, 3641, -76, 50, 428, -184, -76, 428,
      -13, 50, 554, 428, 572, -76, 68, -454, -13, -76, 50, 50}, {68, -76, 104,
      -58, 572, -40, -76, 428, 572, 50, 554, -13, -58, -58, -40, -76, -13, 5,
      50, -76, -13, -13, 5, -58, -454, 50, 68, -13, 50, -76, 68, -76, -76, -13,
      -13, 5, -76, 5, -40, -58, -13, -58, -454, 68, 68, -76, -76, -40, -40,
      -184, -40, -76, 428, -58, -76, 428, -58, 5, -76, -58, -184, 428, 428,
      -76, -76, 50, -76, 5, 50, -13, -58, -13, -40, -40, -184, 428, -76, -58,
      572, -76, -40, -58, 68, 104, 176, -184, -184, -40, -40, 572, 68, -76,
      -454, 68, -40, -76, -454, 68, 50, 50, -13, -76, 3641, -454, -454, 50, 50,
      572, 572, -76, 50, -13, 428, 554, -76, 428, 50, -76, -184, 428}, {104,
      -58, 68, -76, -40, 572, -58, -58, -40, -76, -13, 5, -76, 428, 572, 50,
      554, -13, -76, 50, 5, -58, -13, -13, 68, -13, -454, 50, -76, 50, -76, 5,
      -40, -58, -13, -58, 68, -76, -76, -13, -13, 5, 68, -454, -76, -40, 68,
      -76, -40, -40, -184, 428, -76, -58, 572, -76, -40, -58, 68, 104, 176,
      -184, -184, -40, -40, 572, 68, -76, -454, 68, -40, -76, -40, -184, -40,
      -76, 428, -58, -76, 428, -58, 5, -76, -58, -184, 428, 428, -76, -76, 50,
      -76, 5, 50, -13, -58, -13, 68, -454, -13, -76, 50, 50, -454, 3641, 50,
      572, -454, 50, -76, 428, 50, -76, -184, 428, 572, -76, 50, -13, 428,
      554}, {-76, 428, -58, 5, -76, -58, -40, -184, -40, -76, 428, -58, 428,
      -76, -184, 428, 50, -76, 50, -13, -76, 5, -13, -58, 68, -76, -76, -13,
      -13, 5, -454, 50, 68, -13, 50, -76, -40, -58, -76, 5, -58, -13, 68, -76,
      -454, 68, -40, -76, -76, 428, 572, 50, 554, -13, 68, -76, 104, -58, 572,
      -40, -40, -76, -58, -58, 5, -13, -13, -13, 50, -76, -58, 5, -184, -40,
      176, -184, 572, -40, -40, -58, 572, -76, 104, 68, -184, 428, -40, -40,
      -58, -76, -76, -40, 68, -454, -76, 68, 50, -13, 572, -76, 554, 428, -454,
      50, 3641, -454, 572, 50, 50, 50, -454, 68, -76, -13, 428, -184, -76, 50,
      428, -76}, {572, -76, -40, -58, 68, 104, -40, -40, -184, 428, -76, -58,
      -184, -40, 176, -184, 572, -40, -454, 68, 68, -76, -76, -40, -76, 5, -40,
      -58, -13, -58, 68, -13, -454, 50, -76, 50, -76, -13, 68, -76, 5, -13,
      -76, -40, 68, -454, -76, 68, -58, -58, -40, -76, -13, 5, 104, -58, 68,
      -76, -40, 572, 572, 50, -76, 428, -13, 554, 5, -58, -76, 50, -13, -13,
      428, -76, -184, 428, 50, -76, -58, 5, -76, 428, -58, -76, -40, -76, -40,
      -184, -58, 428, -13, -58, -13, 50, 5, -76, 50, -76, -76, 428, 428, -184,
      50, 572, -454, 3641, 50, -454, -13, -76, 68, -454, 50, 50, 554, 428, -13,
      50, -76, 572}, {-58, 5, -76, 428, -58, -76, 428, -76, -184, 428, 50, -76,
      -40, -184, -40, -76, 428, -58, -13, 50, -13, -58, -76, 5, -76, -13, 68,
      -76, 5, -13, -40, -58, -76, 5, -58, -13, -454, 50, 68, -13, 50, -76, -76,
      68, -40, -76, -454, 68, -184, -40, 176, -184, 572, -40, -40, -58, 572,
      -76, 104, 68, -184, 428, -40, -40, -58, -76, -76, -40, 68, -454, -76, 68,
      -76, 428, 572, 50, 554, -13, 68, -76, 104, -58, 572, -40, -40, -76, -58,
      -58, 5, -13, -13, -13, 50, -76, -58, 5, -13, 50, 554, 428, 572, -76, 50,
      -454, 572, 50, 3641, -454, 428, -184, -76, 50, 428, -76, 50, 50, -454,
      68, -76, -13}, {-40, -58, 572, -76, 104, 68, -184, -40, 176, -184, 572,
      -40, -40, -40, -184, 428, -76, -58, 68, -454, -76, -40, 68, -76, -40,
      -58, -76, 5, -58, -13, -76, -13, 68, -76, 5, -13, 68, -13, -454, 50, -76,
      50, -40, -76, -76, 68, 68, -454, 428, -76, -184, 428, 50, -76, -58, 5,
      -76, 428, -58, -76, -40, -76, -40, -184, -58, 428, -13, -58, -13, 50, 5,
      -76, -58, -58, -40, -76, -13, 5, 104, -58, 68, -76, -40, 572, 572, 50,
      -76, 428, -13, 554, 5, -58, -76, 50, -13, -13, -76, 50, 428, -184, -76,
      428, 572, 50, 50, -454, -454, 3641, 554, 428, -13, 50, -76, 572, -13,
      -76, 68, -454, 50, 50}, {-76, 428, 572, 50, 554, -13, 68, -76, 104, -58,
      572, -40, -40, -76, -58, -58, 5, -13, -13, -13, 50, -76, -58, 5, -40,
      -184, -40, -76, 428, -58, -76, 428, -58, 5, -76, -58, -184, 428, 428,
      -76, -76, 50, -76, 5, 50, -13, -58, -13, -454, 50, 68, -13, 50, -76, 68,
      -76, -76, -13, -13, 5, -76, 5, -40, -58, -13, -58, -454, 68, 68, -76,
      -76, -40, -184, 428, -40, -40, -58, -76, 176, -184, -184, -40, -40, 572,
      572, -76, -40, -58, 68, 104, -76, 68, -40, -76, -454, 68, 50, 50, -454,
      68, -76, -13, 572, -76, 50, -13, 428, 554, 3641, -454, -454, 50, 50, 572,
      428, -76, -184, 428, 50, -76}, {-58, -58, -40, -76, -13, 5, 104, -58, 68,
      -76, -40, 572, 572, 50, -76, 428, -13, 554, 5, -58, -76, 50, -13, -13,
      -40, -40, -184, 428, -76, -58, 572, -76, -40, -58, 68, 104, 176, -184,
      -184, -40, -40, 572, 68, -76, -454, 68, -40, -76, 68, -13, -454, 50, -76,
      50, -76, 5, -40, -58, -13, -58, 68, -76, -76, -13, -13, 5, 68, -454, -76,
      -40, 68, -76, -40, -76, -40, -184, -58, 428, -184, 428, 428, -76, -76,
      50, -76, 428, -58, 5, -76, -58, 5, -76, -58, -13, 50, -13, -13, -76, 68,
      -454, 50, 50, -76, 428, 50, -76, -184, 428, -454, 3641, 50, 572, -454,
      50, -76, 572, 428, 554, 50, -13}, {-40, -184, -40, -76, 428, -58, -76,
      428, -58, 5, -76, -58, -184, 428, 428, -76, -76, 50, -76, 5, 50, -13,
      -58, -13, -76, 428, 572, 50, 554, -13, 68, -76, 104, -58, 572, -40, -40,
      -76, -58, -58, 5, -13, -13, -13, 50, -76, -58, 5, 68, -76, -76, -13, -13,
      5, -454, 50, 68, -13, 50, -76, -40, -58, -76, 5, -58, -13, 68, -76, -454,
      68, -40, -76, 176, -184, -184, -40, -40, 572, -184, 428, -40, -40, -58,
      -76, -40, -58, 572, -76, 104, 68, -40, -76, -76, 68, 68, -454, 572, -76,
      50, -13, 428, 554, 50, 50, -454, 68, -76, -13, -454, 50, 3641, -454, 572,
      50, -184, 428, 428, -76, -76, 50}, {-40, -40, -184, 428, -76, -58, 572,
      -76, -40, -58, 68, 104, 176, -184, -184, -40, -40, 572, 68, -76, -454,
      68, -40, -76, -58, -58, -40, -76, -13, 5, 104, -58, 68, -76, -40, 572,
      572, 50, -76, 428, -13, 554, 5, -58, -76, 50, -13, -13, -76, 5, -40, -58,
      -13, -58, 68, -13, -454, 50, -76, 50, -76, -13, 68, -76, 5, -13, -76,
      -40, 68, -454, -76, 68, -184, 428, 428, -76, -76, 50, -40, -76, -40,
      -184, -58, 428, -58, 5, -76, 428, -58, -76, -58, -13, 5, -76, -13, 50,
      -76, 428, 50, -76, -184, 428, -13, -76, 68, -454, 50, 50, 50, 572, -454,
      3641, 50, -454, 428, 554, -76, 572, -13, 50}, {428, -76, -184, 428, 50,
      -76, -58, 5, -76, 428, -58, -76, -40, -76, -40, -184, -58, 428, -13, -58,
      -13, 50, 5, -76, -184, -40, 176, -184, 572, -40, -40, -58, 572, -76, 104,
      68, -184, 428, -40, -40, -58, -76, -76, -40, 68, -454, -76, 68, -76, -13,
      68, -76, 5, -13, -40, -58, -76, 5, -58, -13, -454, 50, 68, -13, 50, -76,
      -76, 68, -40, -76, -454, 68, 572, 50, -76, 428, -13, 554, -40, -76, -58,
      -58, 5, -13, 68, -76, 104, -58, 572, -40, -13, -13, -58, 5, 50, -76, 554,
      428, -13, 50, -76, 572, 428, -184, -76, 50, 428, -76, 50, -454, 572, 50,
      3641, -454, 50, 50, -76, -13, -454, 68}, {-184, -40, 176, -184, 572, -40,
      -40, -58, 572, -76, 104, 68, -184, 428, -40, -40, -58, -76, -76, -40, 68,
      -454, -76, 68, 428, -76, -184, 428, 50, -76, -58, 5, -76, 428, -58, -76,
      -40, -76, -40, -184, -58, 428, -13, -58, -13, 50, 5, -76, -40, -58, -76,
      5, -58, -13, -76, -13, 68, -76, 5, -13, 68, -13, -454, 50, -76, 50, -40,
      -76, -76, 68, 68, -454, -40, -76, -58, -58, 5, -13, 572, 50, -76, 428,
      -13, 554, 104, -58, 68, -76, -40, 572, -58, 5, -13, -13, -76, 50, 428,
      -184, -76, 50, 428, -76, 554, 428, -13, 50, -76, 572, 572, 50, 50, -454,
      -454, 3641, -76, -13, 50, 50, 68, -454}, {572, 50, -76, 428, -13, 554,
      -40, -76, -58, -58, 5, -13, 68, -76, 104, -58, 572, -40, -13, -13, -58,
      5, 50, -76, -40, -76, -40, -184, -58, 428, -184, 428, 428, -76, -76, 50,
      -76, 428, -58, 5, -76, -58, 5, -76, -58, -13, 50, -13, -184, 428, -40,
      -40, -58, -76, 176, -184, -184, -40, -40, 572, 572, -76, -40, -58, 68,
      104, -76, 68, -40, -76, -454, 68, -454, 50, 68, -13, 50, -76, 68, -76,
      -76, -13, -13, 5, -76, 5, -40, -58, -13, -58, -454, 68, 68, -76, -76,
      -40, 50, 50, -76, -13, -454, 68, -76, 572, 428, 554, 50, -13, 428, -76,
      -184, 428, 50, -76, 3641, -454, -454, 50, 50, 572}, {-40, -76, -58, -58,
      5, -13, 572, 50, -76, 428, -13, 554, 104, -58, 68, -76, -40, 572, -58, 5,
      -13, -13, -76, 50, -184, 428, -40, -40, -58, -76, 176, -184, -184, -40,
      -40, 572, 572, -76, -40, -58, 68, 104, -76, 68, -40, -76, -454, 68, -40,
      -76, -40, -184, -58, 428, -184, 428, 428, -76, -76, 50, -76, 428, -58, 5,
      -76, -58, 5, -76, -58, -13, 50, -13, 68, -13, -454, 50, -76, 50, -76, 5,
      -40, -58, -13, -58, 68, -76, -76, -13, -13, 5, 68, -454, -76, -40, 68,
      -76, -76, -13, 50, 50, 68, -454, 428, -76, -184, 428, 50, -76, -76, 572,
      428, 554, 50, -13, -454, 3641, 50, 572, -454, 50}, {-40, -76, -40, -184,
      -58, 428, -184, 428, 428, -76, -76, 50, -76, 428, -58, 5, -76, -58, 5,
      -76, -58, -13, 50, -13, 572, 50, -76, 428, -13, 554, -40, -76, -58, -58,
      5, -13, 68, -76, 104, -58, 572, -40, -13, -13, -58, 5, 50, -76, 176,
      -184, -184, -40, -40, 572, -184, 428, -40, -40, -58, -76, -40, -58, 572,
      -76, 104, 68, -40, -76, -76, 68, 68, -454, 68, -76, -76, -13, -13, 5,
      -454, 50, 68, -13, 50, -76, -40, -58, -76, 5, -58, -13, 68, -76, -454,
      68, -40, -76, -76, 572, 428, 554, 50, -13, 50, 50, -76, -13, -454, 68,
      -184, 428, 428, -76, -76, 50, -454, 50, 3641, -454, 572, 50}, {-184, 428,
      -40, -40, -58, -76, 176, -184, -184, -40, -40, 572, 572, -76, -40, -58,
      68, 104, -76, 68, -40, -76, -454, 68, -40, -76, -58, -58, 5, -13, 572,
      50, -76, 428, -13, 554, 104, -58, 68, -76, -40, 572, -58, 5, -13, -13,
      -76, 50, -184, 428, 428, -76, -76, 50, -40, -76, -40, -184, -58, 428,
      -58, 5, -76, 428, -58, -76, -58, -13, 5, -76, -13, 50, -76, 5, -40, -58,
      -13, -58, 68, -13, -454, 50, -76, 50, -76, -13, 68, -76, 5, -13, -76,
      -40, 68, -454, -76, 68, 428, -76, -184, 428, 50, -76, -76, -13, 50, 50,
      68, -454, 428, 554, -76, 572, -13, 50, 50, 572, -454, 3641, 50, -454},
      {-184, 428, 428, -76, -76, 50, -40, -76, -40, -184, -58, 428, -58, 5,
      -76, 428, -58, -76, -58, -13, 5, -76, -13, 50, 176, -184, -184, -40, -40,
      572, -184, 428, -40, -40, -58, -76, -40, -58, 572, -76, 104, 68, -40,
      -76, -76, 68, 68, -454, 572, 50, -76, 428, -13, 554, -40, -76, -58, -58,
      5, -13, 68, -76, 104, -58, 572, -40, -13, -13, -58, 5, 50, -76, -76, -13,
      68, -76, 5, -13, -40, -58, -76, 5, -58, -13, -454, 50, 68, -13, 50, -76,
      -76, 68, -40, -76, -454, 68, 428, 554, -76, 572, -13, 50, -184, 428, 428,
      -76, -76, 50, 50, 50, -76, -13, -454, 68, 50, -454, 572, 50, 3641, -454},
      {176, -184, -184, -40, -40, 572, -184, 428, -40, -40, -58, -76, -40, -58,
      572, -76, 104, 68, -40, -76, -76, 68, 68, -454, -184, 428, 428, -76, -76,
      50, -40, -76, -40, -184, -58, 428, -58, 5, -76, 428, -58, -76, -58, -13,
      5, -76, -13, 50, -40, -76, -58, -58, 5, -13, 572, 50, -76, 428, -13, 554,
      104, -58, 68, -76, -40, 572, -58, 5, -13, -13, -76, 50, -40, -58, -76, 5,
      -58, -13, -76, -13, 68, -76, 5, -13, 68, -13, -454, 50, -76, 50, -40,
      -76, -76, 68, 68, -454, -184, 428, 428, -76, -76, 50, 428, 554, -76, 572,
      -13, 50, -76, -13, 50, 50, 68, -454, 572, 50, 50, -454, -454, 3641}};
}
