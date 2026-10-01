//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <cmath>

#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/HelAmps_sm_lepton_masses.h"
#include "Graniitti/Particle/MForm.h"



namespace MG5_PP_ZJ {

using namespace MG5_sm_lepton_masses;

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g d~ > z d~ WEIGHTED<=3 @2
// *   Decay: z > mu+ mu- WEIGHTED<=2
// Process: g s~ > z s~ WEIGHTED<=3 @2
// *   Decay: z > mu+ mu- WEIGHTED<=2

//--------------------------------------------------------------------------
// Initialize process.

void MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::initProc(std::string param_card_name) {
  InitParameters(SLHAReader(param_card_name));
}

// Initialize model parameters before constructing external wavefunctions
void MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::InitParameters(SLHAReader slha)
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm_lepton_masses();
  pars.setIndependentParameters(slha);
  gra::mg5::ValidateModel(slha, pars.Particles());
  pars.setIndependentCouplings();
  // pars.printIndependentParameters();
  // pars.printIndependentCouplings();
  // Reinitialization must replace the external mass table
  mME.clear();
  // Set external particle masses for this matrix element
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  mME.push_back(pars.mdl_MM);
  mME.push_back(pars.mdl_MM);
  mME.push_back(pars.ZERO);
  jamp2[0].assign(1, 0.0);
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::sigmaKin()
{
  // Set the parameters which change event by event
  pars.setDependentParameters(alphaS);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
  if (alphaQEDZero()) { pars.setAlphaQEDZero(); }


  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.;

  // Local variables and constants
  const int ncomb = 32;
  bool goodhel[ncomb] = {};
  int ntry = 0, sum_hel = 0, ngood = 0;
  int igood[ncomb + 1] = {};
  int jhel = 0;
  std::complex<double> * * wfs;
  double t[nprocesses];
  // Helicities for the process
  const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1}, {-1,
      -1, -1, -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1}, {-1, -1, 1, -1,
      -1}, {-1, -1, 1, -1, 1}, {-1, -1, 1, 1, -1}, {-1, -1, 1, 1, 1}, {-1, 1,
      -1, -1, -1}, {-1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1},
      {-1, 1, 1, -1, -1}, {-1, 1, 1, -1, 1}, {-1, 1, 1, 1, -1}, {-1, 1, 1, 1,
      1}, {1, -1, -1, -1, -1}, {1, -1, -1, -1, 1}, {1, -1, -1, 1, -1}, {1, -1,
      -1, 1, 1}, {1, -1, 1, -1, -1}, {1, -1, 1, -1, 1}, {1, -1, 1, 1, -1}, {1,
      -1, 1, 1, 1}, {1, 1, -1, -1, -1}, {1, 1, -1, -1, 1}, {1, 1, -1, 1, -1},
      {1, 1, -1, 1, 1}, {1, 1, 1, -1, -1}, {1, 1, 1, -1, 1}, {1, 1, 1, 1, -1},
      {1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {96, 96};

  ntry = ntry + 1;

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

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]);
        t[0] = matrix_2_gdx_zdx_z_mupmum();
        // Mirror initial state momenta for mirror process
        perm[0] = 1;
        perm[1] = 0;
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]);
        // Mirror back
        perm[0] = 0;
        perm[1] = 1;
        // Calculate matrix elements
        t[1] = matrix_2_gdx_zdx_z_mupmum();
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
      t[0] = matrix_2_gdx_zdx_z_mupmum();
      // Mirror initial state momenta for mirror process
      perm[0] = 1;
      perm[1] = 0;
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]);
      // Mirror back
      perm[0] = 0;
      perm[1] = 1;
      // Calculate matrix elements
      t[1] = matrix_2_gdx_zdx_z_mupmum();
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt;
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i];



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::sigmaHat()
{
  // Select between the different processes
  if(id1 == -3 && id2 == 21)
  {
    // Add matrix elements for processes with beams (-3, 21)
    return matrix_element[1];
  }
  else if(id1 == -1 && id2 == 21)
  {
    // Add matrix elements for processes with beams (-1, 21)
    return matrix_element[1];
  }
  else if(id1 == 21 && id2 == -3)
  {
    // Add matrix elements for processes with beams (21, -3)
    return matrix_element[0];
  }
  else if(id1 == 21 && id2 == -1)
  {
    // Add matrix elements for processes with beams (21, -1)
    return matrix_element[0];
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.;
  }
}

// Compute the generated subprocess selected by the incoming flavours
int MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::selectedProcess() const
{
  // Select between the different processes
  if(id1 == -3 && id2 == 21)
  {
    // Add matrix elements for processes with beams (-3, 21)
    return 1;
  }
  else if(id1 == -1 && id2 == 21)
  {
    // Add matrix elements for processes with beams (-1, 21)
    return 1;
  }
  else if(id1 == 21 && id2 == -3)
  {
    // Add matrix elements for processes with beams (21, -3)
    return 0;
  }
  else if(id1 == 21 && id2 == -1)
  {
    // Add matrix elements for processes with beams (21, -1)
    return 0;
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return -1;
  }
}

// Compute the identical-flavour multiplicity of the selected subprocess
double MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::selectedProcessMultiplicity() const
{
  // Select between the different processes
  if(id1 == -3 && id2 == 21)
  {
    // Add matrix elements for processes with beams (-3, 21)
    return 1.0;
  }
  else if(id1 == -1 && id2 == 21)
  {
    // Add matrix elements for processes with beams (-1, 21)
    return 1.0;
  }
  else if(id1 == 21 && id2 == -3)
  {
    // Add matrix elements for processes with beams (21, -3)
    return 1.0;
  }
  else if(id1 == 21 && id2 == -1)
  {
    // Add matrix elements for processes with beams (21, -1)
    return 1.0;
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.0;
  }
}

// Compute orthogonal complex components of the generated MG5 color sum
std::vector<std::complex<double>> MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::colorAmplitudes() const
{
  constexpr int ncolor = 1;
  std::complex<double> jamp[ncolor];
  jamp[0] = -amp[0] - amp[1];
  const std::vector<std::complex<double>> color_jamp(jamp, jamp + ncolor);
  const std::vector<double> denominators = {1};
  const std::vector<std::vector<double>> color_factors = {{4}};
  return gra::mg5helas::ColorMetricAmplitudes(color_jamp, denominators, color_factors);
}

// Compute raw MG5 leading-color flow amplitudes
std::vector<std::complex<double>> MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::leadingColorAmplitudes() const
{
  constexpr int ncolor = 1;
  std::complex<double> jamp[ncolor];
  jamp[0] = -amp[0] - amp[1];
  return std::vector<std::complex<double>>(jamp, jamp + ncolor);
}

// Compute fixed-basis complex helicity and orthogonal color components
std::vector<gra::mg5helas::HelicityComponent>
MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::helicityAmplitudes()
{
  const int selected = selectedProcess();
  if (selected < 0) { return {}; }
  pars.setDependentParameters(alphaS);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
  if (alphaQEDZero()) { pars.setAlphaQEDZero(); }

  const int ncomb = 32;
  const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1}, {-1,
      -1, -1, -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1}, {-1, -1, 1, -1,
      -1}, {-1, -1, 1, -1, 1}, {-1, -1, 1, 1, -1}, {-1, -1, 1, 1, 1}, {-1, 1,
      -1, -1, -1}, {-1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1},
      {-1, 1, 1, -1, -1}, {-1, 1, 1, -1, 1}, {-1, 1, 1, 1, -1}, {-1, 1, 1, 1,
      1}, {1, -1, -1, -1, -1}, {1, -1, -1, -1, 1}, {1, -1, -1, 1, -1}, {1, -1,
      -1, 1, 1}, {1, -1, 1, -1, -1}, {1, -1, 1, -1, 1}, {1, -1, 1, 1, -1}, {1,
      -1, 1, 1, 1}, {1, 1, -1, -1, -1}, {1, 1, -1, -1, 1}, {1, 1, -1, 1, -1},
      {1, 1, -1, 1, 1}, {1, 1, 1, -1, -1}, {1, 1, 1, -1, 1}, {1, 1, 1, 1, -1},
      {1, 1, 1, 1, 1}};
  const int denominators[nprocesses] = {96, 96};
  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) { perm[i] = i; }
  if (selected == 1) { std::swap(perm[0], perm[1]); }

  // Keep the complete MadGraph denominator until the stable final state is known
  // The subprocess_sum restores its final-state symmetry factor once
  const double average = std::sqrt(selectedProcessMultiplicity() /
                                   static_cast<double>(denominators[selected]));
  std::vector<gra::mg5helas::HelicityComponent> components;
  for (int ihel = 0; ihel < ncomb; ++ihel) {
    calculate_wavefunctions(perm, helicities[ihel]);
    const auto colors = colorAmplitudes();
    const auto leading_flows = leadingColorAmplitudes();
    if (colors.empty() || leading_flows.empty()) { return {}; }
    std::array<int, 2> physical_helicities = {0, 0};
    physical_helicities[perm[0]] = helicities[ihel][0];
    physical_helicities[perm[1]] = helicities[ihel][1];
    std::vector<int> outgoing;
    for (int leg = ninitial; leg < nexternal; ++leg) {
      outgoing.push_back(helicities[ihel][leg]);
    }
    for (std::size_t color = 0; color < colors.size(); ++color) {
      gra::mg5helas::HelicityComponent component{physical_helicities, outgoing, color,
                                                   average * colors[color], {}};
      // Store raw flows once because the orthogonal color components already span this helicity
      if (color == 0) {
        component.flow_values.reserve(leading_flows.size());
        for (const auto flow : leading_flows) {
          component.flow_values.push_back(average * flow);
        }
      }
      components.push_back(std::move(component));
    }
  }
  return components;
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j;

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]);
  oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
  FFV2_4_3(w[2], w[3], pars.GC_50, pars.GC_59, pars.mdl_MZ, pars.mdl_WZ,
      w[4]);
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[5]);
  FFV1_2(w[5], w[0], pars.GC_11, pars.ZERO, pars.ZERO, w[6]);
  FFV1_1(w[1], w[0], pars.GC_11, pars.ZERO, pars.ZERO, w[7]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV2_3_0(w[6], w[1], w[4], pars.GC_50, pars.GC_58, amp[0]);
  FFV2_3_0(w[5], w[7], w[4], pars.GC_50, pars.GC_58, amp[1]);

}
double MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx::matrix_2_gdx_zdx_z_mupmum()
{
  int i, j;
  // Local variables
  const int ngraphs = 2;
  const int ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {1};
  static const double cf[ncolor][ncolor] = {{4}};

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1];

  // Sum and square the color flows to get the matrix element
  double matrix = 0;
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.;
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j];
    matrix = matrix + std::real(ztemp * std::conj(jamp[i]))/denom[i];
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += std::real(jamp[i] * std::conj(jamp[i]));

  return matrix;
}




}  // namespace MG5_PP_ZJ
