// Functional methods for Spherical Harmonic Expansions
// 
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <iostream>
#include <vector>
#include <complex>

// Own
#include "Graniitti/MSpherical.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MRandom.h"
#include "Graniitti/MAux.h"


// Libraries
#include "rang.hpp"


using gra::math::msqrt;
using gra::math::pow2;
using gra::math::PI;
using gra::math::zi;
using gra::aux::indices;


namespace gra {
namespace spherical {


// Monte Carlo integral I ~= V 1/N \sum_{i=1}^N f(\vec{x}_i) = V <f(x)>
// True integral being  I = \int_\Omega f(\vec{x}) d\vec{x}
//
// This inner-product (overlap integral) matrix is identity if the phase space
// is flat (uniform). In a geometrically restricted phase space, these basis
// functions get mixed => lost orthogonality.

MMatrix<double> GetGMixing(const std::vector<Omega>& events, const std::vector<std::size_t>& ind, int LMAX, int mode) {
  
  const int NCOEF = (LMAX+1)*(LMAX+1);

  std::cout << "GetGMixing: mode = " << mode << std::endl;
  std::cout << "Generated reference MC phase space events = " << ind.size() << std::endl;

  // Construct the efficiency coefficients EPSILON_LM with linear
  // indexing
  MMatrix<double>   E(NCOEF, NCOEF, 0.0);
  MMatrix<double>  E2(NCOEF, NCOEF, 0.0); // for the uncertainty

  // Loop over GENERATED MC events and do the integral in effect via
  // uniform MC sampling
  // We evaluate the integral:
  //
  // eta_LM = \int eta(Omega) Re Y_LM(Omega) dOmega, Omega =
  // (costheta,phi)
  //
  int fiducial = 0;
  int selected = 0;

  const double VOL = 4.0*PI; // [costheta] x [phi] plane area

  for (const auto& k : ind) {

    const bool fidtrue = events[k].fiducial;
    const bool seltrue = events[k].selected;

    // Full reference phase space
    if (mode == 0) {
      // all fine
    }
    // Geometric acceptance
    if (mode == 1) {
      if (!fidtrue) { continue; } else { ++fiducial; }
    }
    // Geometric x Efficiency
    if (mode == 2) {
      if (!fidtrue) { continue; } else { ++fiducial; }
      if (!seltrue) { continue; } else { ++selected; }
    }

    const double costheta = events[k].costheta;
    const double phi      = events[k].phi;

    for (int l = 0; l <= LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {
        const int index = LinearInd(l, m);

        const std::complex<double> Y = gra::math::Y_complex_basis(costheta, phi, l, m);
        const double ReY = gra::math::NReY(Y, l, m);

        for (int lprime = 0; lprime <= LMAX; ++lprime) {
          for (int mprime = -lprime; mprime <= lprime; ++mprime) {

            const int indexprime = LinearInd(lprime, mprime);
            const std::complex<double> Yprime = gra::math::Y_complex_basis(costheta, phi, lprime, mprime);
            const double ReYprime = gra::math::NReY(Yprime, lprime, mprime);

            // Sum to the MC integral and
            // its squared version (for
            // uncertainty)
            const double f = VOL * ReY * ReYprime; // Note volume term
            E[index][indexprime]  += f;
            E2[index][indexprime] += f * f;
          }
        }
      }
    }
  }

  if (mode == 1 || mode == 2) {
    printf("Fiducial reference MC phase space events = %d (acceptance %0.3f percent) \n",
        fiducial, fiducial / static_cast<double>(ind.size()) * 100);
  }
  if (mode == 2) {
    printf("Selected referemce MC phase space events = %d (efficiency %0.3f percent) \n",
        selected, selected / static_cast<double>(fiducial) * 100);
  }
  if (mode == 1 || mode == 2) {
    if (fiducial < 50 || selected < 50) {
      printf("GetGMixing:: Too low MC event count in the mass bin! \n");
    }
  }
  std::cout << std::endl;

  // Do the normalization
  const double N_generated = static_cast<double>(ind.size()); // Generated events within this mass interval
  std::vector<double> E_error(E.size_row(), 0.0);
  std::cout << "Symmetric mixing matrix integral coefficients G_{ll'}^{mm'}:" << std::endl;
  std::vector<double> rowsum(NCOEF, 0.0);
  std::cout << "         " << std::endl;

  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      printf("G_%d%d  \t", l, m);
    }
  }
  std::cout << std::endl;

  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = LinearInd(l, m);

      printf("G_%d%d : ", l, m);
      double colsum = 0;

      for (int lprime = 0; lprime <= LMAX; ++lprime) {
        for (int mprime = -lprime; mprime <= lprime; ++mprime) {

          const int indexprime = LinearInd(lprime, mprime);

          E[index][indexprime]  /= N_generated;
          E2[index][indexprime] /= N_generated;
          colsum += E[index][indexprime];

          double error = CalcError(E2[index][indexprime], E[index][indexprime], N_generated);
          printf("%6.3f +- %6.3f \t", E[index][indexprime], error);
        }
      }
      printf(" | %6.3f \n", colsum);
    }
  }
  
  gra::aux::PrintBar("-");
  std::cout << "       " << std::endl;
  for (std::size_t j = 0; j < E.size_row(); ++j) {
    double rowsum = 0;
    for (std::size_t i = 0; i < E.size_col(); ++i) {
      rowsum += E[i][j];
    }
    printf("%6.3f\t", rowsum);
  }
  printf("\n\n\n\n");

  // Calculate matrix condition number via SVD
  Eigen::MatrixXd A = gra::aux::Matrix2Eigen(E);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
  Eigen::MatrixXd sval = svd.singularValues();

  std::cout << "SVD singular values of the matrix:" << std::endl;
  for (std::size_t i = 0; i < E.size_row(); ++i) {
    printf("%0.4f ", sval(i));
  }
  std::cout << std::endl;
  const double conditionnumber = sval(0) / sval(E.size_row() - 1);
  if (conditionnumber < 10) {
    std::cout << rang::fg::green;
  } else {
    std::cout << rang::fg::red;
  }

  printf("Condition number: sMax/sMin = %0.4f \n", conditionnumber);
  std::cout << std::endl << std::endl << std::endl;
  std::cout << rang::fg::reset;
  return E;
}


// Acceptance expansion coefficients using MC events
std::pair<std::vector<double>,std::vector<double>>
  GetELM(const std::vector<Omega>& MC, const std::vector<std::size_t>& ind, int LMAX, int mode) {
  
  const int NCOEF = (LMAX+1)*(LMAX+1);

  std::cout << "GetELM: mode = " << mode << std::endl;
  std::cout << "Generated reference MC phase space events = " << ind.size() << std::endl;

  // Construct the efficiency coefficients E_LM with linear indexing
  std::vector<double>  E(NCOEF, 0.0);
  std::vector<double> E2(NCOEF, 0.0); // for the uncertainty

  // Loop over GENERATED MC events and do the integral in effect via
  // uniform MC sampling
  // We evaluate the integral:
  //
  // E_LM = \int E(Omega) Re Y_LM(Omega) dOmega, Omega = (costheta,phi)

  int fiducial = 0;
  int selected = 0;
  const double V = sqrt(4.0 * PI); // Normalization volume

  for (const auto& k : ind) {

    bool fid = MC[k].fiducial;
    bool sel = MC[k].selected;

    // Full reference phase space
    if (mode == 0) {
      // all fine
    }
    // Geometric acceptance
    if (mode == 1) {
      if (!fid) { continue; } else { ++fiducial; }
    }
    // Geometric x Efficiency
    if (mode == 2) {
      if (!fid) { continue; } else { ++fiducial; }
      if (!sel) { continue; } else { ++selected; }
    }
    
    const double costheta = MC[k].costheta;
    const double phi      = MC[k].phi;

    for (int l = 0; l <= LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {

        const std::complex<double> Y = gra::math::Y_complex_basis(costheta, phi, l, m);

        // Sum to the MC integral and its square (for uncertainty)
        double f = V * gra::math::NReY(Y, l, m); // Note volume V term

        const int index = LinearInd(l, m);
        E[index]  += f;
        E2[index] += pow2(f);
      }
    }
  }
  if (mode == 1 || mode == 2) {
    printf(
        "Fiducial reference MC phase space events = %d (geometric-kinematic acceptance %0.3f percent) \n\n",
        fiducial, fiducial / static_cast<double>(ind.size()) * 100);
  }
  if (mode == 2) {
    printf(
        "Selected reference MC phase space events = %d (fiducial efficiency %0.3f percent) \n\n",
        selected, selected / static_cast<double>(fiducial) * 100);
  }
  if (mode == 1) {
    if (fiducial < 1000) {
      std::cout << "GetELM:: Very low [fiducial] MC event count = "<< fiducial << " in the mass bin!" << std::endl;
    }
  }
  if (mode == 2) {
    if (selected < 1000) {
      std::cout << "GetELM:: Very low [selected] MC event count = "<< selected << " in the mass bin!" << std::endl;
    }
  }

  // Do the normalization
  const double N_generated = (double) ind.size(); // Generated events within this hypercell
  std::vector<double> E_error(E.size(), 0.0);

  std::cout << "Acceptance decomposition coefficients:" << std::endl;
  double sum = 0;
  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = LinearInd(l, m);

      E[index]  /= N_generated;
      E2[index] /= N_generated;

      // 1 sigma MC integration uncertainty
      const double error   = msqrt((E2[index] - std::pow(E[index], 2)) / N_generated);
      E_error[index] = error;
      sum += E[index];

      printf("E_%d%d \t= %12.8f +- %12.8f   (rel.error %9.3f percent) \n",
            l, m, E[index], error, std::abs(error / E[index]) * 100);
    }
  }
  printf("SUM_lm =   %0.8f \n", sum);
  std::cout << "Uncertanties by Monte Carlo errors = sqrt{(<f^2> - <f>^2) / n}" << std::endl;
  std::cout << std::endl;

  return {E, E_error};
}


// Calculate expansion coefficients directly via algebraic expansion
//
std::vector<double> SphericalMoments(const std::vector<Omega>& input,
                                     const std::vector<std::size_t>& ind, int LMAX, int mode) {

  const double V = sqrt(4.0 * PI); // Normalization volume
  const unsigned int NCOEF = (LMAX + 1) * (LMAX + 1);
  std::vector<double> t(NCOEF, 0.0);

  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {

      double sum = 0;
      const unsigned int index = LinearInd(l, m);
      for (const auto& k : ind) {

        // Check is this event selected by cuts
        const bool fidtrue = input[k].fiducial;
        const bool seltrue = input[k].selected;

        // Full reference phase space
        if (mode == 0) {
          // all fine
        }
        // Geometric acceptance
        if (mode == 1) {
          if (!fidtrue) { continue; }
        }
        // Geometric x Efficiency
        if (mode == 2) {
          if (!fidtrue) { continue; }
          if (!seltrue) { continue; }
        }
        const double costheta        = input[k].costheta;
        const double phi             = input[k].phi;
        const std::complex<double> Y = gra::math::Y_complex_basis(costheta, phi, l, m);

        sum += gra::math::NReY(Y, l, m);
      }
      t[index] = V * sum; // Note volume V term
    }
  }
  return t;
}


// Calculate Y_lm values for each event
MMatrix<double> YLM(const std::vector<Omega>& events, int LMAX) {

  std::cout << "YLM:" << std::endl;

  const unsigned int NCOEFF = (LMAX+1)*(LMAX+1);
  MMatrix<double> Y_lm(events.size(), NCOEFF, 0.0);
  
  // Loop over events
  for (const auto& k : indices(events)) {

    // Get event (costheta, phi)
    const double costheta = events[k].costheta;
    const double      phi = events[k].phi;

    // Loop over (l,m) coefficients
    for (int l = 0; l <= LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {

        const std::complex<double> Y = gra::math::Y_complex_basis(costheta, phi, l, m);

        // Save its value
        const int index = LinearInd(l, m);
        Y_lm[k][index]  = gra::math::NReY(Y, l, m);
      }
    }
  }
  return Y_lm;
}


// Spherical harmonic dot product between coefficients gives the integral.
// This is the property of orthonormality.
// \int G(Omega) x(omega) dOmega <=> \sum_i G_i x_i
double HarmDotProd(const std::vector<double>& G, const std::vector<double>& x, const std::vector<bool>& ACTIVE, int LMAX) {

  double sum = 0.0;
  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {

      const int index = LinearInd(l,m);
      if (ACTIVE[index]) {
        sum += G[index] * x[index];
      }
    }
  }
  return sum;
}


// TODO: add correlations in to the sum
double HarmDotProdError() {
  return 0.0;
}


void PrintOutMoments(const std::vector<double>& x, const std::vector<double>& x_error, const std::vector<bool>& ACTIVE, int LMAX) {

	for (int l = 0; l <= LMAX; ++l) {
		for (int m = -l; m <= l; ++m) {

			const unsigned int index = LinearInd(l,m);
			printf("t_%d%d \t= %8.1f +- %5.1f   ", l, m, x[index], x_error[index]);

			if (ACTIVE[index]) {
				std::cout << "[" << rang::fg::green
				          << "active" << rang::fg::reset
				          << "]" << std::endl;
			} else {
				std::cout << "[" << rang::fg::red
				          << "inactive" << rang::fg::reset
				          << "]" << std::endl;
			}
		}
	}
	std::cout << std::endl;
}


// Calculate indices for this interval
//
std::vector<std::size_t> GetIndices(const std::vector<Omega>& events,
				const std::vector<double>& M,
				const std::vector<double>& Pt,
				const std::vector<double>& Y) {

	std::vector<std::size_t> ind;

	for (const auto& i : indices(events)) {
		if ( events[i].M  >  M[0]  && events[i].M  <  M[1]
			&& events[i].Pt > Pt[0]  && events[i].Pt < Pt[1]
			&& events[i].Y  >  Y[0]  && events[i].Y  <  Y[1]) {

			ind.push_back(i);
		}
	}

  gra::aux::PrintBar("-");
  std::cout << rang::fg::green;
  printf("MASS RANGE: [%0.3f, %0.3f] GeV, PT RANGE: [%0.3f, %0.3f] GeV, Y RANGE: [%0.3f, %0.3f] : Events in this hyperbin %lu/%lu \n\n",
          M[0], M[1], Pt[0], Pt[1], Y[0], Y[1], ind.size(), events.size());
  std::cout << rang::fg::reset;

  return ind;
}


// TestIntegrals spherical harmonics
// This is a test function, all integrals should give 1 if the normalization is correct
void TestSphericalIntegrals(int LMAX) {

	MRandom random;
	const unsigned int NCOEF = (LMAX+1)*(LMAX+1);
	std::cout << "TestSphericalIntegrals: " << std::endl << std::endl;

	// Construct the efficiency coefficients EPSILON_LM with linear indexing
	std::vector<double>  E(NCOEF, 0.0);
	std::vector<double> E2(NCOEF, 0.0); // for the uncertainty

	// We MC evaluate the integral:
	// \int |Re Y_LM(Omega)|^2 dOmega, Omega = (costheta,phi)
	const double N = 100000; // Number of samples

	// Integral domain volume: |costheta| x |phi| = [-1,1] x [0,2pi]
	const double V = 4.0 * gra::math::PI;

	for (std::size_t k = 0; k < N; ++k) {

		const double costheta = random.U(-1.0, 1.0);
		const double phi      = random.U(0.0, 2.0*PI);
		
		for (int l = 0; l <= LMAX; ++l) {
			for (int m = -l; m <= l; ++m) {

				const std::complex<double> Y = gra::math::Y_complex_basis(costheta, phi, l, m);

				// Function value
				double f = std::pow(std::abs(gra::math::NReY(Y, l, m)), 2) * V;

				// Sum to the MC integral and its squared version (for uncertainty)
				const int ind = LinearInd(l, m);
				E[ind]  += f;
				E2[ind] += pow2(f);
			}
		}
	}

	for (int l = 0; l <= LMAX; ++l) {
		for (int m = -l; m <= l; ++m) {
			const int ind = LinearInd(l, m);
			E[ind]  /= N;
			E2[ind] /= N;

			// 1 sigma MC integration uncertainty
			double error = CalcError(E2[ind], E[ind], N);
			printf("I_%d%d \t= %12.8f +- %12.8f   (rel.error %9.3f percent) \n",
		    		l, m, E[ind], error, std::abs(error / E[ind]) * 100);
		}
	}
	std::cout << std::endl;
}

// Find out linear index
int LinearInd(int l, int m) {
	return l * (l + 1) + m;
}

// Standard error
double CalcError(double f2, double f, double N) {
	return sqrt((f2 - std::pow(f, 2)) / N);
}

// Print matrix to a file
void PrintMatrix(FILE* fp, const std::vector<std::vector<double>>& A) {

	// Print out coefficients
	for (std::size_t i = 0; i < A.size(); ++i) {
		for (std::size_t j = 0; j < A[i].size(); ++j) {
			fprintf(fp, "%0.1f ", A[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}

} // spherical namespace
} // gra namespace
