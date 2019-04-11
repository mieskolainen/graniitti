// Spherical harmonic expansion and analysis class
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHARMONIC_H
#define MHARMONIC_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MRandom.h"
#include "Graniitti/MDimArray.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MSpherical.h"
#include "Graniitti/MTensor.h"


namespace gra {

class MHarmonic {

public:

	// Parameters
	struct HPARAM {

		std::vector<double> M  = {0,0,0};
		std::vector<double> Y  = {0,0,0};
		std::vector<double> PT = {0,0,0};

		int    LMAX;            // Maximum spherical harmonic truncation degree (e.g. 4)
		bool   REMOVEODD;       // Fix odd moments to zero (due to lacking spesific spin states, for example)
		bool   REMOVENEGATIVEM; // Fix negative m to zero (due to parity conservation, for example)
		bool   EML;             // Use Extended Maximum Likelihood fit
		double LAMBDA;          // L1-norm regularization (put 0 for no regularization)
		
		void Print() {
			std::cout << "LMAX            = " << LMAX << std::endl;
			std::cout << "REMOVEODD       = " << (REMOVEODD ? "true" : "false") << std::endl;
			std::cout << "REMOVENEGATIVEM = " << (REMOVENEGATIVEM ? "true" : "false") << std::endl;
			std::cout << "EML             = " << (EML ? "true" : "false") << std::endl;
			std::cout << "LAMBDA          = " << LAMBDA << std::endl;
		}
	};

	// Constructor, destructor
	 MHarmonic();
	~MHarmonic();

	void   Init(const HPARAM& hp);
	void   HyperLoop(void (*fitfunc)(int&, double*, double&, double*, int),
					 const std::vector<gra::spherical::Omega>& MC,
					 const std::vector<gra::spherical::Omega>& DATA, const HPARAM& hp);
	void   MomentFit(const std::vector<std::size_t>& cell, void (*fitfunc)(int&, double*, double&, double*, int));
	double PrintOutHyperCell(const std::vector<std::size_t>& cell);
	void   logLfunc(int& npar, double* gin, double& f, double* par, int iflag) const;
	

	bool   PrintLoop(const std::string& output) const;
	void   PlotAll() const;
	void   PlotFigures(const MTensor<gra::spherical::SH>& tensor, const std::string& DATATYPE,
					   unsigned int OBSERVABLE,
					   const std::string& outputfile, int barcolor) const;

	void   PlotFigures2D(const MTensor<gra::spherical::SH>& tensor,
				         const std::string& DATATYPE, std::vector<int> OBSERVABLE,
                         const std::string& outputfile,
                         int barcolor) const;

	// Parameters
	HPARAM param;
	int    NCOEF;

private:

	MMatrix<double> Y_lm;                             // Calculate once for speed
	std::vector<gra::spherical::Omega> DATA_events;   // Input data

	// ====================================================================
	// Needed by MINUIT fit/cost function

	// Currently active
	std::vector<std::size_t>  DATA_ind;
	std::vector<std::size_t>  activecell;	// Hypercell indices
	std::vector<double> t_lm;              	// Fitted moments
	std::vector<double> t_lm_error;        	// Their errors

	// Error and covariance matrices
	MMatrix<double> errmat;
	MMatrix<double> covmat;

	// ====================================================================

	// Active moments bookkeeping here
	std::vector<bool> ACTIVE;
	unsigned int ACTIVENDF = 0;

	// Tensor data:
	//
	// ref = reference phase space
	// fid = fiducial phase space
	// det = detector level
	MTensor<gra::spherical::SH> ref;
	MTensor<gra::spherical::SH> fid;
	MTensor<gra::spherical::SH> det;

	// x-axis center points
	std::vector<std::vector<double>> xcenter;
};

} // gra namespace ends


#endif