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

	// Bin edges
	struct EDGE {

		double min = 0;
		double max = 0;

		double center() const {
			return (max+min)/2.0;
		}
	};

	// Parameters
	struct HPARAM {

		std::vector<double> M  = {0,0,0};
		std::vector<double> Y  = {0,0,0};
		std::vector<double> PT = {0,0,0};

		int    LMAX            = 2;      // Maximum spherical harmonic truncation degree (non-negative integer)
		bool   REMOVEODD       = true;   // Fix odd moments to zero (due to lacking spesific spin states, for example)
		bool   REMOVENEGATIVEM = true;   // Fix negative m to zero (due to parity conservation, for example)
		bool   EML             = false;  // Use Extended Maximum Likelihood fit
		double SVDREG          = 0.001;  // SVD regularization strength in algebraic inverse (put 0 for no regularization) 
		double L1REG           = 0.001;  // L1-norm regularization in EML fit (put 0 for no regularization)
		
		void Print() {
			std::cout << "HARMONIC EXPANSION PARAMETERS:" << std::endl << std::endl;

			std::cout << "LMAX            = " << LMAX << std::endl;
			std::cout << "REMOVEODD       = " << (REMOVEODD ? "true" : "false") << std::endl;
			std::cout << "REMOVENEGATIVEM = " << (REMOVENEGATIVEM ? "true" : "false") << std::endl;
			std::cout << "EML             = " << (EML ? "true" : "false") << std::endl;
			std::cout << "SVDREG          = " << SVDREG << std::endl;
			std::cout << "L1REG           = " << L1REG << std::endl;
		}
	};

	// Constructor, destructor
	MHarmonic();
	~MHarmonic();

	void   Init(const HPARAM& hp);

	void   HyperLoop(void (*fitfunc)(int&, double*, double&, double*, int),
					 const std::vector<gra::spherical::Omega>& MC,
					 const std::vector<gra::spherical::Data>& DATA, const HPARAM& hp);
	
	void   MomentFit(const gra::spherical::Meta& META, const std::vector<std::size_t>& cell,
					 void (*fitfunc)(int&, double*, double&, double*, int));

	double PrintOutHyperCell(const gra::spherical::Meta& META, const std::vector<std::size_t>& cell);
	void   logLfunc(int& npar, double* gin, double& f, double* par, int iflag) const;

	bool   PrintLoop(const std::string& output) const;
	void   PlotAll(const std::string& outputpath) const;

	void   PlotFigures(const std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>>& tensor, unsigned int OBSERVABLE,
                       const std::string& TYPESTRING, int barcolor,
                       const std::string& outputpath) const;

	void   PlotFigures2D(const std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>>& tensor, const std::vector<int>& OBSERVABLE2,
                         const std::string& TYPESTRING, int barcolor,
                         const std::string& outputpath) const;

	void   Plot1DEfficiency(unsigned int OBSERVABLE, const std::string& outputpath) const;

	// Parameters
	HPARAM param;
	int    NCOEF;

private:

	MMatrix<double> Y_lm;                             // Calculate once for speed
	std::vector<gra::spherical::Omega> DATA_events;   // Input data

	// ====================================================================
	// Needed by MINUIT fit/cost function

	// Currently active
	std::vector<std::size_t> DATA_ind;
	std::vector<std::size_t> activecell;			  // Hypercell indices
	std::vector<double>      t_lm;					  // Fitted moments
	std::vector<double>      t_lm_error;			  // Their errors

	// Error and covariance matrices
	MMatrix<double> errmat;
	MMatrix<double> covmat;

	// ====================================================================

	// Active moments bookkeeping here
	std::vector<bool> ACTIVE;
	unsigned int ACTIVENDF = 0;

	// Detector expansion tensor data
	MTensor<gra::spherical::SH_DET> det_DET;
	MTensor<gra::spherical::SH_DET> fid_DET;
	MTensor<gra::spherical::SH_DET> fla_DET;
	
	// Expanded data in a map for different data sources
	std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>> det;	
	std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>> fid;
	std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>> fla;

	// Discretization
	std::vector<std::vector<EDGE>> grid;
};

} // gra namespace ends

#endif
