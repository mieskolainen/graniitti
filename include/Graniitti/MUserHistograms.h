// Container class for different type of histograms
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MUSERHISTOGRAMS_H
#define MUSERHISTOGRAMS_H

// C++
#include <complex>
#include <vector>

// Own
#include "Graniitti/MH1.h"
#include "Graniitti/MH2.h"
#include "Graniitti/MKinematics.h"

namespace gra {

class MUserHistograms {
       public:
	// Constructor, destructor
	MUserHistograms() {
	}
	~MUserHistograms() {
	}

	void InitHistograms();
	void FillHistograms(double totalweight, const gra::LORENTZSCALAR &scalar);
	void PrintHistograms();
	void SetHistograms(unsigned int in) {
		HIST = in;
	}
	void FillCosThetaPhi(double totalweight, const gra::LORENTZSCALAR &scalar);

	// Histograms indexed by name std::string
	std::map<std::string, MH1<double>> h1;
	std::map<std::string, MH2> h2;

	unsigned int HIST = 0; // Histogramming level
};

} // gra namespace ends

#endif