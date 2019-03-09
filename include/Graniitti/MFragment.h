// Toy fragmentation class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFRAGMENTATION_H
#define MFRAGMENTATION_H

// C++
#include <complex>
#include <random>
#include <vector>
#include <map>

// Own
#include "Graniitti/MPDG.h"
#include "Graniitti/MRandom.h"


namespace gra {

class MFragment {

public:
	
	MFragment()  {}
	~MFragment() {}

	static double TubeFragment(const M4Vec& mother, double M0, const std::vector<double>& m,
		std::vector<M4Vec>& p, double q, double T, double maxpt, std::mt19937_64& rng);

	static bool SolveAlpha(double& alpha, double M0, const std::vector<double>& m,
			const std::valarray<double>& mt, const std::valarray<double>& y);

	static void ExpPowRND(double q, double T, double maxpt,
		const std::vector<double>& mass, std::vector<double>& x, std::mt19937_64& rng);

	static void GetDecayStatus(const std::vector<int>& pdgcode, std::vector<bool>& isstable);
	
	static void GetForwardMass(double& mass1, double& mass2, bool& excite1, bool& excite2, unsigned int excite, MRandom& random);
	static void GetSingleForwardMass(double& mass, MRandom& random);
	static void NstarDecayTable(double m0, std::vector<int>& pdgcode, std::mt19937_64& rng);

	static int PickParticles(double M, unsigned int N, int B, int S, int Q,
					  std::vector<double>& mass, std::vector<int>& pdgcode, const MPDG& PDG, std::mt19937_64& rng);

	// Uniform random numbers from [a,b)
  	template <typename T>
  	static inline double U(double a, double b, T& rng) {

		// C++11 thread_local is also static
    	thread_local std::uniform_real_distribution<double> flat; // Default [0,1)
    	return a + (b - a) * flat(rng);
  	}
};


} // gra namespace ends


#endif
