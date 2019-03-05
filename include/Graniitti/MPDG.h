// PDG class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MPDG_H
#define MPDG_H

// C++
#include <complex>
#include <random>
#include <vector>
#include <map>

#include "Graniitti/MKinematics.h"


namespace gra {

namespace PDG {
	
	// Common PDG ids, http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
	// For MC internal, use 81-100
	constexpr const int PDG_p     = 2212;               // Proton
	constexpr const int PDG_n     = 2112;               // ...
	constexpr const int PDG_pip   = 211;
	constexpr const int PDG_pim   = -211;
	constexpr const int PDG_pi0   = 111;
	constexpr const int PDG_Kp    = 321;
	constexpr const int PDG_Km    = -321;
	constexpr const int PDG_gamma = 22;
	constexpr const int PDG_gluon = 21;
	constexpr const int PDG_muon  = 13;
	
	constexpr const int PDG_pomeron = 990;
	constexpr const int PDG_reggeon = 110;
	constexpr const int PDG_MSTAR   = 90000;            // Meson  resonance
	constexpr const int PDG_NSTAR   = 90210;            // Baryon resonance
	constexpr const int PDG_system  = 93;               // "Container" / "Jet" / "System"

	// HepPMC STABLE definitions conventions
	constexpr const int PDG_BEAM         = 4;           // Beam particles
	constexpr const int PDG_STABLE       = 1;           // Final state "stable" particles
	constexpr const int PDG_DECAY        = 2;           // Such as pi0
	constexpr const int PDG_INTERMEDIATE = 81;          // Such as Pomerons etc

	// For CPU efficiency [GeV]
	constexpr const double mp   = 0.938272081;
	constexpr const double mK   = 0.493677;
	constexpr const double mpi  = 0.13957018;
	constexpr const double mpi0 = 0.13497660;


	// [hbar] = [M][L]^{2][T]^{-1} and [c] = [L][T]^{-1}
	// set (hbar = c = 1)
	// then [M] = [L]^{-1} = T^{-1}

	// http://pdg.lbl.gov/2018/reviews/rpp2018-rev-phys-constants.pdf
	// Basic constants
	constexpr const double c    = 2.99792458E8;         // c    = [m/s] (EXACT/DEFINITION)
	constexpr const double hbar = 6.582119514E-25;      // hbar = [GeV*s]
	constexpr const double eV   = 1.6021766208E-19;     // e    = [Joule], e~0.303 in nat.u.

	// Basic definition
	constexpr const double barn2m2 = 1E-28;             // 1 barn to [m^2]


	// Standard conversions
	constexpr const double GeV2m = hbar * c;            // 1 GeV^{-1} to [m]
	constexpr const double GeV2fm = GeV2m * 1E15;       // 1 GeV^{-1} to [fermi] (1 fm)
	constexpr const double GeV2m2 = GeV2m*GeV2m;        // 1 GeV^{-2} to [m^2]
	constexpr const double GeV2s = hbar;                // 1 GeV^{-1} to [s]

	constexpr const double GeV2barn = GeV2m2 / barn2m2; // 1 GeV^{-2} to barns
	constexpr const double GeV2mb = GeV2barn * 1E3;     // 1 GeV^{-2} to mbarns
	//constexpr const double GeV2ub = GeV2barn * 1E6;   // 1 GeV^{-2} to microbarn
	//constexpr const double GeV2nb = GeV2barn * 1E9;   // 1 GeV^{-2} to nanobarn
	//constexpr const double GeV2pb = GeV2barn * 1E12;  // 1 GeV^{-2} to picobarns
	//constexpr const double mb2GeV = 1.0 / GeV2mb;     // 1 mb to GeV^{-2}


	// More conversions
	//constexpr const double GeV2J   = eV * 1E9; // 1 GeV^{1} to [kg*m^2/s^2]=[Joule]
	//constexpr const double GeV2kg  = GeV2J / GeV2m2 * pow2(GeV2s); // 1 GeV^{1} to [kg]
	//constexpr const double GeV2N   = GeV2J / GeV2m; // 1 GeV^{2} to [N]=[kg*m/s^2] (Force)
	//constexpr const double GeV2mom = GeV2J / GeV2m * GeV2s; // 1 GeV^{1} to [kg*m/s] (Momentum)

	// cross section: [value] x [GeV^{-2}] = [value] x [hbar x c]^2
	// decay rate:    [value] x [GeV] = [value] / hbar
	// length:        [value] x [GeV^{-1}] = [value] x [hbar x c]

} // Namespace PDG ends


class MPDG {

public:

	MPDG()  {}
	~MPDG() {}

	void ReadParticleData(const std::string& filepath);
	
	void TokenizeProcess(const std::string& str, int depth, std::vector<gra::MDecayBranch>& branches) const;
	bool IsDecay(const std::string& str) const;
	
	void PrintPDGTable() const;
	const gra::MParticle& FindByPDG(int pdgcode) const;
	const gra::MParticle& FindByPDGName(const std::string& pdgname) const;
	
	// PDG tables
	std::map<int, gra::MParticle> PDG_table;
};


} // gra namespace ends


#endif