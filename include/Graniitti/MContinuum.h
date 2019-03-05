// Continuum 2->N phase space class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MCONTINUUM_H
#define MCONTINUUM_H

// C++
#include <complex>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/M4Vec.h"

// HepMC3
#include "HepMC/GenEvent.h"


namespace gra {

class MContinuum : public MProcess {

public:
	MContinuum();
	MContinuum(std::string process);
	virtual ~MContinuum();

	void post_Constructor();

	double operator() (const std::vector<double>& randvec, AuxIntData& aux) {
		return EventWeight(randvec, aux);
	}
	double EventWeight(const std::vector<double>& randvec, AuxIntData& aux);
	bool EventRecord(HepMC::GenEvent& evt);
	void PrintInit(bool silent) const;

private:

	bool LoopKinematics(const std::vector<double>& p1p,
	                      const std::vector<double>& p2p);
	bool FiducialCuts() const;
	void ConstructProcesses();
	
	// 3*N-4 dimensional phase space, 2->N
	bool BNRandomKin(uint Nf, const std::vector<double>& randvec);
	bool BNBuildKin(uint Nf, double pt1, double pt2, double phi1, double phi2,
	              const std::vector<double>& kt, const std::vector<double>& phi,
	              const std::vector<double>& y);

	void BLinearSystem(std::vector<M4Vec>& p,
		const std::vector<M4Vec>& q,
		const M4Vec& p1f,
		const M4Vec& p2f) const;

	double BNIntegralVolume() const;
	double BNPhaseSpaceWeight() const;
	
	// Auxialary (kt) vectors
	std::vector<M4Vec> pkt_;
};

} // gra namespace ends

#endif