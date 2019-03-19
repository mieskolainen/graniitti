// Simple parton model 2->2 phase space class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MPARTON_H
#define MPARTON_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/M4Vec.h"

// HepMC33
#include "HepMC3/GenEvent.h"


namespace gra {

class MParton : public MProcess {

public:
	MParton();
	MParton(std::string process, const std::vector<aux::OneCMD>& syntax);
	virtual ~MParton();

	void   post_Constructor();

	double operator() (const std::vector<double>& randvec, AuxIntData& aux) {
		return EventWeight(randvec, aux);
	}
	double EventWeight(const std::vector<double>& randvec, AuxIntData& aux);
	bool   EventRecord(HepMC3::GenEvent& evt);
	void   PrintInit(bool silent) const;
	
private:

	// Internal
	bool   LoopKinematics(const std::vector<double>& p1p,
	                      const std::vector<double>& p2p);
	bool   FiducialCuts() const;
	void   ConstructProcesses();

	// 2->2 dim phase space
	bool   B2RandomKin(const std::vector<double>& randvec);
	bool   B2BuildKin(double xbj1, double xbj2);
	void   B2RecordEvent(HepMC3::GenEvent& evt);
	
	double B2IntegralVolume() const;
	double B2PhaseSpaceWeight() const;

	void   DecayWidthPS(double& exact) const;
	
};


} // gra namespace ends


#endif
