# Python pre-processor of MadGraph C++ output to GRANIITTI

# Set here the new process name
PROCESS = input("Give the name of the process folder: ")
#PROCESS = 'gg_gg'


# ------------------------------------------------------------------------

CLASSNAME = 'AMP_MG5_' + PROCESS;

message_orig = '// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch'
message_new  = '// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch \n// @@@@ MadGraph to GRANIITTI autoconversion done @@@@'


# ------------------------------------------------------------------------
# Process the header file

fin  = open('./' + PROCESS + '/' + 'CPPProcess.h', 'rt')
fout = open('./output/' + CLASSNAME + '.h', 'wt')

for line in fin:
	
	line = line.replace('virtual','')
	line = line.replace(message_orig, message_new)

	line = line.replace('#include "Parameters_sm.h', 
		'#include "Graniitti/Amplitude/Parameters_sm.h" \n#include "Graniitti/MForm.h" \n#include "Graniitti/MAux.h" \n#include "Graniitti/MKinematics.h" \n#include "Graniitti/MMath.h')
	
	line = line.replace('CPPProcess', CLASSNAME)
	line = line.replace('Parameters_sm * pars;', 'Parameters_sm pars; // GRANIITTI')
	line = line.replace('void sigmaKin()', 'double CalcAmp2(gra::LORENTZSCALAR &lts, double alphas)')
	
	line = line.replace(CLASSNAME + '() {}', CLASSNAME + '() {initProc(gra::aux::GetBasePath(2) + "/MG5cards/param_card_' + PROCESS + '.dat"); }');
	
	fout.write(line)

fin.close()
fout.close()

print('Output to ' + CLASSNAME + '.h created')

# ------------------------------------------------------------------------
# Process the source file

fin  = open('./' + PROCESS + '/' + 'CPPProcess.cc', 'rt')
fout = open('./output/' + CLASSNAME + '.cc', 'wt')


for line in fin:

	line = line.replace(message_orig, message_new)

	line = line.replace('#include "HelAmps_sm.h', '#include "Graniitti/Amplitude/HelAmps_sm.h')
	line = line.replace('#include "CPPProcess.h', '#include "Graniitti/Amplitude/CPPProcess.h')
	
	line = line.replace('void CPPProcess::sigmaKin()', 'double CPPProcess::CalcAmp2(gra::LORENTZSCALAR &lts, double alphas)')
	line = line.replace('CPPProcess', CLASSNAME)
	
	line = line.replace('pars = Parameters_sm::getInstance();', 'pars = Parameters_sm(); // GRANIITTI')
	line = line.replace('pars->', 'pars.')
	line = line.replace('static bool firsttime = true;', 'static bool firsttime = false; // GRANIITTI')
	line = line.replace('ntry = ntry + 1;','ntry = 1; // GRANIITTI')

	line = line.replace('setDependentParameters()', 'setDependentParameters(alphas)')	

	s = line.find('printIndependentParameters()')
	if s != -1:
		line = '// pars.printIndependentParameters(); // GRANIITTI \n'

	s = line.find('printIndependentCouplings()')
	if s != -1:
		line = '// pars.printIndependentCouplings(); // GRANIITTI \n'

	s = line.find('setDependentCouplings()')
	if s != -1:
		extra = 'pars.setAlphaQEDZero(); // GRANIITTI: Set at scale alpha_QED(Q2=0) \n\n'
		fout.write(extra)

	s = line.find('if (sum_hel == 0')
	if s != -1:
		extra = 'goto SKIPLABEL; // GRANIITTI: Skip this block \n'
		fout.write(extra)
	
	fout.write(line)

	# ----------------------------------------------------------------
	# Find the first lines of the main calculation function
	s = line.find('// Local variables and constants')
	blockline = """
	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	// *** MADGRAPH CONVENTION IS [E,px,py,pz] ! ***

	// *** Set masses for HELAS ***
	std::vector<double> masses = {0, 0}; // Massless two initial states
	std::vector<gra::M4Vec> pf;

	for (std::size_t i = 0; i < lts.decaytree.size(); ++i) {
		masses.push_back(lts.decaytree[i].p4.M());
		pf.push_back(lts.decaytree[i].p4);
	}
	mME = masses;
	
	gra::M4Vec p1_ = lts.q1;
	gra::M4Vec p2_ = lts.q2;
		
	
	// Do kinematic transform
	gra::kinematics::OffShell2LightCone(p1_, p2_, pf);

	// Set initial state 4-momentum
	p.clear();
	double p1[] = {p1_.E(), p1_.Px(), p1_.Py(), p1_.Pz()}; 
	p.push_back(&p1[0]);

	double p2[] = {p2_.E(), p2_.Px(), p2_.Py(), p2_.Pz()};
	p.push_back(&p2[0]);

	// Set final state 4-momentum
	double pthis[pf.size()][4];
	for (std::size_t i = 0; i < pf.size(); ++i) {
		pthis[i][0] = pf[i].E();
		pthis[i][1] = pf[i].Px();
		pthis[i][2] = pf[i].Py();
		pthis[i][3] = pf[i].Pz();
		p.push_back(&pthis[i][0]);
	}
  	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	"""
	if s != -1:
		fout.write(blockline)
		continue;

	# ----------------------------------------------------------------
	# Number of combinations
	s = line.find('const int ncomb =')

	if s != -1:
		extra = 'lts.hamp = std::vector<std::complex<double>>(ncomb); // GRANIITTI \n\n'
		fout.write(extra)
		continue
	# ----------------------------------------------------------------


	# Find the last line of the main calculation function

	s = line.find('matrix_element[i] /= denominators[i];')
	blockline = """
	SKIPLABEL:
	
	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  	// Define permutation
  	for (int i = 0; i < nexternal; ++i) { perm[i] = i; }

  	// Loop over helicity combinations
  	for (int ihel = 0; ihel < ncomb; ++ihel) {
    		calculate_wavefunctions(perm, helicities[ihel]);

    		// Sum of subamplitudes (s,t,u,...)
    		for (int k = 0; k < namplitudes; ++k) { lts.hamp[ihel] += amp[k]; }
	}

	// Total amplitude squared over all helicity combinations individually
	double amp2 = 0.0;
	for (int ihel = 0; ihel < ncomb; ++ihel) {
	amp2 += gra::math::abs2(lts.hamp[ihel]);
	}
	amp2 /= denominators[0];  // spin average matrix element squared

	return amp2;  // amplitude squared
	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	"""
		
	if s != -1:
		fout.write(blockline)

fin.close()
fout.close()

print('Output to ' + CLASSNAME + '.cc created')


# ------------------------------------------------------------------------
# Copy parameter card

import os
os.system('cp ./' + PROCESS + '/param_card.dat ./output/' + 'param_card_' + PROCESS +'.dat')
print('Output to param_card_' + PROCESS + '.dat created')


os.system('cp ./output/*.h ../../include/Graniitti/Amplitude/')
os.system('cp ./output/*.cc ../../src/Amplitude/')
os.system('cp ./output/*.dat ../../MG5cards/')





