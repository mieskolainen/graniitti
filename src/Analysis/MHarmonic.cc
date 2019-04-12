// "Meta class" for Spherical Harmonics Expansion - contains
// dataprocessing/Minuit/ROOT plotting functions.
//
// All spherical harmonic processing is separated to MSpherical namespace.
// 
// 
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// ROOT
#include "TLatex.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TProfile.h"
#include "TStyle.h"

// Own
#include "Graniitti/Analysis/MHarmonic.h"
#include "Graniitti/MSpherical.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MPDG.h"

// Eigen
#include <Eigen/Dense>


using gra::aux::indices;
using gra::math::msqrt;
using gra::math::zi;
using gra::math::PI;
using gra::math::pow2;


namespace gra {


// Constructor
MHarmonic::MHarmonic() {
}


// Initialize
void MHarmonic::Init(const HPARAM& hp) {

	param = hp;
	NCOEF = (param.LMAX + 1)*(param.LMAX + 1);

	// ------------------------------------------------------------------
	// SETUP parameters of the fit
	ACTIVE.resize(NCOEF); // Active moments
	ACTIVENDF = 0;

	for (int l = 0; l <= param.LMAX; ++l) {
		for (int m = -l; m <= l; ++m) {
			const int index = gra::spherical::LinearInd(l, m);
			ACTIVE[index] = true;

			// FIX ODD MOMENTS TO ZERO
			if (param.REMOVEODD && ((l % 2) != 0)) {
				ACTIVE[index] = false;
			}
			// FIX NEGATIVE M TO ZERO
			if (param.REMOVENEGATIVEM && (m < 0)) {
				ACTIVE[index] = false;
			}
			if (ACTIVE[index]) {
				++ACTIVENDF;
			}
		}
	}
	
	// ------------------------------------------------------------------
	std::cout << std::endl;
	param.Print();
	std::cout << std::endl;
	// ------------------------------------------------------------------

	// Init arrays used by Minuit function
	t_lm       = std::vector<double>(NCOEF,    0.0);
	t_lm_error = std::vector<double>(NCOEF,    0.0);
	errmat     = MMatrix<double>(NCOEF, NCOEF, 0.0);
	covmat     = MMatrix<double>(NCOEF, NCOEF, 0.0);

	// Test functions
	gra::spherical::TestSphericalIntegrals(param.LMAX);

	std::cout << rang::fg::yellow <<
				"<Spherical Harmonic Based (costheta,phi)_r.f. Decomposition and Efficiency inversion>" << std::endl << std::endl;
	std::cout << "TERMINOLOGY:" << rang::fg::reset << std::endl;
	std::cout << "  {G} Generated == Events in the angular flat phase space (no cuts on final states, only on the system)" << std::endl;
	std::cout << "  {F} Fiducial  == Events within the strict fiducial (geometric-kinematic) final state phase space (cuts on final states)" << std::endl;
	std::cout << "  {S} Selected  == Events after the detector efficiency losses" << std::endl;
	std::cout << std::endl;
	std::cout << "{S} subset of {F} subset of {G} (this strict hierarchy might be violated in some special cases)" << std::endl;
	std::cout << std::endl;
	std::cout << "  The basic idea is to define {G} such that minimal extrapolation " << std::endl;
	std::cout << "  is required from the strict fiducial (geometric) phase space {F}. " << std::endl;
	std::cout << "  Flatness requirement of {G} is strictly required to represent moments in an unmixed basis (non-flat phase space <=> geometric moment mixing)." << std::endl;
	std::cout << std::endl;
	std::cout << rang::fg::yellow << "EXAMPLE OF A FORMALLY VALID DEFINITION:" << rang::fg::reset << std::endl;
	std::cout << "  G = {|Y(system)| < 0.9}" << std::endl;
	std::cout << "  F = {|eta(pi)|   < 0.9 && pt(pi) > 0.1 GeV}" << std::endl;
	std::cout << std::endl << std::endl;
	std::cout << "Note also the rotation between inactive coefficients and active one, due to moment mixing." << std::endl;
	std::cout << std::endl << std::endl;
	
	//pause(5);
}


// Destructor
MHarmonic::~MHarmonic() {
}


bool MHarmonic::PrintLoop(const std::string& output) const {

	/*
	// ------------------------------------------------------------------
	// Print out results for EXTERNAL analysis
	FILE* fp;

	char buff[100];
	snprintf(buff, sizeof(buff), "./fits/%s_%s.m", output.c_str(), FRAME.c_str());
	std::string outputfile = buff;

	fp = fopen(outputfile.c_str(), "w");

	if (fp == NULL) {
		fprintf(stderr, "Cannot open outputfile %s !\n", outputfile.c_str());
		return false;
	}

	fprintf(fp, "FRAME = %s; \n",     FRAME.c_str());
	fprintf(fp, "param.LMAX = %d; \n",      param.LMAX);
	fprintf(fp, "LAMBDA = %0.3f; \n", LAMBDA);
	fprintf(fp, "ACTIVE = [");

	for (std::size_t i = 0; i < ACTIVE.size(); ++i) {
		fprintf(fp, "%d ", (int)ACTIVE[i]);
	}
	fprintf(fp, "]; \n\n");
	fprintf(fp, "M_CELL = [");

	for (std::size_t i = 0; i < masscenter.size(); ++i) {
		fprintf(fp, "%0.3f ", masscenter[i]);
	}
	fprintf(fp, "]; \n\n");

	fprintf(fp, "A{1} = [");
	PrintMatrix(fp, ref.t_lm_EML);
	fprintf(fp, "];\n");
	fprintf(fp, "A{2} = [");
	PrintMatrix(fp, fid.t_lm_EML);
	fprintf(fp, "];\n");
	fprintf(fp, "A{3} = [");
	PrintMatrix(fp, det.t_lm_EML);
	fprintf(fp, "];\n");
	
	fprintf(fp, "A{4} = [");
	PrintMatrix(fp, ref.t_lm_MPP);
	fprintf(fp, "];\n");
	fprintf(fp, "A{5} = [");
	PrintMatrix(fp, fid.t_lm_MPP);
	fprintf(fp, "];\n");
	fprintf(fp, "A{6} = [");
	PrintMatrix(fp, det.t_lm_MPP);
	fprintf(fp, "];\n");

	fclose(fp);
	*/
	return true;
}


void MHarmonic::PlotAll(const std::string& legendstr, const std::string& outputpath) const {

	for (const auto& OBSERVABLE : {0,1,2}) {

	// **** EFFICIENCY DECOMPOSITION ****
	PlotFigures(ref, "E_lm",     OBSERVABLE, "Response[FLAT_REFERENCE]<REF>",  17, legendstr, outputpath);
	PlotFigures(fid, "E_lm",     OBSERVABLE, "Response[FIDUCIAL_ACCEPTANCE]<FID>", 33, legendstr, outputpath);
	PlotFigures(det, "E_lm",     OBSERVABLE, "Response[ACCEPTANCE_x_EFFICIENCY]<DET>", 29, legendstr, outputpath);

	// **** ALGEBRAIC INVERSE MOMENTS ****
	PlotFigures(ref, "t_lm_MPP", OBSERVABLE, "Moments[MPP]<REF>", 17, legendstr, outputpath);
	PlotFigures(fid, "t_lm_MPP", OBSERVABLE, "Moments[MPP]<FID>", 33, legendstr, outputpath);
	PlotFigures(det, "t_lm_MPP", OBSERVABLE, "Moments[MPP]<DET>", 29, legendstr, outputpath);

	// **** EXTENDED MAXIMUM LIKELIHOOD INVERSE MOMENTS ****
	if (param.EML) {
	PlotFigures(ref, "t_lm_EML", OBSERVABLE, "Moments[EML]<REF>", 17, legendstr, outputpath);
	PlotFigures(fid, "t_lm_EML", OBSERVABLE, "Moments[EML]<FID>", 33, legendstr, outputpath);
	PlotFigures(det, "t_lm_EML", OBSERVABLE, "Moments[EML]<DET>", 29, legendstr, outputpath);
	}

	}

	// 2D
	std::vector<std::vector<int>> OBSERVABLES = {{0,1},{0,2},{1,2}}; // Different pairs
	
	for (const auto& OBSERVABLE2 : OBSERVABLES) {

	// **** EFFICIENCY DECOMPOSITION ****
	PlotFigures2D(ref, "E_lm",     OBSERVABLE2, "Response[FLAT_REFERENCE]<REF>",  17, legendstr, outputpath);
	PlotFigures2D(fid, "E_lm",     OBSERVABLE2, "Response[FIDUCIAL_ACCEPTANCE]<FID>", 33, legendstr, outputpath);
	PlotFigures2D(det, "E_lm",     OBSERVABLE2, "Response[ACCEPTANCE_x_EFFICIENCY]<DET>", 29, legendstr, outputpath);

	// **** ALGEBRAIC INVERSE MOMENTS ****
	PlotFigures2D(ref, "t_lm_MPP", OBSERVABLE2, "Moments[MPP]<REF>", 17, legendstr, outputpath);
	PlotFigures2D(fid, "t_lm_MPP", OBSERVABLE2, "Moments[MPP]<FID>", 33, legendstr, outputpath);
	PlotFigures2D(det, "t_lm_MPP", OBSERVABLE2, "Moments[MPP]<DET>", 29, legendstr, outputpath);

	// **** EXTENDED MAXIMUM LIKELIHOOD INVERSE MOMENTS ****
	if (param.EML) {
	PlotFigures2D(ref, "t_lm_EML", OBSERVABLE2, "Moments[EML]<REF>", 17, legendstr, outputpath);
	PlotFigures2D(fid, "t_lm_EML", OBSERVABLE2, "Moments[EML]<FID>", 33, legendstr, outputpath);
	PlotFigures2D(det, "t_lm_EML", OBSERVABLE2, "Moments[EML]<DET>", 29, legendstr, outputpath);
	}

	}
}


// Print 1D-figures
void MHarmonic::PlotFigures(const MTensor<gra::spherical::SH>& tensor,
							const std::string& DATATYPE, unsigned int OBSERVABLE,
                            const std::string& outputfile,
                            int barcolor, const std::string& legendstr, const std::string& outputpath) const {

	TCanvas* c1 = new TCanvas("c1", "c1", 700, 600); // horizontal, vertical
	c1->Divide(sqrt(NCOEF), std::sqrt(NCOEF), 0.002, 0.001); // y, x

	// Loop over moments
	int kk = 0;
	std::vector<TGraphErrors*> gr(NCOEF, NULL);

	std::string xlabel = "";
	if      (OBSERVABLE == 0) {
		xlabel = "M (GeV)";
	}
	else if (OBSERVABLE == 1) {
		xlabel = "P_{T} (GeV)";
	}
	else if (OBSERVABLE == 2) {
		xlabel = "Y";
	} else {
		throw std::invalid_argument("MHarmonic::PlotFigures: Unknown observable " + std::to_string(OBSERVABLE));
	}

	const std::size_t BINS = tensor.size(OBSERVABLE);


	// Set max
	double maxcount = 0.0;

	for (int l = 0; l <= param.LMAX; ++l) {
		for (int m = -l; m <= l; ++m) {
			const int index = gra::spherical::LinearInd(l, m);

			// Set canvas position
			c1->cd(kk + 1);

			// Loop over mass
			double x[BINS]     = {0};
			double y[BINS]     = {0};
			double x_err[BINS] = {0};
			double y_err[BINS] = {0};

			for (std::size_t k = 0; k < BINS; ++k) {

				// Get x-axis point
				x[k] = xcenter[OBSERVABLE][k];

				// Set indices {0,0,0, ..., 0}
				std::vector<std::size_t> cell(xcenter.size(),0);
				cell[OBSERVABLE] = k;

				// CHOOSE DATATYPE
				if      (DATATYPE == "t_lm_MPP") {
					y[k]     = tensor(cell).t_lm_MPP[index];
					y_err[k] = tensor(cell).t_lm_MPP_error[index];
				}
				else if (DATATYPE == "t_lm_EML") {
					y[k]     = tensor(cell).t_lm_EML[index];
					y_err[k] = tensor(cell).t_lm_EML_error[index]; 				
				}
				else if (DATATYPE == "E_lm") {
					y[k]     = tensor(cell).E_lm[index];
					y_err[k] = tensor(cell).E_lm_error[index]; 
				} else {
					throw std::invalid_argument("MHarmonic::PlotFigures: Unknown MODE = " + DATATYPE);
				}

				// Save maximum for visualization
				if (l == 0 && m == 0) {
					maxcount = y[k] > maxcount ? y[k] : maxcount;
				}
			}
			
			// Data displayed using TGraphErrors
			gr[kk] = new TGraphErrors(BINS, x, y, x_err, y_err);

			const std::string titletxt = (DATATYPE == "E_lm") ? "MC" : param.TYPE; // Efficiency always MC 
			if (!(l == 0 && m == 0)) {
				gr[kk]->SetTitle(Form("%s: #it{lm} = <%d,%d>", titletxt.c_str(), l, m));
			} else {
				gr[kk]->SetTitle(Form("%s: %s: #it{lm} = <%d,%d>", titletxt.c_str(), outputfile.c_str(), l, m));
			}

			// Set x-axis label
			gr[kk]->GetXaxis()->SetTitle(xlabel.c_str());

			gr[kk]->SetMarkerStyle(20);
			gr[kk]->SetMarkerSize(0.05);

			gStyle->SetBarWidth(0.5);
			gr[kk]->SetFillColor(barcolor);

			gr[kk]->GetXaxis()->SetTitleSize(0.05);
			gr[kk]->GetXaxis()->SetLabelSize(0.05);

			gr[kk]->GetYaxis()->SetTitleSize(0.05);
			gr[kk]->GetYaxis()->SetLabelSize(0.05);

			gStyle->SetTitleFontSize(0.08);

			// gr->SetMarkerColor(4); // blue
			
			// Y-axis range
			if (DATATYPE != "E_lm") {
		   		if (!(l==0 && m==0)) { gr[kk]->GetHistogram()->SetMaximum(     2.0/3.0 * maxcount); }
		  		gr[kk]->GetHistogram()->SetMinimum((l == 0 && m == 0) ? 0.0 : -2.0/3.0 * maxcount);  // Y  
			} else {
				gr[kk]->GetHistogram()->SetMaximum((l == 0 && m == 0) ? 1.0 :  0.3);
		  		gr[kk]->GetHistogram()->SetMinimum((l == 0 && m == 0) ? 0.0 : -0.3);  // Y  
			}
			
			gr[kk]->Draw("A B");
			// gr[kk]->Draw("ALP");

			// Draw legend string
			TText* t = new TText(0.82, 0.85, legendstr.c_str());
			t->SetNDC();
			t->SetTextAlign(22);
			t->SetTextColor(kRed+2);
			t->SetTextFont(43);
			t->SetTextSize(14);
			//t->SetTextAngle(45);
			t->Draw("same");

			// Draw horizontal line
			TLine* line = new TLine(x[0], 0.0, x[BINS-1], 0.0);

			if (!(l == 0 && m == 0)) {
			    line->SetLineColor(kRed);
			    line->SetLineWidth(1.0);
			    line->Draw("same");
			}
			++kk;
		}
	}

	// Require that we have data
	if (BINS > 1) {
	aux::CreateDirectory("./figs");
	aux::CreateDirectory("./figs/harmonicfit");
	aux::CreateDirectory("./figs/harmonicfit/" + outputpath);
	const std::string subpath = "OBS_" + std::to_string(OBSERVABLE) + "_" + legendstr;
	aux::CreateDirectory("./figs/harmonicfit/" + outputpath + "/" + subpath);
	c1->Print(Form("./figs/harmonicfit/%s/%s/%s.pdf", outputpath.c_str(), subpath.c_str(), outputfile.c_str()));
	}

	for (std::size_t i = 0; i < gr.size(); ++i) {
		delete gr[i];
	}

	delete c1;
}


// Print 2D-figures
void MHarmonic::PlotFigures2D(const MTensor<gra::spherical::SH>& tensor,
							const std::string& DATATYPE, std::vector<int> OBSERVABLE,
                            const std::string& outputfile,
                            int barcolor, const std::string& legendstr, const std::string& outputpath) const {

	TCanvas* c1 = new TCanvas("c1", "c1", 700, 500); // horizontal, vertical
	c1->Divide(sqrt(NCOEF), std::sqrt(NCOEF), 0.002, 0.001);

	// Loop over moments
	std::vector<TH2D*> gr(NCOEF, NULL);

	std::vector<std::string> label(2);
	for (std::size_t i = 0; i < 2; ++i) {
		
		if      (OBSERVABLE[i] == 0) {
			label[i] = "M (GeV)";
		}
		else if (OBSERVABLE[i] == 1) {
			label[i] = "P_{T} (GeV)";
		}
		else if (OBSERVABLE[i] == 2) {
			label[i] = "Y";
		} else {
			throw std::invalid_argument("MHarmonic::PlotFigures2D: Unknown observable " + std::to_string(OBSERVABLE[i]));
		}
	}
	int kk = 0;

	std::vector<std::size_t> BINS = {tensor.size(OBSERVABLE[0]), tensor.size(OBSERVABLE[1])};

	for (int l = 0; l <= param.LMAX; ++l) {
		for (int m = -l; m <= l; ++m) {
			const int index = gra::spherical::LinearInd(l, m);

			// Set canvas position
			c1->cd(kk + 1);

			// Loop over mass
			double x[BINS[0]] = {0};
			double y[BINS[1]] = {0};
			double z[BINS[0]][BINS[1]] = {{0}};
			//double x_err[BINS] = {1e-2};
			//double y_err[BINS] = {1e-2};

			for (std::size_t k1 = 0; k1 < BINS[0]; ++k1) {
				x[k1] = xcenter[OBSERVABLE[0]][k1];

				for (std::size_t k2 = 0; k2 < BINS[1]; ++k2) {
					y[k2] = xcenter[OBSERVABLE[1]][k2];

					// Set indices {0,0,0, ..., 0}
					std::vector<std::size_t> cell(xcenter.size(),0);
					cell[OBSERVABLE[0]] = k1;
					cell[OBSERVABLE[1]] = k2;

					// CHOOSE DATATYPE
					if      (DATATYPE == "t_lm_MPP") {
						z[k1][k2] = tensor(cell).t_lm_MPP[index];
					}
					else if (DATATYPE == "t_lm_EML") {
						z[k1][k2] = tensor(cell).t_lm_EML[index];			
					}
					else if (DATATYPE == "E_lm") {
						z[k1][k2] = tensor(cell).E_lm[index];
					} else {
						throw std::invalid_argument("MHarmonic::PlotFigures: Unknown MODE = " + DATATYPE);
					}
				}
			}
			
			// Boundaries
			double xw = (x[1]-x[0])/2;
			double yw = (y[1]-y[0])/2;			

			gr[kk] = new TH2D(("h2_" + std::to_string(kk)).c_str(),"title",
							   BINS[0], x[0]-xw, x[BINS[0]-1]+xw,
							   BINS[1], y[0]-yw, y[BINS[1]-1]+yw);


			const std::string titletxt = (DATATYPE == "E_lm") ? "MC" : param.TYPE; // Efficiency always MC 
			if (!(l == 0 && m == 0)) {
				gr[kk]->SetTitle(Form("%s: #it{lm} = <%d,%d>", titletxt.c_str(), l, m));
			} else {
				gr[kk]->SetTitle(Form("%s: %s: #it{lm} = <%d,%d>", titletxt.c_str(), outputfile.c_str(), l, m));
			}

			// Set x-axis label
			gr[kk]->GetXaxis()->SetTitle(label[0].c_str());
			gr[kk]->GetYaxis()->SetTitle(label[1].c_str());			

			gr[kk]->GetXaxis()->SetTitleSize(0.05);
			gr[kk]->GetXaxis()->SetLabelSize(0.05);

			gr[kk]->GetYaxis()->SetTitleSize(0.05);
			gr[kk]->GetYaxis()->SetLabelSize(0.05);

			gStyle->SetTitleFontSize(0.08);
			
			// Plot
			for (std::size_t i = 0; i < BINS[0]; ++i) {
			   for (std::size_t j = 0; j < BINS[1]; ++j) {
			      gr[kk]->Fill(x[i], y[j], z[i][j]);
			   }
			}
			gr[kk]->Draw("COLZ");

			//gr[kk]->Draw("A B");

			// gr[kk]->Draw("ALP");
			// c1->Update();
			++kk;
		}
	}

	// Require that we have truly 2D data
	if (BINS[0] > 1 && BINS[1] > 1) {
	aux::CreateDirectory("./figs");
	aux::CreateDirectory("./figs/harmonicfit");
	aux::CreateDirectory("./figs/harmonicfit/" + outputpath);
	const std::string subpath = "2D_OBS_" + std::to_string(OBSERVABLE[0]) + "_OBS_" + std::to_string(OBSERVABLE[1]) + "_" + legendstr;
	aux::CreateDirectory("./figs/harmonicfit/" + outputpath + "/" + subpath);
	c1->Print(Form("./figs/harmonicfit/%s/%s/%s.pdf", outputpath.c_str(), subpath.c_str(), outputfile.c_str()));
	}

	for (std::size_t i = 0; i < gr.size(); ++i) {
		delete gr[i];
	}
	delete c1;
}


// Loop over system mass, pt, rapidity
void MHarmonic::HyperLoop(void (*fitfunc)(int&, double*, double&, double*, int),
		const std::vector<gra::spherical::Omega>& MC,
		const std::vector<gra::spherical::Omega>& DATA, const HPARAM& hp) {

	// --------------------------------------------------------
	// Pre-Calculate once Spherical Harmonics for MINUIT fit
	DATA_events = DATA;
	Y_lm = gra::spherical::YLM(DATA_events, param.LMAX);
	// --------------------------------------------------------

	double chi2 = 0.0;

	const unsigned int MASS_BINS = hp.M[0];
	const unsigned int PT_BINS   = hp.PT[0];
	const unsigned int Y_BINS    = hp.Y[0];

	ref = MTensor<gra::spherical::SH>({MASS_BINS, PT_BINS, Y_BINS});
	fid = MTensor<gra::spherical::SH>({MASS_BINS, PT_BINS, Y_BINS});
	det = MTensor<gra::spherical::SH>({MASS_BINS, PT_BINS, Y_BINS});
	
	// System Mass limits
	const double M_MIN   = hp.M[1];
	const double M_MAX   = hp.M[2];

	// System Pt limits
	const double PT_MIN  = hp.PT[1];
	const double PT_MAX  = hp.PT[2];

	// System rapidity limits
	const double Y_MIN   = hp.Y[1];
	const double Y_MAX   = hp.Y[2];

	const double  M_STEP = (M_MAX  - M_MIN)  / MASS_BINS;
	const double PT_STEP = (PT_MAX - PT_MIN) / PT_BINS;
	const double  Y_STEP = (Y_MAX  - Y_MIN)  / Y_BINS;

	// Center points
	xcenter.clear();
	xcenter.push_back(std::vector<double>(MASS_BINS, 0.0));
	xcenter.push_back(std::vector<double>(PT_BINS,   0.0));
	xcenter.push_back(std::vector<double>(Y_BINS,    0.0));


	// ------------------------------------------------------------------
	// Hyperloop
	for (std::size_t i = 0; i < MASS_BINS; ++i) {
		const double M_min = i * M_STEP + M_MIN;
		const double M_max = M_min + M_STEP;

		xcenter[0][i] = (M_max  + M_min)  / 2.0;

	for (std::size_t j = 0; j < PT_BINS;   ++j) {
		const double Pt_min = j * PT_STEP + PT_MIN;
		const double Pt_max = Pt_min + PT_STEP;

		xcenter[1][j] = (Pt_max + Pt_min) / 2.0;		

	for (std::size_t k = 0; k < Y_BINS;    ++k) {
		const double Y_min = k * Y_STEP + Y_MIN;
		const double Y_max = Y_min + Y_STEP;

		xcenter[2][k] = (Y_max  + Y_min)  / 2.0;

		// --------------------------------------------------------------
		// Get indices for this interval
		const std::vector<std::size_t> MC_ind = 
			gra::spherical::GetIndices(MC, {M_min, M_max}, {Pt_min, Pt_max}, {Y_min, Y_max});

		DATA_ind = gra::spherical::GetIndices(DATA_events, {M_min, M_max}, {Pt_min, Pt_max}, {Y_min, Y_max});

		const unsigned int MINEVENTS = 75;
		if (DATA_ind.size() < MINEVENTS) {
			std::cout << rang::fg::red << "WARNING: Less than "
					  << MINEVENTS << " in the cell!" << rang::fg::reset << std::endl;
		}

		// --------------------------------------------------------------
		// Acceptance mixing matrices
		ref({i,j,k}).MIXlm = gra::spherical::GetGMixing(MC, MC_ind, param.LMAX, 0);
		fid({i,j,k}).MIXlm = gra::spherical::GetGMixing(MC, MC_ind, param.LMAX, 1);
		det({i,j,k}).MIXlm = gra::spherical::GetGMixing(MC, MC_ind, param.LMAX, 2);
		// --------------------------------------------------------------

		// --------------------------------------------------------------
		// Get efficiency decomposition for this interval
		std::pair<std::vector<double>, std::vector<double>> E0 = gra::spherical::GetELM(MC, MC_ind, param.LMAX, 0);
		std::pair<std::vector<double>, std::vector<double>> E1 = gra::spherical::GetELM(MC, MC_ind, param.LMAX, 1);
		std::pair<std::vector<double>, std::vector<double>> E2 = gra::spherical::GetELM(MC, MC_ind, param.LMAX, 2);

		ref({i,j,k}).E_lm       = E0.first;
		ref({i,j,k}).E_lm_error = E0.second;

		fid({i,j,k}).E_lm       = E1.first;
		fid({i,j,k}).E_lm_error = E1.second;

		det({i,j,k}).E_lm       = E2.first;
		det({i,j,k}).E_lm_error = E2.second;
		// --------------------------------------------------------------


		// ==============================================================
		// ALGORITHM 1: DIRECT / OBSERVED / ALGEBRAIC decomposition

		fid({i,j,k}).t_lm_MPP = gra::spherical::SphericalMoments(DATA_events, DATA_ind, param.LMAX, 1);
		det({i,j,k}).t_lm_MPP = gra::spherical::SphericalMoments(DATA_events, DATA_ind, param.LMAX, 2);

		// Forward matrix
		Eigen::MatrixXd M     = gra::aux::Matrix2Eigen(det({i,j,k}).MIXlm);

		// Matrix pseudoinverse
		Eigen::VectorXd b     = gra::aux::Vector2Eigen(det({i,j,k}).t_lm_MPP);
		Eigen::VectorXd x     = gra::math::PseudoInverse(M, param.SVDREG) * b;

		// Collect the inversion result
		ref({i,j,k}).t_lm_MPP = gra::aux::Eigen2Vector(x);

		// ==============================================================


		// Update these to proper errors (TBD. PUT ZERO FOR NOW!!)
		ref({i,j,k}).t_lm_MPP_error = std::vector<double>(ref({i,j,k}).t_lm_MPP.size(), 0.0);
		fid({i,j,k}).t_lm_MPP_error = std::vector<double>(fid({i,j,k}).t_lm_MPP.size(), 0.0);
		det({i,j,k}).t_lm_MPP_error = std::vector<double>(det({i,j,k}).t_lm_MPP.size(), 0.0);

		
		std::cout << "Algebraic Moore-Penrose/SVD inverted (unmixed) moments in the (angular flat) reference phase space:" << std::endl;
		gra::spherical::PrintOutMoments(ref({i,j,k}).t_lm_MPP, ref({i,j,k}).t_lm_MPP_error, ACTIVE, param.LMAX);

		std::cout << "Algebraic (mixed) moments in the fiducial phase space:" << std::endl;
		gra::spherical::PrintOutMoments(fid({i,j,k}).t_lm_MPP, fid({i,j,k}).t_lm_MPP_error, ACTIVE, param.LMAX);

		std::cout << "Algebraic (mixed) moments in the detector space:" << std::endl;
		gra::spherical::PrintOutMoments(det({i,j,k}).t_lm_MPP, det({i,j,k}).t_lm_MPP_error, ACTIVE, param.LMAX);


		// ==============================================================
		// ALGORITHM 2: Extended Maximum Likelihood Fit

		if (param.EML) {

		// *** Get TMinuit based fit decomposition ***
		MomentFit({i,j,k}, fitfunc);

		// Save result
		ref({i,j,k}).t_lm_EML = t_lm;
		fid({i,j,k}).t_lm_EML = fid({i,j,k}).MIXlm * t_lm; // Moment forward rotation: Matrix * Vector
		det({i,j,k}).t_lm_EML = det({i,j,k}).MIXlm * t_lm; // Moment forward rotation: Matrix * Vector

		// Uncertanties
		ref({i,j,k}).t_lm_EML_error = t_lm_error;
		fid({i,j,k}).t_lm_EML_error = fid({i,j,k}).MIXlm * t_lm_error; // Simple rotation for now (should vector Taylor expand)
		det({i,j,k}).t_lm_EML_error = det({i,j,k}).MIXlm * t_lm_error; // Simple rotation for now (should vector Taylor expand)
		// ==============================================================


		std::cout << "Extended Maximum-Likelihood inverse fitted (unmixed) moments in the (angular flat) reference phase space:" << std::endl;
		gra::spherical::PrintOutMoments(ref({i,j,k}).t_lm_EML, ref({i,j,k}).t_lm_EML_error, ACTIVE, param.LMAX);

		std::cout << "Extended Maximum-Likelihood Re-Back-Projected (mixed) moments in the fiducial phase space:" << std::endl;
		gra::spherical::PrintOutMoments(fid({i,j,k}).t_lm_EML, fid({i,j,k}).t_lm_EML_error, ACTIVE, param.LMAX);

		std::cout << "Extended Maximum-Likelihood Re-Back-Projected (mixed) moments in the detector space:" << std::endl;
		gra::spherical::PrintOutMoments(det({i,j,k}).t_lm_EML, det({i,j,k}).t_lm_EML_error, ACTIVE, param.LMAX);
		
		}

		// --------------------------------------------------------------
		// Print results
		chi2 += PrintOutHyperCell({i,j,k});

		// --------------------------------------------------------------
		// Make comparison of synthetic MC data vs estimate
		
		if (param.TYPE == "MC") {

		int fiducial = 0;
		int selected = 0;
		for (const auto& l : DATA_ind) {
			if (DATA_events[l].fiducial){ ++fiducial; }
			if (DATA_events[l].fiducial && DATA_events[l].selected) { ++selected; }
		}
		std::cout << std::endl;

		std::cout << rang::fg::yellow << "MC GROUND TRUTH: " << rang::fg::reset << std::endl;
		printf("   MC 'synthetic data' events generated               = %u \n", (unsigned int) DATA_ind.size());
		printf("   MC 'synthetic data' events fiducial                = %d (acceptance %0.1f percent) \n",
			fiducial, fiducial / (double) DATA_ind.size() * 100);
		printf("   MC 'synthetic data' events fiducial and selected   = %d (efficiency %0.1f percent) \n", 
			selected, selected / (double) fiducial * 100);
		std::cout << std::endl;
		}
	}}}

	if (param.EML) {
	gra::aux::PrintBar("=");
	double reducedchi2 = chi2 / (double)(ACTIVENDF * MASS_BINS * PT_BINS * Y_BINS);
	if (reducedchi2 < 3) {
		std::cout << rang::fg::green;
	} else {
		std::cout << rang::fg::red;
	}
	printf("Total chi2(MPP - EML) / (ACTIVENDF x BINS) = %0.2f / (%d x %d) = %0.2f \n", 
		    chi2, ACTIVENDF, MASS_BINS * PT_BINS * Y_BINS, reducedchi2);
	std::cout << rang::fg::reset;
	gra::aux::PrintBar("=");
	std::cout << std::endl << std::endl;
	}
}


// Print out results for the hypercell
double MHarmonic::PrintOutHyperCell(const std::vector<std::size_t>& cell) {
	
	double chi2 = 0;

	// Extended Maximum Likelihood
	if (param.EML == true) {

	// Loop over moments
	for (int l = 0; l <= param.LMAX; ++l) {
		for (int m = -l; m <= l; ++m) {
			const  int index = gra::spherical::LinearInd(l, m);
			
			const double obs = det(cell).t_lm_MPP[index];
			const double fit = det(cell).t_lm_EML[index];
			
			// Moment active
			if (ACTIVE[index]) {
				chi2 += pow2(obs - fit) / pow2(obs);
			}
		}
	}
	std::cout << std::endl;
	
	const double reducedchi2 = chi2 / (double) ACTIVENDF;
	if (reducedchi2 < 3) {
		std::cout << rang::fg::green;
	} else {
		std::cout << rang::fg::red;
	}
	printf("chi2(MPP - EML) / ndf = %0.3f / %d = %0.3f \n", chi2, ACTIVENDF, reducedchi2);
	std::cout << rang::fg::reset << std::endl;

	const double sum_ref = ref(cell).t_lm_EML[gra::spherical::LinearInd(0,0)];
	printf("EML: Estimate of events in this hyperbin in the reference phase space = %0.1f +- %0.1f \n", sum_ref, 0.0);

	const double sum_FID = gra::spherical::HarmDotProd(fid(cell).E_lm, ref(cell).t_lm_EML, ACTIVE, param.LMAX);
	printf("EML: Estimate of events in this hyperbin in the fiducial  phase space = %0.1f +- %0.1f \n", sum_FID, 0.0);

	const double sum_DET = gra::spherical::HarmDotProd(det(cell).E_lm, ref(cell).t_lm_EML, ACTIVE, param.LMAX);
	printf("EML: Estimate of events in this hyperbin in the detector        space = %0.1f +- %0.1f \n", sum_DET, 0.0);

	}
	std::cout << std::endl;

	// Algebraic inverse
	const double sum_ref = ref(cell).t_lm_MPP[gra::spherical::LinearInd(0,0)];
	printf("MPP: Estimate of events in this hyperbin in the reference phase space = %0.1f +- %0.1f \n", sum_ref, 0.0);

	const double sum_FID = gra::spherical::HarmDotProd(fid(cell).E_lm, ref(cell).t_lm_MPP, ACTIVE, param.LMAX);
	printf("MPP: Estimate of events in this hyperbin in the fiducial  phase space = %0.1f +- %0.1f \n", sum_FID, 0.0);

	const double sum_DET = gra::spherical::HarmDotProd(det(cell).E_lm, ref(cell).t_lm_MPP, ACTIVE, param.LMAX);
	printf("MPP: Estimate of events in this hyperbin in the detector        space = %0.1f +- %0.1f \n", sum_DET, 0.0);

	return chi2;
}


// MINUIT based fit routine for the Extended Maximum Likelihood formalism
void MHarmonic::MomentFit(const std::vector<std::size_t>& cell,
						  void (*fitfunc)(int&, double*, double&, double*, int)) {

	std::cout << "MomentFit: Starting ..." << std::endl << std::endl;

	 // **** This must be set for the loss function ****
	activecell = cell;

	// Init TMinuit
	TMinuit* gMinuit = new TMinuit(NCOEF); // initialize TMinuit with a maximum of N params
	gMinuit->SetFCN(fitfunc);

	// Set Print Level
	// -1 no output
	// 1 standard output
	gMinuit->SetPrintLevel(-1);

	// Set error Definition
	// 1 for Chi square
	// 0.5 for negative log likelihood
	//    gMinuit->SetErrorDef(0.5);
	double arglist[2];
	int ierflg = 0;
	arglist[0] = 0.5; // 0.5 <=> We use negative log likelihood cost function
	gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

	//double fminbest = 1e32;
	// Try different initial values to find out true minimum
	const std::size_t TRIALMAX = 1;

	for (std::size_t trials = 0; trials < TRIALMAX; ++trials) {

		for (int l = 0; l <= param.LMAX; ++l) {
			for (int m = -l; m <= l; ++m) {
				const int index = gra::spherical::LinearInd(l, m);

				const std::string str = "t_" + std::to_string(l) + std::to_string(m);

				const double max = 1e9;  // Physically bounded ultimately by the number of events
				const double min = -1e5;

				// ======================================================
				// *** Use the algebraic inverse solution as a good starting value ***
				const double start_value = ref(activecell).t_lm_MPP[index];
				// ======================================================
				
				const double step_value = 0.1; // in units of Events

				gMinuit->mnparm(index, str, start_value, step_value, min, max, ierflg);

				// After first trial, fix t_00 (error are not estimated with constant parameters)
				//if (trials > 0) {
				//	gMinuit->mnparm(0, "t_00", t_lm[0], 0, 0, 0, ierflg);
				//	gMinuit->FixParameter(0);
				//}
				// FIX ODD MOMENTS TO ZERO
				if (param.REMOVEODD && ((l % 2) != 0)) {
					gMinuit->mnparm(index, str, 0, 0, 0, 0, ierflg);
					gMinuit->FixParameter(index);
				}
				// FIX NEGATIVE M TO ZERO
				if (param.REMOVENEGATIVEM && (m < 0)) {
					gMinuit->mnparm(index, str, 0, 0, 0, 0, ierflg);
					gMinuit->FixParameter(index);
				}
			}
		}
		// Scan main parameter
		// gMinuit->mnscan();
		arglist[0] = 100000;  // Minimum number of function calls
		arglist[1] = 0.00001; // Minimum tolerance

		// First simplex to find approximate answers
		gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);

		// Numerical Hessian (2nd derivatives matrix), inverse of this -> covariance matrix
		gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

		// Confidence intervals based on the profile likelihood ratio
		gMinuit->mnexcm("MINOS", arglist, 2, ierflg);

        // Calculate error & covariance matrix
        gMinuit->mnemat(&errmat[0][0], NCOEF);

        for (int i = 0; i < NCOEF; ++i) {
            for (int j = 0; j < NCOEF; ++j){
                covmat[i][j] = sqrt(errmat[i][i]*errmat[j][j]);
	            if (covmat[i][j] > 1E-80) {
	                covmat[i][j] = errmat[i][j]/covmat[i][j];
	            }
	            else covmat[i][j] = 0.0;
            }
        }
       	errmat.Print("Error matrix");
        covmat.Print("Covariance matrix");

		// Print results
		double fmin, fedm, errdef = 0.0;
		int nvpar, nparx, istat = 0;
		gMinuit->mnstat(fmin, fedm, errdef, nvpar, nparx, istat);

		// Collect fit result
		for (int l = 0; l <= param.LMAX; ++l) {
			for (int m = -l; m <= l; ++m) {
				const int index = gra::spherical::LinearInd(l, m);

				// Set results into t_lm[], t_lm_error[]
				gMinuit->GetParameter(index, t_lm[index], t_lm_error[index]);
			}
		}

	} // TRIALS LOOP END

	// Sum of N random variables Y = X_1 + X_2 + ... X_N
	// sigma_Y^2 = \sum_i^N \sigma_i^2 + 2 \sum_i \sum_i < j
	// covariance(X_i, X_j)

	// Print out results (see MINUIT manual for these parameters)
	double fmin, fedm, errdef = 0.0;
	int nvpar, nparx, istat = 0;
	gMinuit->mnstat(fmin, fedm, errdef, nvpar, nparx, istat);
	gMinuit->mnprin(4, fmin);

	std::cout << "<MINUIT (MIGRAD+MINOS)> done." << std::endl << std::endl;

	delete gMinuit;
}


// Unbinned Extended Maximum Likelihood function
// Extended means that the number of events (itself) is a Poisson distributed
// random variable and that is incorporated to the fit.
void MHarmonic::logLfunc(int& npar, double* gin, double& f, double* par, int iflag) const {

	// Collect fit t_LM coefficients from MINUIT
	std::vector<double> T(NCOEF, 0.0);

	for (const auto& i : indices(T)) {
		T[i] = par[i];
		// printf("T[%d] = %0.5f \n", i);
	}
	// This is the number of events in the current phase space point at detector level
	const int METHOD = 2;
	double nhat = 0.0;
	
	// Equivalent estimator 1
	if (METHOD == 1) {
		const std::vector<double> t_lm_det = det(activecell).MIXlm * T; // Matrix * Vector
		nhat = t_lm_det[0];
	}
	// Equivalent estimator 2
	if (METHOD == 2) {
		nhat = gra::spherical::HarmDotProd(det(activecell).E_lm, T, ACTIVE, param.LMAX);
	}


	// For each event, calculate \sum_{LM} t_{lm}
	// Re[Y_{lm}(costheta,phi_k)], k is the event index
	std::vector<double> I0;
	const double V = msqrt(4.0 * PI); // Normalization volume

	for (const auto& k : DATA_ind) {

		// Event is accepted
		if ( DATA_events[k].fiducial && DATA_events[k].selected ) {
			// fine
		} else { continue; }

		// Loop over (l,m) terms
		double sum = 0.0;
		for (int l = 0; l <= param.LMAX; ++l) {
			for (int m = -l; m <= l; ++m) {
				const int index = gra::spherical::LinearInd(l, m);

				if (ACTIVE[index]) {

					// Calculate here
            		//const std::complex<double> Y =
            		//	gra::math::Y_complex_basis(DATA_events[k].costheta,DATA_events[k].phi, l,m);
            		//const double ReY = gra::math::NReY(Y,l,m);
            		
            		// Pre-calculated for speed
            		const double ReY = Y_lm[k][index];

            		// Add
            		sum += T[index] * ReY;
				}
			}
		}
		I0.push_back(V * sum);
	}
	// Fidelity term
	double logL = 0;

	for (const auto& k : indices(I0)) {

		if (I0[k] > 0) {
			logL += std::log(I0[k]);
			// Positivity is not satisfied, make strong penalty
		} else {
			logL = -1e32;
			// The Re[Y_lm] (purely real) formulation here does
			// not, explicitly, enforce non-negativity.
			// Where as |A|^2 formulation would enforce it by
			// construction,
			// but that would have -l <= m <= l parameters, i.e.,
			// over redundant representation problem
			// on the other hand.
		}
	}
	
	// Note the minus sign on logL term
	f = -logL + nhat;

	// High penalty, total event count estimate has gone out of physical
	if (nhat < 0){
		f = 1e32;
	}

	// L1-norm regularization (Laplace prior), -> + ln(P_laplace)
	double l1term = 0.0;
	const unsigned int START = 1;
	for (std::size_t i = START; i < T.size(); ++i) {
		l1term += std::abs(T[i]);
	}
	f += l1term * param.L1REG * nhat; // nhat for scale normalization

	// printf("MHarmonic:: cost-functional = %0.4E, nhat = %0.1f \n", f, nhat);
}

} // gra namespace ends
