// 2D histogram class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MH2_H
#define MH2_H

// C++
#include <complex>
#include <vector>

// Own
#include "Graniitti/MMatrix.h"


namespace gra {

class MH2 {

public:
	MH2(int xbins, double xmin, double xmax,
		int ybins, double ymin, double ymax, std::string namestr = "noname");
	MH2();
	~MH2();

	void ResetBounds(int xbins, double xmin, double xmax, 
					 int ybins, double ymin, double ymax);
	
	void   Fill(double xvalue, double yvalue);
	void   Fill(double xvalue, double yvalue, double weight);
	void   Clear();
	std::pair<double,double> WeightMeanAndError() const;
	double SumWeights() const;
	double SumWeights2() const;
	long long int SumBinCounts() const;
	long long int FillCount() const {
		return fills;
	}

	double        GetBinWeight(int xbin, int ybin) const;
	long long int GetBinCount(int xbin, int ybin)  const;
	void   GetBinIdx(double xvalue, double yvalue, int& xbin, int& ybin) const;
	double GetMaxWeight() const;
	double GetMinWeight() const;
	void   Print() const;
	double ShannonEntropy() const;
	
	// Set logarithmic binning
	// (user needs to take care that XMIN and XMAX > 0)
	void SetLogX() {
		LOGX = true;
	}
	void SetLogY() {
		LOGX = true;
	}
	void SetLogXY() {
		LOGX = true;
		LOGY = true;
	}

	// Overload + operator to add two histograms
	MH2 operator+(const MH2& rhs) {

		if ((this->XBINS != rhs.XBINS) || (this->YBINS != rhs.YBINS)) {
			throw std::domain_error("MH2 + operator: Histograms with different number of bins");
		}

		MH2 h(this->XBINS, this->XMIN, this->XMAX, this->YBINS, this->YMIN, this->YMAX, this->name);

		h.fills     = this->fills + rhs.fills;
		h.underflow = {this->underflow[0] + rhs.underflow[0], this->underflow[1] + rhs.underflow[1]};
		h.overflow  = {this->overflow[0]  + rhs.overflow[0],  this->overflow[1]  + rhs.overflow[1]};

		// DATA
		h.weights  = this->weights + rhs.weights;
		h.weights2 = this->weights2 + rhs.weights2;
		h.counts   = this->counts + rhs.counts;
		
		return h;
	}

private:

   	std::string name; // Histogram name

	// Boundary conditions
	double XMIN = 0.0;
	double XMAX = 0.0;
	int   XBINS = 0;

	double YMIN = 0.0;
	double YMAX = 0.0;
	int   YBINS = 0;
	
	// Number of underflow and overflow counts
	long long int fills                  = 0;
	std::vector<long long int> overflow  = {0,0};
	std::vector<long long int> underflow = {0,0};
	
	// Logarithmic binning
	bool LOGX = false;
	bool LOGY = false;
	bool ValidBin(int xbin, int ybin) const;
	int  GetIdx(double value, double minval, double maxval, int nbins, bool logbins) const;
	
	// weights (in unweighted case weights = counts)
	MMatrix<double> weights;
	MMatrix<double> weights2;

	// counts
	MMatrix<long long int> counts;
};


} // gra namespace ends

#endif