// Shuvaev PDF and Sudakov suppression class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MSUDAKOV_H
#define MSUDAKOV_H

// C++
#include <complex>
#include <memory>
#include <vector>

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MTimer.h"


namespace gra {


// Interpolation container
class IArray2D {

public:
	IArray2D() {};
	~IArray2D() {};
	
	std::string name[2];
	double    MIN[2] = {0};
	double    MAX[2] = {0};
	uint        N[2] = {0};
	double   STEP[2] = {0};
	bool    islog[2] = {false};
	
	// Setup discretization
	void Set(uint VAR, std::string _name, double _min, double _max, double _N, double _logarithmic) {

		// Out of index
		if (VAR > 1) { throw std::invalid_argument("IArray2D::Set: Error: VAR = " + std::to_string(VAR) + " > 1"); }

		// AT least two intervals
		if (_N < 2) { throw std::invalid_argument("IArray2D::Set: Error: N = " + std::to_string(_N) + " < 2"); }

		// Sanity
		if (_min >= _max) {
			throw std::invalid_argument("IArray2D::Set: Error: Variable " + _name + 
				" MIN = " + std::to_string(_min) + " >= MAX = " + std::to_string(_max));
		}

		// Sanity
		if (_logarithmic && _min < 1e-9) {
			throw std::invalid_argument("IArray2D::Set: Error: Variable " + _name + 
				" is using logarithmic stepping with boundary MIN = " + std::to_string(_min));
		}
		islog[VAR]  = _logarithmic;

		// Logarithm taken here!
		MIN[VAR]    = islog[VAR] ? std::log(_min) : _min;
		MAX[VAR]    = islog[VAR] ? std::log(_max) : _max;
		N[VAR]      = _N;
		
		STEP[VAR]   = (MAX[VAR] - MIN[VAR]) / N[VAR];
		name[VAR]   = _name;
	}
	
	// Call this last
	void InitArray() {
		// Note N+1 !
		F = std::vector<std::vector<std::vector<double>>>(N[0] + 1,
		    std::vector<std::vector<double>>(N[1] + 1,
		    std::vector<double>(4, 0.0)));
	}
	
	// N+1 per dimension!
	std::vector<std::vector<std::vector<double>>> F;

	// CMS energy (needed for boundary conditions)
	double sqrts = 0.0;
	
	std::string GetHashString() const {
		std::string str =
			std::to_string(islog[0]) + std::to_string(islog[1]) +
		    std::to_string(MIN[0])   + std::to_string(MIN[1]) +
		    std::to_string(MAX[0])   + std::to_string(MAX[1]) +
		    std::to_string(N[0])     + std::to_string(N[1])   +
		    std::to_string(sqrts);
		return str;
	}

	bool WriteArray(const std::string& filename, bool overwrite) const;
	bool ReadArray(const std::string& filename);
	std::pair<double,double> Interpolate2D(double A, double B) const;
};


// Sudakov suppression and skewed pdf
class MSudakov {

public:
	MSudakov();
	~MSudakov();

	void Init(double sqrts_in, const std::string& PDFSET, bool init_arrays = true);
	void InitLHAPDF(const std::string& PDFSET);
	std::pair<double,double> Shuvaev_H(double q2, double x);
	std::pair<double,double> Sudakov_T(double qt2, double M);

	double fg_xQ2M(double x, double q2, double M) const;
	double AlphaS_Q2(double q2) const;
	double NumFlavor(double q2) const;
	double xg_xQ2(double x, double Q2) const;
	void TestPDF() const;
	
	bool initialized = false;
	
private:
	
	void   InitArrays();
	double diff_xg_xQ2_wrt_Q2(double x, double q2) const;
	double AP_gg(double delta) const;
	double AP_qg(double delta, double qt2) const;
	void   CalculateArray(IArray2D& arr, std::pair<double,double> (MSudakov::*f)(double,double));

	std::string PDFSETNAME;
	std::unique_ptr<LHAPDF::PDF> PdfPtr;

	IArray2D veto;      // Sudakov veto
	IArray2D spdf;      // Shuvaev pdf
};


} // gra namespace ends


#endif