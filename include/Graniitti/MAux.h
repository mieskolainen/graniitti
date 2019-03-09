// I/O aux functions
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MAUX_H
#define MAUX_H


// C++
#include <regex>
#include <complex>
#include <mutex>
#include <random>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <unistd.h>
#include <limits.h>
//#include <experimental/filesystem>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MSudakov.h"
#include "Graniitti/MMatrix.h"

// HepMC3
#include "HepMC/FourVector.h"

// Eigen
#include <Eigen/Dense>

// Libraries
#include "rang.hpp"


namespace gra {
namespace aux {

// ----------------------------------------------------------------------
// "PROGRAM GLOBALS"

// Model parameters
extern std::string MODELPARAM;

// Sudakov routines
extern MSudakov sudakov;

// GLOBAL for multithreading: HepMC3 outputfile and LHAPDF and mutex lock
extern std::mutex g_mutex;

// For multithreaded VEGAS, to handle the exceptions from forked threads
static std::exception_ptr globalExceptionPtr = nullptr;

// ----------------------------------------------------------------------

// M4Vec to HepMC::FourVector
inline HepMC::FourVector M4Vec2HepMC(const M4Vec& v) {
	return HepMC::FourVector(v.X(),v.Y(),v.Z(), v.E());
}
// HepMC::FourVector to M4Vec
inline M4Vec HepMC2M4Vec(const HepMC::FourVector& v) {
	return M4Vec(v.x(),v.y(),v.z(), v.e());
}

// Eigen to std::vector
inline std::vector<double> Eigen2Vector(const Eigen::VectorXd& x) {

	std::vector<double> y(x.size());
	for (int i = 0; i < x.size(); ++i) {
		y[i] = x[i];
	}
	return y;
}

// std::vector to Eigen
template <typename T>
inline Eigen::VectorXd Vector2Eigen(const std::vector<T>& x) {

	Eigen::VectorXd y(x.size());
	for (std::size_t i = 0; i < x.size(); ++i) {
		y[i] = x[i];
	}
	return y;
}

// MMatrix to Eigen matrix
template <typename T>
inline Eigen::MatrixXd Matrix2Eigen(const MMatrix<T>& M) {

	Eigen::MatrixXd eM(M.size_row(), M.size_col());
	for (std::size_t i = 0; i < M.size_row(); ++i) {
		for (std::size_t j = 0; j < M.size_col(); ++j) {
			eM(i,j) = M[i][j];
		}
	}
	return eM;
}

// System information
std::string GetExecutablePath();
std::string GetBasePath(std::size_t level);

/*
std::string GetCurrentPath();
bool FileExist(const std::experimental::filesystem::path& p,
	std::experimental::filesystem::file_status s = std::experimental::filesystem::file_status{});
*/

std::uintmax_t GetFileSize(const std::string& filename);
void GetProcessMemory(double& peak_use, double& resident_use);
void GetDiskUsage(const std::string& path, int64_t& size, int64_t& free, int64_t& used);
unsigned long long TotalSystemMemory();
std::string SystemName();
std::string HostName();
const std::string DateTime();

// Progress bar
void PrintProgress(double ratio);
void ClearProgress();

// djb2hash function
unsigned long djb2hash(const std::string& s);

// Simple CSV reader
void ReadCSV(const std::string& inputfile,
             std::vector<std::vector<std::string>>& output);

// Input processing
std::string GetInputData(const std::string& inputfile);

bool IsIntegerDigits(const std::string& str);

// Trim extra spaces of a string
void TrimExtraSpace(std::string& value);

// Number to string with formatting
template <typename T>
std::string ToString(const T value, const uint n = 6) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << value;
    return out.str();
}


// Quantum numbers as a string
std::string ParityToString(int value);
std::string Charge3XtoString(int q3);
std::string Spin2XtoString(int J2);


// String splitting
std::vector<std::string> SplitStr2Str(std::string input, const char delim = ',');
std::vector<int>         SplitStr2Int(std::string input, const char delim = ',');
std::vector<std::string> Extract(const std::string& str);


// Split string to int or double
template <class T>
std::vector<T> SplitStr(std::string input, T type, const char delim = ',') {
	std::vector<T> output;
	std::stringstream ss(input);

	// Get inputfiles by comma
	while (ss.good()) {
		std::string substr;
		std::getline(ss, substr, delim);

		// Detect type >>
		// int
		if (std::is_same<T, int>::value) {
			output.push_back(std::stoi(substr));
		}
		// double
		else if (std::is_same<T, double>::value) {
			output.push_back(std::stod(substr));
		}
	}
	return output;
}

// Check if file exists
bool FileExist(const std::string& name);

void PrintWarning();
void PrintGameOver();
void PrintFlashScreen(rang::fg pcolor);
void PrintVersion();

// Bar print
void PrintBar(std::string str, uint N = 74);

// Create directory
void CreateDirectory(std::string fullpath);

double      GetVersion();
std::string GetVersionString();
std::string GetVersionTLatex();
std::string GetWebTLatex();

// Assert functions

// Check cut, lower bound needs to be smaller than upper bound
template <typename T>
bool AssertCut(std::vector<T> cut, const std::string& name = "", bool dothrow = false) {

	if (cut.size() != 2) {
		throw std::invalid_argument("AssertCut: Input '" + name + "' vector size not 2");
	}
	if (cut[1] <= cut[0] && dothrow) {
		std::string message = "AssertCut: Input '" + name + "' with " + std::to_string(cut[1])
							+ " <= " + std::to_string(cut[0]);
		throw std::invalid_argument(message);
	}
	return true;
}


// Assert compare function, threshold 0.01 means 1 percent accuracy
template <typename T>
bool AssertRatio(T value, T reference, T threshold, const std::string& name = "", bool dothrow = false) {

	bool ok = false;
	const T ratio = value / reference;
	if (ratio < (1.0 + threshold) && ratio > (1.0 - threshold)) {
		ok = true;
	}
	if (!ok && dothrow) {
		throw std::invalid_argument("AsserRatio: Input '"
			+ name + "' = " + std::to_string(value) + " not within reference = "
			+ std::to_string(reference) + " under threshold = " + std::to_string(threshold));		
	}
	return ok;
}

// Assert range [a,b]
template <typename T>
bool AssertRange(T value, std::vector<T> range, const std::string& name = "", bool dothrow = false) {

	if (range.size() != 2) {
		throw std::invalid_argument("AssertRange: Input '" + name + "' , range vector size is not 2!");
	}
	bool ok = false;
	if (value >= range[0] && value <= range[1]) {
		ok = true;
	}
	if (!ok && dothrow) {
		throw std::invalid_argument("AssertRange: Input '"
			+ name + "' = " + std::to_string(value) + " out of range ["
			+ std::to_string(range[0]) + "," + std::to_string(range[1]) + "]");
	}
	return ok;
}

// Assert value if found from a set of numbers
template <typename T>
bool AssertSet(T value, std::vector<T> set, const std::string& name = "", bool dothrow = false) {

	bool ok = false;
	for (const auto& i : set) {
		if (value == i) { 
			ok = true;
			break;
		}
	}
	if (!ok && dothrow) {
		throw std::invalid_argument("AssertSet: Input '"
			+ name + "' = " + std::to_string(value) + " not found from the input set");
	}
	return ok;
}

// ----------------------------------------------------------------------
// Templates for enhanced index based looping
// for (const auto& i : indices(vector)) { vector[i] = foo; ... }

template <typename T>
struct index_range {

	struct iterator {
		bool operator!=(iterator x) const { return index != x.index; }
		iterator operator++() { ++index; return *this; }
		T operator*() const { return index; }

		T index;
	};

	iterator begin() const { return {0}; }
	iterator end()   const { return {n}; }

	T n;
};

template <typename T, typename Index = typename T::size_type>
index_range<Index> indices(const T& container) {
  return {container.size()};
}
// ----------------------------------------------------------------------

} // aux namespace
} // gra namespace

#endif