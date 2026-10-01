// I/O aux functions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MAUX_H
#define MAUX_H

// C++
#include <limits.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <mutex>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>
//#include <experimental/filesystem>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MFloat.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// HepMC3
#include "HepMC3/FourVector.h"

// Libraries
#include "rang.hpp"

namespace gra {
namespace aux {

// Sentinel for particles whose spin is not a fixed quantum number
constexpr int kNullSpinX2 = -1;

// Sentinel for particles whose isospin is not defined
constexpr int kNullIsospinX2 = -1;

// Example:
//
// @SOMECOMMAND:value
// @PDG[995]{M:300, W:0}

// Format options:
//
// @id:value is represented with map<"_SINGLET_", value>
//
// @id:value
// @id{arg}
// @id[]{arg}
// @id[target1,target2,...]{arg}
// @j={u,d,s,c,b,g}

struct OneCMD {
  std::string id;
  std::vector<std::string> target;
  std::map<std::string, std::string> arg;
  std::vector<std::string> values;

  void Print() {
    std::cout << "id : " << id << std::endl;

    std::cout << "[target] : [";
    for (std::size_t i = 0; i < target.size(); ++i) {
      std::cout << target[i];
      if (i < target.size() - 1) {
        std::cout << ",";
      }
    }
    std::cout << "] " << std::endl;

    std::cout << "{arg} : " << std::endl;
    for (const auto &x : arg) {
      std::cout << x.first << ":" << x.second << std::endl;
    }
    if (!values.empty()) {
      std::cout << "{values} : [";
      for (std::size_t i = 0; i < values.size(); ++i) {
        std::cout << values[i];
        if (i + 1 < values.size()) {
          std::cout << ",";
        }
      }
      std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
  }
};

// M4Vec to HepMC3::FourVector
inline HepMC3::FourVector M4Vec2HepMC3(const M4Vec &v) {
  return HepMC3::FourVector(v.X(), v.Y(), v.Z(), v.E());
}
// HepMC3::FourVector to M4Vec
inline M4Vec HepMC2M4Vec(const HepMC3::FourVector &v) {
  return M4Vec(v.x(), v.y(), v.z(), v.e());
}

// System information
std::string GetTimeStamp(const std::string format = "%d-%m-%Y_%H-%M-%S");
void PrintArgv(int argc, char *argv[]);
// Validate a plain LHAPDF set identifier before loading or downloading
void ValidateLHAPDFName(const std::string &name);
void AutoDownloadLHAPDF(const std::string pdfname);
std::string ExecCommand(const std::string &cmd);
std::string GetExecutablePath();
std::string GetBasePath(std::size_t level);
std::string ResolveProjectPath(const std::string &relative_path,
                               std::size_t max_levels = 6);

/*
std::string GetCurrentPath();
bool FileExist(const std::experimental::filesystem::path& p,
                std::experimental::filesystem::file_status s =
std::experimental::filesystem::file_status{});
*/

std::uintmax_t GetFileSize(const std::string &filename);
void GetProcessMemory(double &peak_use, double &resident_use);
void GetDiskUsage(const std::string &path, int64_t &size, int64_t &free,
                  int64_t &used);
unsigned long long TotalSystemMemory();
std::string SystemName();
std::string HostName();
std::string DateTime();

std::string execsystem(const char *cmd);

// Progress bar
bool IsTerminal();
void PrintProgress(double ratio);
void ClearProgress();

// djb2hash function
unsigned long djb2hash(const std::string &s);

// Base case: for arithmetic types, use std::to_string
template <typename T>
std::enable_if_t<std::is_arithmetic_v<T>, std::string>
dvec2str(const T &value) {
  return std::to_string(value);
}

// Convert container content to a string (recursion e.g. for vector of vectors)
template <typename T> std::string dvec2str(const std::vector<T> &vec) {
  std::ostringstream oss;
  oss << '[';
  if (!vec.empty()) {
    oss << dvec2str(vec.front());
    for (std::size_t i = 1; i < vec.size(); ++i)
      oss << ", " << dvec2str(vec[i]);
  }
  oss << ']';
  return oss.str();
}

// Simple CSV reader
void ReadCSV(const std::string &inputfile,
             std::vector<std::vector<std::string>> &output);

// Input processing
// Read file contents without interpreting or changing the input
std::string ReadFile(const std::string &inputfile);

// Read JSON input and apply registered model card overrides
std::string GetInputData(const std::string &inputfile);

// Read a JSON card before reference expansion
std::string GetInputDataRaw(const std::string &inputfile, bool overrides = true);

bool IsIntegerDigits(const std::string &str);

// Parse one complete integer token with optional surrounding whitespace
int ParseInt(const std::string &text, const std::string &context);

// Parse one complete finite real token with optional surrounding whitespace
double ParseDouble(const std::string &text, const std::string &context);

// Trim extra spaces of a string
void TrimExtraSpace(std::string &value);
void TrimLeadSpace(std::string &value);
void TrimTrailSpace(std::string &value);
void TrimEmptySpace(std::string &value);
void TrimAllSpace(std::string &value);

// Boolean to string
std::string bool_cast(bool b);
bool ParseBool(std::string value, const std::string &context = "boolean");

// Number to string with formatting
template <typename T>
std::string ToString(const T value, const unsigned int n = 6) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << value;
  return out.str();
}

// Quantum numbers as a string
std::string ParityToString(int value);
std::string Charge3XtoString(int q3);
std::string Spin2XtoString(int J2);
std::string NullableSpin2XtoString(int J2);

// String splitting
std::vector<std::string> SplitStr2Str(std::string input, const char delim = ',',
                                      bool trimextraspace = true);
std::vector<int> SplitStr2Int(std::string input, const char delim = ',');
std::vector<std::string> Extract(const std::string &str);

// Split string to int or double
template <class T>
std::vector<T> SplitStr(std::string input, T, const char delim = ',') {
  static_assert(std::is_same_v<T, int> || std::is_same_v<T, double>);
  std::vector<T> output;
  std::stringstream ss(input);

  // Get inputfiles by comma
  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, delim);

    TrimExtraSpace(substr);

    // Detect type >>
    // int
    if constexpr (std::is_same_v<T, int>) {
      output.push_back(ParseInt(substr, "SplitStr"));
    } else {
      output.push_back(ParseDouble(substr, "SplitStr"));
    }
  }
  return output;
}

// Check if file exists
bool FileExist(const std::string &name);

void PrintNotice();
// Print a warning header with an optional compact layout
void PrintWarning(bool compact = false);
void PrintGameOver();
void PrintFlashScreen(rang::fg pcolor);
void PrintVersion();

// Bar print
void PrintBar(std::string str, unsigned int N = 74);

// ASCII table formatting
std::string FormatTable(const std::vector<std::string> &header,
                        const std::vector<std::vector<std::string>> &rows);
void PrintTable(const std::vector<std::string> &header,
                const std::vector<std::vector<std::string>> &rows,
                std::ostream &os = std::cout);

// Create directory
void CreateDirectory(std::string fullpath, std::error_code *error = nullptr);

// Version information
void CheckUpdate();
void CreateVersionJSON();

double GetVersion();
std::string GetVersionType();
std::string GetVersionDate();
std::string GetVersionUpdate();
std::string GetVersionString();
std::string GetVersionTLatex();
std::string GetWebTLatex();

// Assert functions

// Check cut, lower value needs to be smaller than upper value
template <typename T>
bool AssertCut(std::vector<T> cut, const std::string &name = "",
               bool dothrow = false) {
  if (cut.size() != 2) {
    throw std::invalid_argument("AssertCut: Input '" + name +
                                "' vector size not 2");
  }
  if (!(cut[1] > cut[0])) {
    if (dothrow) {
      std::string message = "AssertCut: Input '" + name + "' with [" +
                            std::to_string(cut[0]) + "," +
                            std::to_string(cut[1]) +
                            "] (bounds must be strictly ordered)";
      throw std::invalid_argument(message);
    }
    return false;
  }
  return true;
}

// Check cut obeys given boundaries
template <typename T>
bool AssertCutRange(std::vector<T> cut, std::vector<T> bounds,
                    const std::string &name = "", bool dothrow = false) {
  if (cut.size() != 2) {
    throw std::invalid_argument("AssertCutRange: Input '" + name +
                                "' vector size not 2");
  }
  if (bounds.size() != 2) {
    throw std::invalid_argument("AssertCutRange: Input '" + name +
                                "' bounds vector size not 2");
  }
  if (!AssertCut(cut, name, dothrow) || !AssertCut(bounds, name + " bounds", dothrow)) {
    return false;
  }

  if (cut[0] < bounds[0] || cut[1] > bounds[1]) {
    if (dothrow) {
      std::string message =
          "AssertCutRange: Input '" + name + "' with [" +
          std::to_string(cut[0]) + "," + std::to_string(cut[1]) + "]" +
          " invalid given bounds: [" + std::to_string(bounds[0]) + "," +
          std::to_string(bounds[1]) + "]";
      throw std::invalid_argument(message);
    }
    return false;
  }
  return true;
}

// Compare one value with a relative threshold and exact zero support
// |value-reference|/|reference| <= threshold
template <typename T>
bool AssertRatio(T value, T reference, T threshold,
                 const std::string &name = "", bool dothrow = false) {
  static_assert(std::is_arithmetic_v<T>,
                "AssertRatio requires an arithmetic scalar type");
  const long double observed = static_cast<long double>(value);
  const long double expected = static_cast<long double>(reference);
  const long double tolerance = static_cast<long double>(threshold);
  if (!std::isfinite(observed) || !std::isfinite(expected) ||
      !std::isfinite(tolerance)) {
    if (dothrow) {
      throw std::invalid_argument("AssertRatio: non-finite input");
    }
    return false;
  }
  if (tolerance < 0.0L) {
    if (dothrow) {
      throw std::invalid_argument("AssertRatio: negative threshold");
    }
    return false;
  }
  const bool ok =
      math::IsZero(expected)
          ? math::IsZero(observed)
          : std::abs(observed - expected) / std::abs(expected) <= tolerance;
  if (!ok && dothrow) {
    throw std::invalid_argument(
        "AssertRatio: Input '" + name + "' = " + std::to_string(value) +
        " not within reference = " + std::to_string(reference) +
        " under threshold = " + std::to_string(threshold));
  }
  return ok;
}

// Assert range [a,b]
template <typename T>
bool AssertRange(T value, std::vector<T> range, const std::string &name = "",
                 bool dothrow = false) {
  if (range.size() != 2) {
    throw std::invalid_argument("AssertRange: Input '" + name +
                                "' , range vector size is not 2!");
  }
  if (range[1] < range[0]) {
    if (dothrow) {
      throw std::invalid_argument("AssertRange: invalid ordered range");
    }
    return false;
  }
  bool ok = false;
  if (value >= range[0] && value <= range[1]) {
    ok = true;
  }
  if (!ok && dothrow) {
    throw std::invalid_argument("AssertRange: Input '" + name +
                                "' = " + std::to_string(value) +
                                " out of range [" + std::to_string(range[0]) +
                                "," + std::to_string(range[1]) + "]");
  }
  return ok;
}

// Assert value if found from a set of numbers
template <typename T>
bool AssertSet(T value, std::vector<T> set, const std::string &name = "",
               bool dothrow = false) {
  bool ok = false;
  for (const auto &i : set) {
    if (value == i) {
      ok = true;
      break;
    }
  }
  if (!ok && dothrow) {
    throw std::invalid_argument("AssertSet: Input '" + name +
                                "' = " + std::to_string(value) +
                                " not found from the input set");
  }
  return ok;
}

// ----------------------------------------------------------------------
// Templates for enhanced index based looping
// for (const auto& i : indices(vector)) { vector[i] = foo; ... }

template <typename T> struct index_range {
  struct iterator {
    constexpr bool operator!=(iterator x) const { return index != x.index; }
    constexpr iterator operator++() {
      ++index;
      return *this;
    }
    constexpr T operator*() const { return index; }

    T index;
  };

  constexpr iterator begin() const { return {0}; }
  constexpr iterator end() const { return {n}; }

  T n;
};

template <typename T,
          typename Index = decltype(std::declval<const T &>().size())>
constexpr index_range<Index> indices(const T &container) {
  return {container.size()};
}
// ----------------------------------------------------------------------

std::vector<OneCMD> SplitCommands(const std::string &fullstr);
std::vector<std::size_t> FindOccurance(const std::string &str,
                                       const std::string &sub);

} // namespace aux
} // namespace gra

#endif
