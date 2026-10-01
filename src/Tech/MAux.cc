// I/O aux functions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

//#include <experimental/filesystem>

// C system functions
#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>

// Networking
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>

// Libraries
#include "HepMC3/FourVector.h"
#include "LHAPDF/LHAPDF.h"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Tech/MJsonOverride.h"
#include "Graniitti/Tech/MJson.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/PDF/MSudakov.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;

namespace gra {
namespace aux {
namespace {

// Close one pipe owned by a unique pointer
struct PipeCloser {
  // Close one nonnull process pipe without throwing from cleanup
  void operator()(FILE *pipe) const noexcept {
    if (pipe != nullptr) {
      const int status = pclose(pipe);
      (void)status;
    }
  }
};

// Strip JSON comments while preserving strings and separating adjacent tokens
std::string StripJsonComments(const std::string &data) {
  std::string output;
  output.reserve(data.size());
  bool quoted = false;
  for (std::size_t i = 0; i < data.size(); ++i) {
    const char c = data[i];
    if (quoted) {
      output += c;
      if (c == '\\' && i + 1 < data.size()) { output += data[++i]; }
      else if (c == '"') { quoted = false; }
    } else if (c == '"') {
      quoted = true;
      output += c;
    } else if (c == '/' && i + 1 < data.size() && data[i + 1] == '/') {
      output += ' ';
      i += 2;
      while (i < data.size() && data[i] != '\n' && data[i] != '\r') { ++i; }
      if (i < data.size()) { output += data[i]; }
    } else if (c == '/' && i + 1 < data.size() && data[i + 1] == '*') {
      output += ' ';
      i += 2;
      while (i + 1 < data.size() && !(data[i] == '*' && data[i + 1] == '/')) {
        if (data[i] == '\n' || data[i] == '\r') { output += data[i]; }
        ++i;
      }
      if (i + 1 >= data.size()) {
        throw std::invalid_argument("GetInputData: unterminated JSON block comment");
      }
      ++i;
    } else {
      output += c;
    }
  }
  return output;
}

// Remove trailing JSON commas only outside quoted strings
std::string StripJsonTrailingCommas(const std::string &data) {
  std::string output;
  output.reserve(data.size());
  bool quoted = false;
  for (std::size_t i = 0; i < data.size(); ++i) {
    const char c = data[i];
    if (quoted) {
      output += c;
      if (c == '\\' && i + 1 < data.size()) { output += data[++i]; }
      else if (c == '"') { quoted = false; }
      continue;
    }
    if (c == '"') { quoted = true; }
    if (c == ',') {
      std::size_t next = i + 1;
      while (next < data.size() && std::isspace(static_cast<unsigned char>(data[next]))) { ++next; }
      std::size_t previous = i;
      while (previous > 0 && std::isspace(static_cast<unsigned char>(data[previous - 1]))) { --previous; }
      const bool has_value = previous > 0 && data[previous - 1] != '[' &&
                             data[previous - 1] != '{' && data[previous - 1] != ',';
      if (has_value && next < data.size() && (data[next] == '}' || data[next] == ']')) { continue; }
    }
    output += c;
  }
  return output;
}

}  // namespace

// -------------------------------------------------------
// FIXED HERE manually

double      GetVersion()     { return 1.50; }
std::string GetVersionType() { return "beta"; }
std::string GetVersionDate() { return "30.09.2026"; }
std::string GetVersionUpdate() {
  return "New end-to-end helicity-amplitudes; New multi-eikonals; New icetune code; Update Conda env, CMake setup, HepMC3, Eigen; Bugfixes";
}

void PrintVersion() {
  std::cout << GetVersionString() << std::endl;
  std::cout << rang::style::bold << "<github.com/mieskolainen/graniitti>" << rang::style::reset
            << std::endl
            << std::endl;
  std::cout << "References: arXiv:1910.06300, arXiv:2304.06010 [hep-ph]" << std::endl;
  std::cout << std::endl;
  std::cout << "(c) 2017-2026 Mikael Mieskolainen" << std::endl;
  std::cout << "<mikael.mieskolainen@cern.ch>" << std::endl;
  std::cout << std::endl;
  std::cout << "<opensource.org/licenses/GPL-3.0>" << std::endl;
  std::cout << "<opensource.org/licenses/MIT>" << std::endl;
}

// -------------------------------------------------------

// Print input arguments
void PrintArgv(int argc, char *argv[]) {
  std::cout << rang::fg::green << "$ ";
  for (int i = 0; i < argc; ++i) {
    const std::string s = std::string(argv[i]);
    if (s.find(' ') != std::string::npos) {  // e.g. "MP[CON]<F> -> pi+ pi-"
      std::cout << "\"" << s << "\""
                << " ";
    } else {
      std::cout << s << " ";
    }
  }
  std::cout << rang::fg::reset << std::endl;
}

// Validate a plain LHAPDF set identifier before loading or downloading
void ValidateLHAPDFName(const std::string &name) {
  const std::string letters =
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
  if (name.empty() || name == "null" ||
      letters.find(name.front()) == std::string::npos ||
      name.find_first_not_of(letters + "_.-") != std::string::npos) {
    throw std::invalid_argument("ValidateLHAPDFName: invalid PDF set name '" + name + "'");
  }
}

// Download a validated LHAPDF set with quoted shell arguments
void AutoDownloadLHAPDF(const std::string pdfname) {
  ValidateLHAPDFName(pdfname);
  std::cout << rang::fg::red
            << "aux::AutoDownloadLHAPDF: Trying automatic download:" << rang::fg::reset
            << std::endl;
  std::cout << std::endl;

  // Get instal path and remove "\n"
  std::string INSTALLPATH = aux::ExecCommand("lhapdf-config --prefix");
  INSTALLPATH.erase(std::remove(INSTALLPATH.begin(), INSTALLPATH.end(), '\n'), INSTALLPATH.end());

  if (INSTALLPATH.find("command not found") != std::string::npos) {
    throw std::invalid_argument("aux::AutoDownloadLHAPDF: Failure: lhapdf-config command missing");
  }
  
  // Quote the installation directory, including any embedded apostrophes
  std::string directory = "'";
  for (const char c : INSTALLPATH + "/share/LHAPDF") {
    directory += (c == '\'') ? "'\\''" : std::string(1, c);
  }
  directory += "'";

  // Download and unpack the selected set
  std::string cmd1 = "wget 'https://lhapdfsets.web.cern.ch/lhapdfsets/current/" + pdfname +
                     ".tar.gz' -O- | tar xz -C " + directory;
  std::cout << cmd1 << std::endl;
  std::string OUTPUT1 = aux::ExecCommand(cmd1);
}

// Compute timestamp in the given format
std::string GetTimeStamp(const std::string format) {
  const std::time_t t = std::time(nullptr);
  std::tm           tm{};
  // Keep calendar fields local because localtime returns a buffer shared by concurrent callers
  if (::localtime_r(&t, &tm) == nullptr) {
    throw std::runtime_error("aux::GetTimeStamp: local time conversion failed");
  }
  std::ostringstream oss;
  oss << std::put_time(&tm, format.c_str());
  std::string timestamp(oss.str());

  return timestamp;
}

// Run terminal command, get output to std::string
std::string ExecCommand(const std::string &cmd) {
  std::array<char, 128>                    buffer;
  std::string                              result;
  std::unique_ptr<FILE, PipeCloser> pipe(popen(cmd.c_str(), "r"));
  if (!pipe) { throw std::runtime_error("aux::ExecCommand:: popen() failed!"); }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) { result += buffer.data(); }
  return result;
}

std::string GetExecutablePath() {
  char    buff[2048];
  ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff) - 1);
  if (len != -1) {
    buff[len] = '\0';
    return std::string(buff);
  }
  return "";
}

// level 0 returns same as GetExecutablePath:
// ~ /home/user/graniitti/bin/gr
//
// level 1 returns
// ~ /home/user/graniitti/bin
//
// etc..
//
std::string GetBasePath(std::size_t level) {
  std::string s   = GetExecutablePath();
  char        sep = '/';

  for (std::size_t k = 0; k < level; ++k) {
    // Search backwards from end of string
    size_t i = s.rfind(sep);
    if (i != std::string::npos) {
      s = s.substr(0, i);
    } else {
      return s;
    }
  }
  return s;
}

// Preserve absolute paths and locate relative inputs from the executable directory
std::string ResolveProjectPath(const std::string &relative_path, std::size_t max_levels) {
  if (std::filesystem::path(relative_path).is_absolute()) { return relative_path; }
  auto parent_path = [](const std::string &path) {
    const std::size_t pos = path.rfind('/');
    if (pos == std::string::npos) { return path; }
    return path.substr(0, pos);
  };

  std::string base = parent_path(GetExecutablePath());
  for (std::size_t level = 0; level <= max_levels; ++level) {
    const std::string candidate = base + "/" + relative_path;
    if (gra::aux::FileExist(candidate)) { return candidate; }
    const std::string parent = parent_path(base);
    if (parent == base) { break; }
    base = parent;
  }

  return GetBasePath(2) + "/" + relative_path;
}

/*
// alias
namespace fs = std::experimental::filesystem;

std::string GetCurrentPath() {
                fs::path cwd = std::experimental::filesystem::current_path();
                return cwd.string();
}


// Folder exists
bool FileExist(const fs::path& p, fs::file_status s) {
        if (fs::status_known(s) ? fs::exists(s) : fs::exists(p)) {
                return true;
        } else {
                return false;
        }
}
*/

// Get filesize in bytes
std::uintmax_t GetFileSize(const std::string &filename) {
  /*
  namespace fs = std::experimental::filesystem;
  fs::path p = fs::current_path() / filename;
  return fs::file_size(p);
  */

  struct stat stat_buf;
  int         rc = stat(filename.c_str(), &stat_buf);
  return rc == 0 ? stat_buf.st_size : 0;
}

// Get Process Memory Usage (linux/BSD/OSX) in bytes
void GetProcessMemory(double &peak_use, double &resident_use) {
  // Peak memory
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
  peak_use = rusage.ru_maxrss * 1024L;

  // Current memory
  long  rss = 0L;
  FILE *fp  = NULL;
  if ((fp = fopen("/proc/self/statm", "r")) == NULL) {
    resident_use = 0.0;
    return;
  }
  if (fscanf(fp, "%*s%ld", &rss) != 1) {  // Reading problem
    resident_use = 0.0;
  } else {  // fine
    resident_use = rss * (size_t)sysconf(_SC_PAGESIZE);
  }
  fclose(fp);
}

// Get disk usage
void GetDiskUsage(const std::string &path, int64_t &size, int64_t &free, int64_t &used) {
  int64_t frsize;
  int64_t blocks;
  int64_t bfree;

  struct statvfs buf;
  int            ret = statvfs(path.c_str(), &buf);

  if (!ret) {
    frsize = buf.f_frsize;  // block size
    blocks = buf.f_blocks;  // blocks
    bfree  = buf.f_bfree;   // free blocks

    size = frsize * blocks;
    free = frsize * bfree;
    used = size - free;
  }
}

// Get available memory in bytes
unsigned long long TotalSystemMemory() {
  long pages     = sysconf(_SC_PHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  return pages * page_size;
}

// System information
std::string SystemName() {
  struct utsname name;
  uname(&name);

  std::string sysname(name.sysname);
  std::string nodename(name.nodename);
  std::string release(name.release);
  std::string version(name.version);
  std::string machine(name.machine);

  std::string s = " ";

  return sysname + s + release + s + version + s + machine;
}

// Get system hostname
std::string HostName() {
  char hostname[2048];
  gethostname(hostname, 2048);
  return std::string(hostname);
}


// Compute current date and time using the same per-call calendar storage as formatted timestamps
std::string DateTime() { return GetTimeStamp("%Y-%m-%d %X"); }

// Compute whether standard output is attached to an interactive terminal
// Compute whether standard output currently points to a terminal
bool IsTerminal() { return isatty(fileno(stdout)) != 0; }

// Print out progress bar visualization
void PrintProgress(double ratio) {
  if (ratio > 1.0) { ratio = 1.0; }

#define BAR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
  const int WIDTH = 62;
  const int pos   = static_cast<int>(ratio * 100);
  const int left  = static_cast<int>(ratio * WIDTH);
  const int right = WIDTH - left;

  if (IsTerminal()) { // if no terminal, do not print!
    std::cout << rang::fg::green;
    printf("\r%3d%% [%.*s%*s]", pos, left, BAR, right, "");
    std::cout << rang::fg::reset;
    std::cout << std::flush;
  }
}

// Clear the line, and move cursor to the left
void ClearProgress() {
  if (IsTerminal()) {
    std::cout << "\33[2K"
              << "\r";
  }
}

// djb2 hash function (used for saving unique filenames)
unsigned long djb2hash(const std::string &s) {
  unsigned long hash = 5381;  // Magic number
  for (auto c : s) {
    hash = (hash << 5) + hash + c;  // same as: hash * 33 + c
  }
  return hash;
}

// Convert a vector of doubles to string
/*
std::string dvec2str(const std::vector<double>& vec) {
    std::ostringstream oss;
    oss << '[' << (vec.empty() ? "" : std::accumulate(std::next(vec.begin()), vec.end(), std::to_string(vec[0]),
        [](std::string a, double b) { return a + ", " + std::to_string(b); })) << ']';
    return oss.str();
}
*/


// Read CSV file
void ReadCSV(const std::string &inputfile, std::vector<std::vector<std::string>> &output) {
  std::ifstream file(inputfile);
  if (!file.is_open()) { throw std::invalid_argument("MAux::ReadCSV: Cannot open inputfile " + inputfile); }
  std::string line;

  // Read every line from the stream
  while (getline(file, line)) {
    std::istringstream       stream(line);
    std::vector<std::string> columns;
    std::string              element;

    // Every line element separated by separator
    while (getline(stream, element, ',')) { columns.push_back(element); }
    output.push_back(columns);
  }
  if (file.bad()) { throw std::ios_base::failure("MAux::ReadCSV: Cannot read inputfile " + inputfile); }
}

// Read file contents without interpreting or changing the input
std::string ReadFile(const std::string &inputfile) {
  std::ifstream ifs(inputfile, std::ios::binary);
  if (!ifs.is_open()) {
    throw std::invalid_argument("MAux::ReadFile: Cannot open inputfile " + inputfile);
  }
  std::string data((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
  if (ifs.bad()) { throw std::ios_base::failure("MAux::ReadFile: Cannot read inputfile " + inputfile); }
  return data;
}

// Read JSON input and apply registered model card overrides
std::string GetInputDataRaw(const std::string &inputfile, bool overrides) {
  std::string data = ReadFile(inputfile);

  // Accept comments and trailing commas without changing JSON string values
  data = StripJsonTrailingCommas(StripJsonComments(data));

  // Apply registered command line JSON overrides for matching card reads
  if (overrides && gra::json_override::HasCardOverrides()) {
    nlohmann::json j = nlohmann::json::parse(data);
    gra::json_override::ApplyRegisteredCardOverrides(inputfile, j);
    data = j.dump();
  }

  return data;
}

// Read a JSON card and resolve references through the common loader
std::string GetInputData(const std::string &inputfile) { return gra::MJson{}.Read(inputfile).dump(); }

// Check for one optional sign followed by decimal integer digits
bool IsIntegerDigits(const std::string &str) {
  if (str.empty()) { return false; }
  const std::size_t first = (str.front() == '-' || str.front() == '+') ? 1 : 0;
  return first < str.size() && str.find_first_not_of("0123456789", first) == std::string::npos;
}

// Parse one complete integer token with optional surrounding whitespace
int ParseInt(const std::string &text, const std::string &context) {
  try {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (text.find_first_not_of(" \t\r\n\f\v", consumed) != std::string::npos) {
      throw std::invalid_argument("trailing characters");
    }
    return value;
  } catch (const std::exception &error) {
    throw std::invalid_argument(context + ": invalid integer '" + text + "': " + error.what());
  }
}

// Parse one complete finite real token with optional surrounding whitespace
double ParseDouble(const std::string &text, const std::string &context) {
  try {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (!std::isfinite(value) || text.find_first_not_of(" \t\r\n\f\v", consumed) != std::string::npos) {
      throw std::invalid_argument("non-finite value or trailing characters");
    }
    return value;
  } catch (const std::exception &error) {
    throw std::invalid_argument(context + ": invalid real value '" + text + "': " + error.what());
  }
}

// Compute particle parity as a string
std::string ParityToString(int value) {
  if (value > 0) {
    return "+";
  } else if (value == 0) {
    return "";
  } else {
    return "-";
  }
}

// Compute particle charge as a string
std::string Charge3XtoString(int q3) {
  const std::string sign  = (q3 < 0) ? "-" : " ";
  const int         absq3 = std::abs(q3);

  if (absq3 == 6)
    return sign + "2";
  else if (absq3 == 3)
    return sign + "1";
  else if (absq3 == 2)
    return sign + "2/3";
  else if (absq3 == 1)
    return sign + "1/3";
  else
    return "0";
}

// Compute spin as a string
std::string Spin2XtoString(int J2) {
  if (J2 < 0)
    return "-" + Spin2XtoString(-J2);
  else if (J2 == 10)
    return "5";
  else if (J2 == 9)
    return "9/2";
  else if (J2 == 8)
    return "4";
  else if (J2 == 7)
    return "7/2";
  else if (J2 == 6)
    return "3";
  else if (J2 == 5)
    return "5/2";
  else if (J2 == 4)
    return "2";
  else if (J2 == 3)
    return "3/2";
  else if (J2 == 2)
    return "1";
  else if (J2 == 1)
    return "1/2";
  else
    return "0";
}

// Compute nullable particle spin as a string
std::string NullableSpin2XtoString(int J2) {
  if (J2 == kNullSpinX2) { return "null"; }
  return Spin2XtoString(J2);
}

// Split a string to strings separated by delimiter
std::vector<std::string> SplitStr2Str(std::string input, const char delim, bool trimextraspace) {
  std::vector<std::string> output;
  std::stringstream        ss(input);

  // String by string
  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, delim);

    if (trimextraspace) { TrimExtraSpace(substr); }
    output.push_back(substr);
  }
  return output;
}

// Split string to ints
std::vector<int> SplitStr2Int(std::string input, const char delim) {
  std::vector<int>  output;
  std::stringstream ss(input);

  // Get inputfiles by comma
  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, delim);

    TrimExtraSpace(substr);
    output.push_back(ParseInt(substr, "SplitStr2Int"));
  }
  return output;
}

// Trim leading, extra and trailing spaces
void TrimExtraSpace(std::string &value) {
  value = std::regex_replace(value, std::regex(R"(^ +| +$|( ) +)"), "$1");
}

void TrimLeadSpace(std::string &value) {
  value = std::regex_replace(value, std::regex(R"(^ +)"), "$1");
}

void TrimTrailSpace(std::string &value) {
  value = std::regex_replace(value, std::regex(R"( +$)"), "$1");
}

void TrimEmptySpace(std::string &value) {
  value = std::regex_replace(value, std::regex(R"( +)"), "$1");
}

void TrimAllSpace(std::string &value) {
  value = std::regex_replace(value, std::regex(R"([^\S\r\n]+)"), "$1");
}


// Extract words from a string
std::vector<std::string> Extract(const std::string &str) {
  std::vector<std::string> words;
  std::stringstream        ss(str);
  std::string              buff;

  while (ss >> buff) { words.push_back(buff); }
  return words;
}

// Check if file exists
bool FileExist(const std::string &name) {
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

void PrintNotice() {
  std::cout << rang::fg::red << "<NOTICE>\n\n"
            << "ZZZ    ZZ  ZZZZZZ  ZZZZZZZZ ZZ ZZZZZZ ZZZZZZZ  \n"
               "ZZZZ   ZZ ZZ    ZZ    ZZ    ZZ ZZ     ZZ       \n"
               "ZZ ZZ  ZZ ZZ    ZZ    ZZ    ZZ ZZ     ZZZZZ    \n"
               "ZZ  ZZ ZZ ZZ    ZZ    ZZ    ZZ ZZ     ZZ       \n"
               "ZZ   ZZZZ  ZZZZZZ     ZZ    ZZ ZZZZZZ ZZZZZZZ  \n"
            << rang::fg::reset << std::endl;
}

// Print a warning header with an optional compact layout
void PrintWarning(bool compact) {
  std::cout << rang::fg::red << "<WARNING>";
  if (compact) {
    std::cout << rang::fg::reset << std::endl;
    return;
  }
  std::cout << "\n\n"
            << "ZZ     ZZ  ZZZZZ  ZZZZZZ  ZZZ    ZZ ZZ ZZZ    ZZ  ZZZZZZ   \n"
               "ZZ     ZZ ZZ   ZZ ZZ   ZZ ZZZZ   ZZ ZZ ZZZZ   ZZ ZZ        \n"
               "ZZ  Z  ZZ ZZZZZZZ ZZZZZZ  ZZ ZZ  ZZ ZZ ZZ ZZ  ZZ ZZ   ZZZ  \n"
               "ZZ ZZZ ZZ ZZ   ZZ ZZ   ZZ ZZ  ZZ ZZ ZZ ZZ  ZZ ZZ ZZ    ZZ  \n"
               " ZZZ ZZZ  ZZ   ZZ ZZ   ZZ ZZ   ZZZZ ZZ ZZ   ZZZZ  ZZZZZZ   \n"
            << rang::fg::reset << std::endl;
}

void PrintGameOver() {
  std::cout << "<GAME OVER>\n\n"
            << " ZZZZZZ   ZZZZZ  ZZZ    ZZZ ZZZZZZZ    ZZZZZZ  ZZ    ZZ ZZZZZZZ ZZZZZZ   \n"
               "ZZ       ZZ   ZZ ZZZZ  ZZZZ ZZ        ZZ    ZZ ZZ    ZZ ZZ      ZZ   ZZ  \n"
               "ZZ   ZZZ ZZZZZZZ ZZ ZZZZ ZZ ZZZZZ     ZZ    ZZ ZZ    ZZ ZZZZZ   ZZZZZZ   \n"
               "ZZ    ZZ ZZ   ZZ ZZ  ZZ  ZZ ZZ        ZZ    ZZ  ZZ  ZZ  ZZ      ZZ   ZZ  \n"
               " ZZZZZZ  ZZ   ZZ ZZ      ZZ ZZZZZZZ    ZZZZZZ    ZZZZ   ZZZZZZZ ZZ   ZZ  \n"
            << std::endl;
}


// Execute terminal command get output
std::string execsystem(const char *cmd) {
  std::array<char, 1024>                   buffer;
  std::string                              result;
  std::unique_ptr<FILE, PipeCloser> pipe(popen(cmd, "r"));
  if (!pipe) { throw std::runtime_error("popen() failed!"); }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) { result += buffer.data(); }
  return result;
}


// Check updates online
void CheckUpdate() {
  auto PrintMessage = []() {
    std::cout << std::endl;
    PrintBar("-", 80);
    std::cout << "Check you have 'curl' and internet access -- could not update version information"
              << std::endl;
    PrintBar("-", 80);
    std::cout << std::endl;
  };

  const std::string cmd =
      "curl -s https://raw.githubusercontent.com/mieskolainen/GRANIITTI/master/VERSION.json &> "
      "/dev/stdout";
  const std::string data = execsystem(cmd.c_str());

  if (data.size() > 0) {
    nlohmann::json j;

    try {
      j = nlohmann::json::parse(data);

      const double      online_version = j.at("version");
      const std::string online_date    = j.at("date");

      if (GetVersion() < online_version) {  // Older version than online
        std::cout << std::endl;
        PrintBar("-", 80);
        std::cout << rang::style::bold << rang::fg::green << "New version " << online_version
                  << " (" << online_date << ") available at <github.com/mieskolainen/graniitti>"
                  << rang::fg::reset << rang::style::reset << std::endl;
        std::cout << std::endl;
        std::cout << "Updates are: " << j.at("update") << std::endl;
        std::cout << std::endl;
        std::cout << "To update, copy-and-run: " << std::endl;
        std::cout
            << "git pull origin master && source ./install/setenv.sh && make superclean && make -j4"
            << std::endl
            << std::endl;
        std::cout << "If compilation fails, check that HepMC3 and LHAPDF6 are as required"
                  << std::endl;
        std::cout << "You can re-install them with cd install && source autoinstall.sh"
                  << std::endl;
        PrintBar("-", 80);
        std::cout << std::endl;
      } else if (std::abs(GetVersion() - online_version) < 1e-6) {  // Same as online
        std::cout << std::endl;
        PrintBar("-", 80);
        std::cout << rang::style::bold << "This version " << GetVersion() << " ("
                  << GetVersionDate() << ") is up to date with the online version"
                  << rang::style::reset << std::endl;
        PrintBar("-", 80);
        std::cout << std::endl;
      } else {  // This version is newer than online
        std::cout << std::endl;
        PrintBar("-", 80);
        std::cout << rang::style::bold << rang::fg::green << "This version " << GetVersion() << " ("
                  << GetVersionDate() << ") is newer than the online version " << online_version
                  << " (" << online_date << ")" << rang::style::reset << std::endl;
        PrintBar("-", 80);
        std::cout << std::endl;
      }

    } catch (...) { PrintMessage(); }
  } else {
    PrintMessage();
  }

  // Create this version JSON
  CreateVersionJSON();
}

void CreateVersionJSON() {
  // Get program path and output this VERSION.json
  const std::string output = GetBasePath(2) + "/VERSION.json";
  const std::string MSG =
      "{\n \"name\": \"%s\",\n \"version\": %0.3f,\n \"type\": \"%s\",\n \"date\": \"%s\",\n "
      "\"update\": \"%s\"\n}\n";

  FILE *file = fopen(output.c_str(), "w");
  fprintf(file, MSG.c_str(), "GRANIITTI", GetVersion(), GetVersionType().c_str(),
          GetVersionDate().c_str(), GetVersionUpdate().c_str());
  fclose(file);
}

std::string bool_cast(bool b) {
  std::ostringstream ss;
  ss << std::boolalpha << b;
  return ss.str();
}

bool ParseBool(std::string value, const std::string &context) {
  TrimExtraSpace(value);
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

  if (value == "true" || value == "1") { return true; }
  if (value == "false" || value == "0") { return false; }

  throw std::invalid_argument("ParseBool: " + context +
                              " expects true/false (case-insensitive) or 1/0, got '" + value + "'");
}

std::string GetVersionString() {
  char buff[100];
  snprintf(buff, sizeof(buff), "Version %0.3f (%s) %s", GetVersion(), GetVersionType().c_str(),
           GetVersionDate().c_str());
  std::string str = buff;
  return str;
}

std::string GetVersionTLatex() {
  char         buff[100];
  const double version = GetVersion();
  snprintf(buff, sizeof(buff), "#color[16]{#scale[0.6]{GRANIITTI #scale[0.8]{%0.3f}}}", version);
  std::string str = buff;
  return str;
}

std::string GetWebTLatex() {
  char buff[100];
  snprintf(buff, sizeof(buff), "#color[16]{#scale[0.5]{#LTgithub.com/mieskolainen#GT}}");
  std::string str = buff;
  return str;
}

void PrintFlashScreen(rang::fg pcolor) {
  std::cout << std::endl;
  gra::aux::PrintBar("-");
  std::cout << pcolor
            << ".``````````````````````     ````   ``                        ``   ```   ``\n"
               ".```..````````````    ` ``          `                             ``    ``\n"
               ".``..:``````````` ```.``````.`.```                                 `    ``\n"
               ".``.`.```````   ````     ``-.```.````                               ``` ``\n"
               "..-.`-.``    `.`   `..``---/--.......``                      `````````````\n"
               "..:..-.`````.`   `.--oys+ohms:..-.`..```                   ```````````````\n"
               ".-:.-...`..     .ooshydddddddh+::...``...               ``````````````````\n"
               "..-......       `yhhdhhmhddddmdy-...-...`                   ``````````````\n"
               "::::::.`    `    `--/ohdhdhddyo-`   `..`..                               `\n"
               "so/-.`  `.`.-````-...++/syhdho`        `..                               `\n"
               "``..--.:o-..--..-/-::-.-:/ohh/           `                            `` `\n"
               "`.-::-.:yyh/`--.`s+ohso:/shhd:                                      ``  `-\n"
               ".:/-s+o+syhs--:/-..-oydyshhhho`                                 ````` ``./\n"
               ".-+:-+shhhys+ysyh//sshyy/:oohyo.```````` `````` ````````````...::..  `.-oh\n"
               "-..-..:/+shhhyhhhyhysso.   ``:ys:----............------:::::++oo-  `.-:ohd\n"
               "--+:--//:::/+oyyshhs+`   ``   .+so+++++++++++/-://+++++++++o+/-```-/yhhys+\n"
               "`/dysyhs/:::-:o+/-.`   `...    `..::+oooooo+:-``.:+oooo+/:-.```.:/ssso/.``\n"
               "``:oossso-.--.-.    ..`...`..`.``...`-/ossys+ssoso+ys:--..::/:/--.-//:```.\n"
               "```````` `   `-.` `....```...:++.`./o:+-:/-:://----...--.`-::::+ys/.`.:/sy\n"
               "```...-...`/...... ` `..`...`yo+:/o.:oshssy/.  ``  `.--.--:+sys+-` -+ossss\n"
               "``.//-/++-sh+-.-..`-`.-.--.- :hhho/+//:/yyyhhs+::+o+/++ssssso:.`-::oyssyyy\n"
               "`/yhoshhhh////++:.....`:yy+/...sddhs/:////:-/++ooo+oyssyoo:.`:/+oosssshddd\n"
               "`dddyhydhdyyddhddo``.+.`+/oyoh-..+yddho/-......--:----://+osyyhho++/--::/+\n"
               "/hddhmhmhmdmddhdds:+s+:+/-:++hhh//.-/++osssso++/-.````.-://::-.....`......\n"
            << std::endl
            << std::endl;
  std::cout << rang::fg::reset;
}


// Print horizontal bar
void PrintBar(std::string str, unsigned int N) {
  for (std::size_t k = 0; k < N; ++k) { std::cout << str; }
  std::cout << std::endl;
}

// Format aligned columns with horizontal separators for empty rows
std::string FormatTable(const std::vector<std::string> &header,
                        const std::vector<std::vector<std::string>> &rows) {
  std::size_t ncols = header.size();
  for (const auto &row : rows) { ncols = std::max(ncols, row.size()); }

  std::vector<std::size_t> width(ncols, 0);
  for (const auto &i : indices(header)) {
    width[i] = std::max(width[i], header[i].size());
  }
  for (const auto &row : rows) {
    for (const auto &i : indices(row)) {
      width[i] = std::max(width[i], row[i].size());
    }
  }

  auto print_row = [&width, ncols](std::ostringstream &out,
                                   const std::vector<std::string> &row) {
    for (std::size_t i = 0; i < ncols; ++i) {
      if (i > 0) { out << " | "; }
      const std::string cell = (i < row.size()) ? row[i] : "";
      out << std::left << std::setw(static_cast<int>(width[i])) << cell;
    }
    out << '\n';
  };

  auto print_separator = [&width](std::ostringstream &out) {
    for (const auto &i : indices(width)) {
      if (i > 0) { out << "-+-"; }
      out << std::string(width[i], '-');
    }
    out << '\n';
  };

  std::ostringstream out;
  print_row(out, header);
  print_separator(out);

  for (const auto &row : rows) {
    if (row.empty()) { print_separator(out); }
    else { print_row(out, row); }
  }

  return out.str();
}

// Print an aligned table to the selected stream
void PrintTable(const std::vector<std::string> &header,
                const std::vector<std::vector<std::string>> &rows,
                std::ostream &os) {
  os << FormatTable(header, rows);
}

// Create missing directories and report filesystem errors
void CreateDirectory(std::string fullpath, std::error_code *error) {
  if (error) { error->clear(); }
  if (fullpath.empty()) { return; }
  if (error) { std::filesystem::create_directories(fullpath, *error); }
  else { std::filesystem::create_directories(fullpath); }
}

// Get commandline arguments which are split by @... @... tagging syntax
std::vector<OneCMD> SplitCommands(const std::string &fullstr) {
  // Find @ blocks [ ... ]
  std::vector<std::size_t> pos = FindOccurance(fullstr, "@");
  std::vector<std::string> subcmd;

  if (!pos.empty()) {
    // Split multiple @id{} @id{} blocks using positions in the original text
    for (const auto &i : indices(pos)) {
      const std::size_t start = pos[i] + 1;
      const std::size_t end = (i + 1 < pos.size()) ? pos[i + 1] : fullstr.size();
      subcmd.push_back(fullstr.substr(start, end - start));
    }
  }

  // Loop over @ blocks
  std::vector<OneCMD> cmd;

  for (const auto &i : indices(subcmd)) {
    // ** Remove whitespace **
    TrimExtraSpace(subcmd[i]);

    std::string                        id;
    std::vector<std::string>           target;
    std::map<std::string, std::string> arg;
    std::vector<std::string>           values;

    // Check we find [] brackets
    std::size_t left  = subcmd[i].find("[");
    std::size_t right = subcmd[i].find("]");

    // Incomplete brackets
    if (left == std::string::npos && right != std::string::npos) {
      throw std::invalid_argument("SplitCommands: Incomplete @R[]{} syntax with missing [");
    }
    // Incomplete brackets
    else if (left != std::string::npos && right == std::string::npos) {
      throw std::invalid_argument("SplitCommands: Incomplete @R[]{} syntax with missing ]");
    }
    // Found [] brackets
    else if (left != std::string::npos && right != std::string::npos) {
      const std::size_t argument = subcmd[i].find_first_of("{:=");
      if (right < left || (argument != std::string::npos && argument < right) ||
          subcmd[i].find('[', left + 1) != std::string::npos ||
          subcmd[i].find(']', right + 1) != std::string::npos) {
        throw std::invalid_argument("SplitCommands: malformed [] targets");
      }
      // Read (multiple) targets out
      const std::string TARGET_str = subcmd[i].substr(left + 1, right - left - 1);
      if (TARGET_str != "") {  // do not process empty
        target = SplitStr2Str(TARGET_str, ',');
      }

      // Remove [...] block including brackets
      subcmd[i].erase(left, right - left + 1);
    }

    // Check we find {} brackets
    std::size_t Lpos = subcmd[i].find("{", 0);
    std::size_t Rpos = subcmd[i].find("}", 0);

    // Incomplete brackets
    if (Lpos == std::string::npos && Rpos != std::string::npos) {
      throw std::invalid_argument("SplitCommands: Incomplete @R[]{} syntax with missing {");
    }
    // Incomplete brackets
    else if (Lpos != std::string::npos && Rpos == std::string::npos) {
      throw std::invalid_argument("SplitCommands: Incomplete @R[]{} syntax with missing }");
    }

    if (Lpos != std::string::npos && Rpos != std::string::npos &&
        (Rpos < Lpos || subcmd[i].find('{', Lpos + 1) != std::string::npos ||
         subcmd[i].find('}', Rpos + 1) != std::string::npos || Rpos + 1 != subcmd[i].size())) {
      throw std::invalid_argument("SplitCommands: malformed {} arguments or trailing text");
    }

    // Final-state parton set command @j={u,d,s,c,b,g}
    const bool value_list =
        Lpos != std::string::npos && Lpos > 0 && subcmd[i][Lpos - 1] == '=';
    if (value_list) {
      id = subcmd[i].substr(0, Lpos - 1);
      TrimExtraSpace(id);
      if (id != "j") {
        throw std::invalid_argument(
            "gra::SplitCommands: value-list syntax is supported only for @j={...}");
      }
      if (!target.empty()) {
        throw std::invalid_argument("gra::SplitCommands: @j does not accept [] targets");
      }
      if (subcmd[i].find('{', Lpos + 1) != std::string::npos ||
          subcmd[i].find('}', Rpos + 1) != std::string::npos || Rpos + 1 != subcmd[i].size()) {
        throw std::invalid_argument("gra::SplitCommands: malformed @j={...} value list");
      }
      values = SplitStr2Str(subcmd[i].substr(Lpos + 1, Rpos - Lpos - 1), ',');
      if (!values.empty() && values.back().empty()) { values.pop_back(); }
      if (values.empty() ||
          std::any_of(values.begin(), values.end(),
                      [](const std::string &value) { return value.empty(); })) {
        throw std::invalid_argument("gra::SplitCommands: @j={...} contains an empty value");
      }
    }

    // Singlet command @ID:VALUE or @ID=VALUE
    if (!value_list && Lpos == std::string::npos && Rpos == std::string::npos) {
      char delimiter = ':';
      if (subcmd[i].find(delimiter) == std::string::npos &&
          subcmd[i].find('=') != std::string::npos) {
        delimiter = '=';
      }
      std::vector<std::string> strip = SplitStr2Str(subcmd[i], delimiter);

      // ID
      id = strip[0];

      if (strip.size() == 1) {
        arg["_SINGLET_"] = "true";  // default true, no : given
      } else if (strip.size() == 2) {
        arg["_SINGLET_"] = strip[1];
      } else {
        throw std::invalid_argument("gra::SplitCommands: @Syntax invalid with '" + subcmd[i] + "'");
      }

      // Block command @blaa{key:val,key:val,...}
    } else if (!value_list) {  // We have {} block

      // ID
      id = subcmd[i].substr(0, Lpos);

      // Content {}

      // Now split all arguments inside {} by comma ','
      std::vector<std::string> keyvals =
          SplitStr2Str(subcmd[i].substr(Lpos + 1, Rpos - Lpos - 1), ',');

      // Add all
      for (const auto &i : indices(keyvals)) {
        // Split using syntax definition: key:val
        std::vector<std::string> strip = SplitStr2Str(keyvals[i], ':');
        if (strip.size() != 2) {
          throw std::invalid_argument("gra::SplitCommands: @Syntax not good with brackets {" +
                                      keyvals[i] + "}");
        }

        // Strip spaces from the key
        TrimExtraSpace(strip[0]);

        arg[strip[0]] = strip[1];  // add to map
      }
    }

    TrimExtraSpace(id);
    if (id.empty()) { throw std::invalid_argument("SplitCommands: missing command identifier"); }

    // Add this command block
    OneCMD o;
    o.id     = id;
    o.target = target;
    o.arg    = arg;
    o.values = values;
    o.Print();  // For debug
    cmd.push_back(o);
  }

  const auto selector_count =
      std::count_if(cmd.begin(), cmd.end(), [](const OneCMD &entry) { return entry.id == "j"; });
  if (selector_count > 1) {
    throw std::invalid_argument("gra::SplitCommands: duplicate @j final-state selector");
  }
  return cmd;
}

// Find string occurances
std::vector<std::size_t> FindOccurance(const std::string &str, const std::string &sub) {
  // Holds all the positions that sub occurs within str
  std::vector<size_t> positions;
  std::size_t         pos = str.find(sub, 0);

  while (pos != std::string::npos) {
    positions.push_back(pos);
    pos = str.find(sub, pos + 1);
  }
  return positions;
}


}  // namespace aux
}  // namespace gra
