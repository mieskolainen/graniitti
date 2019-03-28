// I/O aux functions
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <fstream>
#include <iostream>
#include <regex>
#include <vector>
//#include <experimental/filesystem>

// C system functions
#include <unistd.h>
#include <sys/utsname.h>
#include <sys/statvfs.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>


// LHAPDF
#include "LHAPDF/LHAPDF.h"

// HepMC3
#include "HepMC3/FourVector.h"

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MSudakov.h"

// Libraries
#include "rang.hpp"

namespace gra {
namespace aux {


// ======================================================================
// "PROGRAM GLOBALS"

// Model tune
std::string MODELPARAM = "";

// Sudakov routines
MSudakov* GlobalSudakovPtr = nullptr;

// GLOBAL for multithreading: HepMC3 outputfile and LHAPDF and mutex lock
std::mutex g_mutex;

// Check do we have terminal output (=true), or output to file (=false)
static const bool IS_TERMINAL = isatty(fileno(stdout)) != 0;

// Normal pdfs
LHAPDF::PDF* GlobalPdfPtr = nullptr;
int pdf_trials = 0;

// ======================================================================


// Download LHAPDFset automatically
void AutoDownloadLHAPDF(const std::string pdfname) {

    std::cout << rang::fg::red << "aux::AutoDownloadLHAPDF: Trying automatic download:" << rang::fg::reset << std::endl;
    std::cout << std::endl;

    // Get instal path and remove "\n"
    std::string INSTALLPATH = aux::ExecCommand("lhapdf-config --prefix");
	INSTALLPATH.erase(std::remove(INSTALLPATH.begin(), INSTALLPATH.end(), '\n'), INSTALLPATH.end());

	if (INSTALLPATH.find("command not found") != std::string::npos) {
		throw std::invalid_argument("aux::AutoDownloadLHAPDF: Failure: lhapdf-config command missing");
	}

    // Download, untar
    std::string cmd1 = "wget http://www.hepforge.org/archive/lhapdf/pdfsets/6.2/" + pdfname + ".tar.gz -O- | " +
    				   " tar xz -C " + INSTALLPATH + "/share/LHAPDF";
    std::cout << cmd1 << std::endl;
    std::string OUTPUT1 = aux::ExecCommand(cmd1);
}


// Run terminal command, get output to std::string
std::string ExecCommand(const std::string& cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("aux::ExecCommand:: popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}


std::string GetExecutablePath() {
    char buff[2048];
    ssize_t len = ::readlink("/proc/self/exe", buff, sizeof(buff)-1);
    if (len != -1) {
      buff[len] = '\0';
      return std::string(buff);
    }
    return "";
}

// level 0 returns same as GetExecutablePath:
// ~ /home/user/graniitti/bin/gr
// level 1 returns
// ~ /home/user/graniitti/bin
std::string GetBasePath(std::size_t level) {

	std::string s = GetExecutablePath();
	char sep = '/';

	for (std::size_t k = 0; k < level; ++k) {

		// Search backwards from end of string
		size_t i = s.rfind(sep);
	   	if (i != std::string::npos) {
			s = s.substr(0,i);
	   	} else {
	   		return s;
	   	}
	}
	return s;
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
std::uintmax_t GetFileSize(const std::string& filename) {
	
	/*
	namespace fs = std::experimental::filesystem;
	fs::path p = fs::current_path() / filename;
	return fs::file_size(p);
	*/

    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : 0;
}


// Get Process Memory Usage (linux/BSD/OSX) in bytes
void GetProcessMemory(double& peak_use, double& resident_use) {

	// Peak memory
	struct rusage rusage;
	getrusage(RUSAGE_SELF, &rusage);
	peak_use = rusage.ru_maxrss * 1024L;

	// Current memory
	long rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL ){
		resident_use = 0.0;
		return;
	}
	if ( fscanf(fp, "%*s%ld", &rss ) != 1 ) { // Reading problem
		resident_use = 0.0;
	} else { // fine
		resident_use = rss * (size_t)sysconf(_SC_PAGESIZE);
	}
	fclose(fp);
}


// Get disk usage
void GetDiskUsage(const std::string& path, int64_t& size, int64_t& free, int64_t& used) {

	int64_t frsize;
	int64_t blocks;
	int64_t bfree;

	struct statvfs buf;
	int ret = statvfs(path.c_str(), &buf);

	if (!ret) {
		frsize = buf.f_frsize; // block size
		blocks = buf.f_blocks; // blocks
		bfree  = buf.f_bfree;  // free blocks

		size = frsize * blocks;
		free = frsize * bfree;
		used = size - free;
	}
}


// Get available memory in bytes
unsigned long long TotalSystemMemory() {

    long pages = sysconf(_SC_PHYS_PAGES);
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

	return sysname +s+ release +s+ version +s+ machine;
}


// Get system hostname
std::string HostName() {

	char hostname[2048];
	gethostname(hostname, 2048);
	
	std::string str(hostname);
	return str;
}
	
// Get current date/time, format is YYYY-MM-DD HH:mm:ss
// http://en.cppreference.com/w/cpp/chrono/c/strftime
// for more information about date/time format
const std::string DateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);

	strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
    return buf;
}

// Print out progress bar visualization
void PrintProgress(double ratio) {

	if (ratio > 1.0) { ratio = 1.0; }

#define BAR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
	const int WIDTH = 62;
	const int pos = static_cast<int>(ratio * 100);
	const int left = static_cast<int>(ratio * WIDTH);
	const int right = WIDTH - left;
	
	if (gra::aux::IS_TERMINAL) { // if no terminal, do not print!
		std::cout << rang::fg::green;
		printf("\r%3d%% [%.*s%*s]", pos, left, BAR, right, "");
		std::cout << rang::fg::reset;
		std::cout << std::flush;
	}
}

// Clear the line, and move cursor to the left
void ClearProgress() {
	std::cout << "\33[2K" << "\r";
}


// djb2 hash function (used for saving unique filenames)
unsigned long djb2hash(const std::string& s) {
	unsigned long hash = 5381; // Magic number
	for (auto c : s) {
		hash = (hash << 5) + hash + c; // same as: hash * 33 + c
	}
	return hash;
}

// Read CSV file
void ReadCSV(const std::string& inputfile,
             std::vector<std::vector<std::string>>& output) {
	std::fstream file;
	file.open(inputfile);
	std::string line;

	// Read every line from the stream
	while (getline(file, line)) {
		std::istringstream stream(line);
		std::vector<std::string> columns;
		std::string element;

		// Every line element separated by separator
		while (getline(stream, element, ',')) {
			columns.push_back(element);
		}
		output.push_back(columns);
	}
	file.close();
}


// Get JSON file input
std::string GetInputData(const std::string& inputfile) {

	// Check if exists
	if (gra::aux::FileExist(inputfile) == false) {
		std::string str = "gra::aux::GetInputData: Inputfile '" +
		                   inputfile + "' does not exist";
		throw std::invalid_argument(str);
	}
	// Create a JSON object from file
	std::ifstream ifs(inputfile);
	std::string data((std::istreambuf_iterator<char>(ifs)),
	                 (std::istreambuf_iterator<char>()));

	// https://json5.org/ features:
	// We use C++11 raw string literals here R"(text)"
	// so it is easier to type in expressions without problems with slashes

	// JSON5: <Single and multi-line comments are allowed>
	// First this. Remove C style block comments /* */
	data = std::regex_replace(data, std::regex(R"(\/\*(\*(?!\/)|[^*])*\*\/)"), "");
	// Second this. Remove C++ style single line comments "//"
	data = std::regex_replace(data, std::regex(R"([/]+.*)"), "");
	
	// Remove "\n" or "\t" or "\r", needed for the following regex operations
	data = std::regex_replace(data, std::regex(R"(\n|\t|\r)"), "");

	// JSON5: <Additional white space characters are allowed.>
	// Compress all extra spaces (compress to 1 space)
	// This also processes through the strings "", not problems for use,
	// but (IMPROVE THIS!)
	data = std::regex_replace(data, std::regex(R"(\s{2,})"), "");
	
	// JSON5: <Objects may have a single trailing comma.>
	// Replace ",}" or ", }" with "}"
	data = std::regex_replace(data, std::regex(R"(,}|, })"), "}");
	
	// JSON5: <Arrays may have a single trailing comma.>
	// Replace ",]" or ", ]" with "]"
	data = std::regex_replace(data, std::regex(R"(,]|, ])"), "]");
	
	//std::cout << std::endl << data << std::endl << std::endl;

	return data;
}

// Check is it integer digits
bool IsIntegerDigits(const std::string& str) {
	return str.find_first_not_of("-0123456789") == std::string::npos;
}

// Return particle parity as a string
std::string ParityToString(int value) {
	if (value > 0) {
		return "+";
	}
	else if (value == 0){
		return "";
	}
	else {
		return "-";
	}
}


// Return particle charge as a string
std::string Charge3XtoString(int q3) {
	const std::string sign = (q3 < 0) ? "-" : " ";
	const int absq3 = std::abs(q3);

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


// Return spin as a string
std::string Spin2XtoString(int J2) {
	if (J2 == 10)
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


// Split a string to strings
std::vector<std::string> SplitStr2Str(std::string input, const char delim) {
	std::vector<std::string> output;
	std::stringstream ss(input);

	// String by string
	while (ss.good()) {
		std::string substr;
		std::getline(ss, substr, delim);
		output.push_back(substr);
	}
	return output;
}

// Split string to ints
std::vector<int> SplitStr2Int(std::string input, const char delim) {
	std::vector<int> output;
	std::stringstream ss(input);

	// Get inputfiles by comma
	while (ss.good()) {
		std::string substr;
		std::getline(ss, substr, delim);
		output.push_back(std::stoi(substr));
	}
	return output;
}

// Trim extra space
void TrimExtraSpace(std::string& value) {
	value = std::regex_replace(value, std::regex("^ +| +$|( ) +"), "$1");
}

// Extract words from a string
std::vector<std::string> Extract(const std::string& str) {
	std::vector<std::string> words;
	std::stringstream ss(str);
	std::string buff;

	while (ss >> buff) {
		words.push_back(buff);
	}
	return words;
}

// Check if file exists
bool FileExist(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

void PrintWarning() {

std::cout << rang::fg::red << std::endl;
std::cout << "██╗    ██╗ █████╗ ██████╗ ███╗   ██╗██╗███╗   ██╗ ██████╗ " << std::endl;
std::cout << "██║    ██║██╔══██╗██╔══██╗████╗  ██║██║████╗  ██║██╔════╝ " << std::endl;
std::cout << "██║ █╗ ██║███████║██████╔╝██╔██╗ ██║██║██╔██╗ ██║██║  ███╗" << std::endl;
std::cout << "██║███╗██║██╔══██║██╔══██╗██║╚██╗██║██║██║╚██╗██║██║   ██║" << std::endl;
std::cout << "╚███╔███╔╝██║  ██║██║  ██║██║ ╚████║██║██║ ╚████║╚██████╔╝" << std::endl;
std::cout << " ╚══╝╚══╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝╚═╝  ╚═══╝ ╚═════╝ " << std::endl;
std::cout << rang::fg::reset << std::endl;

}


void PrintGameOver() {
	std::cout << std::endl;
	std::cout <<
	    " ██████╗  █████╗ ███╗   ███╗███████╗     ██████╗ ██╗   "
	    "██╗███████╗██████╗ " << std::endl;
	std::cout <<
	    "██╔════╝ ██╔══██╗████╗ ████║██╔════╝    ██╔═══██╗██║   "
	    "██║██╔════╝██╔══██╗" << std::endl;
	std::cout <<
	    "██║  ███╗███████║██╔████╔██║█████╗      ██║   ██║██║   ██║█████╗ "
	    " ██████╔╝" << std::endl;
	std::cout <<
	    "██║   ██║██╔══██║██║╚██╔╝██║██╔══╝      ██║   ██║╚██╗ ██╔╝██╔══╝ "
	    " ██╔══██╗" << std::endl;
	std::cout <<
	    "╚██████╔╝██║  ██║██║ ╚═╝ ██║███████╗    ╚██████╔╝ ╚████╔╝ "
	    "███████╗██║  ██║" << std::endl;
	std::cout <<
	    " ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝╚══════╝     ╚═════╝   ╚═══╝  "
	    "╚══════╝╚═╝  ╚═╝" << std::endl;
	std::cout << std::endl;
}

double GetVersion() {
	return 0.25;
}

std::string GetVersionString() {
	char buff[100];
	snprintf(buff, sizeof(buff), "Version %0.2f (BETA) 240319", GetVersion());
	std::string str = buff;
	return str;
}
std::string GetVersionTLatex() {
	char buff[100];
	snprintf(buff, sizeof(buff), "#color[15]{#scale[0.6]{GRANIITTI #scale[0.8]{%0.2f}#it{#beta}}}", GetVersion());
	std::string str = buff;
	return str;
}

std::string GetWebTLatex() {
	char buff[100];
	snprintf(buff, sizeof(buff), "#color[15]{#scale[0.5]{#LTgithub.com/mieskolainen#GT}}");
	std::string str = buff;
	return str;
}

void PrintVersion() {

	std::cout << GetVersionString() << std::endl;
	std::cout << rang::style::bold
	          << "Available at <github.com/mieskolainen/graniitti>"
	          << rang::style::reset << std::endl
	          << std::endl;
    std::cout << "References: arXiv:1811.01730 [hep-ph]" << std::endl;
    std::cout << std::endl;
	std::cout << "(c) 2017-2019 Mikael Mieskolainen" << std::endl;
	std::cout << "<mikael.mieskolainen@cern.ch>"     << std::endl;
	std::cout << std::endl;
	std::cout << "<opensource.org/licenses/GPL-3.0>" << std::endl;
	std::cout << "<opensource.org/licenses/MIT>"     << std::endl;
}

void PrintFlashScreen(rang::fg pcolor) {

	std::cout << std::endl;
	gra::aux::PrintBar("-");
	std::cout << pcolor;

	std::cout << ".``````````````````````     ````   ``           "
	             "             ``   ```   ``"
	          << std::endl;
	std::cout << ".```..````````````    ` ``          `           "
	             "                  ``    ``"
	          << std::endl;
	std::cout << ".``..:``````````` ```.``````.`.```              "
	             "                   `    ``"
	          << std::endl;
	std::cout << ".``.`.```````   ````     ``-.```.````           "
	             "                    ``` ``"
	          << std::endl;
	std::cout << "..-.`-.``    `.`   `..``---/--.......``         "
	             "             `````````````"
	          << std::endl;
	std::cout << "..:..-.`````.`   `.--oys+ohms:..-.`..```        "
	             "           ```````````````"
	          << std::endl;
	std::cout << ".-:.-...`..     .ooshydddddddh+::...``...       "
	             "        ``````````````````"
	          << std::endl;
	std::cout << "..-......       `yhhdhhmhddddmdy-...-...`       "
	             "            ``````````````"
	          << std::endl;
	std::cout << "::::::.`    `    `--/ohdhdhddyo-`   `..`..      "
	             "                         `"
	          << std::endl;
	std::cout << "so/-.`  `.`.-````-...++/syhdho`        `..      "
	             "                         `"
	          << std::endl;
	std::cout << "``..--.:o-..--..-/-::-.-:/ohh/           `      "
	             "                      `` `"
	          << std::endl;
	std::cout << "`.-::-.:yyh/`--.`s+ohso:/shhd:                  "
	             "                    ``  `-"
	          << std::endl;
	std::cout << ".:/-s+o+syhs--:/-..-oydyshhhho`                 "
	             "                ````` ``./"
	          << std::endl;
	std::cout << ".-+:-+shhhys+ysyh//sshyy/:oohyo.```````` `````` "
	             "````````````...::..  `.-oh"
	          << std::endl;
	std::cout << "-..-..:/+shhhyhhhyhysso.   "
	             "``:ys:----............------:::::++oo-  `.-:ohd"
	          << std::endl;
	std::cout << "--+:--//:::/+oyyshhs+`   ``   "
	             ".+so+++++++++++/-://+++++++++o+/-```-/yhhys+"
	          << std::endl;
	std::cout << "`/dysyhs/:::-:o+/-.`   `...    "
	             "`..::+oooooo+:-``.:+oooo+/:-.```.:/ssso/.``"
	          << std::endl;
	std::cout << "``:oossso-.--.-.    "
	             "..`...`..`.``...`-/ossys+ssoso+ys:--..::/:/--.-/"
	             "/:```."
	          << std::endl;
	std::cout << "```````` `   `-.` "
	             "`....```...:++.`./o:+-:/-:://"
	             "----...--.`-::::+ys/.`.:/sy"
	          << std::endl;
	std::cout << "```...-...`/...... ` `..`...`yo+:/o.:oshssy/.  "
	             "``  `.--.--:+sys+-` -+ossss"
	          << std::endl;
	std::cout << "``.//-/++-sh+-.-..`-`.-.--.- "
	             ":hhho/+//:/yyyhhs+::+o+/++ssssso:.`-::oyssyyy"
	          << std::endl;
	std::cout << "`/yhoshhhh////++:.....`:yy+/...sddhs/:////:-/"
	             "++ooo+oyssyoo:.`:/+oosssshddd"
	          << std::endl;
	std::cout << "`dddyhydhdyyddhddo``.+.`+/oyoh-..+yddho/"
	             "-......--:----://+osyyhho++/--::/+"
	          << std::endl;
	std::cout << "/hddhmhmhmdmddhdds:+s+:+/-:++hhh//.-/++osssso++/"
	             "-.````.-://::-.....`......"
	          << std::endl;
	std::cout << std::endl << std::endl;
	std::cout << rang::fg::reset;
}

// Print horizontal bar
void PrintBar(std::string str, unsigned int N) {

	for (std::size_t k = 0; k < N; ++k) {
		std::cout << str;
	}
	std::cout << std::endl;
}

// Create output directory if it does not exist
void CreateDirectory(std::string fullpath) {
	struct stat st = {0};
	if (stat(fullpath.c_str(), &st) == -1) {
	    mkdir(fullpath.c_str(), 0700);
	}
}

// Get commandline arguments which are split by @... @... tagging syntax
std::vector<OneCMD> SplitCommands(const std::string& fullstr) {

	// Find @ blocks [ ... ]
	std::vector<std::size_t> pos = FindOccurance(fullstr, "@");
	std::vector<std::string> subcmd;

	if (pos.size() != 0) {

		// Take full string from the first @ till end
		std::string str = fullstr.substr(pos[0]);

		// Split multiple @id{} @id{} ... @id{} occurances
		if (pos.size() >= 2) {
			for (std::size_t i = 0; i < pos.size(); ++i) {
				const int start = pos[i] + 1; // +1 so we skip @
				const int end   = (i < pos.size()-1) ? pos[i+1] : str.size();
				subcmd.push_back(str.substr(start, end - start + 1));
			}
		} else {
			subcmd.push_back(str.substr(1)); // +1 so we skip @
		}
	}

	// Loop over @ blocks
	std::vector<OneCMD> cmd;

	for (std::size_t i = 0; i < subcmd.size(); ++i) {

		// ** Remove whitespace **
		TrimExtraSpace(subcmd[i]);

		std::string id;
		std::map<std::string, std::any> arg;
		
		// Check we find {} brackets
		std::size_t Lpos = subcmd[i].find("{",0);
		std::size_t Rpos = subcmd[i].find("}",0);

		// Singlet command @ID:VALUE
		if (Lpos == std::string::npos && Rpos == std::string::npos) {

			std::vector<std::string> strip = SplitStr2Str(subcmd[i], ':');

			// ID
			id = strip[0];

			if (strip.size() == 1) {
				arg["_SINGLET_"] = "true"; // default true, no : given
			}
			else if (strip.size() == 2) {
				arg["_SINGLET_"] = strip[1];
			} else {
				throw std::invalid_argument("gra::SplitCommands: @Syntax invalid with '" + subcmd[i] + "'");
			}
			
		// Block command @blaa{key:val,key:val,...}
		} else { // We have {} block

			// ID
			id = subcmd[i].substr(0, Lpos);

			// Content {}

			// Now split all arguments inside {} by comma ','
			std::vector<std::string> keyvals = SplitStr2Str(subcmd[i].substr(Lpos+1, Rpos-Lpos-1), ',');

			// Add all
			for (std::size_t i = 0; i < keyvals.size(); ++i) {

				// Split using syntax definition: key:val
				std::vector<std::string> strip = SplitStr2Str(keyvals[i], ':');
				if (strip.size() != 2) {
					throw std::invalid_argument("gra::SplitCommands: @Syntax not good with '"+ keyvals[i] + "'");
				}
				arg[strip[0]] = strip[1]; // add to map
			}
		}

		// Add this command block
		OneCMD o;
		o.id = id;
		o.arg = arg;
		o.Print(); // For debug
		cmd.push_back(o);
	}
	return cmd;
}

// Find string occurances
std::vector<std::size_t> FindOccurance(const std::string& str, const std::string& sub) {

	// Holds all the positions that sub occurs within str
	std::vector<size_t> positions;
	std::size_t pos = str.find(sub, 0);
   	
	while(pos != std::string::npos) {
	    positions.push_back(pos);
	    pos = str.find(sub,pos+1);
	}
	return positions;
}


} // aux namespace ends
} // gra namespace ends
