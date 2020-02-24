// (Sub)-Processes and Amplitudes
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MKinematics.h"
#include "Graniitti/MSubProc.h"

// Libraries
#include "rang.hpp"

namespace gra {

// Add all processes
void MSubProc::CreateProcesses() {
  DeleteProcesses();

  pr.push_back(new PROC_0());
  pr.push_back(new PROC_1());
  pr.push_back(new PROC_2());
  pr.push_back(new PROC_3());
  pr.push_back(new PROC_4());
  pr.push_back(new PROC_5());
  pr.push_back(new PROC_6());
  pr.push_back(new PROC_7());
  pr.push_back(new PROC_8());
  pr.push_back(new PROC_9());
  pr.push_back(new PROC_10());
  pr.push_back(new PROC_11());
  pr.push_back(new PROC_12());
  pr.push_back(new PROC_13());
  pr.push_back(new PROC_14());
  pr.push_back(new PROC_15());
  pr.push_back(new PROC_16());
  pr.push_back(new PROC_17());
  pr.push_back(new PROC_18());
  pr.push_back(new PROC_19());
  pr.push_back(new PROC_20());
  pr.push_back(new PROC_21());
  pr.push_back(new PROC_22());
  pr.push_back(new PROC_23());
  pr.push_back(new PROC_24());
  pr.push_back(new PROC_25());
}

// Generic constructor, used for the first initialization
MSubProc::MSubProc(const std::vector<std::string>& istate, const std::string& mc) {
  for (const auto& i : aux::indices(istate)) { ConstructDescriptions(istate[i], mc); }
}

// Set spesific initial state and channel
void MSubProc::Initialize(const std::string& istate, const std::string& channel) {
  ISTATE  = istate;
  CHANNEL = channel;
  DeleteProcesses();
}

// Destructor
MSubProc::~MSubProc() { DeleteProcesses(); }

// Delete all processes
void MSubProc::DeleteProcesses() {
  for (const auto& i : aux::indices(pr)) { delete pr[i]; }
  pr.clear();  // Finally empty
}

void MSubProc::ConstructDescriptions(const std::string& istate, const std::string& mc) {
  CreateProcesses();

  // Construct labels
  for (const auto& i : aux::indices(pr)) {
    if (pr[i]->ISTATE == istate) {
      Processes.insert(std::make_pair(pr[i]->ISTATE + "[" + pr[i]->CHANNEL + "]<" + mc + ">",
                                      pr[i]->DESCRIPTION));
    }
  }

  DeleteProcesses();  // Important!
}

// Activate spesific process
void MSubProc::ActivateProcess() {
  CreateProcesses();

  while (true) {
    for (const auto& i : aux::indices(pr)) {
      if (pr[i]->ISTATE == ISTATE && pr[i]->CHANNEL == CHANNEL) {
        // Fine
      } else {
        delete pr[i];
        pr.erase(pr.begin() + i);
        break;
      }
    }
    if (pr.size() <= 1) { break; }
  }
  if (pr.size() == 0) {
    throw std::invalid_argument("MSubProc::ActivateProcess: Unknown ISTATE = " + ISTATE +
                                " or CHANNEL = " + CHANNEL);
  }
}

// Print out process lists
std::vector<std::string> MSubProc::PrintProcesses() const {
  // Iterate through processes
  std::vector<std::string>                                        procstr;
  std::map<std::string, std::vector<std::string>>::const_iterator it = Processes.begin();

  while (it != Processes.end()) {
    procstr.push_back(it->first);

    if (it->second.size() == 0) {
      printf("%25s  =      ", it->first.c_str());
    } else if (it->second.size() == 1) {
      printf("%25s  =  %-43s", it->first.c_str(), it->second[0].c_str());
    } else if (it->second.size() == 2) {
      printf("%25s  =  %-43s  |  %-20s", it->first.c_str(), it->second[0].c_str(),
             it->second[1].c_str());
    } else if (it->second.size() >= 3) {
      printf("%25s  =  %-43s  |  %-20s  | %-10s", it->first.c_str(), it->second[0].c_str(),
             it->second[1].c_str(), it->second[2].c_str());
    }
    std::cout << "\n";
    ++it;
  }
  return procstr;
}

// Check if process exists
bool MSubProc::ProcessExist(std::string str) const {
  if (Processes.find(str) != Processes.end()) { return true; }
  return false;
}

// Return process description string
std::vector<std::string> MSubProc::GetProcessDescriptor(std::string str) const {
  if (!ProcessExist(str)) {
    throw std::invalid_argument("MSubProc::GetProcessDescriptor: Process by name " + str +
                                " does not exist");
  }
  return Processes.find(str)->second;
}

// Wrapper function
double MSubProc::GetBareAmplitude2(gra::LORENTZSCALAR& lts) {
  if (pr.size() == 0) { ActivateProcess(); }
  return pr[0]->Amp2(lts);
}

}  // namespace gra
