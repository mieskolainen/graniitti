// (Sub)-Processes and Amplitude containers
//
// (c) 2017-2021 Mikael Mieskolainen
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
std::vector<std::shared_ptr<MProc>> MSubProc::CreateAllProcesses() const {
  std::vector<std::shared_ptr<MProc>> v;

  v.push_back(std::shared_ptr<MProc>{new PROC_0()});
  v.push_back(std::shared_ptr<MProc>{new PROC_1()});
  v.push_back(std::shared_ptr<MProc>{new PROC_2()});
  v.push_back(std::shared_ptr<MProc>{new PROC_3()});
  v.push_back(std::shared_ptr<MProc>{new PROC_4()});
  v.push_back(std::shared_ptr<MProc>{new PROC_5()});
  v.push_back(std::shared_ptr<MProc>{new PROC_6()});
  v.push_back(std::shared_ptr<MProc>{new PROC_7()});
  v.push_back(std::shared_ptr<MProc>{new PROC_8()});
  v.push_back(std::shared_ptr<MProc>{new PROC_9()});
  v.push_back(std::shared_ptr<MProc>{new PROC_10()});
  v.push_back(std::shared_ptr<MProc>{new PROC_11()});
  v.push_back(std::shared_ptr<MProc>{new PROC_12()});
  v.push_back(std::shared_ptr<MProc>{new PROC_13()});
  v.push_back(std::shared_ptr<MProc>{new PROC_14()});
  v.push_back(std::shared_ptr<MProc>{new PROC_15()});
  v.push_back(std::shared_ptr<MProc>{new PROC_16()});
  v.push_back(std::shared_ptr<MProc>{new PROC_17()});
  v.push_back(std::shared_ptr<MProc>{new PROC_18()});
  v.push_back(std::shared_ptr<MProc>{new PROC_19()});
  v.push_back(std::shared_ptr<MProc>{new PROC_20()});
  v.push_back(std::shared_ptr<MProc>{new PROC_21()});
  v.push_back(std::shared_ptr<MProc>{new PROC_22()});
  v.push_back(std::shared_ptr<MProc>{new PROC_23()});
  v.push_back(std::shared_ptr<MProc>{new PROC_24()});
  v.push_back(std::shared_ptr<MProc>{new PROC_25()});
  v.push_back(std::shared_ptr<MProc>{new PROC_26()});
  v.push_back(std::shared_ptr<MProc>{new PROC_27()});

  return v;
}


// Generic constructor, used for the first initialization
MSubProc::MSubProc(const std::vector<std::string>& istate, const std::string& mc) {
  for (const auto& i : aux::indices(istate)) { ConstructDescriptions(istate[i], mc); }
}


// Set spesific initial state and channel
void MSubProc::Initialize(const std::string& istate, const std::string& channel) {
  ISTATE  = istate;
  CHANNEL = channel;
}


// Construct textual descriptions
void MSubProc::ConstructDescriptions(const std::string& istate, const std::string& mc) {
  const std::vector<std::shared_ptr<MProc>> p = CreateAllProcesses();

  for (const auto& i : aux::indices(p)) {
    if (p[i]->ISTATE == istate) {
      Processes.insert(
          std::make_pair(p[i]->ISTATE + "[" + p[i]->CHANNEL + "]<" + mc + ">", p[i]->DESCRIPTION));
    }
  }
}


// Activate spesific process
void MSubProc::ActivateProcess() {
  std::vector<std::shared_ptr<MProc>> p = CreateAllProcesses();

  while (true) {
    for (const auto& i : aux::indices(p)) {
      if (p[i]->ISTATE == ISTATE && p[i]->CHANNEL == CHANNEL) {
        // Fine
      } else {
        p.erase(p.begin() + i);
        break;
      }
    }
    if (p.size() <= 1) { break; }
  }
  if (p.size() == 0) {
    throw std::invalid_argument("MSubProc::ActivateProcess: Unknown ISTATE = " + ISTATE +
                                " or CHANNEL = " + CHANNEL);
  }

  // move it!
  pr = std::move(p[0]);
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
  if (pr == nullptr) { ActivateProcess(); }
  return pr->Amp2(lts);
}

}  // namespace gra
