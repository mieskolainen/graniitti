// (Sub)-Processes and Amplitudes
// 
// (c) 2017-2019 Mikael Mieskolainen
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

// Proper constructor
MSubProc::MSubProc(const std::string &_ISTATE, const std::string &_CHANNEL) {
  ISTATE  = _ISTATE;
  CHANNEL = _CHANNEL;

  ConstructDescriptions(ISTATE);
}

// Constructor
MSubProc::MSubProc(const std::vector<std::string> &first) {
  for (const auto& i : aux::indices(first)) {
    ConstructDescriptions(first[i]);
  }
}

// Destructor
MSubProc::~MSubProc() {
  DeleteProcesses();
}

// Delete all processes
void MSubProc::DeleteProcesses() {
  for (const auto& i : aux::indices(pr)) {
    delete pr[i];
  }
  pr.clear(); // Finally empty
}

void MSubProc::ConstructDescriptions(const std::string &first) {

  CreateProcesses();

  std::map<std::string, std::string> channels;
  for (const auto & i : aux::indices(pr)) {
    if (pr[i]->ISTATE == first) {
      channels.insert(std::pair<std::string, std::string>(pr[i]->CHANNEL, pr[i]->DESCRIPTION));
    }
  }
  descriptions.insert(std::pair<std::string, std::map<std::string, std::string>>(first, channels));

  DeleteProcesses(); // Important!
}

// Construct process
void MSubProc::ActivateProcess() {

  CreateProcesses();

  while (true) {
    for (std::size_t i = 0; i < pr.size(); ++i) {
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
    throw std::invalid_argument("MSubProc::ActivateProcess: Unknown ISTATE = " + ISTATE + " or CHANNEL = " + CHANNEL);
  }
}

// Wrapper function
double MSubProc::GetBareAmplitude2(gra::LORENTZSCALAR &lts) {

  if (pr.size() == 0) {
      ActivateProcess();
  }
  return pr[0]->Amp2(lts);
}

}  // namespace gra
