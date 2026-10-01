// GRANIITTI shared model state and model-data path resolution
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/MGlobals.h"

#include <exception>

#include "Graniitti/Tech/MAux.h"

namespace gra {

// Default model tune for standalone particle and resonance readers
std::string MODELPARAM;

// Compute true when a model tune names an explicit filesystem path
bool IsExplicitModelTunePath(const std::string &modelparam) {
  return !modelparam.empty() && (modelparam[0] == '/' || modelparam[0] == '.' ||
                                 modelparam.find('/') != std::string::npos ||
                                 modelparam.find('\\') != std::string::npos);
}

// Resolve a model tune name or explicit path to its directory
std::string ResolveModelTuneDir(const std::string &modelparam) {
  return IsExplicitModelTunePath(modelparam)
             ? modelparam
             : gra::aux::ResolveProjectPath("modeldata/" + modelparam);
}

// Resolve one file inside a named or explicit model tune directory
std::string ResolveModelDataFile(const std::string &modelparam,
                                 const std::string &filename) {
  return ResolveModelTuneDir(modelparam) + "/" + filename;
}

// Multithreading state
std::mutex g_mutex;

} // namespace gra
