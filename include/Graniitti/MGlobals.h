// GRANIITTI modeldata path resolution and output lock
// 
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGLOBALS_H
#define MGLOBALS_H

// C++
#include <mutex>
#include <string>

namespace gra {
// ======================================================================
// Model-data paths are resolved independently of the generator frontend

// Model tune
extern std::string MODELPARAM;
bool IsExplicitModelTunePath(const std::string &modelparam);
std::string ResolveModelTuneDir(const std::string &modelparam);
std::string ResolveModelDataFile(const std::string &modelparam,
                                 const std::string &filename);

// Multithreading lock
extern std::mutex g_mutex;

// ======================================================================

} // namespace gra

#endif
