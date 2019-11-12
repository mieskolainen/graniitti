// GRANIITTI Monte Carlo all global variables collected here
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGLOBALS_H
#define MGLOBALS_H

// Own
#include "Graniitti/MSudakov.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"


namespace gra {
// ======================================================================
// These variables are initialized by MGraniitti.cc

// Model tune
extern std::string MODELPARAM;

// Multithreading lock
extern std::mutex g_mutex;

// For multithreaded VEGAS, to handle the exceptions from forked threads
extern std::exception_ptr globalExceptionPtr;

// ======================================================================

}  // namespace gra

#endif