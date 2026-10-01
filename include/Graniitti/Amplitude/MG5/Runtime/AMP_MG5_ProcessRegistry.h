// Generated MadGraph amplitude process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PROCESSREGISTRY_H
#define AMP_MG5_PROCESSREGISTRY_H

#include <optional>
#include <string>
#include <vector>

#include "Graniitti/Process/MAmpMatch.h"

namespace gra::amplitude {

// Compute all generated MG5 amplitude processes by value
std::vector<Process> AllProcesses();

// Compute processes in one generated amplitude family by value
std::vector<Process> Processes(
    const std::string &process_family);

// Compute one exact generated process as a vector
std::vector<Process> Processes(
    const std::string &process_family, const std::string &process_name);

// Find one generated process by its family and process names
std::optional<Process> FindProcess(
    const std::string &process_family, const std::string &process_name);

// Find one generated process by typed topology and ordered stable PDGs
std::optional<Process> FindProcess(
    const std::string &process_family, const AmplitudeTopology &topology,
    const std::vector<int> &stable_pdgs);

// Compute the parameter card owned by one exact generated process
std::optional<std::string> ParameterCard(
    const std::string &process_family, const std::string &process_name);

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_ddbar
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_ddbar()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_ddbar")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_ssbar
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_ssbar()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_ssbar")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_ccbar
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_ccbar()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_ccbar")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_bbbar
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_bbbar()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_bbbar")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_uubar
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_uubar()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_uubar")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_uubarg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_uubarg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_uubarg")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_ddbarg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_ddbarg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_ddbarg")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_ssbarg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_ssbarg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_ssbarg")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_ccbarg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_ccbarg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_ccbarg")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_bbbarg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_bbbarg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_bbbarg")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_gg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_gg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_gg")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_ggg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_ggg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_ggg")) {}
};

// Own converter-generated processes for PHOTON
class MG5ProcessRegistry_yy_ll
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_yy_ll()
      : ProcessRegistry(
            ::gra::amplitude::Processes("PHOTON", "yy_ll")) {}
};

// Own converter-generated processes for PHOTON
class MG5ProcessRegistry_yy_uubarg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_yy_uubarg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("PHOTON", "yy_uubarg")) {}
};

// Own converter-generated processes for DURHAM
class MG5ProcessRegistry_gg_gggg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_gg_gggg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("DURHAM", "gg_gggg")) {}
};

// Own converter-generated processes for PHOTON
class MG5ProcessRegistry_yy_uubargg
    : public ProcessRegistry {
 public:
  // Construct the immutable generated process view
  MG5ProcessRegistry_yy_uubargg()
      : ProcessRegistry(
            ::gra::amplitude::Processes("PHOTON", "yy_uubargg")) {}
};

}  // namespace gra::amplitude

#endif
