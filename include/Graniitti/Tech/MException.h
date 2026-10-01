// Custom library exceptions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEXCEPTION_H
#define MEXCEPTION_H

// C++
#include <stdexcept>
#include <string>

namespace gra {

// Signal one event-local amplitude failure that sampling must reject
class AmplitudeFailure final : public std::runtime_error {
 public:
  // Construct an amplitude failure with diagnostic context
  explicit AmplitudeFailure(const std::string &message) : std::runtime_error(message) {}
};

// Signal one exceptional event-local phase-space construction failure
class PhaseSpaceFailure final : public std::runtime_error {
 public:
  // Construct a phase-space failure with diagnostic context
  explicit PhaseSpaceFailure(const std::string &message) : std::runtime_error(message) {}
};

}  // namespace gra

#endif
