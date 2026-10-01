// Immutable GRANIITTI model tune snapshot
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMODELTUNE_H
#define MMODELTUNE_H

// C++
#include <filesystem>
#include <map>
#include <memory>
#include <string>

// Own
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Regge/MSoftModel.h"

// Libraries
#include "json.hpp"

namespace gra {

// Program wide numerical controls shared by all amplitude models
struct MGlobalNumerics {
  double coupling_min;
};

// Store the shared charged Standard Model masses in GeV and inverse electromagnetic coupling
struct SMParam {
  double e, mu, tau;
  double u, d, s, c, b, t;
  double w;
  double alpha_em_inv;
};

// Immutable GENERAL and NUMERICS model tune
class MModelTune {
 public:
  // Load and validate one complete model tune from its GENERAL file
  static std::shared_ptr<const MModelTune> Load(const std::string &general_file);

  // Select another soft model in an independent immutable tune
  std::shared_ptr<const MModelTune> WithSoft(const std::string &name) const;

  // Compute the GENERAL source path
  const std::string &GeneralFile() const noexcept { return general_file_; }

  // Compute the directory containing the model and decay cards
  std::string Directory() const { return std::filesystem::path(general_file_).parent_path().string(); }

  // Compute the NUMERICS source path
  const std::string &NumericsFile() const noexcept { return numerics_file_; }

  // Compute one immutable GENERAL block
  const nlohmann::json &General(const std::string &name) const;

  // Compute one immutable NUMERICS block
  const nlohmann::json &Numerics(const std::string &name) const;

  // Compute the complete immutable GENERAL document
  const nlohmann::json &General() const noexcept { return general_; }

  // Compute the complete immutable NUMERICS document
  const nlohmann::json &Numerics() const noexcept { return numerics_; }

  // Compute one immutable MP, XP, GP or TP continuum steering card
  const nlohmann::json &Continuum(const std::string &model) const;

  // Access the immutable soft physics snapshot
  const SoftModelPtr &Soft() const noexcept { return soft_; }

  // Compute the program wide numerical controls
  const MGlobalNumerics &Global() const noexcept { return global_; }

  // Compute the immutable Standard Model inputs
  const SMParam &SM() const noexcept { return sm_; }

  // Access the immutable proton-structure parameter selection
  const form::ParamStore &Structure() const noexcept { return structure_; }

  // Compute the immutable flat-amplitude parameters
  const form::FlatParam &Flat() const noexcept { return flat_; }

 private:
  // Construct one fully parsed immutable tune
  MModelTune(std::string general_file, std::string numerics_file, nlohmann::json general, nlohmann::json numerics,
             SoftModelPtr soft, MGlobalNumerics global, SMParam sm, form::ParamStore structure, form::FlatParam flat,
             std::map<std::string, nlohmann::json> continuum);

  std::string                           general_file_;
  std::string                           numerics_file_;
  nlohmann::json                        general_;
  nlohmann::json                        numerics_;
  SoftModelPtr                          soft_;
  MGlobalNumerics                       global_;
  SMParam                              sm_;
  form::ParamStore                      structure_;
  form::FlatParam                       flat_;
  std::map<std::string, nlohmann::json> continuum_;
};

using MModelTunePtr = std::shared_ptr<const MModelTune>;

}  // namespace gra

#endif
