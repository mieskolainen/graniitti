// Load complex fixed-helicity MadLoop amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_MadLoop.h"
#include "Graniitti/Kinematics/M4Vec.h"

#include <dlfcn.h>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

#if defined(__GLIBC__)
#include <link.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "json.hpp"

namespace gra::mg5 {
namespace {

// Detect unresolved t or u in massless two body scattering before loop reduction
bool SingularMasslessPoint(const std::vector<double>& momentum) {
  if (momentum.size() != 16) { return false; }
  std::array<M4Vec, 4> p;
  for (std::size_t i = 0; i < p.size(); ++i) {
    const std::size_t k = 4 * i;
    p[i]                = M4Vec(momentum[k + 1], momentum[k + 2], momentum[k + 3], momentum[k]);
  }
  const double s = (p[0] + p[1]).M2();
  if (!std::isfinite(s) || s <= 0.0) { return true; }
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * s;
  for (const auto& leg : p) {
    if (std::abs(leg.M2()) > tolerance) { return false; }
  }
  // CutTools cannot construct its null basis when a massless channel is unresolved
  return std::min(std::abs((p[0] - p[2]).M2()), std::abs((p[0] - p[3]).M2())) <= tolerance;
}

// Compute true when one regular runtime data file can be opened
bool Readable(const std::string &path) {
  return std::ifstream(path, std::ios::binary).good();
}

// Compute one process-local lock path for a serialized loader fallback
std::string LockPath(const std::string &name) {
  return "/tmp/graniitti_madloop_" + name + "_" +
         std::to_string(static_cast<unsigned long>(getuid())) + "_" +
         std::to_string(static_cast<unsigned long>(getpid())) + ".lock";
}

} // namespace

// Load and initialize one generated MadLoop runtime directory
MadLoop::MadLoop(const std::string &directory, const std::string &name,
                 std::size_t nexternal, std::size_t nhelicity,
                 std::size_t ncolor)
    : name_(name), nexternal_(nexternal), nhelicity_(nhelicity),
      ncolor_(ncolor) {
  if (name_.empty() || nexternal_ < 3 || nhelicity_ == 0 || ncolor_ == 0) {
    error_ = "Invalid MadLoop runtime dimensions";
    return;
  }
  Load(directory);
}

// Release the generated runtime and serialization handle
MadLoop::~MadLoop() {
  if (handle_ != nullptr) {
    dlclose(handle_);
  }
  if (lock_fd_ >= 0) {
    close(lock_fd_);
  }
}

// Compute an empty string when the runtime files are readable
std::string MadLoop::Check(const std::string &directory) {
  const std::string library = directory + "/libMadLoopAdapter.so";
  const std::string parameter_card = directory + "/param_card.dat";
  const std::string resource_card =
      directory + "/MadLoop5_resources/MadLoopParams.dat";
  if (!Readable(library) || !Readable(parameter_card) || !Readable(directory + "/runtime.json") ||
      !Readable(resource_card)) {
    return "MadLoop runtime is incomplete below " + directory;
  }
  return {};
}

// Check whether the generated runtime was loaded and initialized
bool MadLoop::Ready() const noexcept { return ready_; }

// Compute the load or initialization diagnostic
std::string MadLoop::Error() const {
  std::lock_guard<std::mutex> guard(mutex_);
  return error_;
}

// Compute the runtime parameter-card strong coupling
double MadLoop::AlphaS() const {
  std::lock_guard<std::mutex> guard(mutex_);
  return alpha_s_;
}

// Compute the runtime parameter-card electromagnetic coupling
double MadLoop::AlphaQED() const {
  std::lock_guard<std::mutex> guard(mutex_);
  return alpha_qed_;
}

// Compute the physical model poles for phase space
ParticleMap MadLoop::Particles() const {
  std::lock_guard<std::mutex> guard(mutex_);
  if (!ready_) { throw std::invalid_argument(error_); }
  return particles_;
}

// Read or update evaluated model parameters under the runtime lock
void MadLoop::Model(int update) {
  std::vector<double> poles(2 * pdgs_.size());
  std::array<double, 2> alpha{};
  model_(values_.data(), poles.data(), alpha.data(), update);
  for (std::size_t i = 0; i < pdgs_.size(); ++i) {
    auto &particle = particles_.at(pdgs_[i]);
    particle.mass = particle.signed_mass ? std::abs(poles[2 * i]) : poles[2 * i];
    particle.width = particle.signed_mass ? std::abs(poles[2 * i + 1]) : poles[2 * i + 1];
  }
  alpha_s_ = alpha[0];
  alpha_qed_ = alpha[1];
}

// Initialize all external UFO parameters from their generated SLHA indices
void MadLoop::InitParameters(const SLHAReader &card) {
  std::lock_guard<std::mutex> guard(mutex_);
  if (!ready_ || !Lock()) { throw std::invalid_argument(error_); }
  for (std::size_t i = 0; i < parameters_.size(); ++i) {
    const auto &[block, indices] = parameters_[i];
    values_[i] = card.get_block_entry(block, indices, defaults_[i]);
  }
  Model(1);
  Unlock();
  ValidateModel(card, particles_);
  if ((scale_qcd_ && (!std::isfinite(alpha_s_) || alpha_s_ <= 0.0)) ||
      (scale_qed_ && (!std::isfinite(alpha_qed_) || alpha_qed_ <= 0.0))) {
    throw std::invalid_argument("MadLoop coupling rescaling requires positive UFO reference couplings");
  }
}

// Load symbols and initialize immutable model and resource data
void MadLoop::Load(const std::string &directory) {
  const std::string library = directory + "/libMadLoopAdapter.so";
  const std::string parameter_card = directory + "/param_card.dat";
  const std::string resource = directory + "/MadLoop5_resources";
  error_ = Check(directory);
  if (!error_.empty()) {
    return;
  }
  try {
    const auto data = nlohmann::json::parse(std::ifstream(directory + "/runtime.json"));
    scale_qcd_ = data.at("qcd_order").get<int>() > 0;
    scale_qed_ = data.at("qed_order").get<int>() > 0 && data.at("alpha_zero").get<bool>();
    for (const auto &parameter : data.at("parameters")) {
      parameters_.emplace_back(parameter.at("block").get<std::string>(),
                               parameter.at("indices").get<std::vector<int>>());
    }
    values_.resize(parameters_.size());
    for (const auto &particle : data.at("particles")) {
      const int pdg = particle.at("pdg").get<int>();
      pdgs_.push_back(pdg);
      particles_.emplace(pdg, ParticleParam{0.0, 0.0, particle.at("signed_mass").get<bool>()});
    }
  } catch (const std::exception &exception) {
    error_ = "Invalid MadLoop model metadata: " + std::string(exception.what());
    return;
  }

#if defined(__GLIBC__)
  dlerror();
  handle_ = dlmopen(LM_ID_NEWLM, library.c_str(), RTLD_NOW | RTLD_LOCAL);
  isolated_ = handle_ != nullptr;
#endif
  if (handle_ == nullptr) {
    dlerror();
    handle_ = dlopen(library.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (handle_ == nullptr) {
      const char *message = dlerror();
      error_ = message == nullptr ? "Could not load MadLoop runtime" : message;
      return;
    }
    lock_fd_ =
        open(LockPath(name_).c_str(), O_CREAT | O_RDWR | O_CLOEXEC, 0600);
    if (lock_fd_ < 0) {
      error_ = "Could not open the MadLoop serialization lock";
      return;
    }
  }

  dlerror();
  init_ = reinterpret_cast<InitFunction>(dlsym(handle_, "gra_madloop_init"));
  eval_ = reinterpret_cast<EvalFunction>(dlsym(handle_, "gra_madloop_eval"));
  model_ = reinterpret_cast<ModelFunction>(dlsym(handle_, "gra_madloop_model"));
  const auto sizes = reinterpret_cast<SizeFunction>(dlsym(handle_, "gra_madloop_sizes"));
  const char *symbol_error = dlerror();
  if (symbol_error != nullptr || init_ == nullptr || eval_ == nullptr || model_ == nullptr || sizes == nullptr) {
    error_ =
        symbol_error == nullptr ? "MadLoop C ABI is incomplete" : symbol_error;
    return;
  }

  std::array<int, 5> dimensions{};
  sizes(dimensions.data());
  const std::array<std::size_t, 5> expected{nexternal_, nhelicity_, ncolor_, values_.size(), pdgs_.size()};
  for (std::size_t i = 0; i < dimensions.size(); ++i) {
    if (dimensions[i] < 0 || static_cast<std::size_t>(dimensions[i]) != expected[i]) {
      error_ = "MadLoop model metadata disagrees with the generated library dimensions";
      return;
    }
  }

  if (!Lock()) {
    return;
  }
  int status = -1;
  init_(parameter_card.data(), static_cast<int>(parameter_card.size()),
        resource.data(), static_cast<int>(resource.size()), &status);
  if (status == 0) {
    Model(0);
    defaults_ = values_;
  }
  Unlock();
  if (status != 0) {
    error_ =
        "MadLoop initialization failed with status " + std::to_string(status);
    return;
  }
  try {
    ready_ = true;
    InitParameters(SLHAReader(parameter_card));
  } catch (const std::exception &exception) {
    ready_ = false;
    error_ = exception.what();
    return;
  }
}

// Acquire the process-wide fallback lock for shared Fortran state
bool MadLoop::Lock() {
  if (isolated_) {
    return true;
  }
  if (lock_fd_ < 0 || flock(lock_fd_, LOCK_EX) != 0) {
    error_ = "Could not lock the shared MadLoop runtime";
    return false;
  }
  return true;
}

// Release the process-wide fallback lock
void MadLoop::Unlock() {
  if (!isolated_ && lock_fd_ >= 0) {
    (void)flock(lock_fd_, LOCK_UN);
  }
}

// Evaluate every complex color and helicity amplitude
bool MadLoop::Evaluate(const std::vector<double> &momentum,
                       MadLoopResult &result) {
  result = {};
  if (!ready_) {
    return false;
  }
  std::lock_guard<std::mutex> guard(mutex_);
  error_.clear();
  if (momentum.size() != 4 * nexternal_ ||
      !std::all_of(momentum.cbegin(), momentum.cend(),
                   [](double value) { return std::isfinite(value); })) {
    error_ = "MadLoop received invalid momentum";
    return false;
  }
  if (SingularMasslessPoint(momentum)) {
    error_ = "MadLoop received numerically singular massless two body kinematics";
    return false;
  }
  if (!Lock()) {
    return false;
  }

  // Restore this model when dlopen shares Fortran COMMON blocks
  if (!isolated_) { Model(1); }
  const std::size_t amplitudes = nhelicity_ * ncolor_;
  std::vector<double> amplitude_re(amplitudes, 0.0);
  std::vector<double> amplitude_im(amplitudes, 0.0);
  std::vector<double> pole_re(2 * amplitudes, 0.0);
  std::vector<double> pole_im(2 * amplitudes, 0.0);
  result.squared.assign(nhelicity_, 0.0);
  result.return_code.assign(nhelicity_, 0);
  int status = -1;
  eval_(momentum.data(), amplitude_re.data(), amplitude_im.data(),
        pole_re.data(), pole_im.data(), result.squared.data(),
        result.return_code.data(), &status);
  Unlock();
  if (status != 0) {
    error_ = "MadLoop evaluation failed with status " + std::to_string(status);
    return false;
  }

  result.amplitude.resize(amplitudes);
  result.pole.resize(2 * amplitudes);
  for (std::size_t index = 0; index < amplitudes; ++index) {
    result.amplitude[index] = {amplitude_re[index], amplitude_im[index]};
    for (std::size_t order = 0; order < 2; ++order) {
      const std::size_t offset = 2 * index + order;
      result.pole[offset] = {pole_re[offset], pole_im[offset]};
    }
  }
  return Validate(result);
}

// Validate finite amplitudes and cancellation of loop poles
bool MadLoop::Validate(const MadLoopResult &result) {
  for (std::size_t index = 0; index < result.amplitude.size(); ++index) {
    const double scale = std::max(1.0, std::abs(result.amplitude[index]));
    if (!std::isfinite(result.amplitude[index].real()) ||
        !std::isfinite(result.amplitude[index].imag()) ||
        !std::isfinite(result.pole[2 * index].real()) ||
        !std::isfinite(result.pole[2 * index].imag()) ||
        !std::isfinite(result.pole[2 * index + 1].real()) ||
        !std::isfinite(result.pole[2 * index + 1].imag()) ||
        std::abs(result.pole[2 * index]) > 1.0e-7 * scale ||
        std::abs(result.pole[2 * index + 1]) > 1.0e-7 * scale) {
      error_ = "MadLoop returned nonfinite data or an uncancelled pole";
      return false;
    }
  }
  return std::all_of(result.squared.cbegin(), result.squared.cend(),
                     [](double value) { return std::isfinite(value); });
}

} // namespace gra::mg5
