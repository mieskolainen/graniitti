// Persistent storage helpers for eikonal numerical tables
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEIKONAL_CACHE_H
#define MEIKONAL_CACHE_H

// C++
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

// Libraries
#include "Graniitti/Tech/MAux.h"
#include "json.hpp"

// POSIX file locking and temporary cache files
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace gra::eikonal {

// Hold one exclusive interprocess lock while a table is built
class MCacheLock {
 public:
  // Acquire one exclusive cache file lock
  explicit MCacheLock(const std::string& filename) {
    std::error_code error;
    aux::CreateDirectory(std::filesystem::path(filename).parent_path().string(), &error);
    if (error) { return; }
    descriptor_ = ::open(filename.c_str(), O_CREAT | O_RDWR, 0600);
    if (descriptor_ < 0) { return; }
    int status = 0;
    do { status = ::flock(descriptor_, LOCK_EX); } while (status != 0 && errno == EINTR);
    if (status != 0) {
      ::close(descriptor_);
      descriptor_ = -1;
    }
  }

  // Release the cache file lock
  ~MCacheLock() {
    if (descriptor_ >= 0) {
      ::flock(descriptor_, LOCK_UN);
      ::close(descriptor_);
    }
  }

  MCacheLock(const MCacheLock&)            = delete;
  MCacheLock& operator=(const MCacheLock&) = delete;

  // Compute whether the exclusive cache lock was acquired
  bool Acquired() const noexcept { return descriptor_ >= 0; }

 private:
  int descriptor_ = -1;
};

// Compute a recoverable name for one replaced cache file
inline std::filesystem::path BackupPath(const std::filesystem::path& filename) {
  std::error_code error;
  for (std::size_t version = 0;; ++version) {
    const std::string     suffix = version == 0 ? "._old" : "._old." + std::to_string(version);
    std::filesystem::path backup(filename.string() + suffix);
    if (!std::filesystem::exists(backup, error)) { return backup; }
  }
}

// Atomically publish one JSON table while preserving the previous table
inline void PublishCache(const std::string& filename, const nlohmann::json& cache, const std::string& label) {
  aux::CreateDirectory(std::filesystem::path(filename).parent_path().string());
  std::string       pattern = filename + ".tmp.XXXXXX";
  std::vector<char> path(pattern.cbegin(), pattern.cend());
  path.push_back('\0');
  const int descriptor = ::mkstemp(path.data());
  if (descriptor < 0) { throw std::runtime_error(label + " temporary file creation failed"); }
  ::close(descriptor);
  const std::filesystem::path temporary(path.data());
  std::ofstream               output(temporary, std::ios::trunc);
  if (!output.is_open()) { throw std::runtime_error(label + " temporary file open failed"); }
  output << cache.dump() << '\n';
  output.close();
  if (output.fail()) { throw std::runtime_error(label + " temporary file close failed"); }

  const std::filesystem::path target(filename);
  std::filesystem::path       backup;
  std::error_code             error;
  if (std::filesystem::exists(target, error)) {
    backup = BackupPath(target);
    std::filesystem::rename(target, backup, error);
    if (error) { throw std::runtime_error(label + " archive failed: " + error.message()); }
  }
  std::filesystem::rename(temporary, target, error);
  if (!error) { return; }
  const std::string publication_error = error.message();
  if (!backup.empty()) {
    error.clear();
    std::filesystem::rename(backup, target, error);
  }
  throw std::runtime_error(label + " publication failed: " + publication_error);
}

}  // namespace gra::eikonal

#endif
