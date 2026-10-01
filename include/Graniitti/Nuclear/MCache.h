// Thread-safe persistent storage for exact nuclear numerical tables
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARCACHE_H
#define MNUCLEARCACHE_H

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "Graniitti/Tech/MAux.h"
#include "json.hpp"

namespace gra::nuclear {

// Hold one exclusive interprocess lock while a nuclear table is built
class MCacheLock {
 public:
  // Acquire one exclusive cache-file lock
  explicit MCacheLock(const std::string &filename) {
    std::error_code error;
    aux::CreateDirectory(std::filesystem::path(filename).parent_path().string(), &error);
    if (error) { return; }
    descriptor_ = ::open(filename.c_str(), O_CREAT | O_RDWR, 0600);
    if (descriptor_ >= 0 && ::flock(descriptor_, LOCK_EX) != 0) {
      ::close(descriptor_);
      descriptor_ = -1;
    }
  }

  // Release the cache-file lock
  ~MCacheLock() {
    if (descriptor_ >= 0) {
      ::flock(descriptor_, LOCK_UN);
      ::close(descriptor_);
    }
  }

  MCacheLock(const MCacheLock &)            = delete;
  MCacheLock &operator=(const MCacheLock &) = delete;

  // Compute whether the exclusive cache lock was acquired
  bool Acquired() const { return descriptor_ >= 0; }

 private:
  int descriptor_ = -1;
};

// Compute the persistent filename for one exact nuclear table
inline std::string NuclearCacheFilename(const std::string &name, const std::string &key) {
  return aux::GetBasePath(2) + "/nuclear/" + name + '_' + std::to_string(aux::djb2hash(key)) + ".json";
}

// Read one JSON cache with an exact type, version and physics key
inline bool ReadNuclearCache(const std::string &filename, const std::size_t version, const std::string &type,
                             const std::string &key, nlohmann::json &cache) {
  try {
    std::ifstream input(filename);
    if (!input.is_open()) { return false; }
    input >> cache;
    return cache.is_object() && cache.at("version").get<std::size_t>() == version &&
           cache.at("type").get<std::string>() == type && cache.at("key").get<std::string>() == key;
  } catch (const std::exception &) { return false; }
}

// Construct the common identity fields of one JSON cache
inline nlohmann::json NuclearCache(const std::size_t version, const std::string &type, const std::string &key) {
  return {{"version", version}, {"type", type}, {"key", key}};
}

// Atomically publish one exact JSON table
inline void PublishNuclearCache(const std::string &filename, const nlohmann::json &cache) {
  std::error_code error;
  aux::CreateDirectory(std::filesystem::path(filename).parent_path().string(), &error);
  if (error) { throw std::runtime_error("Nuclear cache directory creation failed: " + error.message()); }
  std::string       temporary = filename + ".tmp.XXXXXX";
  std::vector<char> path(temporary.begin(), temporary.end());
  path.push_back('\0');
  const int descriptor = ::mkstemp(path.data());
  if (descriptor < 0) { throw std::runtime_error("Nuclear cache temporary file creation failed"); }
  ::close(descriptor);
  const std::string temporary_filename(path.data());
  std::ofstream     output(temporary_filename, std::ios::trunc);
  if (!output.is_open()) { throw std::runtime_error("Nuclear cache temporary file open failed"); }
  output << cache.dump() << '\n';
  output.close();
  if (output.fail()) { throw std::runtime_error("Nuclear cache temporary file close failed"); }

  std::filesystem::rename(temporary_filename, filename, error);
  if (error) { throw std::runtime_error("Nuclear cache publication failed: " + error.message()); }
}

}  // namespace gra::nuclear

#endif
