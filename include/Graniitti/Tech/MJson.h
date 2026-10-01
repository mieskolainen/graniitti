// Generic JSON card loading with local file and JSON Pointer references
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MJSON_H
#define MJSON_H

#include <map>
#include <set>
#include <string>
#include <utility>

#include "json.hpp"

namespace gra {

// Resolve one collection of JSON cards with storage owned by the reader
class MJson {
 public:
  // Select whether registered command line card edits are visible
  explicit MJson(bool overrides = true) : overrides_(overrides) {}

  // Compute the source file and pointer through references at any ancestor
  std::pair<std::string, std::string> Origin(const std::string &path, const std::string &pointer);

  // Access the unexpanded source document for reference-aware edits
  const nlohmann::json &Document(const std::string &path);

  // Seed a source document while invalidating previously expanded values
  void Set(const std::string &path, const nlohmann::json &document);

  // Load one card and expand its references before model initialization
  nlohmann::json Read(const std::string &path);

  // Resolve references in an already parsed document at its source path
  nlohmann::json Resolve(const nlohmann::json &document, const std::string &path);

 private:
  bool overrides_;
  using Key = std::pair<std::string, std::string>;
  std::map<std::string, nlohmann::json> documents_;
  std::map<Key, nlohmann::json> values_;
  std::set<Key> active_;

  // Resolve one node including references encountered along its pointer
  nlohmann::json Node(const std::string &path, const std::string &pointer);
};

}  // namespace gra

#endif
