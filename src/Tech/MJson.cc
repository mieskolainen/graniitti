// Generic JSON card loading with local file and JSON Pointer references
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tech/MJson.h"

#include <filesystem>
#include <stdexcept>

#include "Graniitti/Tech/MAux.h"

namespace gra {
namespace {

// Compute the canonical source path used to detect references across files
std::string CardPath(const std::string &path) {
  return std::filesystem::weakly_canonical(std::filesystem::absolute(path)).string();
}

// Decode percent escapes in a local URI reference
std::string DecodeURI(const std::string &text) {
  std::string result;
  for (std::size_t i = 0; i < text.size(); ++i) {
    if (text[i] != '%') { result += text[i]; continue; }
    if (i + 2 >= text.size()) { throw std::invalid_argument("incomplete URI percent escape"); }
    const std::string digits = text.substr(i + 1, 2);
    if (digits.find_first_not_of("0123456789abcdefABCDEF") != std::string::npos) {
      throw std::invalid_argument("invalid URI percent escape");
    }
    result += static_cast<char>(std::stoi(digits, nullptr, 16));
    i += 2;
  }
  (void)nlohmann::json(result).dump();
  return result;
}

// Decode a reference to a whole card or one JSON Pointer within it
std::pair<std::string, std::string> Reference(const nlohmann::json &node, const std::string &path) {
  if (node.size() != 1 || !node.at("$ref").is_string()) {
    throw std::invalid_argument("a reference must contain only a string $ref");
  }
  const std::string value = node.at("$ref");
  const auto hash = value.find('#');
  const std::string file = DecodeURI(value.substr(0, hash));
  if (file.find(':') != std::string::npos || file.find('?') != std::string::npos || file.find('\0') != std::string::npos || file.starts_with("//")) {
    throw std::invalid_argument("$ref requires a local file path and optional JSON Pointer fragment");
  }
  const std::string pointer = hash == std::string::npos ? "" : DecodeURI(value.substr(hash + 1));
  (void)nlohmann::json::json_pointer(pointer);
  return {file.empty() ? path : CardPath((std::filesystem::path(path).parent_path() / file).string()), pointer};
}

// Encode an object key as one JSON Pointer token
std::string PointerToken(const std::string &key) {
  std::string result;
  for (const char c : key) {
    if (c == '~') { result += "~0"; }
    else if (c == '/') { result += "~1"; }
    else { result += c; }
  }
  return result;
}

}  // namespace

// Access the unexpanded source document for reference-aware edits
const nlohmann::json &MJson::Document(const std::string &path) {
  const auto canonical = documents_.contains(path) ? path : CardPath(path);
  const auto found = documents_.find(canonical);
  if (found != documents_.end()) { return found->second; }
  const auto data = aux::GetInputDataRaw(canonical, overrides_);
  return documents_.emplace(canonical, nlohmann::json::parse(data)).first->second;
}

// Seed a source document while invalidating previously expanded values
void MJson::Set(const std::string &path, const nlohmann::json &document) {
  documents_[CardPath(path)] = document;
  values_.clear();
}

// Compute the source file and pointer through references at any ancestor
std::pair<std::string, std::string> MJson::Origin(const std::string &path, const std::string &pointer) {
  std::string file = documents_.contains(path) ? path : CardPath(path);
  std::string suffix = pointer;
  (void)nlohmann::json::json_pointer(suffix);
  std::set<Key> references;
  while (true) {
    const nlohmann::json *node = &Document(file);
    std::size_t position = 0;
    while (true) {
      if (node->is_object() && node->contains("$ref")) {
        const Key location{file, suffix.substr(0, position)};
        if (!references.insert(location).second) {
          throw std::invalid_argument("circular JSON reference at " + file + "#" + location.second);
        }
        const auto [target, prefix] = Reference(*node, file);
        file = target;
        suffix = prefix + suffix.substr(position);
        break;
      }
      if (position == suffix.size()) { return {file, suffix}; }
      const auto end = suffix.find('/', position + 1);
      const auto token = suffix.substr(position, end == std::string::npos ? end : end - position);
      node = &node->at(nlohmann::json::json_pointer(token));
      position = end == std::string::npos ? suffix.size() : end;
    }
  }
}

// Resolve one node including references encountered along its pointer
nlohmann::json MJson::Node(const std::string &path, const std::string &pointer) {
  const Key key{path, pointer};
  if (const auto found = values_.find(key); found != values_.end()) { return found->second; }
  if (!active_.insert(key).second) { throw std::invalid_argument("circular JSON reference at " + path + "#" + pointer); }
  try {
    const auto [target, source] = Origin(path, pointer);
    if (Key{target, source} != key) {
      auto value = Node(target, source);
      values_.emplace(key, value);
      active_.erase(key);
      return value;
    }
    const nlohmann::json *node = &Document(path).at(nlohmann::json::json_pointer(pointer));
    nlohmann::json value = *node;
    if (node->is_object()) {
      for (const auto &[name, item] : node->items()) {
        (void)item;
        value[name] = Node(path, pointer + "/" + PointerToken(name));
      }
    } else if (node->is_array()) {
      using gra::aux::indices;
      for (const auto i : indices(*node)) { value[i] = Node(path, pointer + "/" + std::to_string(i)); }
    }
    values_.emplace(key, value);
    active_.erase(key);
    return value;
  } catch (const std::exception &error) {
    active_.erase(key);
    throw std::invalid_argument(path + "#" + pointer + ": " + error.what());
  }
}

// Load one card and expand its references before model initialization
nlohmann::json MJson::Read(const std::string &path) { return Node(CardPath(path), ""); }

// Resolve references in an already parsed document at its source path
nlohmann::json MJson::Resolve(const nlohmann::json &document, const std::string &path) {
  const std::string canonical = CardPath(path);
  Set(canonical, document);
  return Node(canonical, "");
}

}  // namespace gra
