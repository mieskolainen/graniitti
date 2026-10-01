// JSON steering card command line overrides
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MJSONOVERRIDE_H
#define MJSONOVERRIDE_H

// C++
#include <cctype>
#include <filesystem>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Libraries
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJson.h"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;

namespace gra {
namespace json_override {

using json = nlohmann::json;

struct PathToken {
  bool is_index = false;
  std::string key;
  std::size_t index = 0;
};

struct OverrideSpec {
  std::string selector;
  std::vector<PathToken> path;
  json value;
  std::string raw;
  std::size_t applied = 0;
  std::set<std::string> files;
};

struct OverrideRegistryState {
  std::mutex mutex;
  std::vector<OverrideSpec> overrides;
  std::map<std::string, std::string> resolved_cards;
  bool frozen = false;
};

// Compute synchronized process-wide JSON card override state
inline OverrideRegistryState &OverrideRegistry() {
  static OverrideRegistryState registry;
  return registry;
}

// Compute true for an escaped character in a CLI override string
inline bool IsEscaped(const std::string &s, std::size_t pos) {
  std::size_t n = 0;
  while (pos > n && s[pos - n - 1] == '\\') {
    ++n;
  }
  return n % 2 == 1;
}

// Find a separator outside quoted and bracketed path sections
inline std::size_t FindTopLevel(const std::string &s, char needle) {
  bool quoted = false;
  char quote = '\0';
  std::size_t depth = 0;

  for (const auto &i : indices(s)) {
    const char c = s[i];
    if (quoted) {
      if (c == quote && !IsEscaped(s, i)) {
        quoted = false;
      }
      continue;
    }
    if (c == '"' || c == '\'') {
      quoted = true;
      quote = c;
      continue;
    }
    if (c == '[') {
      ++depth;
    } else if (c == ']') {
      if (depth > 0) {
        --depth;
      }
    } else if (c == needle && depth == 0) {
      return i;
    }
  }
  return std::string::npos;
}

// Trim leading and trailing ASCII whitespace
inline std::string Trim(const std::string &s) {
  std::size_t first = 0;
  while (first < s.size() &&
         std::isspace(static_cast<unsigned char>(s[first]))) {
    ++first;
  }
  std::size_t last = s.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(s[last - 1]))) {
    --last;
  }
  return s.substr(first, last - first);
}

// Parse a JSON value from an override right hand side
inline json ParseJsonValue(const std::string &value,
                           const std::string &context) {
  try {
    return json::parse(value);
  } catch (const json::exception &e) {
    throw std::invalid_argument("JSON override '" + context +
                                "': invalid JSON value: " + e.what());
  }
}

// Parse a quoted key inside bracket path notation
inline std::pair<std::string, std::size_t>
ParseQuotedKey(const std::string &path, std::size_t pos) {
  const char quote = path[pos];
  std::string out;
  ++pos;
  while (pos < path.size()) {
    const char c = path[pos];
    if (c == quote && !IsEscaped(path, pos)) {
      if (quote == '"') {
        return {ParseJsonValue("\"" + out + "\"", path).get<std::string>(), pos + 1};
      }
      return {out, pos + 1};
    }
    if (c == '\\' && quote == '\'' && pos + 1 < path.size() &&
        (path[pos + 1] == '\\' || path[pos + 1] == '\'')) {
      out.push_back(path[++pos]);
    } else {
      out.push_back(c);
    }
    ++pos;
  }
  throw std::invalid_argument("JSON override path '" + path +
                              "': unterminated quoted key");
}

// Parse a non-negative array index inside bracket path notation
inline std::size_t ParseArrayIndex(const std::string &text,
                                   const std::string &path) {
  if (text.empty()) {
    throw std::invalid_argument("JSON override path '" + path +
                                "': empty array index");
  }
  for (const char c : text) {
    if (!std::isdigit(static_cast<unsigned char>(c))) {
      throw std::invalid_argument("JSON override path '" + path +
                                  "': invalid array index '" + text + "'");
    }
  }
  return static_cast<std::size_t>(std::stoull(text));
}

// Parse comma-separated array indices inside one bracket pair
inline std::vector<std::size_t> ParseArrayIndices(const std::string &text,
                                                  const std::string &path) {
  std::vector<std::size_t> indices;
  std::size_t start = 0;
  while (start <= text.size()) {
    const std::size_t comma = text.find(',', start);
    const std::size_t end = comma == std::string::npos ? text.size() : comma;
    indices.push_back(
        ParseArrayIndex(Trim(text.substr(start, end - start)), path));
    if (comma == std::string::npos) {
      break;
    }
    start = comma + 1;
  }
  return indices;
}

// Parse dot and bracket path notation into object key and array index tokens
inline std::vector<PathToken> ParsePath(const std::string &path) {
  if (path.empty()) {
    throw std::invalid_argument("JSON override path is empty");
  }

  std::vector<PathToken> tokens;
  std::size_t pos = 0;
  while (pos < path.size()) {
    if (path[pos] == '.') {
      if (pos == 0 || path[pos - 1] == '.') {
        throw std::invalid_argument("JSON override path '" + path +
                                    "': empty key");
      }
      ++pos;
      if (pos == path.size()) {
        throw std::invalid_argument("JSON override path '" + path +
                                    "': trailing dot");
      }
      continue;
    }

    if (path[pos] == '[') {
      ++pos;
      if (pos >= path.size()) {
        throw std::invalid_argument("JSON override path '" + path +
                                    "': unterminated bracket");
      }
      if (path[pos] == '"' || path[pos] == '\'') {
        auto parsed = ParseQuotedKey(path, pos);
        pos = parsed.second;
        if (pos >= path.size() || path[pos] != ']') {
          throw std::invalid_argument("JSON override path '" + path +
                                      "': expected closing bracket");
        }
        tokens.push_back({false, parsed.first, 0});
        ++pos;
      } else {
        const std::size_t start = pos;
        while (pos < path.size() && path[pos] != ']') {
          ++pos;
        }
        if (pos >= path.size()) {
          throw std::invalid_argument("JSON override path '" + path +
                                      "': unterminated bracket");
        }
        const auto indices =
            ParseArrayIndices(path.substr(start, pos - start), path);
        for (const std::size_t index : indices) {
          tokens.push_back({true, "", index});
        }
        ++pos;
      }
      if (pos < path.size() && path[pos] != '.' && path[pos] != '[') {
        throw std::invalid_argument("JSON override path '" + path + "': expected dot or bracket after index or key");
      }
      continue;
    }

    const std::size_t start = pos;
    while (pos < path.size() && path[pos] != '.' && path[pos] != '[') {
      ++pos;
    }
    std::string key = Trim(path.substr(start, pos - start));
    if (key.empty()) {
      throw std::invalid_argument("JSON override path '" + path +
                                  "': empty key");
    }
    tokens.push_back({false, key, 0});
  }

  if (tokens.empty()) {
    throw std::invalid_argument("JSON override path '" + path + "': no tokens");
  }
  return tokens;
}

// Parse one command line override specification
inline OverrideSpec ParseSpec(const std::string &spec) {
  const std::size_t eq = FindTopLevel(spec, '=');
  if (eq == std::string::npos) {
    throw std::invalid_argument("JSON override '" + spec +
                                "': expected [card:]path=json_value");
  }

  std::string lhs = Trim(spec.substr(0, eq));
  std::string rhs = Trim(spec.substr(eq + 1));
  if (lhs.empty() || rhs.empty()) {
    throw std::invalid_argument("JSON override '" + spec +
                                "': empty path or value");
  }

  OverrideSpec out;
  out.raw = spec;
  const std::size_t colon = FindTopLevel(lhs, ':');
  if (colon != std::string::npos) {
    out.selector = Trim(lhs.substr(0, colon));
    lhs = Trim(lhs.substr(colon + 1));
  }
  out.path = ParsePath(lhs);
  out.value = ParseJsonValue(rhs, spec);
  return out;
}

// Parse all command line override specifications
inline std::vector<OverrideSpec>
ParseSpecs(const std::vector<std::string> &specs) {
  std::vector<OverrideSpec> out;
  out.reserve(specs.size());
  for (const auto &spec : specs) {
    out.push_back(ParseSpec(spec));
  }
  return out;
}

// Compute true if the override selector targets the main input card
inline bool IsInputOverride(const OverrideSpec &spec) {
  return spec.selector.empty() || spec.selector == "input";
}

// Compute true if the override selector targets a model card
inline bool IsCardOverride(const OverrideSpec &spec) {
  return !IsInputOverride(spec);
}

// Compute true when a selector contains an explicit path separator
inline bool HasPathSeparator(const std::string &selector) {
  return selector.find('/') != std::string::npos ||
         selector.find('\\') != std::string::npos;
}

// Compute a filesystem path in generic string form
inline std::string GenericPath(const std::filesystem::path &path) {
  return path.lexically_normal().generic_string();
}

// Resolve one model-card selector to all matching regular files
inline std::vector<std::filesystem::path>
ResolveCardSelectorTargets(const std::string &selector,
                           const std::filesystem::path &model_tune_dir) {
  std::set<std::string> targets;
  const std::filesystem::path selected_path(selector);
  if (selected_path.is_absolute() || selector.rfind(".", 0) == 0) {
    if (std::filesystem::is_regular_file(selected_path)) {
      targets.insert(GenericPath(std::filesystem::absolute(selected_path)));
    }
  } else if (HasPathSeparator(selector)) {
    const std::filesystem::path target = model_tune_dir / selected_path;
    if (std::filesystem::is_regular_file(target)) {
      targets.insert(GenericPath(std::filesystem::absolute(target)));
    }
  } else {
    for (const auto &entry :
         std::filesystem::recursive_directory_iterator(model_tune_dir)) {
      if (entry.is_regular_file() && entry.path().filename() == selected_path) {
        targets.insert(GenericPath(std::filesystem::absolute(entry.path())));
      }
    }
  }

  std::vector<std::filesystem::path> resolved;
  resolved.reserve(targets.size());
  for (const auto &target : targets) {
    resolved.emplace_back(target);
  }
  return resolved;
}

// Compute true if a model card selector resolves to an existing file
inline bool
SelectorResolvesToExistingFile(const std::string &selector,
                               const std::filesystem::path &model_tune_dir) {
  return !ResolveCardSelectorTargets(selector, model_tune_dir).empty();
}

// Fail early when a model card selector cannot resolve inside the chosen tune
inline void
ValidateCardOverrideSelectors(const std::vector<OverrideSpec> &specs,
                              const std::filesystem::path &model_tune_dir) {
  if (!std::filesystem::exists(model_tune_dir)) {
    throw std::invalid_argument("JSON override: model tune directory '" +
                                GenericPath(model_tune_dir) +
                                "' does not exist");
  }
  for (const auto &spec : specs) {
    if (IsCardOverride(spec) &&
        !SelectorResolvesToExistingFile(spec.selector, model_tune_dir)) {
      throw std::invalid_argument(
          "JSON override '" + spec.raw + "': selector '" + spec.selector +
          "' does not resolve to an existing JSON card under '" +
          GenericPath(model_tune_dir) + "'");
    }
  }
}

// Resolve and validate the unique model cards targeted by registered overrides
inline std::vector<std::string>
ResolveCardOverrideTargets(const std::vector<OverrideSpec> &specs,
                           const std::filesystem::path &model_tune_dir) {
  ValidateCardOverrideSelectors(specs, model_tune_dir);

  std::set<std::string> targets;
  for (const auto &spec : specs) {
    if (!IsCardOverride(spec)) {
      continue;
    }
    for (const auto &target :
         ResolveCardSelectorTargets(spec.selector, model_tune_dir)) {
      targets.insert(GenericPath(target));
    }
  }
  return {targets.begin(), targets.end()};
}

// Compute true if a selector matches the card file being read
inline bool SelectorMatchesFile(const std::string &selector,
                                const std::string &inputfile) {
  if (selector.empty() || selector == "input") {
    return false;
  }

  const std::filesystem::path file_path(inputfile);
  const std::string file = GenericPath(std::filesystem::absolute(file_path));
  const std::string file_name = file_path.filename().generic_string();
  const std::string sel = GenericPath(std::filesystem::path(selector));

  if (sel == file_name) {
    return true;
  }
  if (sel == file) {
    return true;
  }
  if (sel.size() < file.size() &&
      file.compare(file.size() - sel.size(), sel.size(), sel) == 0) {
    return file[file.size() - sel.size() - 1] == '/';
  }
  return false;
}

// Compute a readable path token for diagnostics
inline std::string TokenToString(const PathToken &token) {
  return token.is_index ? "[" + std::to_string(token.index) + "]" : token.key;
}

// Compute the left hand side target name from the original override
inline std::string TargetLabel(const OverrideSpec &spec) {
  const std::size_t eq = FindTopLevel(spec.raw, '=');
  return eq == std::string::npos ? spec.raw : Trim(spec.raw.substr(0, eq));
}

// Print a terminal notice for a JSON card override
inline void PrintOverrideNotice(const OverrideSpec &spec,
                                const json &old_value) {
  std::cout << rang::fg::yellow << "JSON override: " << TargetLabel(spec)
            << " old=" << old_value.dump() << " new=" << spec.value.dump()
            << rang::fg::reset << std::endl;
}

// Compute true when the final path token requests an object key rename
inline bool IsKeyRenameOverride(const OverrideSpec &spec) {
  return !spec.path.empty() && !spec.path.back().is_index &&
         spec.path.back().key == "@key";
}

// Print a terminal notice for a JSON object key rename
inline void PrintKeyRenameNotice(const OverrideSpec &spec,
                                 const std::string &old_key,
                                 const std::string &new_key) {
  std::cout << rang::fg::yellow << "JSON override: " << TargetLabel(spec)
            << " old=" << json(old_key).dump()
            << " new=" << json(new_key).dump() << rang::fg::reset << std::endl;
}

// Compute a mutable JSON child selected by one path token
inline json &SelectChild(json &node, const PathToken &token,
                         const std::string &context,
                         bool create_missing_key = false) {
  if (token.is_index) {
    if (!node.is_array()) {
      throw std::invalid_argument("JSON override '" + context + "': token " +
                                  TokenToString(token) + " requires an array");
    }
    if (token.index >= node.size()) {
      throw std::invalid_argument(
          "JSON override '" + context + "': array index " +
          std::to_string(token.index) + " is out of range");
    }
    return node.at(token.index);
  }

  if (!node.is_object()) {
    throw std::invalid_argument("JSON override '" + context + "': key '" +
                                token.key + "' requires an object");
  }
  if (!node.contains(token.key) && create_missing_key) {
    node[token.key] = nullptr;
  }
  if (!node.contains(token.key)) {
    throw std::invalid_argument("JSON override '" + context +
                                "': missing key '" + token.key + "'");
  }
  return node.at(token.key);
}

// Rename one existing JSON object key selected with a terminal @key token
inline void ApplyKeyRenameOverride(json &root, const OverrideSpec &spec) {
  if (spec.path.size() < 2 || spec.path[spec.path.size() - 2].is_index) {
    throw std::invalid_argument("JSON override '" + spec.raw +
                                "': @key requires an object key target");
  }
  if (!spec.value.is_string() || spec.value.get<std::string>().empty()) {
    throw std::invalid_argument("JSON override '" + spec.raw +
                                "': @key value must be a non-empty string");
  }

  json *node = &root;
  for (std::size_t i = 0; i + 2 < spec.path.size(); ++i) {
    node = &SelectChild(*node, spec.path[i], spec.raw);
  }
  if (!node->is_object()) {
    throw std::invalid_argument("JSON override '" + spec.raw +
                                "': @key target parent must be an object");
  }

  const std::string old_key = spec.path[spec.path.size() - 2].key;
  const std::string new_key = spec.value.get<std::string>();
  if (!node->contains(old_key)) {
    throw std::invalid_argument("JSON override '" + spec.raw +
                                "': missing key '" + old_key + "'");
  }
  if (old_key != new_key && node->contains(new_key)) {
    throw std::invalid_argument("JSON override '" + spec.raw +
                                "': target key '" + new_key +
                                "' already exists");
  }
  if (old_key != new_key) {
    json value = std::move(node->at(old_key));
    node->erase(old_key);
    (*node)[new_key] = std::move(value);
  }
  PrintKeyRenameNotice(spec, old_key, new_key);
}

// Apply one override to a parsed JSON object
inline void ApplyOverride(json &root, const OverrideSpec &spec) {
  if (spec.path.empty()) { throw std::invalid_argument("JSON override: empty path"); }
  if (IsKeyRenameOverride(spec)) {
    ApplyKeyRenameOverride(root, spec);
    return;
  }

  json *node = &root;
  for (std::size_t i = 0; i + 1 < spec.path.size(); ++i) {
    node = &SelectChild(*node, spec.path[i], spec.raw);
  }

  json &target = SelectChild(*node, spec.path.back(), spec.raw, true);
  const json old_value = target;
  target = spec.value;
  PrintOverrideNotice(spec, old_value);
}

// Apply input card overrides directly to the already parsed input card
inline void ApplyInputOverrides(json &root,
                                const std::vector<OverrideSpec> &specs) {
  for (const auto &spec : specs) {
    if (IsInputOverride(spec)) {
      ApplyOverride(root, spec);
    }
  }
}

// Replace registered model card overrides with a fresh set
inline void RegisterCardOverrides(const std::vector<OverrideSpec> &specs) {
  auto &registry = OverrideRegistry();
  std::lock_guard<std::mutex> lock(registry.mutex);
  if (registry.frozen) {
    throw std::logic_error(
        "JSON card overrides are frozen after initialization");
  }
  registry.overrides.clear();
  registry.resolved_cards.clear();
  for (const auto &spec : specs) {
    if (IsCardOverride(spec)) {
      OverrideSpec fresh = spec;
      fresh.applied = 0;
      fresh.files.clear();
      registry.overrides.push_back(std::move(fresh));
    }
  }
}

// Clear registered model card overrides
inline void ClearCardOverrides() {
  auto &registry = OverrideRegistry();
  std::lock_guard<std::mutex> lock(registry.mutex);
  registry.overrides.clear();
  registry.resolved_cards.clear();
  registry.frozen = false;
}

// Compute true when model card overrides are currently registered
inline bool HasCardOverrides() {
  auto &registry = OverrideRegistry();
  std::lock_guard<std::mutex> lock(registry.mutex);
  return !registry.overrides.empty();
}

// Encode command line path tokens with the standard JSON Pointer implementation
inline std::string OverridePointer(const std::vector<PathToken> &tokens, std::size_t count) {
  json::json_pointer pointer;
  for (std::size_t i = 0; i < count; ++i) {
    pointer /= tokens[i].is_index ? std::to_string(tokens[i].index) : tokens[i].key;
  }
  return pointer.to_string();
}

// Edit the referenced source value so every alias observes the command line change
inline std::string ApplyReferencedOverride(MJson &reader, const std::string &inputfile, const OverrideSpec &spec) {
  if (spec.path.empty()) { throw std::invalid_argument("JSON override: empty path"); }
  const bool rename = IsKeyRenameOverride(spec);
  if (rename && spec.path.size() < 2) { throw std::invalid_argument("JSON override: @key requires an object key target"); }
  const std::size_t count = spec.path.size() - (rename ? 2 : 1);
  const auto [parent_file, parent_path] = reader.Origin(inputfile, OverridePointer(spec.path, count));
  auto document = reader.Document(parent_file);
  auto &parent = document.at(json::json_pointer(parent_path));
  if (rename) {
    OverrideSpec local = spec;
    local.path = {spec.path[spec.path.size() - 2], spec.path.back()};
    ApplyKeyRenameOverride(parent, local);
    reader.Set(parent_file, document);
    return parent_file;
  }
  auto &child = SelectChild(parent, spec.path.back(), spec.raw, true);
  if (child.is_object() && child.contains("$ref")) {
    const auto [file, path] = reader.Origin(inputfile, OverridePointer(spec.path, spec.path.size()));
    auto source = reader.Document(file);
    auto &target = source.at(json::json_pointer(path));
    const auto old_value = target;
    target = spec.value;
    reader.Set(file, source);
    PrintOverrideNotice(spec, old_value);
    return file;
  }
  const auto old_value = child;
  child = spec.value;
  reader.Set(parent_file, document);
  PrintOverrideNotice(spec, old_value);
  return parent_file;
}

// Apply matching overrides once and store the resolved immutable card snapshot
inline void ApplyRegisteredCardOverrides(const std::string &inputfile,
                                         json &root) {
  auto &registry = OverrideRegistry();
  std::lock_guard<std::mutex> lock(registry.mutex);

  const std::string key = GenericPath(std::filesystem::absolute(inputfile));
  MJson reader(false);
  for (const auto &[file, source] : registry.resolved_cards) { reader.Set(file, json::parse(source)); }
  if (!registry.resolved_cards.contains(key)) { reader.Set(inputfile, root); }
  std::set<std::string> changed;
  for (auto &spec : registry.overrides) {
    if (SelectorMatchesFile(spec.selector, inputfile) && !spec.files.contains(key)) {
      if (registry.frozen) {
        throw std::logic_error("JSON card overrides are frozen before resolving '" + key + "'");
      }
      changed.insert(ApplyReferencedOverride(reader, inputfile, spec));
      ++spec.applied;
      spec.files.insert(key);
    }
  }
  for (const auto &file : changed) { registry.resolved_cards[file] = reader.Document(file).dump(); }
  root = reader.Document(inputfile);
}

// Fail if any registered model card override never matched a card read
inline void RequireAllCardOverridesApplied() {
  auto &registry = OverrideRegistry();
  std::lock_guard<std::mutex> lock(registry.mutex);
  for (const auto &spec : registry.overrides) {
    if (spec.applied == 0) {
      throw std::invalid_argument(
          "JSON override '" + spec.raw + "': selector '" + spec.selector +
          "' did not match any JSON card read during initialization");
    }
  }
  registry.frozen = true;
}

} // namespace json_override
} // namespace gra

#endif
