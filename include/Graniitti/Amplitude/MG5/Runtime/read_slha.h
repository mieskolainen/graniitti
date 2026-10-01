// Read SLHA parameter cards for generated MG5 amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef READ_SLHA_H
#define READ_SLHA_H

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

class SLHABlock {
 public:
  // Construct one named SLHA block
  explicit SLHABlock(std::string name = "");

  // Store one finite value under a fixed-rank integer index
  void set_entry(const std::vector<int> &indices, double value);

  // Compute one block entry or the requested default
  double get_entry(const std::vector<int> &indices,
                   double default_value = 0.0) const;

  // Compute one block entry without emitting a missing-entry warning
  std::optional<double> find_entry(const std::vector<int> &indices) const;

  // Replace the normalized block name
  void set_name(std::string name);

  // Compute the normalized block name
  const std::string &get_name() const;

  // Compute the index rank after the first stored entry
  std::size_t get_indices() const;

 private:
  std::string _name;
  std::map<std::vector<int>, double> _entries;
  std::size_t _indices = 0;
  bool _has_indices = false;
};

class SLHAReader {
 public:
  // Construct an empty reader or load one parameter card
  explicit SLHAReader(const std::string &file_name = "");

  // Parse one complete parameter card atomically
  void read_slha_file(const std::string &file_name);

  // Compute one indexed block entry or the requested default
  double get_block_entry(const std::string &block_name,
                         const std::vector<int> &indices,
                         double default_value = 0.0) const;

  // Compute one single-index block entry or the requested default
  double get_block_entry(const std::string &block_name, int index,
                         double default_value = 0.0) const;

  // Compute one indexed block entry without emitting a missing-entry warning
  std::optional<double>
  find_block_entry(const std::string &block_name,
                   const std::vector<int> &indices) const;

  // Compute one single-index block entry without a missing-entry warning
  std::optional<double> find_block_entry(const std::string &block_name,
                                         int index) const;

  // Store one indexed block entry
  void set_block_entry(const std::string &block_name,
                       const std::vector<int> &indices, double value);

  // Store one single-index block entry
  void set_block_entry(const std::string &block_name, int index,
                       double value);

 private:
  std::map<std::string, SLHABlock> _blocks;
};

#endif
