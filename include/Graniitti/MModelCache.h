// Run owned caches for expensive immutable model tables
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMODELCACHE_H
#define MMODELCACHE_H

// C++
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <utility>

// Own
#include "Graniitti/MModelTune.h"
#include "Graniitti/PDF/MHardPomeronPDF.h"
#include "Graniitti/PDF/MLHAPDF.h"
#include "Graniitti/PDF/MSudakov.h"
#include "Graniitti/Tech/MFixedStore.h"

namespace gra {

namespace model {

// Identify one immutable parameter block by its C++ type and physics key
struct CacheKey {
  std::type_index type;
  std::string physics;

  // Order heterogeneous parameter blocks for the fixed store
  bool operator<(const CacheKey &other) const {
    if (type != other.type) {
      return type < other.type;
    }
    return physics < other.physics;
  }
};

// Type erased base for immutable parameter blocks
struct CacheValue {
  // Destroy one type erased immutable parameter block
  virtual ~CacheValue() = default;
};

// Hold one typed immutable parameter block in the unified store
template <typename T> struct CacheValueT final : CacheValue {
  // Construct one typed immutable parameter block
  explicit CacheValueT(std::shared_ptr<const T> value_in)
      : value(std::move(value_in)) {}

  std::shared_ptr<const T> value;
};

} // namespace model

// Run owned stores shared by every worker process
class MModelCache {
public:
  // Construct empty stores bound to one immutable model tune
  explicit MModelCache(MModelTunePtr tune = nullptr)
      : sudakov(pdf), hard_pomeron(pdf), tune_(std::move(tune)) {}

  // Compute the immutable model tune owning every parameter block
  const MModelTune &Tune() const { return *TunePtr(); }

  // Compute the shared immutable model tune handle
  MModelTunePtr TunePtr() const {
    std::lock_guard<std::mutex> lock(tune_mutex);
    if (tune_ == nullptr) {
      throw std::logic_error("MModelCache::TunePtr: no model tune is bound");
    }
    return tune_;
  }

  // Bind an unconfigured cache or validate one immutable tune snapshot
  void Bind(const MModelTunePtr &tune) {
    if (tune == nullptr) {
      throw std::invalid_argument("MModelCache::Bind: null model tune");
    }
    std::lock_guard<std::mutex> lock(tune_mutex);
    if (tune_ == nullptr) {
      tune_ = tune;
      return;
    }
    if (tune_ == tune) {
      return;
    }
    if (tune_->GeneralFile() != tune->GeneralFile() ||
        tune_->NumericsFile() != tune->NumericsFile() ||
        tune_->General().dump() != tune->General().dump() ||
        tune_->Numerics().dump() != tune->Numerics().dump()) {
      throw std::invalid_argument(
          "MModelCache::Bind: model tune snapshot changed");
    }
  }

  // Compute one typed immutable parameter block and load a cold key once
  template <typename T, typename Loader>
  std::shared_ptr<const T> Get(const std::string &physics, Loader &&loader) {
    using Value = model::CacheValueT<T>;
    const model::CacheKey key{std::type_index(typeid(T)), physics};
    const auto erased = param.GetOrLoad(key, [&loader] {
      std::shared_ptr<const T> value = std::forward<Loader>(loader)();
      if (value == nullptr) {
        throw std::logic_error(
            "MModelCache::Get: parameter loader returned null");
      }
      return std::shared_ptr<const model::CacheValue>(
          std::make_shared<const Value>(std::move(value)));
    });
    const auto typed = std::dynamic_pointer_cast<const Value>(erased);
    if (typed == nullptr) {
      throw std::logic_error("MModelCache::Get: parameter type mismatch");
    }
    return typed->value;
  }

  MLHAPDFStore pdf;
  MSudakovStore sudakov;
  MHardPomeronPDFStore hard_pomeron;

private:
  MModelTunePtr tune_;
  mutable std::mutex tune_mutex;
  MFixedStore<model::CacheKey, model::CacheValue> param;
};

using MModelCachePtr = std::shared_ptr<MModelCache>;

// Compute a run owned cache bound to one immutable tune snapshot
inline MModelCache &RequireModelCache(MModelCachePtr &cache,
                                      const MModelTunePtr &tune,
                                      const std::string &context) {
  if (tune == nullptr) {
    throw std::invalid_argument(context + ": missing model tune");
  }
  if (cache == nullptr) {
    cache = std::make_shared<MModelCache>(tune);
  } else {
    cache->Bind(tune);
  }
  return *cache;
}

} // namespace gra

#endif
