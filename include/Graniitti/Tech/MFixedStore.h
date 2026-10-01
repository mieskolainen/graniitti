// Thread safe cache for immutable shared physics / model parameters
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFIXEDSTORE_H
#define MFIXEDSTORE_H

// C++
#include <chrono>
#include <exception>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>

namespace gra {

// Coordinate one immutable load per key without serializing independent keys
template <typename Key, typename Value>
class MFixedStore {
public:
  using Ptr = std::shared_ptr<const Value>;

  // Compute the cached value or run the loader once for one cold key
  template <typename Loader>
  Ptr GetOrLoad(const Key &key, Loader &&loader) {
    std::shared_ptr<std::promise<Ptr>> load_promise;
    std::shared_future<Ptr> value_future;

    {
      std::lock_guard<std::mutex> lock(mutex);
      const auto found = cache.find(key);
      if (found != cache.end()) {
        value_future = found->second.value;
        if (found->second.loader_thread == std::this_thread::get_id() &&
            value_future.wait_for(std::chrono::seconds(0)) !=
                std::future_status::ready) {
          throw std::logic_error(
              "MFixedStore::GetOrLoad: recursive load for the same key");
        }
      } else {
        load_promise = std::make_shared<std::promise<Ptr>>();
        value_future = load_promise->get_future().share();
        cache.emplace(key, Entry{value_future, std::this_thread::get_id()});
      }
    }

    if (!load_promise) {
      return value_future.get();
    }

    try {
      Ptr value = std::forward<Loader>(loader)();
      if (!value) {
        throw std::logic_error(
            "MFixedStore::GetOrLoad: loader returned a null value");
      }
      load_promise->set_value(value);
      return value;
    } catch (...) {
      load_promise->set_exception(std::current_exception());
      std::lock_guard<std::mutex> lock(mutex);
      cache.erase(key);
      throw;
    }
  }

private:
  struct Entry {
    std::shared_future<Ptr> value;
    std::thread::id loader_thread;
  };

  std::mutex mutex;
  std::map<Key, Entry> cache;
};

} // namespace gra

#endif
