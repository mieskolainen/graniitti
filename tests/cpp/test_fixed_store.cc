// Thread-safety tests for immutable shared physics stores
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// Own
#include "Graniitti/Tech/MFixedStore.h"

// Libraries
#include <catch.hpp>

namespace {

struct StoredValue {
  int value = 0;
};

// Observe actual keyed lookups before allowing the shared loader to fail
struct StoreRequests {
  std::mutex mutex;
  std::condition_variable condition;
  std::set<std::thread::id> callers;

  // Record each caller once even when the ordered lookup compares twice
  void Mark() {
    std::lock_guard<std::mutex> lock(mutex);
    callers.insert(std::this_thread::get_id());
    condition.notify_all();
  }

  // Wait until every request has entered the store
  bool Wait(std::size_t count) {
    std::unique_lock<std::mutex> lock(mutex);
    return condition.wait_for(lock, std::chrono::seconds(10), [&] { return callers.size() == count; });
  }
};

// Retain ordinary integer ordering while observing calls through the public key API
struct StoreKey {
  int value;
  StoreRequests *requests;

  // Mark lookup entry while the store protects its shared future
  bool operator<(const StoreKey &other) const {
    requests->Mark();
    return value < other.value;
  }
};

} // namespace

TEST_CASE("MFixedStore runs one cold load for concurrent same-key calls",
          "[MFixedStore][threading]") {
  gra::MFixedStore<int, StoredValue> store;
  constexpr std::size_t thread_count = 12;
  std::atomic<std::size_t> ready{0};
  std::atomic<bool> start{false};
  std::atomic<std::size_t> load_count{0};
  std::vector<std::shared_ptr<const StoredValue>> values(thread_count);
  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  // Start every request against the same genuinely cold key
  for (std::size_t i = 0; i < thread_count; ++i) {
    threads.emplace_back([&, i] {
      ready.fetch_add(1, std::memory_order_release);
      while (!start.load(std::memory_order_acquire)) {
        std::this_thread::yield();
      }
      values[i] = store.GetOrLoad(7, [&] {
        load_count.fetch_add(1, std::memory_order_relaxed);
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        return std::make_shared<const StoredValue>(StoredValue{91});
      });
    });
  }
  while (ready.load(std::memory_order_acquire) != thread_count) {
    std::this_thread::yield();
  }
  start.store(true, std::memory_order_release);
  for (auto &thread : threads) {
    thread.join();
  }

  REQUIRE(load_count.load(std::memory_order_relaxed) == 1);
  for (const auto &value : values) {
    REQUIRE(value == values.front());
    REQUIRE(value->value == 91);
  }
}

TEST_CASE("MFixedStore loads independent keys concurrently",
          "[MFixedStore][threading]") {
  gra::MFixedStore<int, StoredValue> store;
  std::mutex mutex;
  std::condition_variable entered_condition;
  std::size_t entered = 0;
  bool concurrent = true;

  // Wait until both independent loaders are active at the same time
  const auto loader = [&](int value) {
    std::unique_lock<std::mutex> lock(mutex);
    ++entered;
    entered_condition.notify_all();
    concurrent = concurrent && entered_condition.wait_for(
                                   lock, std::chrono::seconds(2),
                                   [&] { return entered == 2; });
    return std::make_shared<const StoredValue>(StoredValue{value});
  };

  std::shared_ptr<const StoredValue> first;
  std::shared_ptr<const StoredValue> second;
  std::thread first_thread([&] {
    first = store.GetOrLoad(1, [&] { return loader(11); });
  });
  std::thread second_thread([&] {
    second = store.GetOrLoad(2, [&] { return loader(22); });
  });
  first_thread.join();
  second_thread.join();

  REQUIRE(concurrent);
  REQUIRE(first->value == 11);
  REQUIRE(second->value == 22);
}

TEST_CASE("MFixedStore publishes one failure and permits a later retry",
          "[MFixedStore][threading]") {
  gra::MFixedStore<StoreKey, StoredValue> store;
  StoreRequests requests;
  const StoreKey key{3, &requests};
  constexpr std::size_t thread_count = 8;
  std::atomic<bool> all_waiting{false};
  std::atomic<std::size_t> load_count{0};
  std::vector<std::string> failures(thread_count);
  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  // Make all callers observe the same in-flight exception
  for (std::size_t i = 0; i < thread_count; ++i) {
    threads.emplace_back([&, i] {
      try {
        static_cast<void>(store.GetOrLoad(key, [&]() -> std::shared_ptr<const StoredValue> {
          load_count.fetch_add(1, std::memory_order_relaxed);
          requests.Mark();
          all_waiting.store(requests.Wait(thread_count), std::memory_order_relaxed);
          throw std::invalid_argument("invalid immutable input");
        }));
      } catch (const std::invalid_argument &error) {
        failures[i] = error.what();
      }
    });
  }
  for (auto &thread : threads) {
    thread.join();
  }

  REQUIRE(all_waiting.load(std::memory_order_relaxed));
  REQUIRE(load_count.load(std::memory_order_relaxed) == 1);
  for (const std::string &failure : failures) {
    REQUIRE(failure == "invalid immutable input");
  }

  const auto recovered = store.GetOrLoad(key, [&] {
    load_count.fetch_add(1, std::memory_order_relaxed);
    return std::make_shared<const StoredValue>(StoredValue{33});
  });
  REQUIRE(recovered->value == 33);
  REQUIRE(load_count.load(std::memory_order_relaxed) == 2);
}

TEST_CASE("MFixedStore rejects recursive same-key loading",
          "[MFixedStore][threading]") {
  gra::MFixedStore<int, StoredValue> store;

  REQUIRE_THROWS(
      store.GetOrLoad(5, [&] {
        return store.GetOrLoad(5, [] {
          return std::make_shared<const StoredValue>(StoredValue{5});
        });
      }));
}
