// Timer class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTIMER_H
#define MTIMER_H

// C++
#include <chrono>

namespace gra {
class MTimer {
  public:
	// Declare explicit, means here no conversion to bool allowed
	explicit MTimer(bool reset = true) {
		if(reset) {
			Reset();
		}
	}
	// << operator
	template <typename T, typename Traits>
	friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out,
													 const MTimer& timer) {
		return out << timer.Elapsed().count();
	}
	// Time in sec
	double ElapsedSec() const {
		return Elapsed().count() / 1000.0;
	}
	// Time in msec
	std::chrono::milliseconds Elapsed() const {
		return std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::high_resolution_clock::now() - start);
	}
	// Reset timer
	void Reset() {
		start = std::chrono::high_resolution_clock::now();
	}

  private:
	std::chrono::high_resolution_clock::time_point start;
};

} // gra namespace ends

#endif