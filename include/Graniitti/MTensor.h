// Minimal tensor class
//
// Example: (rank-6 tensor with dim-4 per dimension)
//
// MTensor<double> tensor = MTensor({4,4,4,4,4,4}, 0.0);
// tensor({0,3,2,0,1,2}) = 1.0;
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSOR_H
#define MTENSOR_H

#include <initializer_list>
#include <iomanip>
#include <iostream>

namespace gra {
template <typename T>
class MTensor {
  public:
	MTensor() {
		data = nullptr;
	}
	MTensor(const std::vector<std::size_t>& newdim) {
		dim = newdim;
		data = new T[Prod(dim)];
		std::fill(data, data + Prod(dim), T()); // No initialization
	}
	MTensor(const std::vector<std::size_t>& newdim, T value) {
		dim = newdim;
		data = new T[Prod(dim)];
		std::fill(data, data + Prod(dim), T(value)); // Initialization
	}
	~MTensor() {
		// delete dynamically allocated memory
		delete[] data;
	}
	// Copy constructor
	MTensor(const MTensor& a) {
		dim = a.dim;
		data = new T[Prod(dim)];

		// Copy all elements
		Copy(a);
	}
	// Assignment operator
	MTensor& operator=(const MTensor& rhs) {
		if(data != rhs.data && rhs.data != nullptr) {
			ReSize(rhs.dim);
			Copy(rhs);
		}
		return *this;
	}
	T& operator()(const std::vector<size_t>& ind) {
		const std::size_t i = Index(ind);
		return data[i];
	}
	T& operator()(const std::vector<size_t>& ind) const {
		const std::size_t i = Index(ind);
		return data[i];
	}
	// Size operators
	std::size_t size(std::size_t ind) const {
		return dim[ind];
	}

  private:
	// Multidimensional indexing algorithm
	// Row-major order
	std::size_t Index(const std::vector<size_t>& ind) const {
		if(ind.size() != dim.size()) {
			throw std::invalid_argument("MTensor:: Error: Index vector with rank = " +
										std::to_string(ind.size()) + " c.f. Tensor has rank " +
										std::to_string(dim.size()));
		}
		std::size_t sum = 0;
		for(std::size_t i = 0; i < ind.size(); ++i) {
			if(ind[i] >= dim[i]) {
				throw std::invalid_argument("MTensor:: Error: Input index " + std::to_string(i) +
											" over bounds: " + std::to_string(ind[i]) + " >= " +
											std::to_string(dim[i]));
			}
			std::size_t product = 1;
			for(std::size_t j = i + 1; j < ind.size(); ++j) {
				product *= dim[j];
			}
			sum += product * ind[i];
		}
		return sum;
	}
	// Copy data from a to *this (after ReSize)
	void Copy(const MTensor& a) {
		T* p = data + Prod(dim);
		T* q = a.data + Prod(dim);
		while(p > data) {
			*--p = *--q;
		}
	}
	// Re-Allocate
	void ReSize(std::vector<std::size_t> newdim) {
		if(data != nullptr) {
			delete[] data;
		}
		dim = newdim;
		data = new T[Prod(newdim)];
	}
	// Product to get array volume (memory)
	std::size_t Prod(const std::vector<std::size_t>& x) const {
		std::size_t product = 1;
		for(std::size_t i = 0; i < x.size(); ++i) {
			product *= x[i];
		}
		return product;
	}

	// Dimensions
	std::vector<std::size_t> dim;

	T* data;
};

} // gra namespace ends

#endif
