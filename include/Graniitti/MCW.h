// Monte Carlo Weight Objects [HEADER ONLY file]
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MCW_H
#define MCW_H

// C++
#include <complex>
#include <random>
#include <valarray>
#include <vector>

#include "Graniitti/MMath.h"


namespace gra {
namespace kinematics {

// A simple Monte Carlo weight (container)
class MCW {
 public:
  MCW() {
    W  = 0.0;
    W2 = 0.0;
    N  = 0.0;
  }
  MCW(double w, double w2, double n) {
    W  = w;
    W2 = w2;
    N  = n;
  }
  // + operator
  MCW operator+(const MCW &obj) {
    MCW res;
    res.W  = W + obj.W;
    res.W2 = W2 + obj.W2;
    res.N  = N + obj.N;
    return res;
  }
  // += operator
  MCW &operator+=(const MCW &obj) {
    this->W += obj.W;
    this->W2 += obj.W2;
    this->N += obj.N;
    return *this;  // return by reference
  }
  // * operator (does scale W and W^2 but not N)
  MCW operator*(double scale) { return MCW(GetW() * scale, GetW2() * scale * scale, GetN()); }

  // Push new weight
  void Push(double w) {
    W += w;
    W2 += w * w;
    N += 1.0;
  }
  // MC estimate of the integral
  double Integral() const { return (N > 0.0) ? W / N : 0.0; }
  // MC estimate of the integral error squared (standard error of the mean)
  double IntegralError2() const {
    if (N > 0.0) {
      return (W2 / N - gra::math::pow2(W / N)) / N;
    } else {
      return 0.0;
    }
  }
  double IntegralError() const { return gra::math::msqrt(IntegralError2()); }
  double GetW() const { return W; }
  double GetW2() const { return W2; }
  double GetN() const { return N; }
  void   SetW(double w) { W = w; }
  void   SetW2(double w2) { W2 = w2; }
  void   SetN(double n) { N = n; }

 private:
  double W;   // sum of weights
  double W2;  // sum of weights^2
  double N;   // trials
};


// A weighted sum of MCW container integral values (use with VEGAS, for example)
class MCWSUM {
 public:
  MCWSUM() {}

  void Add(const MCW &x, double weight) {
    wsum += weight;
    wintsum += weight * x.Integral();
    error2sum += (weight * weight) * x.IntegralError2();
  }
  double Integral() const {
    if (wsum > 0.0) {
      return wintsum / wsum;
    } else {
      return 0.0;
    }
  }
  double IntegralError2() const {
    if (wsum > 0.0) {
      return error2sum / gra::math::pow2(wsum);
    } else {
      return 0.0;
    }
  }
  double IntegralError() const { return gra::math::msqrt(IntegralError2()); }

 private:
  // \sum_i weight_i
  double wsum = 0.0;

  // \sum_i weight_i integral_i
  double wintsum = 0.0;

  // \sum_i weight_i^2 error_i^2
  double error2sum = 0.0;
};

}  // namespace kinematics
}  // namespace gra

#endif