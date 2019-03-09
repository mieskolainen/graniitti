// 4-vectors [HEADER ONLY class]
// with standard metric (+,-,-,-) and MC initialization convention (px,py,pz,e)
//
// However, note that ^,% operators index in "normal textbook convention",
// to be more streamlined with Lorent index contractions of amplitudes
//
//  p^0 =  E
//  p^1 =  px
//  p^2 =  py
//  p^3 =  pz
// 
//  p%0 =  E
//  p%1 = -px
//  p%2 = -py
//  p%3 = -pz
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef M4VEC_H
#define M4VEC_H

#include <iostream>
#include <vector>


namespace gra {

// ----------------------------------------------------------------------
// C++17 allows class template argument default deduction without brackets
//
// M4Vec a;
// M4Vec<> b;
// M4Vec<double> c;
//
// a,b,c are all valid.
//
// template <typename T = double>
//
// Complex 4-vectors could be implemented this way, TBD.
// ----------------------------------------------------------------------


class M4Vec {

public:

// All default to zero
    M4Vec() {
        k = {0.0, 0.0, 0.0, 0.0};
    }
    // Initialize
    M4Vec(double x, double y, double z, double t) {
        k = {t, x,y,z};
    }
    // Copy constructor
    M4Vec(const M4Vec& rhs) {
        k = rhs.k;
    }
    
// SET methods
    void SetPx(double v) { k[X_] = v; }
    void SetPy(double v) { k[Y_] = v; }
    void SetPz(double v) { k[Z_] = v; }
    void SetE (double v) { k[E_] = v; }

    void SetX(double v)  { k[X_] = v; }
    void SetY(double v)  { k[Y_] = v; }
    void SetZ(double v)  { k[Z_] = v; }
    void SetT(double v)  { k[E_] = v; }
    
    void Set(double x, double y, double z, double t) {
        k[X_] = x;
        k[Y_] = y;
        k[Z_] = z;
        k[E_] = t;
    }
    
    void SetPxPyPz(double x, double y, double z) {
        k[X_] = x;
        k[Y_] = y;
        k[Z_] = z;
    }

    void SetPxPy(double x, double y) {
        k[X_] = x;
        k[Y_] = y;
    }

    void SetPzE(double z, double e) {
        k[Z_] = z;
        k[E_] = e;
    }

    // Apply metric tensor: eta_{\mu\nu} k^\nu = k_\mu
    void Flip3() {
        k[X_] = -k[X_];
        k[Y_] = -k[Y_];
        k[Z_] = -k[Z_];
    }


// GET methods
    double Px() const { return k[X_]; }
    double Py() const { return k[Y_]; }
    double Pz() const { return k[Z_]; }
    double E()  const { return k[E_]; }

    double X()  const { return k[X_]; }
    double Y()  const { return k[Y_]; }
    double Z()  const { return k[Z_]; }
    double T()  const { return k[E_]; }


// ALGEBRA methods

    // gamma = E/m = 1/sqrt(1-v^2/c^2) = 1/sqrt(1-beta^2)
    double Gamma() const     { return E() / M(); }
    double Beta()  const     { return P3mod() / E(); }

    // Space-time invariants
    double Invariant() const { return E()*E() - P3mod2(); }
    double M()         const { return (M2() > 0.0) ? std::sqrt(M2()) : -std::sqrt(-M2()); }
    double M2()        const { return Invariant(); }
    
    // Transverse 2-vector norm and norm^2
    double Perp()   const { return Pt(); }
    double Perp2()  const { return Pt2(); }
    double Pt()     const { return msqrt(Pt2()); }
    double Pt2()    const { return k[X_]*k[X_] + k[Y_]*k[Y_];}
    
    // Total 3-vector norm and norm^2
    double P3mod()  const { return msqrt(P3mod2()); }
    double P3mod2() const { return k[X_]*k[X_] + k[Y_]*k[Y_] + k[Z_]*k[Z_]; }
    
    // Transverse mass (invariant under boost in z-direction)
    double Mt()     const { return std::sqrt(Mt2()); }
    double Mt2()    const { return M2() + Pt2(); }
    
    // Coincides with transverse mass for single particle
    double Et()     const { return Mt(); }
    double Et2()    const { return Mt2(); }
    
    // Angles
    double Phi()      const { return std::atan2(Py(), Px()); } // y / x
    double Theta()    const { return std::atan2(Pt(), Pz()); } // |Pt| / z
    double CosTheta() const { return std::cos(Theta()); }
    
    
    // Pseudorapidity and rapidity (boost) in z-direction
    double Eta() const {
        return 0.5*std::log( (P3mod() + Pz()) / (P3mod() - Pz()) );
    }
    double Rap() const {
        return 0.5*std::log( (E() + Pz())   / (E() - Pz()) );
    }
    
    // Lightcone variable: k_+ = E + p_z
    double LightconePos() const { return E() + Pz(); }

    // Lightcone variable: k_- = E - p_z
    double LightconeNeg() const { return E() - Pz(); }


// SPINOR-HELICITY CO-VARIABLES

    std::complex<double> ComplexPt() const {
        return Px() + std::complex<double>(0,1)*Py();
    }
    std::complex<double> ExpCPhi() const { 
        return ComplexPt() / msqrt(LightconePos()*LightconeNeg());
    }


// 2-BODY ALGEBRA

    // Minkowski 4-product
    double DotM(const M4Vec& rhs) const { return E()*rhs.E() - (Px()*rhs.Px() + Py()*rhs.Py() + Pz()*rhs.Pz()); }

    // 3-vector dot product
    double Dot3(const M4Vec& rhs) const { return Px()*rhs.Px() + Py()*rhs.Py() + Pz()*rhs.Pz(); }

    // 3-vector cross product (return vector with 0 energy/time)
    M4Vec Cross3(const M4Vec& rhs) const { 
        M4Vec a(Py() * rhs.Pz() - Pz() * rhs.Py(),
                Pz() * rhs.Px() - Px() * rhs.Pz(),
                Px() * rhs.Py() - Py() * rhs.Px(),
                0.0);
        return a;
    }

    // Azimuth angle difference between [-PI,PI]
    double DeltaPhi(const M4Vec& v) const {
        double D = Phi() - v.Phi();
        while (D >= PI) { D -= 2.0*PI; }
        while (D < -PI) { D += 2.0*PI; }
        return D;
    }

    // Azimuth angle between [0,PI]
    double DeltaPhiAbs(const M4Vec& v) const {
        return std::abs(DeltaPhi(v));
    }

// OPERATORS
    double  operator [] (size_t mu) const { return k[mu]; } // for reading only
    double& operator [] (size_t mu)       { return k[mu]; } // for substituting

    // Access operator in normal contravariant (upper index) indexing
    double operator ^ (size_t mu) const {
        return k[mu];
    }
    
    // Access operator with simultaneous lowering with metric (covariant index)
    double operator % (size_t mu) const {
        if (mu == E_) {
            return k[E_];
        } else {
            return -k[mu];
        }
        throw std::out_of_range("M4Vec::operator %% Index out of bounds!");        
    }

    // Minkowski scalar product
    double operator * (const M4Vec& rhs) const {
        return DotM(rhs);
    }

    // 4-vector + 4-vector
    M4Vec operator + (const M4Vec& rhs) const {
        return M4Vec(k[X_] + rhs.k[X_], k[Y_] + rhs.k[Y_], k[Z_] + rhs.k[Z_], k[E_] + rhs.k[E_]);
    }
    M4Vec operator - (const M4Vec& rhs) const {
        return M4Vec(k[X_] - rhs.k[X_], k[Y_] - rhs.k[Y_], k[Z_] - rhs.k[Z_], k[E_] - rhs.k[E_]);
    }

    // Flip sign of all components
    M4Vec operator - () const {
        return M4Vec(-k[X_], -k[Y_], -k[Z_], -k[E_]);
    }

    // 4-vector */ scalar
    M4Vec operator * (const double rhs) const {
        return M4Vec(k[X_] * rhs, k[Y_] * rhs, k[Z_] * rhs, k[E_] * rhs);
    }
    M4Vec operator / (const double rhs) const {
        return M4Vec(k[X_] / rhs, k[Y_] / rhs, k[Z_] / rhs, k[E_] / rhs);
    }

    // Comparison
    bool operator==(const M4Vec& rhs) const {
        const double EPS = 1e-10;
        return std::abs(k[E_] - rhs.k[E_]) < EPS &&
               std::abs(k[X_] - rhs.k[X_]) < EPS &&
               std::abs(k[Y_] - rhs.k[Y_]) < EPS &&
               std::abs(k[Z_] - rhs.k[Z_]) < EPS;
    }
    bool operator!=(const M4Vec& rhs) const { 
        return !(*this == rhs);
    }

    // 4-vector +-= 4-vector
    void operator += (const M4Vec& rhs) {
        k[E_] += rhs.k[E_];
        k[X_] += rhs.k[X_];
        k[Y_] += rhs.k[Y_];
        k[Z_] += rhs.k[Z_];
    }
    void operator -= (const M4Vec& rhs) {
        k[E_] -= rhs.k[E_];
        k[X_] -= rhs.k[X_];
        k[Y_] -= rhs.k[Y_];
        k[Z_] -= rhs.k[Z_];
    }

    // 4-vector */= scalar
    void operator *= (const double rhs) {
        k[E_] *= rhs;
        k[X_] *= rhs;
        k[Y_] *= rhs;
        k[Z_] *= rhs;
    }
    void operator /= (const double rhs) {
        k[E_] /= rhs;
        k[X_] /= rhs;
        k[Y_] /= rhs;
        k[Z_] /= rhs;
    }

    void Print(const std::string name = "") const {
        std::cout << "M4Vec::" << name << " Px / X: "
                  << Px() << ", Py / Y: " << Py() << ", Pz / Z: " << Pz()
                  << ", E / T: " << E() << ", M / dS: " << M() << std::endl;
    }
    
    // Particle 4-position starting starting propagation from (0,0,0,0)
    // - p is the particle 4-momentum in the lab frame
    // - tau0 is the particle flight time in its rest frame
    M4Vec PropagatePosition(double tau0, double scale) const {

        // Flight time in the lab frame
        const double gamma = Gamma();
        const double tau   = gamma * tau0; 

        // Velocity
        const double beta  = Beta();

        // End point 4-position
        return M4Vec(tau * c * beta * Px() / P3mod() * scale,
                     tau * c * beta * Py() / P3mod() * scale,
                     tau * c * beta * Pz() / P3mod() * scale,
                     tau * c * scale);
    }

private:
    
    static constexpr const double PI = 3.141592653589793238462643383279502884197169399375105820974944L;

    // speed of light, c = [m/s] (EXACT/DEFINITION)
    static constexpr const double c = 2.99792458E8; 

    // Indices
    static const int E_ = 0;
    static const int X_ = 1;
    static const int Y_ = 2;
    static const int Z_ = 3;

    // Safe sqrt
    double msqrt(double x) const {
        return std::sqrt(std::max(0.0, x));
    }
    // 4-vector
    std::vector<double> k;
};

} // gra namespace ends

#endif