# Simple Lorentz vectors with HEP metric (+,---)
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np


def hepmc2vec4(p):
    """HepMC3 python binding FourVector to vec4"""
    return vec4(x=p.px(), y=p.py(), z=p.pz(), t=p.e())


class vec4:
    """Lorentz vectors"""

    # Initialize one explicit or zero Lorentz vector
    def __init__(self, x=None, y=None, z=None, t=None):
        if (x is not None) and (y is not None) and (z is not None) and (t is not None):
            self._x, self._y, self._z, self._t = x, y, z, t
        elif (x is None) and (y is None) and (z is None) and (t is None):
            self._x, self._y, self._z, self._t = 0, 0, 0, 0
        else:
            raise Exception("vec4: Unknown initialization.")
    # Needed for python built-in sum()

    # Iterate over Cartesian momentum and energy components
    def __iter__(self):
        return iter((self.x, self.y, self.z, self.t))

    # Support the zero initializer used by built-in sum
    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other)
    # Addition
    def __add__(self, rhs):
        return vec4(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z, self.t + rhs.t)

    # Subtract
    def __sub__(self, other):
        return vec4(self.x - other.x, self.y - other.y, self.z - other.z, self.t - other.t)

    def __rsub__(self, other):  # not commutative operation
        return vec4(other.x - self.x, other.y - self.y, other.z - self.z, other.t - self.t)

    # Multiply
    def __mul__(self, other):
        if hasattr(other, "x"):
            return self.dot4(other)
        return vec4(other * self.x, other * self.y, other * self.z, other * self.t)

    __rmul__ = __mul__  # commutative operation

    # Print
    def __str__(self):
        return f"[x = {self.x}, y = {self.y}, z = {self.z}, t = {self.t}]"

    # Copy momentum components and their shared references without pickle reconstruction
    def __deepcopy__(self, memo):
        result = type(self).__new__(type(self))
        memo[id(self)] = result
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        return result

    # Compute an independent vector copy
    def copy(self):
        return vec4(self.x, self.y, self.z, self.t)

    # Scale every vector component in place
    def scale(self, a):
        self.setXYZT(a * self.x, a * self.y, a * self.z, a * self.t)

    # Compute the Minkowski inner product
    def dot4(self, other):
        return self.t * other.t - (self.x * other.x + self.y * other.y + self.z * other.z)

    # Compute the Euclidean spatial inner product
    def dot3(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z

    # Compute the opening angle between spatial momenta
    def angle(self, other):
        norm = self.p3mod * other.p3mod
        if norm <= 0:
            return 0.0
        return np.arccos(np.clip(self.dot3(other) / norm, -1.0, 1.0))

    # Set the x momentum component
    def setX(self, x):
        self._x = x

    # Set the y momentum component
    def setY(self, y):
        self._y = y

    # Set the z momentum component
    def setZ(self, z):
        self._z = z

    # Set the energy component
    def setE(self, e):
        self._t = e

    # Set all spatial momentum components
    def setXYZ(self, x, y, z):
        self._x, self._y, self._z = x, y, z

    # Set a vector from transverse momentum squared, rapidity, phi and mass squared
    def setPt2RapPhiM2(self, pt2, rap, phi, m2):
        mT = np.sqrt(m2 + pt2)
        e = mT * np.cosh(rap)
        pz = mT * np.sinh(rap)
        px = np.sqrt(pt2) * np.cos(phi)
        py = np.sqrt(pt2) * np.sin(phi)
        self._x, self._y, self._z, self._t = px, py, pz, e

    # Set spatial momentum from transverse momentum, pseudorapidity and phi
    def setPtEtaPhi(self, pt, eta, phi):
        self.setXYZ(pt * np.cos(phi), pt * np.sin(phi), pt * np.sinh(eta))

    # Set spatial momentum from magnitude and spherical angles
    def setMagThetaPhi(self, mag, theta, phi):
        self.setXYZ(mag * np.sin(theta) * np.cos(phi), mag * np.sin(theta) * np.sin(phi), mag * np.cos(theta))

    # Set Cartesian momentum and energy components
    def setXYZT(self, x, y, z, t):
        self._x, self._y, self._z, self._t = x, y, z, t

    # Set momentum and energy using HEP component names
    def setPxPyPzE(self, px, py, pz, e):
        self.setXYZT(px, py, pz, e)

    # Set Cartesian momentum and derive energy from invariant mass
    def setXYZM(self, x, y, z, m):
        self.setXYZT(x, y, z, np.sqrt(m**2 + x**2 + y**2 + z**2))

    # Set spatial momentum from one three-component sequence
    def setP3(self, p3):
        self.setXYZ(p3[0], p3[1], p3[2])

    # Set a massive vector from transverse momentum, pseudorapidity and phi
    def setPtEtaPhiM(self, pt, eta, phi, m):
        self.setPtEtaPhi(pt=pt, eta=eta, phi=phi)
        self.setE(np.sqrt(m**2 + self.x**2 + self.y**2 + self.z**2))

    # Canonicalize a transverse angle to the interval [-pi, pi)
    def phi_PIPI(self, x):
        return (x + np.pi) % (2.0 * np.pi) - np.pi

    # Compute the signed canonical azimuthal separation
    def deltaphi(self, v):
        return self.phi_PIPI(self.phi - v.phi)

    # Compute the absolute canonical azimuthal separation
    def abs_delta_phi(self, v):
        return np.abs(self.deltaphi(v))

    # Compute the pseudorapidity-azimuth distance
    def deltaR(self, v):
        deta = self.eta - v.eta
        return np.sqrt(deta**2 + self.deltaphi(v) ** 2)


    @property
    # 3-vector
    def p3(self):
        return np.array([self.x, self.y, self.z])

    @property
    # Transverse mass
    def mt(self):
        return np.sqrt(self.m2 + self.pt2)

    @property
    # Mass
    def m(self):
        M2 = self.m2
        return np.sqrt(M2) if M2 > 0 else 0.0  # Compute 0 also if q^2 < 0

    @property
    # Mass squared
    def m2(self):
        return self.t**2 - (self.x**2 + self.y**2 + self.z**2)

    @property
    # 3-vector norm squared
    def p3mod2(self):
        return self.x**2 + self.y**2 + self.z**2

    @property
    # 3-vector norm
    def p3mod(self):
        return np.sqrt(self.p3mod2)

    @property
    # Lorentz beta
    def beta(self):
        return self.p3mod / self.e

    @property
    # Lorentz gamma
    def gamma(self):
        return self.e / self.m

    @property
    # Transverse momentum
    def pt(self):
        return np.sqrt(self.pt2)

    @property
    # Transverse momentum squared
    def pt2(self):
        return self.x**2 + self.y**2

    @property
    # Transverse plane angle
    def phi(self):
        return np.arctan2(self.y, self.x)

    @property
    # Longitudinal angle cosine
    def costheta(self):
        return np.cos(self.theta)

    @property
    # Longitudinal angle
    def theta(self):
        return np.arctan2(self.pt, self.z)

    @property
    # Compute longitudinal rapidity with signed lightlike limits
    def rapidity(self):
        
        energy, longitudinal = float(self.e), float(self.pz)
        if not np.isfinite(energy) or not np.isfinite(longitudinal) or energy <= 0.0:
            return np.nan
        gap = energy - abs(longitudinal)
        if gap < 0.0:
            return np.nan
        if gap <= 0.0:
            return np.copysign(np.inf, longitudinal)
        beta = longitudinal / energy
        if abs(beta) < 0.5:
            return np.arctanh(beta)
        magnitude = 0.5 * (np.log(energy) - np.log(gap) + np.log1p(abs(beta)))
        return np.copysign(magnitude, longitudinal)

    @property
    # Pseudorapidity
    def eta(self):
        pt = np.hypot(self.x, self.y)
        if pt > 0.0:
            return np.arcsinh(self.z / pt)
        if abs(self.z) > 0.0:
            return np.copysign(np.inf, self.z)
        return 0.0

    @property
    def abseta(self):
        return np.abs(self.eta)

    @property
    def x(self):
        return self._x

    # Compute the Cartesian y component
    @property
    def y(self):
        return self._y

    # Compute the Cartesian z component
    @property
    def z(self):
        return self._z

    # Compute the energy component
    @property
    def t(self):
        return self._t

    px, py, pz, e = x, y, z, t

    # Rotate the spatial momentum with one three-dimensional matrix
    def rotateSO3(self, R):
        self.setP3(np.asarray(R) @ self.p3)

    # Rotate the spatial momentum around the x axis
    def rotateX(self, angle):
        s = np.sin(angle)
        c = np.cos(angle)
        self.rotateSO3(((1.0, 0.0, 0.0), (0.0, c, -s), (0.0, s, c)))

    # Rotate the spatial momentum around the y axis
    def rotateY(self, angle):
        s = np.sin(angle)
        c = np.cos(angle)
        self.rotateSO3(((c, 0.0, s), (0.0, 1.0, 0.0), (-s, 0.0, c)))

    # Rotate the spatial momentum around the z axis
    def rotateZ(self, angle):
        s = np.sin(angle)
        c = np.cos(angle)
        self.rotateSO3(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0)))

    # Apply one Lorentz boost in or out of a system rest frame
    def boost(self, b, sign=-1):
        if not np.isfinite(b.e) or b.e <= 0.0 or sign not in {-1, 1}:
            raise ValueError(__name__ + ".boost: Require positive finite energy and sign +/-1")

        beta = np.array((b.x, b.y, b.z), dtype=np.longdouble)
        beta *= sign / np.longdouble(b.e)
        beta2 = beta @ beta
        if not np.isfinite(beta2) or beta2 >= 1.0:
            raise ValueError(__name__ + ".boost: Boost momentum must be timelike")
        momentum = np.asarray(tuple(self), dtype=np.longdouble)
        if not np.all(np.isfinite(momentum)):
            raise ValueError(__name__ + ".boost: Require finite four-vector components")
        gamma = 1.0 / np.sqrt(1.0 - beta2)
        dot = beta @ momentum[:3]
        result = np.empty(4, dtype=np.longdouble)
        result[:3] = momentum[:3] + gamma * (gamma * dot / (1.0 + gamma) + momentum[3]) * beta
        result[3] = gamma * (momentum[3] + dot)
        if not np.all(np.abs(result) <= np.finfo(float).max):
            raise ValueError(__name__ + ".boost: Boosted components exceed finite double precision")
        self.setXYZT(*result.astype(float))
