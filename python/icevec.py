# Simple Lorentz vectors with HEP metric (+,---)
#
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np


def hepmc2vec4(p):
    """ HepMC3 python binding FourVector to vec4 """
    return vec4(x=p.px(), y=p.py(), z=p.pz(), t=p.e())


class vec4:
    """ Lorentz vectors """

    def __init__(self, x=None, y=None, z=None, t=None):

        if   (x is not None) and (y is not None) and (z is not None) and (t is not None):
            self._x, self._y, self._z, self._t = x,y,z,t
        elif (x is None) and (y is None) and (z is None) and (t is None):
            """ Empty initialization """
            self._x, self._y, self._z, self._t = 0,0,0,0
        else:
            raise Exception('vec4: Unknown initialization.')

    # -----------------------------------
    # Needed for python built-in sum()

    def __iter__(self):
        for value in self.values():
            yield value

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    # -----------------------------------

    # Addition
    def __add__(self, rhs):
        return vec4(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z, self.t + rhs.t)

    # Subtract
    def __sub__(self, other):
        return vec4(self.x - other.x, self.y - other.y, self.z - other.z, self.t - other.t)
    
    def __rsub__(self, other):   # not commutative operation
        return vec4(other.x - self.x, other.y - self.y, other.z - self.z, other.t - self.t)

    # Multiply
    def __mul__(self, other):

        # 4-dot
        if hasattr(other, 'x'):
            return self.dot4(other)

        # Scalar multiply
        else:
            return vec4(other*self.x, other*self.y, other*self.z, other*self.t)

    __rmul__ = __mul__   # commutative operation
    
    # Print
    def __str__(self):
        return f"[x = {self.x}, y = {self.y}, z = {self.z}, t = {self.t}]"

    def scale(self, a):
        # Scale all components
        self.x, self.y, self.z, self.t = a*self.x, a*self.y, a*self.y, a*self.t

    def dot4(self, other):
        # 4-vector dot product
        return self.t*other.t - (self.x*other.x + self.y*other.y + self.z*other.z)

    def dot3(self, other):
        # 3-vector dot product
        return self.x*other.x + self.y*other.y + self.z*other.z


    def setX(self, x):
        self._x = x
    def setY(self, y):
        self._y = y
    def setZ(self, z):
        self._z = z
    def setE(self, e):
        self._t = e
    def setXYZ(self, x, y, z):
        self._x = x
        self._y = y
        self._z = z
    

    def setPt2RapPhiM2(self, pt2, rap, phi, m2):

        mT = np.sqrt(m2 + pt2)

        e  = mT * np.cosh(rap)
        pz = mT * np.sinh(rap)
        px = np.sqrt(pt2) * np.cos(phi)
        py = np.sqrt(pt2) * np.sin(phi)

        self._x, self._y, self._z, self._t = px, py, pz, e


    def setPtEtaPhi(self, pt, eta, phi):
        self.setXYZ(pt*np.cos(phi), pt*np.sin(phi), pt/np.tan(2.0*np.arctan(np.exp(-eta))))


    def setMagThetaPhi(self, mag, theta, phi):
        self._x = mag * np.sin(theta) * np.cos(phi)
        self._y = mag * np.sin(theta) * np.sin(phi)
        self._z = mag * np.cos(theta)


    def setPxPyPzE(self, px, py, pz, e):
        self._x = px
        self._y = py
        self._z = pz
        self._t = e

    def setXYZT(self, x, y, z, t):
        self._x = x
        self._y = y
        self._z = z
        self._t = t

    def setXYZM(self, x, y, z, m):
        self._x = x
        self._y = y
        self._z = z
        self._t = np.sqrt(m**2 + x**2 + y**2 + z**2)

    def setP3(self, p3):
        self._x = p3[0]
        self._y = p3[1]
        self._z = p3[2]        

    def setPtEtaPhiM(self, pt, eta, phi, m):
        self.setPtEtaPhi(pt, eta, phi)
        self.setE(np.sqrt(m**2 + self.x**2 + self.y**2 + self.z**2))
    

    def phi_PIPI(self, x):
        # Transverse angle between [-pi,pi]
        while (x >= np.pi):
            x -= 2*np.pi
        while (x < -np.pi):
            x += 2*np.pi
        return x

    def deltaphi(self, v):
        return self.phi_PIPI(self.phi - v.phi)

    def abs_delta_phi(self, v):
        return np.abs(self.deltaphi(v))

    def deltaR(self, v):
        # (eta,phi)-plane separation
        deta = self.eta - v.eta
        dphi = self.phi_PIPI(self.phi - v.phi)
        return np.sqrt(deta**2 + dphi**2)
    
    
    ## Return properties
    
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
        return np.sqrt(M2) if M2 > 0 else 0.0 # Return 0 also if q^2 < 0

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
        if (self.x == 0) and (self.y == 0) and (self.z == 0):
            return 0.0
        else:
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
        if (self.x == 0.0) and (self.y == 0.0):
            return 0.0
        else:
            return np.arctan2(self.y, self.x)

    @property
    # Longitudinal angle cosine
    def costheta(self):
        return np.cos(self.theta)

    @property
    # Longitudinal angle
    def theta(self):
        if (self.x == 0.0) and (self.y == 0.0) and (self.z == 0.0):
            return 0.0
        else:
            return np.arctan2(self.pt, self.z)

    @property
    # Longitudinal rapidity
    def rapidity(self):
        return 0.5*np.log((self.e + self.pz)/(self.e - self.pz))

    @property
    # Pseudorapidity
    def eta(self):
        cosTheta = self.costheta
        if (cosTheta**2 < 1):
            return -0.5*np.log((1.0 - cosTheta) / (1.0 + cosTheta))
        if (self.z == 0):
            return 0
        if (self.z > 0):
            return 10e10
        else:
            return -10e10;

    @property
    def abseta(self):
        return np.abs(self.eta)

    @property
    def x(self):
        return self._x
    @property
    def y(self):
        return self._y
    @property
    def z(self):
        return self._z
    @property
    def t(self):
        return self._t


    @property
    def px(self):
        return self._x
    @property
    def py(self):
        return self._y
    @property
    def pz(self):
        return self._z
    @property
    def e(self):
        return self._t


    def rotateSO3(self, R):
        v = np.dot(R, self.p3)
        self._x = v[0]
        self._y = v[1]
        self._z = v[2]

    def rotateX(self, angle):
        s = np.sin(angle)
        c = np.cos(angle)
        y = self.y
        
        self._y = c*y - s*self.z
        self._z = s*y + c*self.z
    
    def rotateY(self, angle):
        s = np.sin(angle)
        c = np.cos(angle)
        z = self.z
        
        self._z = c*z - s*self.x
        self._x = s*z + c*self.x

    def rotateZ(self, angle):
        s = np.sin(angle)
        c = np.cos(angle)
        x = self.x
        
        self._x = c*x - s*self.y
        self._y = s*x + c*self.y


    def boost(self, b, sign=-1):
        """
        Lorentz boost
        Args:
                   b : Boost 4-momentum (e.g. system)
                sign : 1 or -1 (direction of the boost, into the rest (-1) or out (1))
        """
        if np.isclose(b.e, 0.0):
            raise Exception(__name__ + '.boost: Input boost 4-momentum with e = 0')
        
        # Beta and gamma factors    
        betaX = sign*b.px / b.e  # px / E
        betaY = sign*b.py / b.e  # py / E
        betaZ = sign*b.pz / b.e  # pz / E
        gamma = b.gamma          # E  / m

        # Momentum and energy product
        aux1  = betaX*self.px + betaY*self.py + betaZ*self.pz
        aux2  = gamma*(gamma*aux1 / (1.0 + gamma) + self.e)

        # Lorentz boost
        self._x = self.px + aux2*betaX
        self._y = self.py + aux2*betaY
        self._z = self.pz + aux2*betaZ
        self._e = gamma*(self.e + aux1)

