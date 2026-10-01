# Direct screened equivalent-photon cross sections for the HERA derivation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from scipy.constants import alpha, physical_constants
from scipy.integrate import quad
from scipy.special import j0, j1, jv

from . import model as hera_model


# Map Gauss-Legendre nodes and weights onto a finite interval
def gauss(count, low, high):
    x, w = np.polynomial.legendre.leggauss(count)
    return low + (x + 1) * (high - low) / 2, w * (high - low) / 2


# Read particle masses from the same PDG table as the generator
def masses():
    return {int(row[0]): float(row[1]) for line in (hera_model.REPOSITORY_ROOT / 'modeldata/mass_width_2026.mcd').read_text().splitlines()
            if (row := line.split()) and row[0] in {'2212', '211', '321', '11'}}


# Compute the configured Sachs form factors and their Dirac and Pauli combinations
# [REFERENCE: Kelly, Phys. Rev. C 70 (2004) 068202]
def electromagnetic(q2, mass, model, settings):
    tau = q2 / (4 * mass**2)
    moment = physical_constants['proton mag. mom. to nuclear magneton ratio'][0]
    if model == 'DIPOLE':
        electric = (1 + q2 / settings['dipole_scale2'])**-2
        magnetic = moment * electric
    elif model == 'KELLY':
        electric, magnetic = [(1 + row[0] * tau) / (1 + tau * (row[1] + tau * (row[2] + tau * row[3])))
                              for row in settings['kelly']]
        magnetic *= moment
    else:
        raise ValueError(f'Unsupported proton electromagnetic form: {model}')
    return (electric + tau * magnetic) / (1 + tau), (magnetic - electric) / (1 + tau)


# Compute the single-channel SOFT S matrix using the generator Fourier and phase conventions
def smatrix(b, energy, proton, general, numerics):
    model = general['PARAM_SOFT']['MODEL'][general['PARAM_SOFT']['active_model']]
    settings = model['EIKONAL']
    if len(proton.parameters) != 1 or settings['helicity'] or settings['screening_exchanges'] != ['P']:
        raise ValueError('The direct HERA pp calculation requires one SOFT channel with P screening and no proton spin flip')
    pomeron = model['EXCHANGE']['P']
    if pomeron['eta_mode'] not in ('rotating', 'rotating_t0'):
        raise ValueError('The direct HERA pp calculation requires a rotating SOFT phase')
    q, weight = gauss(numerics['momentum'], 0, numerics['momentum_max'])
    trajectory = pomeron['alpha'][0] - pomeron['alpha'][1] * q**2
    if general['PARAM_SOFT']['EXCHANGE_DEF']['P']['trajectory_mode'] == 'pion_loop':
        pion2 = masses()[211]**2
        ratio = 4 * pion2 / q**2
        root = np.sqrt(1 + ratio)
        loop = 4 / ratio / (1 + q**2 / model['pion_loop_scale2'])**2 * (
            2 * ratio - (1 + ratio)**1.5 * np.log1p(2 / (root - 1)) + np.log(1 / pion2))
        trajectory -= (2 * proton.beam_residue_per_gev / 3)**2 * pion2 / (32 * np.pi**3) * loop
    phase = trajectory if pomeron['eta_mode'] == 'rotating' else pomeron['alpha'][0]
    born = (-pomeron['sign'] * np.exp(-0.5j * np.pi * phase) * energy**(2 * (trajectory - 1))
            * (proton.beam_residue_per_gev * hera_model.proton_profile(q*q, proton))**2)
    # chi(b) = integral q dq J0(qb) A_Born(s,-q^2)/(4 pi s)
    chi = j0(np.asarray(b)[..., None] * q) @ (weight * q * born) / (4 * np.pi)
    if settings['unitarization'] == 'exp' or np.isclose(settings['q'], 1):
        return np.exp(1j * chi)
    if settings['unitarization'] != 'q_exp':
        raise ValueError('Unsupported SOFT unitarization')
    return (1 + (1 - settings['q']) * 1j * chi)**(1 / (1 - settings['q']))


# Fold the transverse Dirac-Pauli photon current with a screened gamma-proton amplitude
# [REFERENCE: arXiv:0705.2887, Eqs. (2.2)-(2.8)]
class Spectrum:
    """Direct small-x photoproduction quadrature"""

    # Prepare only fixed quadrature factors, photon currents and the SOFT absorption profile
    def __init__(self, energy, rapidity, mass, proton, general, numerics, physics):
        self.energy, self.mass, self.proton = energy, mass, proton
        self.q, self.qw = gauss(numerics['momentum'], 0, numerics['momentum_max'])
        self.r, rw = gauss(numerics['impact'], 0, numerics['impact_max'])
        self.rw = self.r * rw
        self.hankel = j0(self.q[:, None] * self.r) * (self.q * self.qw)[:, None]
        self.profile = hera_model.proton_profile(self.q**2, proton)
        self.rotating = general['PARAM_REGGE']['photoprod_eta_mode'] == 'rotating'
        self.x = mass / energy * np.exp(np.array([rapidity, -np.asarray(rapidity)]))
        if np.any(self.x <= 0) or np.any(self.x >= physics['xi_max']):
            raise ValueError('Photon momentum fractions exceed the configured small-x domain')
        self.w = energy * np.sqrt(self.x)
        self.flux, electric, flip = [], [], []
        mp = masses()[2212]
        form = general['PARAM_STRUCTURE']['EM']
        for x in self.x.ravel():
            a2 = (x * mp)**2
            f1, f2 = electromagnetic((self.q**2 + a2) / (1 - x), mp, form, physics)
            coefficient = np.sqrt(alpha * (1 - x)) / np.pi
            electric.append(coefficient * ((self.qw * self.q**2 / (self.q**2 + a2) * f1) @ j1(self.q[:, None] * self.r)))
            flips = [coefficient / (2 * mp) * ((self.qw * self.q**3 / (self.q**2 + a2) * f2)
                     @ jv(order, self.q[:, None] * self.r)) for order in (0, 2)]
            flip.append((flips[0]**2 + flips[1]**2) / 2)

            # Integrate the unscreened photon flux in log(q^2+a^2) to resolve its collinear logarithm
            def density(logq, x=x, a2=a2):
                den = np.exp(logq)
                q2 = den - a2
                f1, f2 = electromagnetic(den / (1 - x), mp, form, physics)
                return alpha / np.pi * (1 - x) * q2 / den * (f1*f1 + q2 / (4 * mp**2) * f2*f2)

            self.flux.append(quad(density, np.log(a2), np.log(numerics['momentum_max']**2 + a2), epsabs=numerics['flux_atol'])[0])
        self.flux = np.asarray(self.flux).reshape(self.x.shape)
        self.electric = np.asarray(electric).reshape(*self.x.shape, len(self.r))
        self.density = self.electric**2 + np.asarray(flip).reshape(self.electric.shape)
        phi, pw = gauss(numerics['angle'], 0, 2 * np.pi)
        separation = np.sqrt(np.maximum(0, self.r[:, None, None]**2 + self.r[None, :, None]**2
                                       - 2 * self.r[:, None, None] * self.r[None, :, None] * np.cos(phi)))
        # Resolve the smooth SOFT profile independently of photon and meson fit parameters
        b = np.linspace(0, 2 * numerics['impact_max'], numerics['soft_nodes'])
        absorption = 1 - abs(smatrix(b, energy, proton, general, numerics))**2
        kernel = np.interp(separation, b, absorption)
        measure = 4 * np.pi * np.outer(self.rw, self.rw)
        self.k0 = (kernel @ (pw / (2 * np.pi))) * measure
        self.k1 = (kernel @ (pw * np.cos(phi) / (2 * np.pi))) * measure

    # Integrate both photon directions and their coherent absorption interference in microbarn
    def predict(self, forward, slope, delta, shrinkage, w0, curve=None, screened=True):
        energy = np.log(self.w / w0)
        amplitude = np.sqrt(forward) * np.exp(delta * energy / 2)
        phase = np.full_like(energy, 1 + delta / 4)
        if curve is not None:
            amplitude *= np.sqrt(hera_model.dlog_factor(self.w, w0, *curve))
            phase += hera_model.dlog_power(self.w, w0, *curve) / 4
        amplitude = amplitude * (-np.exp(-0.5j * np.pi * phase))
        b = slope + 4 * shrinkage * energy
        profile = self.profile * np.exp(-0.5 * b[..., None] * self.q**2)
        integral = np.sum(abs(profile)**2 * (2 * self.q * self.qw), axis=-1)
        born = np.sum(self.flux * abs(amplitude)**2 * integral, axis=0)
        if not screened:
            return born
        if self.rotating:
            profile = profile * np.exp(0.5j * np.pi * shrinkage * self.q**2)
        hard = (profile @ self.hankel) * amplitude[..., None]
        diagonal = np.einsum('dyr,rs,dys->y', self.density, self.k0, abs(hard)**2, optimize=True)
        cross = self.electric[0] * hard[1].conj()
        other = self.electric[1] * hard[0].conj()
        interference = 2 * np.real(np.einsum('yr,rs,ys->y', cross, self.k1, other.conj(), optimize=True))
        return born - diagonal - interference
