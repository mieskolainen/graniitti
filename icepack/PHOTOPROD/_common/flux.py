# HERA equivalent photon flux for HERA pion pair measurements
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pyjson5
from core.io.cache import cache
from scipy.constants import alpha
from scipy.integrate import quad


# Integrate EPA with optional rho VMD and sigma_L/sigma_T = Q2/m_rho^2
# [REFERENCE: ZEUS, arXiv:hep-ex/9712020, Eqs. (10)-(13), xi^2 = 1]
@cache
def rho_flux(gencard, cuts, vmd=True):
    root = Path(__file__).resolve().parents[3]
    masses = {}
    for line in (root / 'modeldata/mass_width_2026.mcd').read_text().splitlines():
        fields = line.split()
        if fields and fields[0] in {'11', '113'}:
            masses[int(fields[0])] = float(fields[1])
    energy = pyjson5.decode(Path(gencard).read_text())['SCATTERING']['ENERGY']
    y_min, y_max = np.asarray(cuts['W'])**2 / (4.0 * np.prod(energy))
    mass2 = masses[113]**2

    # Integrate logarithmically to retain the electron mass term near the collinear limit
    def flux_y(y):
        qmin = masses[11]**2 * y*y / (1.0-y)

        # Compute the photon density per unit log Q2
        def flux_q(logq):
            q2 = np.exp(logq)
            transverse = ((1.0 + (1.0-y)**2) - 2.0*(1.0-y)*qmin/q2) / y
            if not vmd:
                return transverse
            return (transverse + 2.0*(1.0-y)/y*q2/mass2) / (1.0+q2/mass2)**2

        return quad(flux_q, np.log(qmin), np.log(cuts['Q2_MAX']))[0]

    return alpha / (2.0*np.pi) * quad(flux_y, y_min, y_max)[0]
