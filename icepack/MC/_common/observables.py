# Shared observables for phase space and cascade proposal comparisons
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import importlib
from functools import partial

import numpy as np
from core.analysis import obs

COMMON = importlib.import_module('icepack._common.obs')
FOUR_BODY = importlib.import_module('core.analysis.observables.four_body')


# Copy a physical projector with matching histogram and display bounds
def observable(name, lower, upper, count=20):
    source = FOUR_BODY if name.startswith('4body_') else COMMON
    output = copy.deepcopy(getattr(source, f'obs_{name}'))
    output['bins'] = np.linspace(lower, upper, count + 1)
    output['xlim'] = (lower, upper)
    return output


# Build the common absolute rate and production observables for one mass region
def production(lower, upper):
    return {
        'obs_cross_section': copy.deepcopy(COMMON.obs_cross_section),
        'obs_M': observable('M', lower, upper),
        'obs_Rap': observable('Rap', -2.5, 2.5),
        'obs_dPhi_pp': observable('dPhi_pp', 0.0, 180.0),
    }


# Build ordered pair masses and Jacob-Wick decay angle distributions
def decay_pair(mass_a, mass_b):
    return {
        'obs_4body_M_A': observable('4body_M_A', *mass_a),
        'obs_4body_M_B': observable('4body_M_B', *mass_b),
        'obs_4body_cos1': observable('4body_cos1', -1.0, 1.0),
        'obs_4body_cos2': observable('4body_cos2', -1.0, 1.0),
        'obs_4body_phi12': observable('4body_phi12', -180.0, 180.0),
    }


# Reconstruct flavour-labelled W or Z decay systems without a pion pairing hypothesis
def flavour_projection(event, *, pairs, field):
    records = obs.proj_central_particle_records(event)
    particles = []
    for pair in pairs:
        daughters = []
        for pid in pair:
            matches = [record['p4'] for record in records if record['pid'] == pid]
            if len(matches) != 1:
                raise ValueError(f'Sampling decay requires exactly one particle with PDG {pid}')
            daughters.append(matches[0])
        particles.append(daughters)
    mothers = [pair[0] + pair[1] for pair in particles]
    if field == '4body_M_A':
        return mothers[0].m
    if field == '4body_M_B':
        return mothers[1].m
    central = mothers[0] + mothers[1]
    helicities = []
    for daughters, mother in zip(particles, mothers, strict=True):
        mother_cm = mother.copy()
        mother_cm.boost(b=central, sign=-1)
        daughters_cm = [particle.copy() for particle in daughters]
        for particle in daughters_cm:
            particle.boost(b=central, sign=-1)
        helicities.append(obs.HXFrame(p=daughters_cm, X=mother_cm)[0])
    if field == '4body_cos1':
        return helicities[0].costheta
    if field == '4body_cos2':
        return helicities[1].costheta
    if field == '4body_phi12':
        phi = helicities[0].phi + helicities[1].phi
        return np.rad2deg(np.arctan2(np.sin(phi), np.cos(phi)))
    raise ValueError(f'Unknown sampling decay observable {field}')


# Use the physical flavour assignments in generated electroweak cascades
def flavour_pair(pairs, mass_a, mass_b):
    output = decay_pair(mass_a, mass_b)
    for item in output.values():
        item['func'] = partial(flavour_projection, pairs=pairs, field=item['tag'])
    return output
