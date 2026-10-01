# Validate flavour-resolved cascade observables used in F and C comparisons
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pytest
from core.kinematics.vec4 import vec4

from icepack.MC._common.observables import flavour_projection
from tests.technical.support.hepmc import decay_event


# Construct two massive parents with different rest-frame decay angles
def flavour_records(pairs):
    parents = [vec4(12.0, 5.0, 30.0, np.sqrt(80.0**2 + 1069.0)),
               vec4(-12.0, -5.0, -30.0, np.sqrt(90.0**2 + 1069.0))]
    records = []
    for pair, parent, mass, theta, phi in zip(
        pairs, parents, (80.0, 90.0), (0.7, 1.1), (0.3, -0.9), strict=True
    ):
        momentum = 0.5 * mass * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
        for pid, sign in zip(pair, (1.0, -1.0), strict=True):
            daughter = vec4(*(sign * momentum), 0.5 * mass)
            daughter.boost(b=parent, sign=1)
            records.append({'pid': pid, 'p4': daughter})
    return records


# Intermediate masses and decay angles must survive record ordering, boosts and rotations
@pytest.mark.parametrize('pairs', [((-13, 14), (11, -12)), ((-13, 13), (2, -2))])
def test_flavour_decay_obs_follow_physical_pairs(pairs):
    records = flavour_records(pairs)
    event = decay_event([(r['pid'], r['p4'], 1) for r in records])
    fields = ('4body_M_A', '4body_M_B', '4body_cos1', '4body_cos2', '4body_phi12')
    before = {field: flavour_projection(event, pairs=pairs, field=field) for field in fields}
    assert before['4body_M_A'] == pytest.approx(80.0)
    assert before['4body_M_B'] == pytest.approx(90.0)
    # Project rest-frame directions onto independent Cartesian helicity axes
    phases = []
    for index, (sign, theta, phi) in enumerate(((1, 0.7, 0.3), (-1, 1.1, -0.9)), start=1):
        z = sign * np.array([12.0, 5.0, 30.0]) / np.sqrt(1069.0)
        y = np.cross(z, [0.0, 0.0, 1.0])
        y /= np.linalg.norm(y)
        x = np.cross(y, z)
        direction = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
        assert before[f'4body_cos{index}'] == pytest.approx(direction @ z, rel=0.0, abs=1e-12)
        phases.append(np.arctan2(direction @ y, direction @ x))
    expected_phi = (np.rad2deg(sum(phases)) + 180.0) % 360.0 - 180.0
    assert before['4body_phi12'] == pytest.approx(expected_phi, rel=0.0, abs=1e-12)
    records.reverse()
    boost = vec4(0.0, 0.0, 3.0, 5.0)
    for record in records:
        record['p4'].rotateZ(0.63)
        record['p4'].boost(b=boost, sign=1)
    event = decay_event([(r['pid'], r['p4'], 1) for r in records])
    after = {field: flavour_projection(event, pairs=pairs, field=field) for field in fields}
    assert after == pytest.approx(before, rel=1e-10, abs=1e-10)


# Missing or repeated flavours must not silently select another decay assignment
@pytest.mark.parametrize('duplicate', [False, True])
def test_flavour_ambiguous_decay(duplicate):
    pairs = ((-13, 14), (11, -12))
    records = flavour_records(pairs)
    if duplicate:
        records.append(records[0])
    else:
        records = records[1:]
    event = decay_event([(r['pid'], r['p4'], 1) for r in records])
    with pytest.raises(ValueError, match='exactly one particle'):
        flavour_projection(event, pairs=pairs, field='4body_M_A')
