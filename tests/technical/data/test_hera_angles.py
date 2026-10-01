# HERA angular units, helicity covariance and differential photon-flux normalization
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from core.analysis import obs
from core.io import steering
from core.kinematics.vec4 import vec4
from core.stats.uncertainty import source_covariance
from scipy.constants import alpha
from scipy.integrate import quad
from scipy.spatial.transform import Rotation

from icepack._common.hepdata import read_table
from icepack.PHOTOPROD._common import angles, exclusive_jpsi, hera_reader
from icepack.PHOTOPROD._common.flux import rho_flux
from icepack.PHOTOPROD.H1_1798511 import reader as h1_reader
from tests.technical.data.test_hera_dissociation import photo_event

ROOT = Path(__file__).resolve().parents[3]


# Preserve original densities and errors while converting degree coordinates to radians
@pytest.mark.parametrize('record,table', [(397423,5),(442537,9),(452353,10),(582237,6),(582237,7)])
def test_azimuth_units(record, table):
    path = ROOT / f'HEPData/PHOTOPROD/HEPData-ins{record}-v1-json/Table{table}.json'
    original = read_table(path)
    data = hera_reader.read(str(path))
    np.testing.assert_allclose(data['y'], [float(row['y'][0]['value']) for row in original['values']])
    expected_error = np.array([float(row['y'][0]['errors'][0]['symerror']) for row in original['values']])
    np.testing.assert_allclose(data['y_err'], expected_error)
    if 'DEGREES' in original['headers'][0]['name']:
        bounds = [[float(row['x'][0][key]) for key in ('low','high')] for row in original['values']]
        np.testing.assert_allclose(data['binedges'], np.deg2rad(bounds))
    else:
        np.testing.assert_allclose(data['bins'], np.linspace(0.0,2*np.pi,len(original['values'])+1))
    assert np.sum(data['binwidth']) == pytest.approx(2*np.pi)


# Use the original row flux, transfer width and correlated H1 statistical submatrix
@pytest.mark.parametrize('observable', ['M', 'M_t'])
def test_h1_mass_slices(observable):
    card, _ = steering.load_dataset(str(ROOT/'icepack/PHOTOPROD/H1_1798511/dataset.json'),cdir=ROOT)
    for entry in [s for s in card['sets'] if 'rows' in s['hist'][0] and s['hist'][0]['obs'] == observable]:
        hist = entry['hist'][0]
        path = ROOT/card['datapath']/hist['file']
        table = read_table(path)
        rows = [table['values'][i] for i in hist['rows']]
        headers = [item['name'] for item in table['headers'][:table['x_count']]]
        data = h1_reader.read(str(path), file_filter='elastic', hist=hist)
        scale = 1.0 / hera_reader.photon_flux(table,rows)
        if observable == 'M_t':
            widths = [float(row['x'][0]['high'])-float(row['x'][0]['low']) for row in rows]
            scale /= widths
        np.testing.assert_allclose(data['mc_scale'],scale)
        np.testing.assert_allclose(data['y'], [float(row['y'][0]['value']) for row in rows])
        stat = np.array([float(row['y'][0]['errors'][0]['symerror']) for row in rows])
        correlation = np.eye(len(rows))
        if 'covariance' in hist:
            column = next(i for i,name in enumerate(headers) if name.startswith('globalBinNumber'))
            identifiers = [int(row['x'][column]['value']) for row in rows]
            lookup = {tuple(int(item['value']) for item in row['x']):float(row['y'][0]['value'])
                      for row in read_table(path.parent/hist['covariance'])['values']}
            correlation = np.array([[lookup.get((i,j),lookup.get((j,i))) for j in identifiers] for i in identifiers])
        np.testing.assert_allclose(source_covariance(data['uncertainties'][0]), correlation*np.outer(stat,stat))


# Check the numerical transverse flux against the analytic Q2 integral for the actual H1 acceptance
def test_transverse_flux():
    from icepack.PHOTOPROD.H1_415281.cuts import cut_param

    card = ROOT/'icepack/PHOTOPROD/H1_415281/gencard.json'
    import pyjson5
    energy = pyjson5.decode(card.read_text())['SCATTERING']['ENERGY']
    mass = next(float(line.split()[1]) for line in (ROOT/'modeldata/mass_width_2026.mcd').read_text().splitlines()
                if line.split() and line.split()[0]=='11')
    ymin,ymax = np.array(cut_param['W'])**2/(4*np.prod(energy))

    # Integrate the closed-form transverse photon spectrum over the measured W range
    def spectrum(y):
        qmin = mass**2*y*y/(1-y)
        qmax = cut_param['Q2_MAX']
        return ((1+(1-y)**2)*np.log(qmax/qmin)-2*(1-y)*(1-qmin/qmax))/y

    expected = alpha/(2*np.pi)*quad(spectrum,ymin,ymax)[0]
    assert rho_flux(str(card),cut_param,vmd=False) == pytest.approx(expected,rel=1e-10)
    assert 0.0 < rho_flux(str(card),cut_param) < expected


# Check helicity angles against recoil-axis projections after a common rotation and boost
@pytest.mark.parametrize('transform', [False,True])
def test_helicity_covariance(transform):
    event,_ = photo_event(1,False,0.938272)
    particles = event.evt.particles()
    for particle in particles:
        if particle.pid()==90210:
            particle.set_pid(2212)
    parent = obs.proj_central_system(event)
    mass = next(p.momentum().m() for p in particles if p.pid()==211)
    momentum = np.sqrt((parent.m/2)**2-mass**2)
    direction = np.array([0.3,0.4,np.sqrt(0.75)])
    for particle in particles:
        if abs(particle.pid())==211:
            sign = 1 if particle.pid()>0 else -1
            p = vec4(*(sign*momentum*direction),parent.m/2)
            p.boost(b=parent,sign=1)
            particle.set_momentum(type(particle.momentum())(p.x,p.y,p.z,p.t))
    event = SimpleNamespace(evt=event.evt,pid=event.pid,cut_param=event.cut_param)
    reference = (angles.costheta(event),angles.phi(event))
    if transform:
        rotation = Rotation.from_rotvec([0.4,-0.8,0.2]).as_matrix()
        boost = vec4(0.2,-0.1,0.3,np.sqrt(1.14))
        for particle in particles:
            p = particle.momentum()
            xyz = rotation@np.array([p.px(),p.py(),p.pz()])
            p4 = vec4(*xyz,p.e())
            p4.boost(b=boost,sign=1)
            particle.set_momentum(type(particle.momentum())(p4.x,p4.y,p4.z,p4.t))
        event = SimpleNamespace(evt=event.evt,pid=event.pid,cut_param=event.cut_param)
    np.testing.assert_allclose([angles.costheta(event),angles.phi(event)],reference,atol=1e-10)
    initial = exclusive_jpsi.unique_momentum(event,2212,True)
    photon = exclusive_jpsi.unique_momentum(event,11,True)-exclusive_jpsi.unique_momentum(event,11,False)
    daughter = next(r['p4'] for r in obs.proj_central_particle_records(event) if r['pid']==211)
    beams = obs.LorentFramePrepare(initial,photon,[daughter],obs.proj_central_system(event))
    recoil = beams[0]+beams[1]
    expected = -np.dot(beams[2][0].p3,recoil.p3)/(np.linalg.norm(beams[2][0].p3)*np.linalg.norm(recoil.p3))
    assert angles.costheta(event)==pytest.approx(expected,abs=1e-10)


# Keep different fiducial regions independent when they reuse one reconstructed event
@pytest.mark.parametrize('field', ['W','M','ABS_T_MIN','ABS_T_MAX'])
def test_shared_cut_regions(field):
    event,_ = photo_event(1,False,0.938272)
    for particle in event.evt.particles():
        if particle.pid()==90210:
            particle.set_pid(2212)
    event = SimpleNamespace(evt=event.evt,pid=event.pid,cut_param={})
    values = exclusive_jpsi.kinematics(event)
    mass = obs.proj_central_system(event).m
    cuts = {'W':[0.5*values['w'],2*values['w']], 'M':[0.5*mass,2*mass],
            'Q2_MAX':values['q2']+1.0, 'ABS_T_MAX':2*values['abs_t']}
    rejected = {**cuts, field:{'W':[1.1*values['w'],2*values['w']], 'M':[1.1*mass,2*mass],
                              'ABS_T_MIN':1.1*values['abs_t'], 'ABS_T_MAX':0.5*values['abs_t']}[field]}
    assert exclusive_jpsi.accepted(event,cuts)
    assert not exclusive_jpsi.accepted(event,rejected)
    assert exclusive_jpsi.accepted(event,cuts)
    separate = SimpleNamespace(evt=event.evt,pid=event.pid,cut_param=rejected)
    assert not exclusive_jpsi.accepted(separate)
    assert not event.cut_param
