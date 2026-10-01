# Checks for the CMS 2752118 observable specific event selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path

from core.io import steering
from core.kinematics.vec4 import vec4

from tests.technical.support.hepmc import decay_event, make_pair_event

ROOT = Path(__file__).resolve().parents[3]
CUT_ROOT = ROOT / "icepack" / "SOFTCEP" / "CMS_2752118"


# Load one CMS event selection from its source file
def load_cuts(relative_path):
    return steering.load_python_module(
        str(CUT_ROOT / relative_path),
    )


COMMON_CUTS = load_cuts("_common/cuts.py")
PHI_CUTS = load_cuts("pipi_0p7/cuts_0p35--0p65.py")
FULL_MASS_CUTS = load_cuts("pipi/cuts.py")
LOW_MASS_CUTS = load_cuts("pipi_0p7/cuts.py")
HIGH_MASS_CUTS = load_cuts("pipi_1p8--2p2/cuts.py")


# Construct exclusive pp pion-pair events with measured forward proton kinematics
def make_event(dphi=1.0, mass=0.5, **kwargs):
    return make_pair_event(pid=211, mass=mass, dphi=math.degrees(dphi), forward_pt=0.775, beam_energy=6500.0, **kwargs)


# Apply missing Delta phi bins only to the Delta phi measurement
def test_dphi_hole_excludes_mass_measurements():
    event = make_event()

    assert not PHI_CUTS.cut_func(event)
    assert FULL_MASS_CUTS.cut_func(event)
    assert LOW_MASS_CUTS.cut_func(event)

    high_mass_event = make_event(mass=2.0)
    assert FULL_MASS_CUTS.cut_func(high_mass_event)
    assert HIGH_MASS_CUTS.cut_func(high_mass_event)


# Keep the published pion pair mass interval on the Delta phi measurement
def test_dphi_selection_keeps_mass_window():

    assert PHI_CUTS.cut_func(make_event(dphi=2.0, mass=0.5))
    assert not PHI_CUTS.cut_func(make_event(dphi=2.0, mass=0.34))
    assert not PHI_CUTS.cut_func(make_event(dphi=2.0, mass=0.66))


# Use exact ten degree bin boundaries for the Delta phi veto
def test_dphi_gap_uses_exact_angle_edges():
    degree = math.pi / 180.0

    assert COMMON_CUTS.dPhi_cuts(0.775, 0.775, 30.0 * degree)
    assert not COMMON_CUTS.dPhi_cuts(0.775, 0.775, 89.9 * degree)
    assert COMMON_CUTS.dPhi_cuts(0.775, 0.775, 90.0 * degree)
    assert COMMON_CUTS.dPhi_cuts(0.775, 0.775, 90.1 * degree)


# Reject events outside the shared pion and proton fiducial selection
def test_mass_fiducial_selection():

    assert not FULL_MASS_CUTS.cut_func(make_pair_event(pid=211, mass=0.5, forward_pt=0.2, beam_energy=6500.0))
    pions = [(pid, vec4(0.0, 0.0, 0.0, 0.13957039), 1) for pid in (211, -211, 211)]
    assert not FULL_MASS_CUTS.cut_func(decay_event(pions))
    assert not FULL_MASS_CUTS.cut_func(make_event(dphi=math.pi, mass=2.0, costheta=0.999))
