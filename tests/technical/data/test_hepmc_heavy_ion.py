# Preserve nuclear event information through HepMC copies and files
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from pyHepMC3 import HepMC3 as hepmc


# Keep impact parameter and both spectator populations across serialization
@pytest.mark.parametrize("transfer", ["copy", "ascii"])
def test_heavy_ion_roundtrip(transfer, tmp_path):
    event = hepmc.GenEvent()
    vertex = hepmc.GenVertex()
    for pid, pz, mass in ((1000822080, 500000.0, 193.7), (1000791970, -470000.0, 183.5)):
        momentum = hepmc.FourVector(0, 0, pz, math.hypot(pz, mass))
        vertex.add_particle_in(hepmc.GenParticle(momentum, pid, 4))
        vertex.add_particle_out(hepmc.GenParticle(momentum, pid, 1))
    event.add_vertex(vertex)
    heavy = hepmc.GenHeavyIon()
    expected = {
        "impact_parameter": 23.4567,
        "sigma_inel_NN": 70.1234,
        "Ncoll": 0,
        "Npart_proj": 0,
        "Npart_targ": 0,
        "Nspec_proj_n": 126,
        "Nspec_proj_p": 82,
        "Nspec_targ_n": 118,
        "Nspec_targ_p": 79,
    }
    for name, value in expected.items():
        setattr(heavy, name, value)
    event.set_heavy_ion(heavy)
    if transfer == "copy":
        result = hepmc.GenEvent(event)
    else:
        path = str(tmp_path / "nuclear.hepmc3")
        writer = hepmc.WriterAscii(path)
        writer.write_event(event)
        writer.close()
        reader = hepmc.ReaderAscii(path)
        result = hepmc.GenEvent()
        reader.read_event(result)
        assert not reader.failed()
        reader.close()
    assert result.heavy_ion() is not None
    for name, value in expected.items():
        assert getattr(result.heavy_ion(), name) == pytest.approx(value)
