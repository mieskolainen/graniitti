# Test HepMC3 attributes used by event visualization
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import subprocess
import sys

from core import iceviz
from pyHepMC3 import HepMC3 as h


# Preserve attributes attached to event, particle and vertex owners through the official reader
def test_iceviz_hepmc3_attributes(tmp_path):
    event = h.GenEvent()
    vertex = h.GenVertex()
    particle = h.GenParticle(h.FourVector(1.0, 0.0, 0.0, 1.0), 22, 1)
    vertex.add_particle_in(h.GenParticle(h.FourVector(1.0, 0.0, 0.0, 1.0), 22, 4))
    vertex.add_particle_out(particle)
    event.add_vertex(vertex)
    expected = {
        0: {"label": "event text"},
        particle.id(): {"tag": "particle text"},
        vertex.id(): {"name": "vertex text"},
    }
    for owner, attributes in expected.items():
        for name, text in attributes.items():
            event.add_attribute(name, h.StringAttribute(text), owner)
    path = tmp_path / "attributes.hepmc3"
    writer = h.WriterAscii(str(path))
    writer.write_event(event)
    writer.close()
    loaded = iceviz.read_event_pyhepmc3(h, path, 0)
    assert iceviz.event_attributes(loaded) == expected


# Exercise the alternative official binding in a separate interpreter to avoid duplicate C++ type registration
def test_iceviz_pyhepmc_attributes(tmp_path):
    script = """import pyhepmc
from core import iceviz
from pathlib import Path
path = Path('attributes.hepmc3')
event = pyhepmc.GenEvent()
vertex = pyhepmc.GenVertex()
particle = pyhepmc.GenParticle((1, 0, 0, 1), 22, 1)
vertex.add_particle_in(pyhepmc.GenParticle((1, 0, 0, 1), 22, 4))
vertex.add_particle_out(particle)
event.add_vertex(vertex)
event.attributes['label'] = 'event text'
particle.attributes['tag'] = 'particle text'
vertex.attributes['name'] = 'vertex text'
with pyhepmc.open(path, 'w') as writer:
    writer.write(event)
loaded = iceviz.read_event_pyhepmc(pyhepmc, path, 0)
assert iceviz.event_attributes(loaded) == {0: {'label': 'event text'}, 2: {'tag': 'particle text'}, -1: {'name': 'vertex text'}}
"""
    subprocess.run([sys.executable, "-c", script], cwd=tmp_path, check=True)
