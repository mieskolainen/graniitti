# Forward particle acceptance on HepMC3 events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from collections.abc import Callable
from dataclasses import dataclass
from math import atan2, hypot, isfinite, isnan, pi

from pyHepMC3 import HepMC3 as hepmc3


# Check a closed interval, allowing an azimuth interval to cross pi
def _inside(value, interval, azimuth=False):
    if interval is None:
        return True
    if isnan(value):
        return False
    low, high = interval
    if azimuth and low > high:
        return value >= low or value <= high
    return low <= value <= high


# Validate a physical acceptance interval at construction
def _interval(value, name, azimuth=False):
    if value is None:
        return
    if len(value) != 2 or any(isnan(x) for x in value):
        raise ValueError(f"Invalid {name} interval")
    if not azimuth and value[0] > value[1]:
        raise ValueError(f"Reversed {name} interval")
    if azimuth and any(not -pi <= x <= pi for x in value):
        raise ValueError("Azimuth limits must lie in [-pi, pi]")


# Compute GeV and mm conversion factors without changing the event
def _units(event):
    energy = {hepmc3.Units.GEV: 1.0, hepmc3.Units.MEV: 0.001}
    length = {hepmc3.Units.MM: 1.0, hepmc3.Units.CM: 10.0}
    return energy[event.momentum_unit()], length[event.length_unit()]


# Require explicit physical particles before interpreting zero multiplicity
def _resolved(event):
    for particle in event.particles():
        if particle.end_vertex() is not None:
            continue
        code = abs(particle.pid())
        if code == 91 or (code >= 1000000000 and particle.status() != 1):
            raise ValueError("Forward nuclear final state is unresolved")


# Compute beam ancestry from HepMC3 vertices, with beam 1 along positive z
def _beam_ancestry(event):
    beams = sorted(event.beams(), key=lambda p: p.momentum().pz(), reverse=True)
    if len(beams) != 2 or not beams[0].momentum().pz() > beams[1].momentum().pz():
        raise ValueError("Beam ancestry requires two distinct longitudinal beams")
    ancestry = {particle.id(): frozenset((i,)) for i, particle in enumerate(beams, 1)}
    active = set()

    # Follow incoming particles while detecting malformed cyclic graphs
    def visit(particle):
        key = particle.id()
        if key in ancestry:
            return ancestry[key]
        if key in active:
            raise ValueError("Cyclic HepMC3 ancestry")
        active.add(key)
        vertex = particle.production_vertex()
        parents = () if vertex is None else vertex.particles_in()
        result = frozenset().union(*(visit(parent) for parent in parents))
        active.remove(key)
        ancestry[key] = result
        return result

    for particle in event.particles():
        visit(particle)
    return ancestry


@dataclass(frozen=True)
class Plane:
    """Acceptance at z in mm, with x, y, radius and ct intervals in mm."""

    z: float
    x: tuple | None = None
    y: tuple | None = None
    radius: tuple | None = None
    ct: tuple | None = None
    transport: Callable | None = None

    # Validate detector geometry independently of event sampling
    def __post_init__(self):
        if not isfinite(self.z):
            raise ValueError("Detector z must be finite")
        for name in ("x", "y", "radius", "ct"):
            _interval(getattr(self, name), name)

    # Compute an arrival position, using external beam optics for charged particles
    def position(self, particle, event):
        if self.transport is not None:
            # A transport callable supplies a HepMC3 FourVector (x,y,z,ct) in mm
            return self.transport(particle, self.z, event)
        if abs(particle.pid()) not in (22, 2112):
            raise ValueError("This particle requires a beam transport callable for plane cuts")
        _, length = _units(event)
        vertex = particle.production_vertex()
        if vertex is None:
            raise ValueError("Plane cuts require a particle production vertex")
        start = vertex.position()
        p = particle.momentum()
        if not abs(p.pz()) > 0.0:
            return None
        flight = (self.z - length * start.z()) / p.pz()
        if flight < 0.0:
            return None
        return hepmc3.FourVector(
            length * start.x() + flight * p.px(),
            length * start.y() + flight * p.py(),
            self.z,
            length * start.t() + flight * p.e(),
        )

    # Apply geometric and arrival-time acceptance at the detector plane
    def accepts(self, particle, event):
        hit = self.position(particle, event)
        if hit is None:
            return False
        if not all(isfinite(v) for v in (hit.x(), hit.y(), hit.z(), hit.t())):
            raise ValueError("Non-finite transported position")
        if abs(hit.z() - self.z) > 1e-9 * max(1.0, abs(self.z)):
            raise ValueError("Transport did not reach the requested plane")
        return all(
            _inside(value, interval)
            for value, interval in (
                (hit.x(), self.x),
                (hit.y(), self.y),
                (hypot(hit.x(), hit.y()), self.radius),
                (hit.t(), self.ct),
            )
        )


@dataclass(frozen=True)
class Acceptance:
    """Final particle cuts in GeV and radians, optionally restricted to one beam."""

    pid: tuple = ()
    side: int = 0
    beam: int | None = None
    energy: tuple | None = None
    pt: tuple | None = None
    pz: tuple | None = None
    eta: tuple | None = None
    abs_eta: tuple | None = None
    rapidity: tuple | None = None
    theta: tuple | None = None
    phi: tuple | None = None
    plane: Plane | None = None

    # Validate particle and kinematic selections once per icepack
    def __post_init__(self):
        if self.side not in (-1, 0, 1) or self.beam not in (None, 1, 2):
            raise ValueError("Use side -1, 0 or 1 and beam 1 or 2")
        if any(not isinstance(pid, int) or pid == 0 for pid in self.pid):
            raise ValueError("Particle selections require nonzero PDG integers")
        for name in ("energy", "pt", "pz", "eta", "abs_eta", "rapidity", "theta", "phi"):
            _interval(getattr(self, name), name, azimuth=name == "phi")
        if self.theta is not None and not 0 <= self.theta[0] <= self.theta[1] <= pi:
            raise ValueError("Polar angle limits must lie in [0, pi]")

    # Select actual HepMC3 final particles without changing their momenta or graph
    def particles(self, event):
        record = getattr(event, "evt", event)
        _resolved(record)
        scale, _ = _units(record)
        ancestry = _beam_ancestry(record) if self.beam is not None else None
        selected = []
        for particle in record.particles():
            if particle.status() != 1 or particle.end_vertex() is not None:
                continue
            if self.pid and particle.pid() not in self.pid:
                continue
            if ancestry is not None and ancestry[particle.id()] != frozenset((self.beam,)):
                continue
            p = particle.momentum()
            if not all(isfinite(v) for v in (p.px(), p.py(), p.pz(), p.e())) or p.e() <= 0:
                raise ValueError("Invalid final particle momentum")
            if self.side and self.side * p.pz() <= 0:
                continue
            # Fold theta toward the selected hemisphere, leaving eta and pz signed
            theta = atan2(p.pt(), (-1 if self.side < 0 else 1) * p.pz())
            values = (
                (p.e() * scale, self.energy),
                (p.pt() * scale, self.pt),
                (p.pz() * scale, self.pz),
                (p.eta(), self.eta),
                (abs(p.eta()), self.abs_eta),
                (p.rap(), self.rapidity),
                (theta, self.theta),
            )
            if not all(_inside(value, interval) for value, interval in values):
                continue
            if not _inside(p.phi(), self.phi, azimuth=True):
                continue
            if self.plane is None or self.plane.accepts(particle, record):
                selected.append(particle)
        return tuple(selected)

    # Sum selected final particle energies in GeV
    def energy_sum(self, event):
        record = getattr(event, "evt", event)
        scale, _ = _units(record)
        return scale * sum(p.momentum().e() for p in self.particles(record))


# Classify final beam-remnant neutrons, with the positive-z beam first
def neutron_class(event):
    record = getattr(event, "evt", event)
    neutrons = Acceptance(pid=(2112,)).particles(record)
    ancestry = _beam_ancestry(record)
    emitting = set().union(*(ancestry[p.id()] for p in neutrons if len(ancestry[p.id()]) == 1))
    return "".join("Xn" if beam in emitting else "0n" for beam in (1, 2))
