# MC observables constructions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numba
import numpy as np

from core.analysis import pdg
from core.io.cache import cache
from core.kinematics.vec4 import hepmc2vec4, vec4


# Normalize a spatial direction independently of momentum units
@numba.njit
def unitvec(v):
    scale = np.max(np.abs(v))
    if scale <= 0.0:
        return np.zeros_like(v)
    scaled = v / scale
    return scaled / np.linalg.norm(scaled)


# Compile one representative kernel before any expensive simulator work starts
def preflight_numba_runtime() -> None:
    try:
        # Load the lazy typed registries before a Ray trial forks histogram workers
        from numba.typed import Dict, dictimpl, typeddict

        from core.stats.uncertainty import stable_l2_norm

        typed_runtime = (Dict, dictimpl, typeddict)
        result = stable_l2_norm(np.asarray([3.0, 4.0], dtype=np.float64))
    except TimeoutError:
        raise
    except Exception as exc:
        raise RuntimeError(f"Numba runtime preflight failed: {exc}") from exc
    if any(symbol is None for symbol in typed_runtime):
        raise RuntimeError("Numba runtime preflight did not load the typed registries")
    if not np.isclose(result, 5.0, rtol=1e-12, atol=1e-12):
        raise RuntimeError(f"Numba runtime preflight returned an invalid result: {result!r}")


# Convert radians to degrees without invoking the Numba compiler
def rad2deg(x):
    return x / (2 * np.pi) * 360


# Convert degrees to radians without invoking the Numba compiler
def deg2rad(x):
    return x / 360 * (2 * np.pi)
@cache
def proj_event_particles(event):
    records = []
    for index, p in enumerate(event.evt.particles()):
        mom = p.momentum()
        is_final = p.status() == pdg.FINAL_STATE and p.end_vertex() is None
        records.append(
            {
                "id": p.id() if hasattr(p, "id") else index,
                "pid": p.pid(),
                "status": p.status(),
                "is_final": is_final,
                "eta": mom.eta(),
                "pt": mom.pt(),
                "p4": hepmc2vec4(mom),
            }
        )
    return records


# Flatten a dataset particle identity selection while retaining repeated entries
def flatten_pid_selection(selection):
    if isinstance(selection, (list, tuple)):
        return [pid_value for item in selection for pid_value in flatten_pid_selection(item)]
    return [int(selection)]


# Compute every particle at or below a set of HepMC roots
def _descendant_particle_ids(roots):
    descendants = set()
    pending = list(roots)
    while pending:
        particle = pending.pop()
        if particle.id() in descendants:
            continue
        descendants.add(particle.id())
        vertex = particle.end_vertex()
        if vertex is not None:
            pending.extend(vertex.particles_out())
    return descendants


# Compute particle identities below the common central production record
@cache
def _central_particle_ids(event):
    if not hasattr(event.evt, "vertices"):
        return {record["id"] for record in proj_event_particles(event)}

    central_roots = []
    for vertex in event.evt.vertices():
        incoming = list(vertex.particles_in())
        propagators = {
            particle.id()
            for particle in incoming
            if particle.pid() == pdg.PDG_PROPAGATOR and particle.status() == pdg.INTERMEDIATE_STATE
        }
        if len(propagators) >= 2:
            central_roots.extend(vertex.particles_out())
            continue

        beams = {particle.id() for particle in incoming if particle.status() == pdg.INITIAL_STATE}
        if len(beams) >= 2:
            central_roots.extend(particle for particle in vertex.particles_out() if particle.pid() == pdg.PDG_SYSTEM)

    if central_roots:
        return _descendant_particle_ids(central_roots)

    # Preserve support for external records without the Graniitti central tags
    beam_branch_roots = []
    for beam in event.evt.particles():
        if beam.status() != pdg.INITIAL_STATE or beam.end_vertex() is None:
            continue
        beam_branch_roots.extend(
            particle
            for particle in beam.end_vertex().particles_out()
            if particle.pid() == beam.pid() or abs(particle.pid()) == pdg.PDG_NSTAR
        )
    beam_branch_ids = _descendant_particle_ids(beam_branch_roots)
    return {particle.id() for particle in event.evt.particles()} - beam_branch_ids


# Compute topology selected stable central particle records
@cache
def proj_central_particle_records(event):
    requested = set(flatten_pid_selection(event.pid))
    central_ids = _central_particle_ids(event)
    return [
        record
        for record in proj_event_particles(event)
        if record["id"] in central_ids and record["pid"] in requested and record["is_final"]
    ]


## Common projectors
#
# N.B: cache functions needs to be called as follows:
# e.g. proj_central_particles(event), but not proj_central_particles(event=event)
#
#
@cache
def proj_init_final_protons(event):
    beams = [
        particle
        for particle in event.evt.particles()
        if abs(particle.pid()) == pdg.PDG_PROTON and particle.status() == pdg.INITIAL_STATE
    ]
    if len(beams) != 2:
        raise ValueError(f"proj_init_final_protons: expected two proton beams, found {len(beams)}")
    beams.sort(key=lambda particle: particle.momentum().pz(), reverse=True)

    final_protons = []
    for beam in beams:
        vertex = beam.end_vertex()
        daughters = (
            []
            if vertex is None
            else [
                particle
                for particle in vertex.particles_out()
                if particle.pid() == beam.pid()
                and particle.status() == pdg.FINAL_STATE
                and particle.end_vertex() is None
            ]
        )
        if len(daughters) != 1:
            raise ValueError(
                f"proj_init_final_protons: expected one elastic proton daughter for beam PID {beam.pid()}, found {len(daughters)}"
            )
        final_protons.append(daughters[0])

    initial = [hepmc2vec4(particle.momentum()) for particle in beams]
    final = [hepmc2vec4(particle.momentum()) for particle in final_protons]
    return initial, final


# Compute one named HepMC attribute as text
def _attribute_text(owner, name):
    attribute_names = getattr(owner, "attribute_names", None)
    if not callable(attribute_names) or name not in set(attribute_names()):
        return None
    attribute_as_string = getattr(owner, "attribute_as_string", None)
    if callable(attribute_as_string):
        return str(attribute_as_string(name))
    attribute = owner.attribute(name)
    value = attribute.value() if callable(getattr(attribute, "value", None)) else attribute
    return str(value)


# Compute one named HepMC event attribute as text
def event_attribute_text(event, name):
    return _attribute_text(event.evt, name)


# Compute one named HepMC particle attribute as text
def _particle_attribute_text(particle, name):
    return _attribute_text(particle, name)


# Compute UPC forward transfers ordered by their explicit beam-leg tags
@cache
def proj_upc_forward_transfers(event):
    records = [None, None]
    beams = [particle for particle in event.evt.particles() if particle.status() == pdg.INITIAL_STATE]
    if len(beams) != 2:
        raise ValueError(f"proj_upc_forward_transfers: expected two beams, found {len(beams)}")

    for beam in beams:
        vertex = beam.end_vertex()
        daughters = [] if vertex is None else list(vertex.particles_out())
        tagged = []
        for daughter in daughters:
            leg_text = _particle_attribute_text(daughter, "graniitti_upc_leg")
            if leg_text is None:
                continue
            try:
                leg = int(leg_text)
            except ValueError as exc:
                raise ValueError("proj_upc_forward_transfers: invalid UPC leg tag") from exc
            if leg not in (1, 2):
                raise ValueError("proj_upc_forward_transfers: UPC leg tag must be 1 or 2")
            tagged.append((leg, daughter))
        if len(tagged) != 1:
            raise ValueError(
                f"proj_upc_forward_transfers: expected one tagged forward daughter per beam, found {len(tagged)}"
            )

        leg, daughter = tagged[0]
        index = leg - 1
        if records[index] is not None:
            raise ValueError(f"proj_upc_forward_transfers: duplicate UPC leg {leg}")
        sector = _particle_attribute_text(daughter, "graniitti_upc_final_sector")
        if sector not in {"coherent", "incoherent"}:
            raise ValueError("proj_upc_forward_transfers: invalid UPC final-sector tag")
        initial = hepmc2vec4(beam.momentum())
        final = hepmc2vec4(daughter.momentum())
        records[index] = {
            "leg": leg,
            "sector": sector,
            "pid": daughter.pid(),
            "initial": initial,
            "final": final,
            "t": (initial - final).m2,
        }

    if any(record is None for record in records):
        raise ValueError("proj_upc_forward_transfers: incomplete UPC leg tags")
    return tuple(records)


# Compute the exact target transfer for a one-ion-incoherent UPC final state
@cache
def proj_upc_incoherent_abs_t(event):
    selected = [record for record in proj_upc_forward_transfers(event) if record["sector"] == "incoherent"]
    if len(selected) != 1:
        raise ValueError(f"proj_upc_incoherent_abs_t: expected one incoherent ion, found {len(selected)}")
    return abs(selected[0]["t"])


@cache
def proj_central_particles(event):
    """Stable particles matching the requested central-system PDG identities"""
    return [record["p4"] for record in proj_central_particle_records(event)]


# Compute the full central system four-momentum
@cache
def proj_central_system(event):
    return sum(proj_central_particles(event), vec4())


# Compute the ordered signed proton momentum transfers
@cache
def proj_t_pair(event):
    initial, final = proj_init_final_protons(event)
    return tuple((initial[i] - final[i]).m2 for i in range(2))


# Transform the selected central particles to the Collins Soper frame
@cache
def proj_1D_CS(event, daughter=None):
    # Central system 4-momentum
    records = proj_central_particle_records(event)
    requested = flatten_pid_selection(event.pid if daughter is None else daughter)
    identical = len(records) == len(requested) and len(set(requested)) == 1
    if daughter is not None:
        if not requested or len(set(requested)) != len(requested) or any(
            sum(r["pid"] == pid for r in records) != 1 for pid in requested
        ):
            raise ValueError("proj_1D_CS: Daughter PIDs must identify unique central particles")
        particles = [sum((r["p4"] for r in records if r["pid"] in requested), vec4())]
    elif not requested or (not identical and sum(r["pid"] == requested[0] for r in records) != 1):
        raise ValueError("proj_1D_CS: First requested PID must identify one central particle")
    else:
        # Use the HepMC particle order for identical daughters without a momentum dependent selection
        particles = [r["p4"] for r in sorted(records, key=lambda r: (r["pid"] != requested[0], r["id"]))]
    X = proj_central_system(event)

    # Beams
    beam = [r["p4"] for r in proj_event_particles(event) if r["status"] == pdg.INITIAL_STATE]

    if len(beam) != 2:
        raise Exception("proj_1D_CS: Did not find two initial states")

    # Do the frame transform
    pb1boost, pb2boost, pfboost = LorentFramePrepare(pbeam1=beam[0], pbeam2=beam[1], particles=particles, X=X)
    pfout = LorentzFrame(pb1boost=pb1boost, pb2boost=pb2boost, pfboost=pfboost, frametype="CS")

    return pfout


# Pair opposite-sign daughters using a configurable resonance mass
@cache
def proj_1D_4body_helicity(event, target_m=0.775):
    """Sequential 2->2 body projector, e.g. resonance -> rho > {pi+ pi-} rho > {pi+ pi-}"""

    records = proj_central_particle_records(event)

    # Pions (muons)
    pi_pos_p4 = [r["p4"] for r in records if np.sign(r["pid"]) == 1]

    if len(pi_pos_p4) != 2:
        raise Exception(
            "proj_1D_4body_helicity: Did not find two (+) central particles -- check generator and python cuts"
        )

    # Anti-pions (anti-muons)
    pi_neg_p4 = [r["p4"] for r in records if np.sign(r["pid"]) == -1]

    if len(pi_neg_p4) != 2:
        raise Exception(
            "proj_1D_4body_helicity: Did not find two (-) central particles -- check generator and python cuts"
        )
    # Find the closest combination (two combinations)
    Pe = [
        [
            [0, 0],  # Combination 0: system A
            [1, 1],
        ],  # Combination 0: system B
        [
            [0, 1],  # Combination 1: system A
            [1, 0],
        ],
    ]  # Combination 1: system B

    loss = np.zeros(len(Pe))

    A, B = 0, 1
    for n in range(len(Pe)):
        ind_A, ind_B = Pe[n][A], Pe[n][B]

        mA = (pi_pos_p4[ind_A[0]] + pi_neg_p4[ind_A[1]]).m
        mB = (pi_pos_p4[ind_B[0]] + pi_neg_p4[ind_B[1]]).m
        loss[n] = (target_m - mA) ** 2 + (target_m - mB) ** 2

    # Pick the better solution
    n = np.argmin(loss)

    ind_A, ind_B = Pe[n][A], Pe[n][B]

    particles = [[pi_pos_p4[ind_A[0]], pi_neg_p4[ind_A[1]]], [pi_pos_p4[ind_B[0]], pi_neg_p4[ind_B[1]]]]
    X = particles[0][0] + particles[0][1] + particles[1][0] + particles[1][1]

    # Over both intermediate systems
    pions_in_HX = []
    system_in_lab = []
    for i in range(len(particles)):
        # Intermediate system Xi in the lab
        Xi_in_lab = particles[i][0] + particles[i][1]
        system_in_lab.append(Xi_in_lab)

        # Intermediate system Xi in the central system resonance X rest frame
        Xi_in_X_frame = Xi_in_lab.copy()
        Xi_in_X_frame.boost(b=X, sign=-1)

        # Put daughters in the same grandmother frame as the helicity axis
        particles_in_X_frame = [q.copy() for q in particles[i]]
        for q in particles_in_X_frame:
            q.boost(b=X, sign=-1)

        # Transform to the intermediate mother Helicity Frame
        pfout = HXFrame(p=particles_in_X_frame, X=Xi_in_X_frame)

        # Get the pions
        pions_in_HX.append(pfout)

    return pions_in_HX, system_in_lab, particles


## Observables
#
#


# Compute the paired resonance mass difference
@cache
def proj_1D_4body_DeltaM_AB(event, target_m=0.775):
    _, system_in_lab, _ = proj_1D_4body_helicity(event, target_m=target_m)
    return system_in_lab[0].m - system_in_lab[1].m


# Compute the first reconstructed resonance mass
@cache
def proj_1D_4body_M_A(event, target_m=0.775):
    """X -> A > {A1 + A2} B > {B1 + B2}"""
    _, system_in_lab, _ = proj_1D_4body_helicity(event, target_m=target_m)
    return system_in_lab[0].m


# Compute the second reconstructed resonance mass
@cache
def proj_1D_4body_M_B(event, target_m=0.775):
    """X -> A > {A1 + A2} B > {B1 + B2}"""
    _, system_in_lab, _ = proj_1D_4body_helicity(event, target_m=target_m)
    return system_in_lab[1].m


# Compute the helicity cosine difference for the selected pairs
@cache
def proj_1D_4body_Deltacos_12(event, target_m=0.775):
    """X -> A > {A1 + A2} B > {B1 + B2}"""
    pions_in_HX, _, _ = proj_1D_4body_helicity(event, target_m=target_m)
    return pions_in_HX[0][0].costheta - pions_in_HX[1][0].costheta


# Compute the first helicity cosine using a spatial projection
@cache
def proj_1D_4body_cos1_dotprod(event, target_m=0.775):
    """X -> A > {A1 + A2} B > {B1 + B2}"""
    pions_in_HX, _, _ = proj_1D_4body_helicity(event, target_m=target_m)

    A, A1 = 0, 0
    z_axis = vec4(0.0, 0.0, 1.0, 0.0)

    # Take the dot product in the same intermediate helicity frame
    costheta = np.cos(pions_in_HX[A][A1].angle(z_axis))

    return costheta


# Compute the first paired daughter helicity cosine
@cache
def proj_1D_4body_cos1(event, target_m=0.775):
    """X -> A > {A1 + A2} B > {B1 + B2}"""
    pions_in_HX, _, _ = proj_1D_4body_helicity(event, target_m=target_m)
    return pions_in_HX[0][0].costheta


# Compute the second paired daughter helicity cosine
@cache
def proj_1D_4body_cos2(event, target_m=0.775):
    """X -> A > {A1 + A2} B > {B1 + B2}"""
    pions_in_HX, _, _ = proj_1D_4body_helicity(event, target_m=target_m)
    return pions_in_HX[1][0].costheta


# Compute the sum of the paired daughter helicity azimuths
@cache
def proj_1D_4body_phi12(event, target_m=0.775):
    """X -> A > {A1 + A2} B > {B1 + B2}"""
    pions_in_HX, _, _ = proj_1D_4body_helicity(event, target_m=target_m)
    angle = pions_in_HX[0][0].phi + pions_in_HX[1][0].phi

    return rad2deg(np.arctan2(np.sin(angle), np.cos(angle)))  # Wrap to [-pi,pi]


# Compute the central system invariant mass
@cache
def proj_1D_M(event):
    return proj_central_system(event).m


# Compute the central system rapidity
@cache
def proj_1D_Rap(event):
    return proj_central_system(event).rapidity


# Compute the central system transverse momentum
@cache
def proj_1D_Pt(event):
    return proj_central_system(event).pt


# Compute the absolute momentum transfer on the positive longitudinal momentum side
@cache
def proj_1D_t1(event):
    return np.abs(proj_t_pair(event)[0])


# Compute the absolute momentum transfer on the negative longitudinal momentum side
@cache
def proj_1D_t2(event):
    return np.abs(proj_t_pair(event)[1])


# Compute the absolute sum of the signed proton momentum transfers
@cache
def proj_1D_Abs_t1t2(event):
    t1, t2 = proj_t_pair(event)
    return np.abs(t1 + t2)


@cache
def proj_1D_dPhi_pp(event):
    """Forward proton pair deltaphi in the lab frame (in deg [0,180])"""
    init_proton, final_proton = proj_init_final_protons(event)

    return rad2deg(final_proton[0].abs_delta_phi(final_proton[1]))


# Compute a selected daughter polar angle in the Collins Soper frame
@cache
def proj_1D_costheta_CS(event, daughter=None):
    """Compute cos(theta) of a selected daughter in the Collins Soper frame [-1,1]"""
    pfout = proj_1D_CS(event, daughter=daughter)

    return pfout[0].costheta


# Compute a selected daughter azimuth in the Collins Soper frame
@cache
def proj_1D_phi_CS(event, in_deg=True, daughter=None):
    """Compute azimuth of a selected daughter in the Collins Soper frame [-180,180] degrees"""
    pfout = proj_1D_CS(event, daughter=daughter)

    # Phi angle in the new rest frame
    phi = pfout[0].phi

    return rad2deg(phi) if in_deg else phi


@cache
def proj_1D_mandelstam_t(event):
    """Mandelstam |t| for elastic scattering"""
    init_proton, final_proton = proj_init_final_protons(event)

    t = (init_proton[0] - final_proton[0]).m2

    return np.abs(t)


@cache
def proj_3D_fp1pt_fp2pt_max_hat_tu(event):
    r"""Forward pt1, pt2 and central sub invariant max(\hat{t}, \hat{u})"""
    init_proton, final_proton = proj_init_final_protons(event)

    # Final forward pt
    pt1 = final_proton[0].pt
    pt2 = final_proton[1].pt

    # Propagator vectors
    q1 = init_proton[0] - final_proton[0]
    # q2 = init_proton[1] - final_proton[1]

    # Central particles
    p = proj_central_particles(event)

    # Central sub-Mandelstam
    that = (q1 - p[0]).m2  # note q1 on both
    uhat = (q1 - p[1]).m2  # in t and u!
    max_hat_tu = max(that, uhat)

    return np.array([pt1, pt2, max_hat_tu])


@cache
def proj_3D_fp1pt_fp2pt_M(event):
    r"""Forward pt1, pt2 and central sub invariant \hat{s}"""
    init_proton, final_proton = proj_init_final_protons(event)

    # Final proton pt
    pt1 = final_proton[0].pt
    pt2 = final_proton[1].pt

    # Central mass
    p = proj_central_particles(event)
    M = sum(p).m

    return np.array([pt1, pt2, M])


@cache
def proj_3D_fp1pt_fp2pt_dPhi_pp(event):
    """Forward pt1, pt2, deltaPhi"""
    init_proton, final_proton = proj_init_final_protons(event)

    # Final proton pt
    pt1 = final_proton[0].pt
    pt2 = final_proton[1].pt

    # DeltaPhi in rad
    dPhi = final_proton[0].abs_delta_phi(final_proton[1])

    return np.array([pt1, pt2, dPhi])
def HXFrame(p, X):
    # ZYZ-sequence rotation angles as defined by the system X direction, note the minus
    Z_angle = -X.phi
    Y_angle = -X.theta

    def rotate(particle):
        particle.rotateZ(Z_angle)
        particle.rotateY(Y_angle)
        particle.rotateZ(np.pi)  # Reflection

    # Rotate particles
    pout = [particle.copy() for particle in p]
    for i in range(len(pout)):
        rotate(pout[i])

    # Boost direction defined as a sum over the rotated particles
    p_boost = vec4()
    for i in range(len(pout)):
        p_boost += pout[i]

    # Boost each particle
    for i in range(len(pout)):
        pout[i].boost(b=p_boost, sign=-1)

    return pout


# Boost the decay and beam vectors into the central system rest frame
def LorentFramePrepare(pbeam1, pbeam2, particles, X):
    """Lorentz Transform preparation function for LorentzFrame()"""

    pfboost = [particle.copy() for particle in particles]
    pbeam1b = pbeam1.copy()
    pbeam2b = pbeam2.copy()

    # Boost each particle to the system X rest frame
    for i in range(len(pfboost)):
        pfboost[i].boost(b=X, sign=-1)

    # Boost the initial state beam vectors
    pbeam1b.boost(b=X, sign=-1)
    pbeam2b.boost(b=X, sign=-1)

    return pbeam1b, pbeam2b, pfboost


# Construct decay axes from the two initial momenta in the central rest frame
def LorentzFrame(pb1boost, pb2boost, pfboost, frametype="CS", direction=1):
    """Rotated Lorentz frame transformation"""

    pb1boost3 = pb1boost.p3
    pb2boost3 = pb2boost.p3

    # Frame rotation x-y-z-axes

    ## NON-ROTATED FRAME AXIS DEFINITION
    if frametype == "CM":
        zaxis = np.array([0.0, 0.0, 1.0], dtype=float)
        yaxis = np.array([0.0, 1.0, 0.0], dtype=float)

    ## COLLINS-SOPER FRAME POLARIZATION AXIS DEFINITION
    elif frametype == "CS" or frametype == "G1":
        zaxis = unitvec(unitvec(pb1boost3) - unitvec(pb2boost3))

    ## ANTI-HELICITY FRAME POLARIZATION AXIS DEFINITION
    elif frametype == "AH" or frametype == "G2":
        zaxis = unitvec(unitvec(pb1boost3) + unitvec(pb2boost3))

    ## HELICITY FRAME POLARIZATION AXIS DEFINITION
    elif frametype == "HX" or frametype == "G3":
        # The negative total incoming momentum points along X in the incoming CM frame
        zaxis = unitvec(-(pb1boost3 + pb2boost3))

    ## PSEUDO-GOTTFRIED-JACKSON AXIS DEFINITION: [1] or [2]
    elif frametype == "PG" or frametype == "G4":
        if direction == -1:
            zaxis = unitvec(pb1boost3)
        elif direction == 1:
            zaxis = unitvec(pb2boost3)
        else:
            raise Exception(f"LorentzFrame: Invalid direction {direction}")
    else:
        raise Exception(f"LorentzFrame: Unknown frame <{frametype}>")

    # y-axis
    if frametype != "CM":
        normal = np.cross(unitvec(pb1boost3), unitvec(pb2boost3))
        yaxis = unitvec(normal) if np.linalg.norm(normal) > 1e-12 else np.zeros(3)

    if np.linalg.norm(zaxis) < 1e-12:
        raise ValueError(f"LorentzFrame: Undefined polar axis for frame <{frametype}>")
    if np.linalg.norm(yaxis) < 1e-12:
        # Collinear beams leave azimuth arbitrary but the physical polar axis remains fixed
        reference = np.eye(3)[np.argmin(np.abs(zaxis))]
        yaxis = unitvec(np.cross(zaxis, reference))

    # x-axis
    xaxis = unitvec(np.cross(yaxis, zaxis))  # x = y [cross product] z

    # Create SO(3) rotation matrix for the new coordinate axes
    R = np.array([xaxis, yaxis, zaxis])  # Axes as rows

    # Rotate all vectors
    for k in range(len(pfboost)):
        pfboost[k].rotateSO3(R)

    return pfboost
