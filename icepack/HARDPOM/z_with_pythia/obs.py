# Observable and histogram definitions for Pythia showered diffractive Z events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numba
import numpy as np
from core.analysis import obs, pdg
from core.io.cache import cache

from . import cuts as DIFFZ_CUTS


# Compute the selected dimuon object prepared by the cut module
def selected_pair(event):
    pair = DIFFZ_CUTS.select_mumu_pair(event)
    if pair is None:
        raise Exception("selected_pair: event did not pass diffz_mumu selection")
    return pair


# Compute the selected dimuon invariant mass
def proj_m_mumu(event):
    return selected_pair(event)["pair"].m


# Compute the selected dimuon transverse momentum
def proj_pt_mumu(event):
    return selected_pair(event)["pair"].pt


# Compute the selected dimuon rapidity
def proj_y_mumu(event):
    return selected_pair(event)["pair"].rapidity


# Compute the selected leading muon transverse momentum
def proj_pt_mu1(event):
    return selected_pair(event)["leading"].pt


# Compute the selected subleading muon transverse momentum
def proj_pt_mu2(event):
    return selected_pair(event)["subleading"].pt


# Compute the selected leading muon pseudorapidity
def proj_eta_mu1(event):
    return selected_pair(event)["leading"].eta


# Compute the selected subleading muon pseudorapidity
def proj_eta_mu2(event):
    return selected_pair(event)["subleading"].eta


# Compute the selected absolute dimuon pseudorapidity separation
def proj_deta_mumu(event):
    pair = selected_pair(event)
    return DIFFZ_CUTS.deta_mumu_from_p4(pair["leading"], pair["subleading"])


# Compute the selected dimuon acoplanarity 1 - |Delta phi| / pi
def proj_acoplanarity(event):
    pair = selected_pair(event)
    return DIFFZ_CUTS.acoplanarity_from_p4(pair["leading"], pair["subleading"])


# --------------------------------------------------------------
# Rapidity gap observables


# Evaluate the sorted boundary integral from per-particle gap penalties
@numba.njit
def sorted_gap_flow_score(
    u: np.ndarray, penalty: np.ndarray, eta_min: float, eta_max: float, beta: float
):
    if u.size == 0:
        return 1.0
    if u.size != penalty.size:
        raise ValueError("sorted_gap_flow_score: u and penalty must have equal size")

    order = np.argsort(u)
    u = u[order]
    penalty = penalty[order]

    tail_penalty = 0.0
    previous_u = eta_max
    numerator = 0.0

    for i in range(u.size - 1, -1, -1):
        numerator += np.exp(-tail_penalty + beta * (previous_u - eta_max)) * (
            -np.expm1(beta * (u[i] - previous_u))
        )
        tail_penalty += penalty[i]
        previous_u = u[i]

    numerator += np.exp(-tail_penalty + beta * (previous_u - eta_max)) * (
        -np.expm1(beta * (eta_min - previous_u))
    )
    normalization = -np.expm1(-beta * (eta_max - eta_min))
    score = numerator / normalization
    return min(1.0, max(0.0, score))


# Evaluate the additive-energy gap flow score on one oriented side
@numba.njit
def energy_gap_flow_score(
    eta: np.ndarray,
    pt: np.ndarray,
    side: int,
    eta_min: float,
    eta_max: float,
    pt_min: float,
    beta: float,
    q0: float,
):
    """
    Energy Gap Flow 'G' or 'cumulant expansion',
    as a new diffraction / rapidity gap observable with pt-dependence.

    This version is by construction in between 0 < G <= 1,
    where values close to 1 means a 'complete detector gap'
    and 0 means essentially no gap ~ non-diffractive.

    Args:
        eta:     pseudorapidity of (charged) particles
        pt:      transverse momentum of (charged) particles (GeV)
        side:    (-1, 1) backward / forward
        eta_min: minimum oriented fiducial eta (e.g. -2.5)
        eta_max: maximum absolute fiducial eta (e.g.  2.5)
        pt_min:  minimum fiducial pt threshold (e.g.  0.2) (GeV)
        beta:    algorithm 'correlation length' parameter (e.g. 0.1)
        q0:      algorithm 'energy scale' parameter (e.g. 2.0) (GeV)

    Combine both hemispheres per event by calling twice:
    max(energy_gap_flow_score(side=-1), energy_gap_flow_score(side=1))
    """

    if side != -1 and side != 1:
        raise ValueError("energy_gap_flow_score: side must be -1 or +1")
    if eta_max <= eta_min:
        raise ValueError("energy_gap_flow_score: eta_max must exceed eta_min")
    if pt_min < 0.0:
        raise ValueError("energy_gap_flow_score: pt_min must be non-negative")
    if beta <= 0.0:
        raise ValueError("energy_gap_flow_score: beta must be positive")
    if q0 <= 0.0:
        raise ValueError("energy_gap_flow_score: q0 (GeV) must be positive")

    oriented = side * eta
    mask = (eta_min < oriented) & (oriented < eta_max) & (pt > pt_min)
    u = oriented[mask]
    penalty = pt[mask] / q0
    return sorted_gap_flow_score(u, penalty, eta_min, eta_max, beta)


# Evaluate the multiplicity-sensitive gap flow score on one oriented side
@numba.njit
def multiplicity_gap_flow_score(
    eta: np.ndarray,
    pt: np.ndarray,
    side: int,
    eta_min: float,
    eta_max: float,
    pt_min: float,
    beta: float,
    p0: float,
):
    r"""
    Multiplicity-sensitive Gap Flow based on the Laplace-void construction.

    Each accepted particle contributes the nonlinear penalty
    $\log(1 + p_T/p_0)$, retaining multiplicity information at fixed total
    transverse momentum while using the same normalized boundary scan.
    """
    if side != -1 and side != 1:
        raise ValueError("multiplicity_gap_flow_score: side must be -1 or +1")
    if eta_max <= eta_min:
        raise ValueError("multiplicity_gap_flow_score: eta_max must exceed eta_min")
    if pt_min < 0.0:
        raise ValueError("multiplicity_gap_flow_score: pt_min must be non-negative")
    if beta <= 0.0:
        raise ValueError("multiplicity_gap_flow_score: beta must be positive")
    if p0 <= 0.0:
        raise ValueError("multiplicity_gap_flow_score: p0 (GeV) must be positive")

    oriented = side * eta
    mask = (eta_min < oriented) & (oriented < eta_max) & (pt > pt_min)
    u = oriented[mask]
    penalty = np.log1p(pt[mask] / p0)
    return sorted_gap_flow_score(u, penalty, eta_min, eta_max, beta)


@numba.njit
def energy_gap_flow_effective_gap(score: float, eta_min: float, eta_max: float, beta: float):
    r"""
    Map the raw Energy Gap Flow (energy_gap_flow_score)
    score to $\Delta\eta^{\mathrm{EGF}}$.

    The transformation is monotonic, maps [0,1] to the configured rapidity
    span, and reduces to the ordinary forward gap in the hard-veto limit.
    """
    if not np.isfinite(score):
        raise ValueError("energy_gap_flow_effective_gap: score must be finite")
    if score < 0.0 or score > 1.0:
        raise ValueError("energy_gap_flow_effective_gap: score must be in [0, 1]")
    if eta_max <= eta_min:
        raise ValueError("energy_gap_flow_effective_gap: eta_max must exceed eta_min")
    if beta <= 0.0:
        raise ValueError("energy_gap_flow_effective_gap: beta must be positive")

    length = eta_max - eta_min
    if score == 0.0:
        return 0.0
    if score == 1.0:
        return length

    normalization = -np.expm1(-beta * length)
    effective_gap = -np.log1p(-normalization * score) / beta
    return min(length, max(0.0, effective_gap))


@numba.njit
def forward_gap_size(
    eta: np.ndarray, pt: np.ndarray, side: int, eta_min: float, eta_max: float, pt_min: float
):
    r"""
    Traditional forward rapidity gap $\Delta\eta^F$.

    Args:
        eta:     pseudorapidity of (charged) particles
        pt:      transverse momentum of (charged) particles
        side:    (-1, 1) backward / forward
        eta_min: minimum oriented fiducial eta (e.g. -2.5)
        eta_max: maximum absolute fiducial eta (e.g.  2.5)
        pt_min:  minimum fiducial pt threshold (e.g. 0.2) (GeV)

    Combine both hemispheres per event by calling twice:
    max(forward_gap_size(side=-1), forward_gap_size(side=1))
    """
    if side != -1 and side != 1:
        raise ValueError("forward_gap_size: side must be -1 or +1")
    if eta_max <= eta_min:
        raise ValueError("forward_gap_size: eta_max must exceed eta_min")
    if pt_min < 0.0:
        raise ValueError("forward_gap_size: pt_min must be non-negative")
    if eta.size != pt.size:
        raise ValueError("forward_gap_size: eta and pt must have equal size")

    outermost_u = eta_min
    for i, value in enumerate(eta):
        oriented = side * value
        if eta_min < oriented < eta_max and pt[i] > pt_min:
            outermost_u = max(outermost_u, oriented)

    return eta_max - outermost_u


# Compute the largest empty pseudorapidity interval including fiducial edges
@numba.njit
def largest_pseudorapidity_gap(
    eta: np.ndarray, pt: np.ndarray, eta_min: float, eta_max: float, pt_min: float
):
    r"""
    Largest conventional pseudorapidity gap inside the fiducial interval.

    The accepted particle pseudorapidities and both fiducial boundaries define
    the interval boundaries. The result is therefore the largest adjacent
    spacing after sorting, and an empty event spans the full acceptance.
    """
    if eta_max <= eta_min:
        raise ValueError("largest_pseudorapidity_gap: eta_max must exceed eta_min")
    if pt_min < 0.0:
        raise ValueError("largest_pseudorapidity_gap: pt_min must be non-negative")
    if eta.size != pt.size:
        raise ValueError("largest_pseudorapidity_gap: eta and pt must have equal size")

    mask = (eta_min < eta) & (eta < eta_max) & (pt > pt_min)
    accepted_eta = np.sort(eta[mask])
    if accepted_eta.size == 0:
        return eta_max - eta_min

    largest = accepted_eta[0] - eta_min
    for index in range(1, accepted_eta.size):
        largest = max(largest, accepted_eta[index] - accepted_eta[index - 1])
    largest = max(largest, eta_max - accepted_eta[-1])
    return largest


@numba.njit
def two_boundary_energy_gap_score(
    eta: np.ndarray, pt: np.ndarray, eta_min: float, eta_max: float, pt_min: float, q0: float
):
    r"""
    Two-boundary Energy Gap Flow.

    Two uniformly sampled boundaries x and y are connected with survival
    weight exp[-sum_i pT_i / Q0], where the sum includes accepted particles
    between the boundaries. Averaging this weight over both boundaries gives

        G_2^EF = L^-2 integral_I dx dy
                 exp[-sum_{min(x,y)<eta_i<max(x,y)} pT_i / Q0].

    The sorted recurrence evaluates the exact integral in O(N log N).
    """
    if eta_max <= eta_min:
        raise ValueError("two_boundary_energy_gap_score: eta_max must exceed eta_min")
    if pt_min < 0.0:
        raise ValueError("two_boundary_energy_gap_score: pt_min must be non-negative")
    if q0 <= 0.0:
        raise ValueError("two_boundary_energy_gap_score: q0 must be positive")
    if eta.size != pt.size:
        raise ValueError("two_boundary_energy_gap_score: eta and pt must have equal size")

    mask = (eta_min < eta) & (eta < eta_max) & (pt > pt_min)
    accepted_eta = eta[mask]
    accepted_pt = pt[mask]
    order = np.argsort(accepted_eta)
    accepted_eta = accepted_eta[order]
    accepted_pt = accepted_pt[order]
    length = eta_max - eta_min
    if accepted_eta.size == 0:
        return 1.0

    previous_eta = eta_min
    diagonal_sum = 0.0
    cross_sum = 0.0
    recurrence = 0.0

    for index in range(accepted_eta.size):
        left_interval = accepted_eta[index] - previous_eta
        diagonal_sum += left_interval * left_interval
        recurrence = np.exp(-accepted_pt[index] / q0) * (left_interval + recurrence)

        if index + 1 < accepted_eta.size:
            right_interval = accepted_eta[index + 1] - accepted_eta[index]
        else:
            right_interval = eta_max - accepted_eta[index]
        cross_sum += right_interval * recurrence
        previous_eta = accepted_eta[index]

    last_interval = eta_max - accepted_eta[-1]
    diagonal_sum += last_interval * last_interval
    score = (diagonal_sum + 2.0 * cross_sum) / (length * length)
    return min(1.0, max(0.0, score))


# Compute the leading-proton light-cone momentum loss on one beam side
@numba.njit
def leading_proton_lightcone_loss(
    eta: np.ndarray,
    pt: np.ndarray,
    pz: np.ndarray,
    energy: np.ndarray,
    side: int,
    beam_lightcone: float,
    eta_min: float,
    pt_max: float,
):
    r"""
    Evaluate the smallest forward-proton light-cone momentum loss

        xi_p,s = 1 - (E'_p + s p'_{z,p}) / (E_beam + |p_{z,beam}|)

    among stable final protons accepted on side s. A side without an accepted
    proton has unit loss, which is the least diffractive score.
    """
    if side != -1 and side != 1:
        raise ValueError("leading_proton_lightcone_loss: side must be -1 or +1")
    if beam_lightcone <= 0.0:
        raise ValueError("leading_proton_lightcone_loss: beam light-cone momentum must be positive")
    if eta_min < 0.0:
        raise ValueError("leading_proton_lightcone_loss: eta_min must be non-negative")
    if pt_max <= 0.0:
        raise ValueError("leading_proton_lightcone_loss: pt_max must be positive")
    if eta.size != pt.size or eta.size != pz.size or eta.size != energy.size:
        raise ValueError("leading_proton_lightcone_loss: particle arrays must have equal size")

    retained_fraction = 0.0
    for index in range(eta.size):
        if side * eta[index] > eta_min and pt[index] < pt_max:
            candidate = (energy[index] + side * pz[index]) / beam_lightcone
            retained_fraction = max(retained_fraction, candidate)

    retained_fraction = min(1.0, max(0.0, retained_fraction))
    return 1.0 - retained_fraction


# Compute true when eta contributes to at least one oriented gap flow side
def in_gap_flow_acceptance(eta, eta_min, eta_max):
    return (eta_min < eta < eta_max) or (eta_min < -eta < eta_max)


# Compute validated charged-particle fiducial parameters from event steering
def charged_fiducial_parameters(event):
    eta_min = event.cut_param["EGF_CH_ETA_MIN"]
    eta_max = event.cut_param["EGF_CH_ETA_MAX"]
    pt_min  = event.cut_param["EGF_CH_PT_MIN"]

    if eta_min < -eta_max:
        raise ValueError(
            "charged_fiducial_parameters: EGF_CH_ETA_MIN must be at least -EGF_CH_ETA_MAX"
        )
    if eta_max <= eta_min:
        raise ValueError("charged_fiducial_parameters: EGF_CH_ETA_MAX must exceed EGF_CH_ETA_MIN")
    if pt_min < 0.0:
        raise ValueError("charged_fiducial_parameters: EGF_CH_PT_MIN must be nonnegative")

    return eta_min, eta_max, pt_min


# Compute charged records accepted by the configured gap flow cuts
@cache
def charged_gap_flow_records(event):
    pair = selected_pair(event)
    selected_ids = {record["id"] for record in pair["records"]}
    records = obs.proj_event_particles(event)
    eta_min, eta_max, pt_min = charged_fiducial_parameters(event)

    return [
        record
        for record in records
        if record["is_final"]
        and record["pt"] > pt_min
        and in_gap_flow_acceptance(record["eta"], eta_min, eta_max)
        and DIFFZ_CUTS.is_charged_final_state_pid(record["pid"])
        and record["id"] not in selected_ids
    ]


# Compute the configured energy gap flow charged final-state multiplicity
def proj_n_charged_central(event):
    return len(charged_gap_flow_records(event))


# Compute both charged energy gap flow scores after removing the selected Z muons
@cache
def charged_energy_gap_flow_scores(event):
    charged_records = charged_gap_flow_records(event)
    eta = np.asarray([record["eta"] for record in charged_records], dtype=np.float64)
    pt = np.asarray([record["pt"] for record in charged_records], dtype=np.float64)
    eta_min, eta_max, pt_min = charged_fiducial_parameters(event)
    beta = event.cut_param["EGF_BETA"]
    q0 = event.cut_param["EGF_Q0"]

    left = energy_gap_flow_score(
        eta=eta, pt=pt, side=-1, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min, beta=beta, q0=q0
    )
    right = energy_gap_flow_score(
        eta=eta, pt=pt, side=1, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min, beta=beta, q0=q0
    )
    return left, right


# Compute both charged multiplicity gap flow scores after removing the selected Z muons
@cache
def charged_multiplicity_gap_flow_scores(event):
    charged_records = charged_gap_flow_records(event)
    eta = np.asarray([record["eta"] for record in charged_records], dtype=np.float64)
    pt = np.asarray([record["pt"] for record in charged_records], dtype=np.float64)
    eta_min, eta_max, pt_min = charged_fiducial_parameters(event)
    beta = event.cut_param["EGF_BETA"]
    p0 = event.cut_param["MGF_P0"]

    left = multiplicity_gap_flow_score(
        eta=eta, pt=pt, side=-1, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min, beta=beta, p0=p0
    )
    right = multiplicity_gap_flow_score(
        eta=eta, pt=pt, side=1, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min, beta=beta, p0=p0
    )
    return left, right


# Compute both charged energy gap flow scores in effective-gap coordinates
@cache
def charged_energy_gap_flow_effective_gaps(event):
    left_score, right_score = charged_energy_gap_flow_scores(event)
    eta_min, eta_max, _ = charged_fiducial_parameters(event)
    beta = event.cut_param["EGF_BETA"]

    left = energy_gap_flow_effective_gap(left_score, eta_min, eta_max, beta)
    right = energy_gap_flow_effective_gap(right_score, eta_min, eta_max, beta)
    return left, right


# Compute both classic rapidity gaps
@cache
def charged_forward_gap_sizes(event):
    charged_records = charged_gap_flow_records(event)
    eta = np.asarray([record["eta"] for record in charged_records], dtype=np.float64)
    pt = np.asarray([record["pt"] for record in charged_records], dtype=np.float64)
    eta_min, eta_max, pt_min = charged_fiducial_parameters(event)

    left = forward_gap_size(
        eta=eta, pt=pt, side=-1, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min
    )
    right = forward_gap_size(
        eta=eta, pt=pt, side=1, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min
    )

    return left, right


# Compute the largest charged-particle pseudorapidity gap after Z-muon removal
@cache
def charged_largest_pseudorapidity_gap(event):
    charged_records = charged_gap_flow_records(event)
    eta = np.asarray([record["eta"] for record in charged_records], dtype=np.float64)
    pt = np.asarray([record["pt"] for record in charged_records], dtype=np.float64)
    eta_min, eta_max, pt_min = charged_fiducial_parameters(event)
    return largest_pseudorapidity_gap(
        eta=eta, pt=pt, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min
    )


# Compute the charged two-boundary energy score after selected Z-muon removal
@cache
def charged_two_boundary_energy_gap_score(event):
    charged_records = charged_gap_flow_records(event)
    eta = np.asarray([record["eta"] for record in charged_records], dtype=np.float64)
    pt = np.asarray([record["pt"] for record in charged_records], dtype=np.float64)
    eta_min, eta_max, pt_min = charged_fiducial_parameters(event)
    q0 = event.cut_param["G2_EF_Q0"]
    return two_boundary_energy_gap_score(
        eta=eta, pt=pt, eta_min=eta_min, eta_max=eta_max, pt_min=pt_min, q0=q0
    )


# Compute incoming proton light-cone momenta on the negative and positive sides
@cache
def beam_proton_lightcone_momenta(event):
    records = obs.proj_event_particles(event)
    beam_minus = [
        record["p4"]
        for record in records
        if abs(record["pid"]) == pdg.PDG_PROTON
        and record["status"] == pdg.INITIAL_STATE
        and record["p4"].pz < 0.0
    ]
    beam_plus = [
        record["p4"]
        for record in records
        if abs(record["pid"]) == pdg.PDG_PROTON
        and record["status"] == pdg.INITIAL_STATE
        and record["p4"].pz > 0.0
    ]
    if not beam_minus or not beam_plus:
        raise Exception("beam_proton_lightcone_momenta: missing incoming proton beam")

    minus = max(beam_minus, key=lambda momentum: -momentum.pz)
    plus = max(beam_plus, key=lambda momentum: momentum.pz)
    return minus.e - minus.pz, plus.e + plus.pz


# Compute fiducial forward-proton light-cone losses on both beam sides
@cache
def forward_proton_lightcone_losses(event):
    records = obs.proj_event_particles(event)
    protons = [
        record for record in records if record["pid"] == pdg.PDG_PROTON and record["is_final"]
    ]
    eta = np.asarray([record["eta"] for record in protons], dtype=np.float64)
    pt = np.asarray([record["pt"] for record in protons], dtype=np.float64)
    pz = np.asarray([record["p4"].pz for record in protons], dtype=np.float64)
    energy = np.asarray([record["p4"].e for record in protons], dtype=np.float64)
    beam_minus, beam_plus = beam_proton_lightcone_momenta(event)
    eta_min = event.cut_param["FORWARD_PROTON_ETA_MIN"]
    pt_max = event.cut_param["FORWARD_PROTON_PT_MAX"]

    left = leading_proton_lightcone_loss(
        eta=eta,
        pt=pt,
        pz=pz,
        energy=energy,
        side=-1,
        beam_lightcone=beam_minus,
        eta_min=eta_min,
        pt_max=pt_max,
    )
    right = leading_proton_lightcone_loss(
        eta=eta,
        pt=pt,
        pz=pz,
        energy=energy,
        side=1,
        beam_lightcone=beam_plus,
        eta_min=eta_min,
        pt_max=pt_max,
    )
    return left, right


# Compute the charged energy gap flow score on the left (negative rapidity)
def proj_g_egf_left(event):
    return charged_energy_gap_flow_scores(event)[0]


# Compute the charged energy gap flow score on the right (positive rapidity)
def proj_g_egf_right(event):
    return charged_energy_gap_flow_scores(event)[1]


# Compute the maximum charged energy gap flow score
def proj_g_egf_max(event):
    left, right = charged_energy_gap_flow_scores(event)
    return max(left, right)


# Compute the charged multiplicity gap flow score on the left or backward side
def proj_g_mgf_left(event):
    return charged_multiplicity_gap_flow_scores(event)[0]


# Compute the charged multiplicity gap flow score on the right or forward side
def proj_g_mgf_right(event):
    return charged_multiplicity_gap_flow_scores(event)[1]


# Compute the maximum charged multiplicity gap flow score
def proj_g_mgf_max(event):
    left, right = charged_multiplicity_gap_flow_scores(event)
    return max(left, right)


# Compute the effective energy gap flow gap on the left or backward side
def proj_deta_egf_left(event):
    return charged_energy_gap_flow_effective_gaps(event)[0]


# Compute the effective energy gap flow gap on the right or forward side
def proj_deta_egf_right(event):
    return charged_energy_gap_flow_effective_gaps(event)[1]


# Compute the larger of the two effective energy gap flow gaps
def proj_deta_egf_max(event):
    left, right = charged_energy_gap_flow_effective_gaps(event)
    return max(left, right)


# Compute the traditional forward rapidity gap on the left or backward side
def proj_deta_f_left(event):
    return charged_forward_gap_sizes(event)[0]


# Compute the traditional forward rapidity gap on the right or forward side
def proj_deta_f_right(event):
    return charged_forward_gap_sizes(event)[1]


# Compute the larger of the two forward rapidity gaps
def proj_deta_f_max(event):
    left, right = charged_forward_gap_sizes(event)
    return max(left, right)


# Compute the largest charged-particle pseudorapidity gap in the event
def proj_deta_largest(event):
    return charged_largest_pseudorapidity_gap(event)


# Compute the charged two-boundary transverse-energy Laplace score
def proj_G2_ef(event):
    return charged_two_boundary_energy_gap_score(event)


# Compute the minimum fiducial forward-proton light-cone momentum loss
def proj_xi_forward_proton_min(event):
    left, right = forward_proton_lightcone_losses(event)
    return min(left, right)


# Compute the logarithmic forward-proton diffraction score
def proj_d_forward_proton(event):
    xi_floor = event.cut_param["FORWARD_PROTON_XI_FLOOR"]
    if xi_floor <= 0.0 or xi_floor > 1.0:
        raise ValueError("proj_d_forward_proton: FORWARD_PROTON_XI_FLOOR must be in (0, 1]")
    xi = max(xi_floor, proj_xi_forward_proton_min(event))
    return -np.log10(xi)


# Compute all diffraction scores used by the ROC plots
def proj_gap_roc_scores(event):
    return (
        proj_g_egf_max(event),
        proj_g_mgf_max(event),
        proj_deta_f_max(event),
        proj_deta_largest(event),
        proj_G2_ef(event),
        proj_n_charged_central(event),
        proj_d_forward_proton(event),
    )


# ----------------------------------------------------------
# Histograms
# ----------------------------------------------------------

obs_m_mumu = {
    "tag": "m_mumu",
    "xlim": (66, 116),
    "ylim": None,
    "xlabel": r"$m_{\mu\mu}$",
    "ylabel": r"$d\sigma/dm_{\mu\mu}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Dimuon mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(66, 116, 51),
    "density": False,
    "func": proj_m_mumu,
}

obs_pt_mumu = {
    "tag": "pt_mumu",
    "xlim": (DIFFZ_CUTS.cut_param["PT_MUMU_MIN"], DIFFZ_CUTS.cut_param["PT_MUMU_MAX"]),
    "ylim": None,
    "xlabel": r"$p_{T}^{\mu\mu}$",
    "ylabel": r"$d\sigma/dp_{T}^{\mu\mu}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Dimuon transverse momentum",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        DIFFZ_CUTS.cut_param["PT_MUMU_MIN"], DIFFZ_CUTS.cut_param["PT_MUMU_MAX"], 61
    ),
    "density": False,
    "func": proj_pt_mumu,
}

obs_y_mumu = {
    "tag": "y_mumu",
    "xlim": (-3, 3),
    "ylim": None,
    "xlabel": r"$y_{\mu\mu}$",
    "ylabel": r"$d\sigma/dy_{\mu\mu}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Dimuon rapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(-3, 3, 49),
    "density": False,
    "func": proj_y_mumu,
}


obs_pt_mu1 = {
    "tag": "pt_mu1",
    "xlim": (10, 120),
    "ylim": None,
    "xlabel": r"$p_{T}^{\mu,1}$",
    "ylabel": r"$d\sigma/dp_{T}^{\mu,1}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Leading muon transverse momentum",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(10, 120, 45),
    "density": False,
    "func": proj_pt_mu1,
}

obs_pt_mu2 = {
    "tag": "pt_mu2",
    "xlim": (10, 100),
    "ylim": None,
    "xlabel": r"$p_{T}^{\mu,2}$",
    "ylabel": r"$d\sigma/dp_{T}^{\mu,2}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Subleading muon transverse momentum",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(10, 100, 45),
    "density": False,
    "func": proj_pt_mu2,
}

obs_eta_mu1 = {
    "tag": "eta_mu1",
    "xlim": (-DIFFZ_CUTS.cut_param["MUON_ETA_MAX"], DIFFZ_CUTS.cut_param["MUON_ETA_MAX"]),
    "ylim": None,
    "xlabel": r"$\eta_{\mu,1}$",
    "ylabel": r"$d\sigma/d\eta_{\mu,1}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Leading muon pseudorapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        -DIFFZ_CUTS.cut_param["MUON_ETA_MAX"], DIFFZ_CUTS.cut_param["MUON_ETA_MAX"], 41
    ),
    "density": False,
    "func": proj_eta_mu1,
}

obs_eta_mu2 = {
    "tag": "eta_mu2",
    "xlim": (-DIFFZ_CUTS.cut_param["MUON_ETA_MAX"], DIFFZ_CUTS.cut_param["MUON_ETA_MAX"]),
    "ylim": None,
    "xlabel": r"$\eta_{\mu,2}$",
    "ylabel": r"$d\sigma/d\eta_{\mu,2}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Subleading muon pseudorapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        -DIFFZ_CUTS.cut_param["MUON_ETA_MAX"], DIFFZ_CUTS.cut_param["MUON_ETA_MAX"], 41
    ),
    "density": False,
    "func": proj_eta_mu2,
}

obs_deta_mumu = {
    "tag": "deta_mumu",
    "xlim": (0, 1.4),
    "ylim": None,
    "xlabel": r"$|\Delta\eta_{\mu\mu}|$",
    "ylabel": r"$d\sigma/d|\Delta\eta_{\mu\mu}|$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Dimuon pseudorapidity separation",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0, 1.4, 29),
    "density": False,
    "func": proj_deta_mumu,
}

obs_acoplanarity = {
    "tag": "acoplanarity",
    "xlim": (0, 0.2),
    "ylim": None,
    "xlabel": r"$1 - |\Delta\phi_{\mu\mu}|/\pi$",
    "ylabel": r"$d\sigma/dA_{\phi}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Dimuon acoplanarity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0, 0.2, 41),
    "density": False,
    "func": proj_acoplanarity,
}

obs_n_charged_central = {
    "tag": "n_charged_central",
    "xlim": (-0.5, 250.5),
    "ylim": None,
    "xlabel": rf"$N_{{\mathrm{{ch}}}}("
    rf"{DIFFZ_CUTS.cut_param['EGF_CH_ETA_MIN']:g}<\eta<"
    rf"{DIFFZ_CUTS.cut_param['EGF_CH_ETA_MAX']:g},\ p_{{T}}>"
    rf"{DIFFZ_CUTS.cut_param['EGF_CH_PT_MIN']:g}\,\mathrm{{GeV}})$",
    "ylabel": r"$d\sigma/dN_{\mathrm{ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged final-state multiplicity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.arange(-0.5, 251.5, 5.0),
    "density": False,
    "func": proj_n_charged_central,
}

obs_g_egf_left = {
    "tag": "g_egf_left",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"$\mathcal{G}_{\mathrm{L}}^{\mathrm{EF,ch}}$",
    "ylabel": r"$d\sigma/d\mathcal{G}_{\mathrm{L}}^{\mathrm{EF,ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged energy gap flow, left/backward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 1.0, 51),
    "density": False,
    "func": proj_g_egf_left,
}

obs_g_egf_right = {
    "tag": "g_egf_right",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"$\mathcal{G}_{\mathrm{R}}^{\mathrm{EF,ch}}$",
    "ylabel": r"$d\sigma/d\mathcal{G}_{\mathrm{R}}^{\mathrm{EF,ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged energy gap flow, right/forward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 1.0, 51),
    "density": False,
    "func": proj_g_egf_right,
}

obs_g_egf_max = {
    "tag": "g_egf_max",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"$\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{EF,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{EF,ch}})$",
    "ylabel": r"$d\sigma/d\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{EF,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{EF,ch}})$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Maximum charged energy gap flow",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 1.0, 51),
    "density": False,
    "func": proj_g_egf_max,
}

obs_g_mgf_left = {
    "tag": "g_mgf_left",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"$\mathcal{G}_{\mathrm{L}}^{\mathrm{mult,ch}}$",
    "ylabel": r"$d\sigma/d\mathcal{G}_{\mathrm{L}}^{\mathrm{mult,ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged multiplicity gap flow, left/backward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 1.0, 51),
    "density": False,
    "func": proj_g_mgf_left,
}

obs_g_mgf_right = {
    "tag": "g_mgf_right",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"$\mathcal{G}_{\mathrm{R}}^{\mathrm{mult,ch}}$",
    "ylabel": r"$d\sigma/d\mathcal{G}_{\mathrm{R}}^{\mathrm{mult,ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged multiplicity gap flow, right/forward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 1.0, 51),
    "density": False,
    "func": proj_g_mgf_right,
}

obs_g_mgf_max = {
    "tag": "g_mgf_max",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"$\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{mult,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{mult,ch}})$",
    "ylabel": r"$d\sigma/d\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{mult,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{mult,ch}})$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Maximum charged multiplicity gap flow",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 1.0, 51),
    "density": False,
    "func": proj_g_mgf_max,
}

obs_deta_egf_left = {
    "tag": "deta_egf_left",
    "xlim": (0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"]),
    "ylim": None,
    "xlabel": r"$\Delta\eta_{\mathrm{L}}^{\mathrm{EGF,ch}}$",
    "ylabel": r"$d\sigma/d\Delta\eta_{\mathrm{L}}^{\mathrm{EGF,ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Effective charged energy gap flow, left/backward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"], 51
    ),
    "density": False,
    "func": proj_deta_egf_left,
}

obs_deta_egf_right = {
    "tag": "deta_egf_right",
    "xlim": (0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"]),
    "ylim": None,
    "xlabel": r"$\Delta\eta_{\mathrm{R}}^{\mathrm{EGF,ch}}$",
    "ylabel": r"$d\sigma/d\Delta\eta_{\mathrm{R}}^{\mathrm{EGF,ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Effective charged energy gap flow, right/forward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"], 51
    ),
    "density": False,
    "func": proj_deta_egf_right,
}

obs_deta_egf_max = {
    "tag": "deta_egf_max",
    "xlim": (0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"]),
    "ylim": None,
    "xlabel": r"$\max(\Delta\eta_{\mathrm{L}}^{\mathrm{EGF,ch}},\Delta\eta_{\mathrm{R}}^{\mathrm{EGF,ch}})$",
    "ylabel": r"$d\sigma/d\max(\Delta\eta_{\mathrm{L}}^{\mathrm{EGF,ch}},\Delta\eta_{\mathrm{R}}^{\mathrm{EGF,ch}})$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Maximum effective charged energy gap flow",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"], 51
    ),
    "density": False,
    "func": proj_deta_egf_max,
}

obs_deta_f_left = {
    "tag": "deta_f_left",
    "xlim": (0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"]),
    "ylim": None,
    "xlabel": r"$\Delta\eta_{\mathrm{L}}^{F,\mathrm{ch}}$",
    "ylabel": r"$d\sigma/d\Delta\eta_{\mathrm{L}}^{F,\mathrm{ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged forward gap, left/backward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"], 51
    ),
    "density": False,
    "func": proj_deta_f_left,
}

obs_deta_f_right = {
    "tag": "deta_f_right",
    "xlim": (0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"]),
    "ylim": None,
    "xlabel": r"$\Delta\eta_{\mathrm{R}}^{F,\mathrm{ch}}$",
    "ylabel": r"$d\sigma/d\Delta\eta_{\mathrm{R}}^{F,\mathrm{ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged forward gap, right/forward",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"], 51
    ),
    "density": False,
    "func": proj_deta_f_right,
}

obs_deta_f_max = {
    "tag": "deta_f_max",
    "xlim": (0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"]),
    "ylim": None,
    "xlabel": r"$\max(\Delta\eta_{\mathrm{L}}^{F,\mathrm{ch}},\Delta\eta_{\mathrm{R}}^{F,\mathrm{ch}})$",
    "ylabel": r"$d\sigma/d\max(\Delta\eta_{\mathrm{L}}^{F,\mathrm{ch}},\Delta\eta_{\mathrm{R}}^{F,\mathrm{ch}})$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Maximum charged forward gap",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"], 51
    ),
    "density": False,
    "func": proj_deta_f_max,
}

obs_deta_largest = {
    "tag": "deta_largest",
    "xlim": (0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"]),
    "ylim": None,
    "xlabel": r"$\Delta\eta_{\max}^{\mathrm{ch}}$",
    "ylabel": r"$d\sigma/d\Delta\eta_{\max}^{\mathrm{ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Largest charged-particle pseudorapidity gap",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(
        0.0, DIFFZ_CUTS.cut_param["EGF_CH_ETA_MAX"] - DIFFZ_CUTS.cut_param["EGF_CH_ETA_MIN"], 51
    ),
    "density": False,
    "func": proj_deta_largest,
}

obs_G2_ef = {
    "tag": "G2_ef",
    "xlim": (0.0, 1.0),
    "ylim": None,
    "xlabel": r"$\mathcal{G}_{2}^{\mathrm{EF,ch}}$",
    "ylabel": r"$d\sigma/d\mathcal{G}_{2}^{\mathrm{EF,ch}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Charged two-boundary energy Laplace score",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 1.0, 51),
    "density": False,
    "func": proj_G2_ef,
}

obs_d_forward_proton = {
    "tag": "d_forward_proton",
    "xlim": (0.0, 4.0),
    "ylim": None,
    "xlabel": r"$D_p=-\log_{10}(\xi_{p,\min})$",
    "ylabel": r"$d\sigma/dD_p$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Forward leading-proton light-cone loss score",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 4.0, 81),
    "density": False,
    "func": proj_d_forward_proton,
}

obs_gap_roc = {
    "tag": "gap_roc",
    "kind": "roc",
    "xlim": (1.0e-3, 1.0),
    "ylim": (0.0, 1.0),
    "xscale": "log",
    "yscale": "linear",
    "xlabel": r"FPR (DY efficiency)",
    "ylabel": r"TPR (hard-SD efficiency)",
    "units": {"x": r"unit", "y": r"unit"},
    "label": r"Diffractive-observable ROC comparison",
    "func": proj_gap_roc_scores,
    "roc": {
        "score_labels": [
            r"$\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{EF,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{EF,ch}})$",
            r"$\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{mult,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{mult,ch}})$",
            r"$\max(\Delta\eta_{\mathrm{L}}^{F,\mathrm{ch}},\Delta\eta_{\mathrm{R}}^{F,\mathrm{ch}})$",
            r"$\Delta\eta_{\max}^{\mathrm{ch}}$",
            r"$\mathcal{G}_{2}^{\mathrm{EF,ch}}$",
            r"$N_{\mathrm{ch}}$",
            r"$D_p=-\log_{10}(\xi_{p,\min})$",
        ],
        "directions": [
            "higher",
            "higher",
            "higher",
            "higher",
            "higher",
            "lower",
            "higher",
        ],
        "score_colors": [
            "#0072B2",
            "#D55E00",
            "#009E73",
            "#000000",
            "#E6AB02",
            "#CC79A7",
            "#332288",
        ],
        "comparisons": [
            {"signal": 4, "background": 2, "linestyle": "-", "label": "MPI off for SD"},
            {"signal": 3, "background": 2, "linestyle": "--", "label": "MPI on for SD"},
        ],
        "label_template": "{score} (AUC = {auc:.3f}) [{comparison}]",
        "show_auc": True,
        "show_diagonal": True,
        "linewidth": 1.4,
    },
}

obs_gap_roc_reversed = {
    "tag": "gap_roc_reversed",
    "kind": "roc",
    "xlim": (1.0e-3, 1.0),
    "ylim": (0.0, 1.0),
    "xscale": "log",
    "yscale": "linear",
    "xlabel": r"FPR (hard-SD efficiency)",
    "ylabel": r"TPR (DY efficiency)",
    "units": {"x": r"unit", "y": r"unit"},
    "label": r"Reversed diffractive-observable ROC comparison",
    "func": proj_gap_roc_scores,
    "roc": {
        "score_labels": [
            r"$\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{EF,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{EF,ch}})$",
            r"$\max(\mathcal{G}_{\mathrm{L}}^{\mathrm{mult,ch}},\mathcal{G}_{\mathrm{R}}^{\mathrm{mult,ch}})$",
            r"$\max(\Delta\eta_{\mathrm{L}}^{F,\mathrm{ch}},\Delta\eta_{\mathrm{R}}^{F,\mathrm{ch}})$",
            r"$\Delta\eta_{\max}^{\mathrm{ch}}$",
            r"$\mathcal{G}_{2}^{\mathrm{EF,ch}}$",
            r"$N_{\mathrm{ch}}$",
            r"$D_p=-\log_{10}(\xi_{p,\min})$",
        ],
        "directions": ["lower", "lower", "lower", "lower", "lower", "higher", "lower"],
        "score_colors": [
            "#0072B2",
            "#D55E00",
            "#009E73",
            "#000000",
            "#E6AB02",
            "#CC79A7",
            "#332288",
        ],
        "comparisons": [
            {"signal": 2, "background": 4, "linestyle": "-", "label": "MPI off for SD"},
            {"signal": 2, "background": 3, "linestyle": "--", "label": "MPI on for SD"},
        ],
        "label_template": "{score} (AUC = {auc:.3f}) [{comparison}]",
        "show_auc": True,
        "show_diagonal": True,
        "linewidth": 1.4,
    },
}

# Reuse the forward ROC steering with a logarithmic hard-SD efficiency axis
obs_gap_roc_logtpr = copy.deepcopy(obs_gap_roc)
obs_gap_roc_logtpr.update(
    {
        "tag": "gap_roc_logtpr",
        "ylim": (1.0e-2, 1.0),
        "yscale": "log",
        "label": r"Diffractive-observable ROC comparison with logarithmic TPR",
    }
)

# Reuse the reversed ROC steering with a logarithmic DY efficiency axis
obs_gap_roc_reversed_logtpr = copy.deepcopy(obs_gap_roc_reversed)
obs_gap_roc_reversed_logtpr.update(
    {
        "tag": "gap_roc_reversed_logtpr",
        "ylim": (1.0e-2, 1.0),
        "yscale": "log",
        "label": r"Reversed diffractive-observable ROC comparison with logarithmic TPR",
    }
)
