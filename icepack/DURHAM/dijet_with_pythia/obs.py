# Durham quark and gluon particle level jet observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import numpy as np
from core.kinematics import fastjet

from . import cuts as DURHAM_CUTS


# Compute cached particle and jet arrays for one accepted event
def jet_data(event):
    return DURHAM_CUTS.particle_level_jets(event)

# Leading dijet invariant mass
def proj_m_jj(event):
    return DURHAM_CUTS.dijet_mass(jet_data(event)["jets"])

# Leading jet transverse momentum
def proj_pt_j1(event):
    jet = jet_data(event)["jets"][0]
    return fastjet.jet_pt(jet[0], jet[1])

# Leading jet absolute pseudorapidity
def proj_abs_eta_j1(event):
    jet = jet_data(event)["jets"][0]
    return abs(fastjet.jet_eta(jet[0], jet[1], jet[2]))

# Azimuthal separation of the leading dijet
def proj_dphi_jj(event):
    jets = jet_data(event)["jets"]
    return fastjet.abs_delta_phi(
        fastjet.jet_phi(jets[0, 0], jets[0, 1]),
        fastjet.jet_phi(jets[1, 0], jets[1, 1]),
    )

# Constituent transverse momenta and angular distances for one jet
def jet_constituent_kinematics(data, jet_index):
    jet = data["jets"][jet_index]
    jet_eta = fastjet.jet_eta(jet[0], jet[1], jet[2])
    jet_phi = fastjet.jet_phi(jet[0], jet[1])
    transverse_momenta = []
    angular_distances = []
    for particle_index in range(data["labels"].size):
        if data["labels"][particle_index] != jet_index:
            continue
        particle_pt = fastjet.jet_pt(data["px"][particle_index], data["py"][particle_index])
        particle_eta = fastjet.jet_eta(
            data["px"][particle_index],
            data["py"][particle_index],
            data["pz"][particle_index],
        )
        particle_phi = fastjet.jet_phi(data["px"][particle_index], data["py"][particle_index])
        transverse_momenta.append(particle_pt)
        angular_distances.append(
            math.sqrt(fastjet.angular_distance2(jet_eta, jet_phi, particle_eta, particle_phi))
        )
    return np.asarray(transverse_momenta), np.asarray(angular_distances)

# Total constituent multiplicity of the leading dijet
def proj_n_constituents_jj(event):
    data = jet_data(event)
    return sum(np.count_nonzero(data["labels"] == index) for index in (0, 1))

# Mean invariant mass of the two leading jets
def proj_mean_jet_mass(event):
    jets = jet_data(event)["jets"]
    masses = [fastjet.jet_mass(jet[0], jet[1], jet[2], jet[3]) for jet in jets[:2]]
    return 0.5 * sum(masses)

# Mean transverse momentum weighted radial width of the leading dijet
def proj_mean_jet_width(event):
    data = jet_data(event)
    widths = []
    for jet_index in (0, 1):
        transverse_momenta, angular_distances = jet_constituent_kinematics(data, jet_index)
        momentum_sum = np.sum(transverse_momenta)
        widths.append(
            np.sum(transverse_momenta * angular_distances) / momentum_sum
            if momentum_sum > 0.0
            else 0.0
        )
    return 0.5 * sum(widths)

# Mean pTD fragmentation observable of the leading dijet
def proj_mean_ptd(event):
    data = jet_data(event)
    values = []
    for jet_index in (0, 1):
        transverse_momenta, _ = jet_constituent_kinematics(data, jet_index)
        momentum_sum = np.sum(transverse_momenta)
        values.append(
            math.sqrt(np.sum(transverse_momenta * transverse_momenta)) / momentum_sum
            if momentum_sum > 0.0
            else 0.0
        )
    return 0.5 * sum(values)

# Mean leading constituent momentum fraction of the leading dijet
def proj_mean_z_leading(event):
    data = jet_data(event)
    values = []
    for jet_index in (0, 1):
        transverse_momenta, _ = jet_constituent_kinematics(data, jet_index)
        momentum_sum = np.sum(transverse_momenta)
        values.append(
            np.max(transverse_momenta) / momentum_sum
            if transverse_momenta.size > 0 and momentum_sum > 0.0
            else 0.0
        )
    return 0.5 * sum(values)

# Construct one standard cross section histogram definition
def histogram(tag, xlim, xlabel, ylabel, units, label, bins, function):
    return {
        "tag": tag,
        "xlim": xlim,
        "ylim": None,
        "xlabel": xlabel,
        "ylabel": ylabel,
        "units": {"x": units, "y": r"pb"},
        "label": label,
        "ylim_ratio": (0.0, 2.0),
        "ytick_ratio_step": 0.5,
        "bins": bins,
        "density": False,
        "func": function,
    }


obs_m_jj = histogram(
    "m_jj",
    (40.0, 200.0),
    r"$m_{jj}$",
    r"$d\sigma/dm_{jj}$",
    r"GeV",
    r"Dijet invariant mass",
    np.linspace(40.0, 200.0, 33),
    proj_m_jj,
)
obs_pt_j1 = histogram(
    "pt_j1",
    (20.0, 100.0),
    r"$p_T^{j_1}$",
    r"$d\sigma/dp_T^{j_1}$",
    r"GeV",
    r"Leading jet transverse momentum",
    np.linspace(20.0, 100.0, 33),
    proj_pt_j1,
)
obs_abs_eta_j1 = histogram(
    "abs_eta_j1",
    (0.0, 2.5),
    r"$|\eta_{j_1}|$",
    r"$d\sigma/d|\eta_{j_1}|$",
    r"unit",
    r"Leading jet absolute pseudorapidity",
    np.linspace(0.0, 2.5, 26),
    proj_abs_eta_j1,
)
obs_dphi_jj = histogram(
    "dphi_jj",
    (0.0, math.pi),
    r"$|\Delta\phi_{jj}|$",
    r"$d\sigma/d|\Delta\phi_{jj}|$",
    r"rad",
    r"Dijet azimuthal separation",
    np.linspace(0.0, math.pi, 33),
    proj_dphi_jj,
)
obs_n_constituents_jj = histogram(
    "n_constituents_jj",
    (0.0, 100.0),
    r"$N_{\mathrm{const}}^{j_1}+N_{\mathrm{const}}^{j_2}$",
    r"$d\sigma/dN_{\mathrm{const}}$",
    r"unit",
    r"Leading dijet constituent multiplicity",
    np.arange(-0.5, 100.5, 1.0),
    proj_n_constituents_jj,
)
obs_mean_jet_mass = histogram(
    "mean_jet_mass",
    (0.0, 40.0),
    r"$\langle m_j\rangle$",
    r"$d\sigma/d\langle m_j\rangle$",
    r"GeV",
    r"Mean leading dijet mass",
    np.linspace(0.0, 40.0, 41),
    proj_mean_jet_mass,
)
obs_mean_jet_width = histogram(
    "mean_jet_width",
    (0.0, 0.6),
    r"$\langle\sum_i p_{T,i}\Delta R_{i,j}/p_{T,j}\rangle$",
    r"$d\sigma/d\langle w_j\rangle$",
    r"unit",
    r"Mean leading dijet width",
    np.linspace(0.0, 0.6, 31),
    proj_mean_jet_width,
)
obs_mean_ptd = histogram(
    "mean_ptd",
    (0.0, 1.0),
    r"$\langle p_T^D\rangle$",
    r"$d\sigma/d\langle p_T^D\rangle$",
    r"unit",
    r"Mean leading dijet transverse momentum dispersion",
    np.linspace(0.0, 1.0, 41),
    proj_mean_ptd,
)
obs_mean_z_leading = histogram(
    "mean_z_leading",
    (0.0, 1.0),
    r"$\langle p_{T,\mathrm{lead}}^{\mathrm{const}}/\sum_i p_{T,i}\rangle$",
    r"$d\sigma/d\langle z_{\mathrm{lead}}\rangle$",
    r"unit",
    r"Mean leading constituent momentum fraction",
    np.linspace(0.0, 1.0, 41),
    proj_mean_z_leading,
)
