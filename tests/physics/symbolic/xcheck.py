# Symbolic normalization cross-checks for process families
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import sympy as sp

# Analytic measures for MFactorized, MCentral, MQuasiElastic and MTensorPomeron
# Photon and Durham derivations live in epa.py, gamma_amp.py and durham.py


# Require an exact symbolic expression to vanish
def require_zero(name, expr):
    value = sp.simplify(expr)
    if value != 0:
        raise AssertionError(f"{name}: expected 0, got {value}")


# Require two exact symbolic expressions to be equal
def require_equal(name, expr, expected):
    require_zero(name, sp.simplify(expr - expected))


# Prove the Moller flux used by the process classes
def prove_moller_flux_conventions():
    s, m1sq, m2sq = sp.symbols("s m1sq m2sq", positive=True)
    lam = s**2 + m1sq**2 + m2sq**2 - 2 * s * (m1sq + m2sq) - 2 * m1sq * m2sq
    beta = sp.sqrt(lam) / s

    # The invariant flux is 4 sqrt[(p1.p2)^2 - m1^2 m2^2]
    # [REFERENCE: PDG 2025, Kinematics review, scattering cross section]
    p1_dot_p2 = (s - m1sq - m2sq) / 2
    generic_two_body_flux = sp.simplify(4 * sp.sqrt(p1_dot_p2**2 - m1sq * m2sq))
    high_energy_flux = generic_two_body_flux.subs({m1sq: 0, m2sq: 0})
    quasielastic_dt_norm = 1 / (16 * sp.pi * lam)
    code_quasielastic_dt_norm = 1 / (16 * sp.pi * (generic_two_body_flux / 2) ** 2)

    require_equal("exact invariant Moller flux", generic_two_body_flux, 2 * s * beta)
    require_equal("high-energy Moller flux", high_energy_flux, 2 * s)
    require_equal("quasielastic dt normalization", code_quasielastic_dt_norm, quasielastic_dt_norm)

    return {
        "generic_flux": generic_two_body_flux,
        "high_energy_flux": high_energy_flux,
        "dt_prefactor": code_quasielastic_dt_norm,
    }


# Prove the factorized 2->3 x 1->M phase-space density
def prove_factorized_phase_space():
    s, beta, e1, e2, jac = sp.symbols("s beta E1 E2 J", positive=True)

    # Reduce (2pi)^4 delta^4(P-sum p) product d^3p/[(2pi)^3 2E]
    # The transverse and longitudinal momentum deltas remove pX_T and p2z
    # dpXz/(2 EX) = dy/2, the energy delta contributes J = 1/|v1z-v2z|
    # [REFERENCE: PDG 2025, Kinematics review, Eqs. (49.11)-(49.12)]
    dphi3_density = (2 * sp.pi) ** 4 / (2 * sp.pi) ** 9 * jac / (2 * e1 * 2 * e2 * 2)
    dsigma_density = sp.simplify(dphi3_density / (2 * s * beta))
    dsigma_with_decay = sp.simplify(dsigma_density / (2 * sp.pi))
    high_energy = {beta: 1, e1: sp.sqrt(s) / 2, e2: sp.sqrt(s) / 2, jac: sp.Rational(1, 2)}

    require_equal("factorized dPhi3 density", dphi3_density, jac / (8 * e1 * e2 * (2 * sp.pi) ** 5))
    require_equal("factorized dSigma density", dsigma_density, jac / (16 * beta * s * e1 * e2 * (2 * sp.pi) ** 5))
    require_equal(
        "factorized decay bridge", dsigma_with_decay, jac / (32 * sp.pi * beta * s * e1 * e2 * (2 * sp.pi) ** 5)
    )
    require_equal(
        "factorized high-energy dSigma",
        sp.simplify(dsigma_density.subs(high_energy)),
        1 / (8 * s**2 * (2 * sp.pi) ** 5),
    )

    return {
        "dPhi3": dphi3_density,
        "dSigma": dsigma_density,
        "dSigma_decay": dsigma_with_decay,
        "high_energy_dSigma": sp.simplify(dsigma_density.subs(high_energy)),
    }


# Prove the direct 2->N continuum phase-space density
def prove_continuum_phase_space():
    s, beta, e1, e2, jac = sp.symbols("s beta E1 E2 J", positive=True)
    out = {}

    for kf in (2, 4, 6):
        nf = kf + 2
        central_transverse_map_det = prove_central_transverse_map_det(kf)
        central_rapidity_measure = sp.Rational(1, 2) ** kf
        forward_pt_measure = 1 / (4 * e1 * e2)
        dphi = sp.simplify(
            central_transverse_map_det
            * central_rapidity_measure
            * forward_pt_measure
            * jac
            / (2 * sp.pi) ** (3 * nf - 4)
        )
        dsigma = sp.simplify(dphi / (2 * s * beta))
        require_equal(
            f"continuum K={kf} dPhi density",
            dphi,
            jac / (kf**2 * 2 ** (kf + 2) * e1 * e2 * (2 * sp.pi) ** (3 * nf - 4)),
        )
        require_equal(
            f"continuum K={kf} dSigma density",
            dsigma,
            jac / (kf**2 * 2 ** (kf + 3) * beta * s * e1 * e2 * (2 * sp.pi) ** (3 * nf - 4)),
        )
        out[f"K={kf}"] = {
            "Nf": nf,
            "central_transverse_map_det": central_transverse_map_det,
            "central_rapidity_measure": central_rapidity_measure,
            "forward_pt_measure": forward_pt_measure,
            "dPhi": dphi,
            "dSigma": dsigma,
        }

    return out


# Prove adjacent differences with fixed total transverse momentum give determinant 1/K
def prove_central_transverse_map_det(kf):
    # Solve sum p_i = P and p_i-p_{i+1} = q_i directly on one axis
    momenta = sp.Matrix(sp.symbols(f"p0:{kf}"))
    constraints = sp.Matrix([sum(momenta), *(momenta[i] - momenta[i + 1] for i in range(kf - 1))])
    inverse = constraints.jacobian(momenta).inv()
    coeff = inverse[: kf - 1, 1:]

    determinant_one_axis = sp.simplify(coeff.det())
    transverse_map_det = sp.simplify(determinant_one_axis**2)

    require_equal(f"central transverse determinant K={kf}", abs(determinant_one_axis), sp.Rational(1, kf))
    require_equal(f"central transverse map determinant K={kf}", transverse_map_det, sp.Rational(1, kf**2))

    return transverse_map_det


# Prove the elastic and diffractive quasielastic phase-space factors
def prove_quasielastic_phase_space():
    # MQuasiElastic::B3PhaseSpaceWeight uses the same factor for EL, SD and DD
    # MProcess::DissociationCrossSectionFactor separately counts both SD sides
    s, pin, pout = sp.symbols("s p_in p_out", positive=True)
    costheta = sp.symbols("cos_theta", real=True)
    t0 = sp.symbols("t0", real=True)
    dphi_dcos = pout / (8 * sp.pi * sp.sqrt(s))
    flux = 4 * pin * sp.sqrt(s)
    t = t0 + 2 * pin * pout * costheta
    density = sp.simplify(dphi_dcos / sp.diff(t, costheta) / flux)
    lam = 4 * s * pin**2
    require_equal("quasielastic dt factor", density, 1 / (16 * sp.pi * lam))
    # Unequal final masses cancel from dt/dcos, there is no SD factor in LIPS
    return {"EL_SD_DD_per_side": density, "SD_side_sum": 2 * density}


# Prove the Breit-Wigner normalization and its physical mass threshold
def prove_delta_bw_normalization():
    mass, gamma = sp.symbols("M Gamma", positive=True)
    z = sp.symbols("z", real=True)
    a = mass * gamma

    delta_bw = a / (sp.pi * (z**2 + a**2))
    integral = sp.integrate(delta_bw, (z, -sp.oo, sp.oo))

    require_equal("delta BW normalization", integral, 1)
    # The whole real axis is the narrow-width extension, not the physical domain
    physical_integral = sp.integrate(delta_bw, (z, -(mass**2), sp.oo))
    require_equal("delta BW physical mass range", physical_integral, sp.Rational(1, 2) + sp.atan(mass / gamma) / sp.pi)
    require_equal("delta BW narrow-width limit", sp.limit(physical_integral, gamma, 0, dir="+"), 1)

    return {"integral_real_axis": integral, "integral_shat_positive": physical_integral, "kernel": delta_bw}


# Prove coherent continuum-resonance summation uses amplitude-level addition
def prove_coherent_amplitude_addition():
    ar, ai, br, bi = sp.symbols("a_r a_i b_r b_i", real=True)
    a = ar + sp.I * ai
    b = br + sp.I * bi

    coherent = sp.expand_complex(abs(a + b) ** 2)
    incoherent = sp.expand_complex(abs(a) ** 2 + abs(b) ** 2)
    interference = sp.simplify(coherent - incoherent)

    require_equal("coherent interference", interference, 2 * ar * br + 2 * ai * bi)

    return {"coherent": coherent, "interference": interference}


# Derive GDecay from the scalar-pair vertices and spin completeness
def prove_tensor_decay_coupling_widths():
    # [REFERENCE: Ewerz, Maniatis and Nachtmann, arXiv:1309.3478, Eqs. (3.35), (3.37), (5.6)]
    # [REFERENCE: PDG Kinematics, https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf, Eq. (49.18)]
    mass, gamma, br, beta, scale, symmetry = sp.symbols("M Gamma BR beta S0 symmetry", positive=True)
    partial_width = gamma * br
    momentum = mass * beta / 2
    # At the pole the form factors are one, choose the daughter momentum along z
    # iG_SSS = i g S0, iG_VPP = -i g (p1-p2)/2, iG_TPP = -i g (p1-p2)(p1-p2)/(2 S0)
    # Contract scalar daughter momenta with vector and spin-two completeness
    # src/Tensor/MTensorPomeron.cc: GDecay and the scalar-pair decay vertices
    theta, phi = sp.symbols("theta phi", real=True)
    p = momentum * sp.Matrix([sp.sin(theta) * sp.cos(phi), sp.sin(theta) * sp.sin(phi), sp.cos(theta)])
    v = -p
    tensor = -2 * p * p.T / scale
    tensor_norm = sum(tensor[i, j] ** 2 for i in range(3) for j in range(3)) - sp.trace(tensor) ** 2 / 3
    spin_sum = {0: scale**2, 1: (v.T * v)[0] / 3, 2: tensor_norm / 5}
    code_factors = {
        0: beta * scale**2 / (16 * sp.pi * mass),
        1: mass * beta**3 / (192 * sp.pi),
        2: mass**3 * beta**5 / (480 * sp.pi * scale**2),
    }
    couplings = {}
    for spin, amp2 in spin_sum.items():
        width_per_g2 = sp.simplify(momentum * amp2 / (8 * sp.pi * mass**2 * symmetry))
        require_equal(f"tensor J{spin} vertex width", width_per_g2, code_factors[spin] / symmetry)
        coupling2 = partial_width * symmetry / code_factors[spin]
        require_equal(f"tensor J{spin} width inversion", coupling2 * width_per_g2, partial_width)
        couplings[f"J{spin}"] = sp.sqrt(coupling2)
    # Identical scalar daughters permit even J only, with a 2! phase-space divisor
    for spin in (0, 2):
        require_equal(
            f"tensor J{spin} identical-pair coupling",
            couplings[f"J{spin}"].subs(symmetry, 2) ** 2,
            2 * couplings[f"J{spin}"].subs(symmetry, 1) ** 2,
        )
    return couplings


# Run all symbolic normalization checks
def main():
    checks = {
        "moller_flux": prove_moller_flux_conventions(),
        "factorized_F": prove_factorized_phase_space(),
        "continuum_C": prove_continuum_phase_space(),
        "quasielastic_Q": prove_quasielastic_phase_space(),
        "delta_bw": prove_delta_bw_normalization(),
        "coherent_addition": prove_coherent_amplitude_addition(),
        "tensor_decay_couplings": prove_tensor_decay_coupling_widths(),
    }
    for name, value in checks.items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
