# Symbolic normalization checks for the kT-EPA photon flux
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import sympy as sp

# The Jacobians are exact, while the transverse EPA and massless hard flux are approximations
# Physical emitters have 0 < xi < 1 and MX >= mp
# Corresponding implementation: src/Photon/MFlux.cc and src/Particle/MForm.cc


# Require two symbolic expressions to be exactly equal
def require_equal(name, left, right):
    diff = left - right
    if isinstance(diff, sp.MatrixBase):
        diff = diff.applyfunc(sp.simplify)
        is_zero = diff == sp.zeros(*diff.shape)
    else:
        diff = sp.simplify(diff)
        is_zero = diff == 0
    if not is_zero:
        raise AssertionError(f"{name}: {left} != {right}; diff={diff}")


# Prove the exact light-cone Q2 and Bjorken-x emitter identities
def prove_exact_emitter_kinematics():
    pplus, xi, pt2, mp2, mx2 = sp.symbols("pplus xi pt2 mp2 mx2", positive=True)

    pminus = mp2 / pplus
    pfinal_plus = (1 - xi) * pplus
    pfinal_minus = (mx2 + pt2) / pfinal_plus
    qplus = pplus - pfinal_plus
    qminus = pminus - pfinal_minus
    q2 = sp.simplify(qplus * qminus - pt2)
    photon_q2 = sp.simplify(-q2)
    two_p_dot_q = sp.simplify(pplus * qminus + pminus * qplus)
    xbj = sp.simplify(q2 / two_p_dot_q)

    expected_q2 = (pt2 + xi * (mx2 - mp2) + xi**2 * mp2) / (1 - xi)
    expected_xbj = photon_q2 / (photon_q2 + mx2 - mp2)
    dz_q2_min = xi**2 * mp2 / (1 - xi)

    require_equal("exact emitter Q2", photon_q2, expected_q2)
    require_equal("exact emitter Bjorken-x", xbj, expected_xbj)
    require_equal("elastic emitter Bjorken-x", xbj.subs(mx2, mp2), 1)
    require_equal("Drees-Zeppenfeld elastic Q2min", photon_q2.subs({mx2: mp2, pt2: 0}), dz_q2_min)
    return {"Q2": photon_q2, "xBj": xbj, "DZ_Q2min": dz_q2_min}


# Prove the exact MFactorized 2->3 skeleton density before photon fluxes
def prove_mfactorized_exact_density():
    s, beta, e1, e2, jac = sp.symbols("s beta E1 E2 J", positive=True)

    # Lorentz measure after momentum deltas, dpXz/(2 EX) = dy/2
    # src/Kinematics/MFactorized.cc: B51PhaseSpaceWeight
    dphi3 = sp.simplify((2 * sp.pi) ** 4 / (2 * sp.pi) ** 9 * jac / (2 * e1 * 2 * e2 * 2))
    flux = 2 * s * beta
    factorized = sp.simplify(dphi3 / flux / (2 * sp.pi))

    require_equal("exact MFactorized dPhi3", dphi3, jac / (256 * sp.pi**5 * e1 * e2))
    require_equal("exact MFactorized dSigma skeleton", factorized, jac / (1024 * sp.pi**6 * beta * e1 * e2 * s))
    return {"dPhi3": dphi3, "dSigma_skeleton": factorized}


# Prove that the high-energy kT-EPA flux factor per proton is 16 pi^2
def prove_kt_epa_required_flux_factor():
    s, x1, x2, c = sp.symbols("s x1 x2 C", positive=True)
    f1, f2, amp2 = sp.symbols("F1 F2 A2", positive=True)

    # F_i is the Budnev density per dx_i dqt_i^2, so the d^2qt density is F_i/pi
    # Lorentz measure, central dM^2/(2pi), incident 2s, E1 E2 -> s/4, J -> 1/2
    high_energy_skeleton = sp.Rational(1, 2) / (8 * (s / 4) * (2 * sp.pi) ** 5) / (2 * sp.pi) / (2 * s)
    code_density = sp.simplify(high_energy_skeleton * c**2 * f1 * f2 * amp2 / (x1 * x2))

    shat = x1 * x2 * s
    jac_x = 1 / s
    subprocess_density = sp.simplify((f1 / sp.pi) * (f2 / sp.pi) * jac_x * amp2 / (2 * shat))

    solved = sp.solve(sp.Eq(code_density, subprocess_density), c)[0]
    require_equal("kT-EPA required per-leg code factor", solved, 16 * sp.pi**2)
    require_equal("kT-EPA high-energy density match", code_density.subs(c, 16 * sp.pi**2), subprocess_density)
    return solved


# Prove the exact light-cone Jacobian and EPA phase-space conversion
def prove_exact_epa_phase_space_factor():
    s, shat, beta, e1, e2, velocity = sp.symbols("s shat beta E1 E2 Delta_v", positive=True)
    x1, x2 = sp.symbols("x1 x2", positive=True)
    beam_plus, beam_minus = sp.symbols("P1plus P2minus", positive=True)
    mt1sq, mt2sq = sp.symbols("mT1sq mT2sq", positive=True)
    u0, v0 = sp.symbols("u0 v0", positive=True)

    xplus = u0 + beam_plus * x1 - mt2sq / (beam_minus * (1 - x2))
    xminus = v0 - mt1sq / (beam_plus * (1 - x1)) + beam_minus * x2
    a, b, a1, a2, b1, b2 = sp.symbols("A B A1 A2 B1 B2", nonzero=True)
    generic_determinant = sp.det(
        sp.Matrix([[a1 * b + a * b1, a2 * b + a * b2], [(a1 / a - b1 / b) / 2, (a2 / a - b2 / b) / 2]])
    )
    require_equal("mass-rapidity chain-rule determinant", generic_determinant, a2 * b1 - a1 * b2)

    determinant = sp.simplify(sp.diff(xplus, x2) * sp.diff(xminus, x1) - sp.diff(xplus, x1) * sp.diff(xminus, x2))
    expected_determinant = beam_plus * beam_minus - mt1sq * mt2sq / (
        beam_plus * beam_minus * (1 - x1) ** 2 * (1 - x2) ** 2
    )
    # MFlux::ExactktEPAPhaseSpaceFactor uses the positive measure on either branch
    dx_jacobian = 1 / sp.Abs(expected_determinant)

    exact_pp_skeleton = 1 / (1024 * sp.pi**6 * velocity * beta * e1 * e2 * s)
    phase_space_factor = sp.simplify(2 * s * beta * e1 * e2 * velocity * dx_jacobian / shat)
    photon_skeleton = dx_jacobian / (512 * sp.pi**6 * shat)
    small_x_high_energy = sp.simplify(
        phase_space_factor.subs(
            {
                beta: 1,
                e1: sp.sqrt(s) / 2,
                e2: sp.sqrt(s) / 2,
                velocity: 2,
                beam_plus: sp.sqrt(s),
                beam_minus: sp.sqrt(s),
                mt1sq: 0,
                mt2sq: 0,
                shat: x1 * x2 * s,
            },
            simultaneous=True,
        )
    )
    finite_x_high_energy = sp.simplify(
        phase_space_factor.subs(
            {
                beta: 1,
                e1: (1 - x1) * sp.sqrt(s) / 2,
                e2: (1 - x2) * sp.sqrt(s) / 2,
                velocity: 2,
                beam_plus: sp.sqrt(s),
                beam_minus: sp.sqrt(s),
                mt1sq: 0,
                mt2sq: 0,
                shat: x1 * x2 * s,
            },
            simultaneous=True,
        )
    )

    require_equal("exact oriented light-cone Jacobian determinant", determinant, -expected_determinant)
    require_equal("absolute light-cone Jacobian", dx_jacobian, 1 / sp.Abs(determinant))
    for transverse_mass, expected in ((sp.Rational(1, 2), sp.Rational(1, 3)), (sp.Rational(3, 2), sp.Rational(1, 5))):
        require_equal(
            "light-cone Jacobian branch",
            dx_jacobian.subs(
                {
                    beam_plus: 2,
                    beam_minus: 2,
                    x1: sp.Rational(1, 2),
                    x2: sp.Rational(1, 2),
                    mt1sq: transverse_mass,
                    mt2sq: transverse_mass,
                }
            ),
            expected,
        )
    require_equal("exact EPA phase-space density", exact_pp_skeleton * phase_space_factor, photon_skeleton)
    require_equal("exact EPA small-x high-energy limit", small_x_high_energy, 1 / (x1 * x2))
    require_equal("exact EPA finite-x high-energy limit", finite_x_high_energy, (1 - x1) * (1 - x2) / (x1 * x2))
    return {
        "dx1dx2_dM2dY": dx_jacobian,
        "factor": phase_space_factor,
        "small_x_high_energy": small_x_high_energy,
        "finite_x_high_energy": finite_x_high_energy,
    }


# Prove the common exact DZ and LUX hard-flux conversion
def prove_collinear_photon_flux_convention():
    s, shat, beta, amp2 = sp.symbols("s shat beta A2", positive=True)
    mass2 = sp.symbols("M2", positive=True)
    rapidity = sp.symbols("Y", real=True)
    photon_pdf1, photon_pdf2 = sp.symbols("f1 f2", positive=True)
    x1 = sp.sqrt(mass2 / s) * sp.exp(rapidity)
    x2 = sp.sqrt(mass2 / s) * sp.exp(-rapidity)

    dx_jacobian = -sp.det(
        sp.Matrix([[sp.diff(x1, mass2), sp.diff(x1, rapidity)], [sp.diff(x2, mass2), sp.diff(x2, rapidity)]])
    )
    proton_moller_flux = 2 * s * beta
    hard_moller_flux = 2 * shat
    hard_flux_correction = s * beta / shat
    code_density = sp.simplify(
        dx_jacobian * photon_pdf1 * photon_pdf2 * amp2 * hard_flux_correction / proton_moller_flux
    )
    factorized_density = sp.simplify(dx_jacobian * photon_pdf1 * photon_pdf2 * amp2 / hard_moller_flux)
    collinear_limit = sp.simplify(hard_flux_correction.subs(shat, x1 * x2 * s))

    require_equal("collinear dx1 dx2 Jacobian", dx_jacobian, 1 / s)
    require_equal("collinear DZ/LUX hard-flux conversion", code_density, factorized_density)
    require_equal("collinear finite-energy correction", collinear_limit, beta / (x1 * x2))
    require_equal("collinear high-energy correction", collinear_limit.subs(beta, 1), 1 / (x1 * x2))
    return {"dx1dx2_dM2dY": dx_jacobian, "hard_flux_correction": hard_flux_correction, "finite_energy": collinear_limit}


# Prove the cancellation-free Drees-Zeppenfeld integral bracket
def prove_drees_zeppenfeld_stable_bracket():
    a, w, u = sp.symbols("A w u", positive=True)
    published = sp.log(a) - sp.Rational(11, 6) + 3 / a - sp.Rational(3, 2) / a**2 + 1 / (3 * a**3)
    stable = -sp.log(1 - w) - w - w**2 / 2 - w**3 / 3
    published_mapped = sp.simplify(published.subs(a, sp.exp(u)))
    stable_mapped = sp.simplify(stable.subs(w, 1 - sp.exp(-u)))

    require_equal("DZ stable integral bracket", published_mapped, stable_mapped)
    require_equal("DZ boundary leading order", sp.limit(stable / w**4, w, 0), sp.Rational(1, 4))
    return {"published": published, "stable": stable}


# Transform the elastic Budnev density from dQ2 to dqt2
def prove_coherent_flux_against_standard_elastic_formula():
    # [REFERENCE: Budnev et al., Phys. Rept. 15 (1975) 181, elastic equivalent-photon spectrum]
    alpha, xi, pt2, mp2, fel, fm = sp.symbols("alpha xi pt2 mp2 F_E F_M", positive=True)
    q2 = sp.symbols("Q2", positive=True)
    den = pt2 + xi**2 * mp2
    delta = pt2 / den

    # MFlux::ElasticSpinHalfFluxTransverse after removing the 16 pi^2 bridge
    # Each magnetic transverse eigenstate carries xi^2/4, the trace carries xi^2/2
    common = alpha / (sp.pi * xi * den)
    code_electric = common * (1 - xi) * delta * fel
    code_magnetic = 2 * common * xi**2 * fm / 4
    virtuality = den / (1 - xi)
    qmin2 = xi**2 * mp2 / (1 - xi)
    electric_dq2 = alpha / (sp.pi * xi * q2) * (1 - xi) * (1 - qmin2 / q2) * fel
    magnetic_dq2 = alpha / (sp.pi * xi * q2) * xi**2 * fm / 2
    standard_electric = electric_dq2.subs(q2, virtuality) * sp.diff(virtuality, pt2)
    standard_magnetic = magnetic_dq2.subs(q2, virtuality) * sp.diff(virtuality, pt2)

    code = sp.simplify(code_electric + code_magnetic)
    standard = sp.simplify(standard_electric + standard_magnetic)
    difference = sp.simplify(code - standard)

    require_equal("CohFlux electric term", code_electric, standard_electric)
    require_equal("CohFlux magnetic term", code_magnetic, standard_magnetic)
    require_equal("CohFlux standard formula", code, standard)
    require_equal("elastic magnetic forward limit", code.subs(pt2, 0), alpha * fm / (2 * sp.pi * xi * mp2))
    return {"code": code, "standard": standard, "code_minus_standard": difference, "magnetic_ratio": 1}


# Prove the full transverse density and generalized-EPA qT projection
def prove_transverse_flux_density_status():
    electric, magnetic = sp.symbols("E M", positive=True)
    theta, delta = sp.symbols("theta delta", real=True)

    density = sp.Matrix([[electric + magnetic, 0], [0, magnetic]])
    parallel_projector = sp.Matrix([[1, 0], [0, 0]])
    perpendicular_projector = sp.Matrix([[0, 0], [0, 1]])
    parallel_flux = sp.trace(density * parallel_projector)
    perpendicular_flux = sp.trace(density * perpendicular_projector)

    require_equal("transverse scalar flux trace", sp.trace(density), electric + 2 * magnetic)
    require_equal("generalized-EPA qT projection", parallel_flux, electric + magnetic)
    require_equal("perpendicular magnetic projection", perpendicular_flux, magnetic)
    require_equal(
        "transverse eigenmode reconstruction",
        parallel_flux * parallel_projector + perpendicular_flux * perpendicular_projector,
        density,
    )
    parallel_helas = sp.Matrix([sp.cos(theta), sp.exp(sp.I * delta) * sp.sin(theta)])
    perpendicular_helas = sp.Matrix([-sp.exp(-sp.I * delta) * sp.sin(theta), sp.cos(theta)])
    helas_basis = sp.Matrix.hstack(parallel_helas, perpendicular_helas)
    mg5_sources = helas_basis * sp.diag(sp.sqrt(2 * parallel_flux), sp.sqrt(2 * perpendicular_flux))
    helas_density = helas_basis * sp.diag(parallel_flux, perpendicular_flux) * helas_basis.conjugate().T

    require_equal("MG5 HELAS parallel source norm", (parallel_helas.conjugate().T * parallel_helas)[0], 1)
    require_equal(
        "MG5 HELAS perpendicular source norm", (perpendicular_helas.conjugate().T * perpendicular_helas)[0], 1
    )
    require_equal("MG5 HELAS source orthogonality", (parallel_helas.conjugate().T * perpendicular_helas)[0], 0)
    require_equal("MG5 photon spin-average density", mg5_sources * mg5_sources.conjugate().T / 2, helas_density)
    return {
        "density": density,
        "parallel": parallel_flux,
        "perpendicular": perpendicular_flux,
        "trace": sp.trace(density),
        "helas_density": helas_density,
    }


# Prove the exact finite-target F1, F2, FL and R identities
def prove_inelastic_callan_gross_status():
    f2, ratio, target_mass = sp.symbols("F2 R rho", positive=True)

    # MForm::F1xQ2 uses rho = 4 mp^2 xBj^2/Q^2 and R = sigma_L/sigma_T
    longitudinal, transverse = sp.symbols("FL FT")
    solution = sp.solve(
        [longitudinal + transverse - (1 + target_mass) * f2, longitudinal - ratio * transverse],
        (longitudinal, transverse),
    )
    fl, two_xf1 = solution[longitudinal], solution[transverse]

    require_equal("inelastic longitudinal relation", fl, (1 + target_mass) * f2 - two_xf1)
    require_equal("inelastic R relation", fl / two_xf1, ratio)
    require_equal("inelastic Bjorken Callan-Gross limit", two_xf1.subs({ratio: 0, target_mass: 0}), f2)
    return {"FL": fl, "2xF1": two_xf1, "R": sp.simplify(fl / two_xf1)}


# Prove IncohFlux reduces to CohFlux with elastic structure functions
def prove_incoherent_flux_elastic_limit():
    xi, pt2, mp2, q2, mx2, coeff = sp.symbols("xi pt2 mp2 Q2 MX2 c", positive=True)
    fel, fm = sp.symbols("F_E F_M", positive=True)

    denom = q2 + mx2 - mp2
    xbj = q2 / denom
    delta = pt2 / (pt2 + xi * (mx2 - mp2) + xi**2 * mp2)

    elastic_subs = {mx2: mp2}
    denom_el = sp.simplify(denom.subs(elastic_subs))
    xbj_el = sp.simplify(xbj.subs(elastic_subs))
    delta_el = sp.simplify(delta.subs(elastic_subs))

    # d(1 - xbj) / dMX2 at MX2 = mp2 gives delta(1 - xbj) = Q2 delta(MX2 - mp2)
    delta_jacobian = 1 / sp.diff(1 - xbj, mx2).subs(mx2, mp2)
    f2_over_denom_el = sp.simplify(fel * delta_jacobian / denom_el)
    two_xbj_f1_over_denom_el = sp.simplify(2 * xbj_el * (fm / 2) * delta_jacobian / denom_el)

    incoh_bracket = sp.simplify(
        (1 - xi) * delta_el**2 * f2_over_denom_el + coeff * xi**2 / xbj_el**2 * delta_el * two_xbj_f1_over_denom_el
    )
    coherent_bracket = sp.simplify((1 - xi) * delta_el**2 * fel + xi**2 * delta_el * fm / 2)

    solved_coeff = sp.solve(sp.Eq(incoh_bracket, coherent_bracket), coeff)[0]
    code_bracket = sp.simplify(incoh_bracket.subs(coeff, sp.Rational(1, 2)))

    require_equal("IncohFlux elastic denominator", denom_el, q2)
    require_equal("IncohFlux elastic Bjorken-x", xbj_el, 1)
    require_equal("IncohFlux elastic F2 projection", f2_over_denom_el, fel)
    require_equal("IncohFlux elastic F1 projection", two_xbj_f1_over_denom_el, fm)
    require_equal("IncohFlux required F1 coefficient", solved_coeff, sp.Rational(1, 2))
    require_equal("IncohFlux elastic limit", code_bracket, coherent_bracket)
    return {
        "denom_elastic": denom_el,
        "xbj_elastic": xbj_el,
        "required_f1_coefficient": solved_coeff,
        "code_minus_coherent": sp.simplify(code_bracket - coherent_bracket),
    }


# Run all exact symbolic kT-EPA checks
def main():
    checks = {
        "exact_emitter_kinematics": prove_exact_emitter_kinematics(),
        "mfactorized_exact": prove_mfactorized_exact_density(),
        "required_flux_factor": prove_kt_epa_required_flux_factor(),
        "exact_phase_space": prove_exact_epa_phase_space_factor(),
        "collinear_flux": prove_collinear_photon_flux_convention(),
        "drees_zeppenfeld": prove_drees_zeppenfeld_stable_bracket(),
        "coherent_flux_standard": prove_coherent_flux_against_standard_elastic_formula(),
        "transverse_flux_density": prove_transverse_flux_density_status(),
        "inelastic_callan_gross": prove_inelastic_callan_gross_status(),
        "incoherent_elastic_limit": prove_incoherent_flux_elastic_limit(),
    }
    for name, value in checks.items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
