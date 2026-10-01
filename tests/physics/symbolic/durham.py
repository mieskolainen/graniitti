# Symbolic checks for MDurham normalization, color factors, and projectors
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from functools import cache
from itertools import permutations
from pathlib import Path
from runpy import run_path

import numpy as np
import sympy as sp

# C++ conventions: MDurham, MSudakov, MDirac and generated MG5 hard amplitudes
# References determine the independent kernels, tests/cpp/test_durham.cc tests their implementation


# Require an exact symbolic expression to vanish
def require_zero(name, expr):
    value = sp.simplify(expr)
    if value != 0:
        raise AssertionError(f"{name}: expected 0, got {value}")


# Require two exact symbolic expressions to be equal
def require_equal(name, expr, expected):
    require_zero(name, sp.simplify(expr - expected))


# Require equality after exact trigonometric and exponential rewriting
def require_equal_trig(name, expr, expected):
    require_zero(name, sp.simplify(sp.expand_trig(expr - expected).rewrite(sp.exp)))


# Require an exact symbolic expression to be non-zero
def require_nonzero(name, expr):
    value = sp.simplify(expr)
    if value == 0:
        raise AssertionError(f"{name}: expected non-zero, got 0")


# Compare the production MG2GRA factor with an independently contracted SU(3) Gram matrix
def check_color_factor(gram, expected=None):
    generator = run_path(str(Path(__file__).resolve().parents[3] / "develop/MG2GRA/modules/mg5_color.py"))
    factor = np.asarray(generator["pivoted_cholesky"](np.array(gram, dtype=complex).tolist()), dtype=complex)
    if factor.size == 0:
        factor = factor.reshape(0, gram.cols)
    tolerance = 64 * np.finfo(float).eps
    np.testing.assert_allclose(factor.conj().T @ factor, np.array(gram, dtype=complex), rtol=tolerance, atol=tolerance)
    if expected is not None:
        np.testing.assert_allclose(factor, np.array(expected, dtype=complex), rtol=tolerance, atol=tolerance)
    return factor


# Compute the sign of epsilon^{0123} with all indices contravariant
def levi_civita4(a, b, c, d):
    index = (a, b, c, d)
    if len(set(index)) != 4:
        return 0
    inversions = sum(1 for i in range(len(index)) for j in range(i + 1, len(index)) if index[i] > index[j])
    return -1 if inversions % 2 else 1


# Compute one vector with the (+,-,-,-) metric applied
def lower_lorentz_components(v):
    return sp.Matrix([v[0], -v[1], -v[2], -v[3]])


# Compute the Minkowski dot product in the code Lorentz order
def lorentz_dot(left, right):
    return sp.simplify(left[0] * right[0] - left[1] * right[1] - left[2] * right[2] - left[3] * right[3])


# Compute one purely transverse four-vector in upper-index convention
def transverse_four_vector(qx, qy):
    return sp.Matrix([0, qx, qy, 0])


# Compute the MDirac rest-frame massive spin-1 basis
def spin1_rest_polarization(m):
    if m == -1:
        return sp.Matrix([0, 1 / sp.sqrt(2), -sp.I / sp.sqrt(2), 0])
    if m == 0:
        return sp.Matrix([0, 0, 0, 1])
    if m == 1:
        return sp.Matrix([0, -1 / sp.sqrt(2), -sp.I / sp.sqrt(2), 0])
    raise ValueError("spin1_rest_polarization expects m = -1, 0, +1")


# Compute the MDirac rest-frame massive spin-2 basis
def spin2_rest_polarization(m):
    from sympy.physics.wigner import clebsch_gordan

    if m not in (-2, -1, 0, 1, 2):
        raise ValueError("spin2_rest_polarization expects m = -2, -1, 0, +1, +2")

    basis = {h: spin1_rest_polarization(h) for h in (-1, 0, 1)}
    eps = sp.zeros(4, 4)
    for m1 in (-1, 0, 1):
        for m2 in (-1, 0, 1):
            coeff = clebsch_gordan(1, 1, 2, m1, m2, m)
            eps += coeff * basis[m1] * basis[m2].T
    return eps


# Compute the exact chi_c1 beam and polarization projector used by MDurham
def axial_beam_projector(p1, p2, eps):
    p1_lo = lower_lorentz_components(p1)
    p2_lo = lower_lorentz_components(p2)
    eps_lo = lower_lorentz_components(eps)
    projector = []
    for mu in range(4):
        value = 0
        for nu in range(4):
            for alpha in range(4):
                for beta in range(4):
                    value += levi_civita4(mu, nu, alpha, beta) * p1_lo[nu] * p2_lo[alpha] * sp.conjugate(eps_lo[beta])
        projector.append(sp.simplify(value))
    return sp.Matrix(projector)


# Compute the exact chi_c2 spin and beam projectors used by MDurham
def tensor_beam_projector(p1, p2, eps):
    p1_lo = lower_lorentz_components(p1)
    p2_lo = lower_lorentz_components(p2)
    spin = eps.conjugate()
    beam = 0
    for mu in range(4):
        for alpha in range(4):
            beam += p1_lo[mu] * p2_lo[alpha] * spin[mu, alpha]
    return spin, sp.simplify(beam)


# Compute a rank-2 tensor with both Lorentz indices lowered
def lower_tensor2(eps):
    metric_sign = (1, -1, -1, -1)
    out = sp.zeros(4, 4)
    for mu in range(4):
        for nu in range(4):
            out[mu, nu] = metric_sign[mu] * metric_sign[nu] * eps[mu, nu]
    return out


# Require two symbolic matrices to be exactly equal
def require_matrix_equal(name, left, right):
    if left.shape != right.shape:
        raise AssertionError(f"{name}: shape mismatch {left.shape} vs {right.shape}")
    for i in range(left.rows):
        for j in range(left.cols):
            require_zero(f"{name}[{i},{j}]", sp.simplify(left[i, j] - right[i, j]))


# Compute SU(3) generators T^a = lambda^a/2 in the fundamental representation
def su3_generators():
    imaginary_unit = sp.I
    half = sp.Rational(1, 2)
    inv_2sqrt3 = sp.Rational(1, 2) / sp.sqrt(3)
    inv_sqrt3 = 1 / sp.sqrt(3)

    return (
        sp.Matrix([[0, half, 0], [half, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, -imaginary_unit * half, 0], [imaginary_unit * half, 0, 0], [0, 0, 0]]),
        sp.Matrix([[half, 0, 0], [0, -half, 0], [0, 0, 0]]),
        sp.Matrix([[0, 0, half], [0, 0, 0], [half, 0, 0]]),
        sp.Matrix([[0, 0, -imaginary_unit * half], [0, 0, 0], [imaginary_unit * half, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, half], [0, half, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, -imaginary_unit * half], [0, imaginary_unit * half, 0]]),
        sp.Matrix([[inv_2sqrt3, 0, 0], [0, inv_2sqrt3, 0], [0, 0, -inv_sqrt3]]),
    )


T = su3_generators()


# Compute Tr(T^{i0} T^{i1} ...) with exact SU(3) matrices
def tr_word(indices):
    prod = sp.eye(3)
    for idx in indices:
        prod *= T[idx]
    return sp.simplify(sp.trace(prod))


# Compute a cached exact trace for repeated color words
@cache
def tr_word_cached(indices):
    """Return a cached exact trace for repeated color words"""
    return tr_word(indices)


# Compute exact SU(3) f^{abc} from [T^b,T^c] = i f^{abc} T^a
def fabc(a, b, c):
    comm = T[b] * T[c] - T[c] * T[b]
    return sp.simplify((2 / sp.I) * sp.trace(T[a] * comm))


# Compute exact SU(3) d^{abc} from {T^b,T^c} = d^{abc} T^a + delta^{bc}/N_C
def dabc(a, b, c):
    anti = T[b] * T[c] + T[c] * T[b]
    return sp.simplify(2 * sp.trace(T[a] * anti))


# Compute cached exact SU(3) f^{abc}
@cache
def fabc_cached(a, b, c):
    """Return cached exact SU(3) f^{abc}"""
    return fabc(a, b, c)


# Compute cached exact SU(3) d^{abc}
@cache
def dabc_cached(a, b, c):
    """Return cached exact SU(3) d^{abc}"""
    return dabc(a, b, c)


# Prove the Durham loop measure pi^2 int d^2Q has no extra factor two
def prove_loop_measure():
    q, qmin, qmax, phi = sp.symbols("q qmin qmax phi", positive=True)
    measure = sp.pi**2 * sp.integrate(q, (q, qmin, qmax)) * sp.integrate(1, (phi, 0, 2 * sp.pi))
    expected = sp.pi**3 * (qmax**2 - qmin**2)
    require_equal("Durham loop measure", measure, expected)
    return expected


# Prove the GL_IR radial map and periodic azimuthal weights used by MDurham
def prove_ir_resolved_loop_quadrature_map():
    unit, q, qmax = sp.symbols("u q Qmax", positive=True)
    k = sp.symbols("k", integer=True, nonnegative=True)
    nphi = sp.symbols("Nphi", integer=True, positive=True)

    qt = qmax * unit**2
    radial_weight = sp.diff(qt, unit)
    polar_jacobian = qt
    mapped_measure = sp.simplify(radial_weight * polar_jacobian)
    mapped_integral = sp.integrate(mapped_measure, (unit, 0, 1))
    direct_integral = sp.integrate(q, (q, 0, qmax))

    phi_node = 2 * sp.pi * k / nphi
    phi_weight = 2 * sp.pi / nphi
    periodic_weight_sum = sp.summation(phi_weight, (k, 0, nphi - 1))

    require_equal("GL_IR transformed radial measure", mapped_measure, 2 * qmax**2 * unit**3)
    require_equal("GL_IR radial measure integral", mapped_integral, direct_integral)
    require_equal("periodic phi weight sum", periodic_weight_sum, 2 * sp.pi)
    return {
        "qt": qt,
        "weight_times_jacobian": mapped_measure,
        "radial_integral": mapped_integral,
        "phi_node": phi_node,
        "phi_weight_sum": periodic_weight_sum,
    }


# Prove normalized CEP singlet contraction of delta^{AB}/N_C gives 1/N_C
def prove_meson_color_factor():
    NC = sp.symbols("NC", integer=True, positive=True)
    nadj = NC**2 - 1
    normalized = sp.simplify(nadj / (nadj * NC))
    unnormalized_sum = sp.simplify(nadj / NC)

    require_equal("normalized meson color factor", normalized, 1 / NC)
    require_equal("SU(3) normalized meson color factor", normalized.subs(NC, 3), sp.Rational(1, 3))
    require_equal("SU(3) unnormalized color sum", unnormalized_sum.subs(NC, 3), sp.Rational(8, 3))
    require_equal("old/new meson color ratio", unnormalized_sum / normalized, nadj)
    return normalized


# Prove the meson-pair prefactor (1/N_C) 64 pi^2 alpha_s^2 / shat
def prove_meson_pair_hard_kernel_normalization():
    NC, alpha_s, shat = sp.symbols("NC alpha_s shat", positive=True)
    gs2 = 4 * sp.pi * alpha_s

    # [REFERENCE: Harland-Lang et al., arXiv:1105.1626, Eqs. (3.12)-(3.13)]
    # This checks Dgg2MMbar's convention for the published prefactor, not its diagram derivation
    # The dimensionless T kernels carry the numerator algebra, with common 4 g_s^4 / shat
    perturbative_prefactor = 4 * gs2**2 / shat
    color_factor = 1 / NC
    norm = sp.simplify(color_factor * perturbative_prefactor)
    expected = (1 / NC) * 64 * sp.pi**2 * alpha_s**2 / shat

    require_equal("meson-pair perturbative prefactor", perturbative_prefactor, 64 * sp.pi**2 * alpha_s**2 / shat)
    require_equal("meson-pair full hard-kernel norm", norm, expected)
    return norm


# Prove the root decay coupling matches a unit-normalized delta-BW density
def prove_charmonium_root_decay_normalization():
    mass, width, branching = sp.symbols("M Gamma BR", positive=True)
    mass2 = sp.symbols("m2", real=True)

    delta_bw_density = mass * width / sp.pi / ((mass2 - mass**2) ** 2 + (mass * width) ** 2)
    delta_bw_integral = sp.integrate(delta_bw_density, (mass2, -sp.oo, sp.oo))
    decay_coupling2 = 2 * mass * width * branching
    decay_phase_space_factor = 1 / (2 * sp.pi)

    root_factor2 = sp.pi / (mass * width)
    normalized_branching = sp.simplify(delta_bw_integral * decay_phase_space_factor * decay_coupling2 * root_factor2)

    require_equal("delta-BW mass-squared normalization", delta_bw_integral, 1)
    require_equal("Durham root decay branching normalization", normalized_branching, branching)
    physical_integral = sp.Rational(1, 2) + sp.atan(mass / width) / sp.pi
    require_equal("root decay narrow-width physical normalization", sp.limit(physical_integral, width, 0, dir="+"), 1)
    return {"root_factor2": root_factor2, "physical_branching_constant_width": branching * physical_integral}


# Match the leading-power Lonnblad-Zlebcik Eq. 14 cross section
def prove_durham_bridge_from_lonnblad_cross_section():
    s, beta, m2, e1, e2, jac, bridge = sp.symbols("s beta M2 E1 E2 J C", positive=True)

    # [REFERENCE: Lonnblad and Zlebcik, arXiv:1608.03765, Eq. (14)]
    # Eq. (14) is a high-energy approximation, not an exact finite-energy cross section
    # Solving with exact phase space below only measures the kinematic difference
    # DQtloop stores T_qV = pi^2 int d^2q q_i V_ij ... .
    # Eq. 14 has coefficient (1/(64 pi^2))*(1/(2M^2)) per dlnM^2
    # for |int d^2q q_i V_ij ...|^2 dSigma_w
    loop_pi2 = sp.pi**2
    graniitti_per_dlnm2 = sp.simplify(
        m2 * jac * bridge**2 * loop_pi2**2 / (32 * sp.pi * beta * s * e1 * e2 * (2 * sp.pi) ** 5)
    )
    target_per_dlnm2 = 1 / (128 * sp.pi**2 * m2)
    solved_bridge = sp.solve(sp.Eq(graniitti_per_dlnm2, target_per_dlnm2), bridge)[0]
    solved_bridge_high_energy = sp.simplify(
        solved_bridge.subs({beta: 1, e1: sp.sqrt(s) / 2, e2: sp.sqrt(s) / 2, jac: sp.Rational(1, 2)})
    )

    require_equal("Lonnblad bridge squared equation", graniitti_per_dlnm2.subs(bridge, solved_bridge), target_per_dlnm2)
    require_equal("Lonnblad formal phase-space match", solved_bridge, sp.sqrt(8 * beta * s * e1 * e2 / jac) / m2)
    require_equal("Lonnblad solved positive high-energy bridge", solved_bridge_high_energy, 2 * s / m2)
    return {"formal_phase_space_match": solved_bridge, "high_energy": solved_bridge_high_energy}


# Prove chi_c kernels are stored as qV and converted by the common loop factor
def prove_charmonium_kernel_loop_convention():
    alpha_s, phi_prime2, mass, qt2 = sp.symbols("alpha_s phi_prime2 M Q_t2", positive=True)
    NC = sp.Integer(3)
    ncolor = NC**2 - 1

    c_chi = sp.simplify(
        16
        * sp.pi
        * alpha_s
        / (2 * sp.sqrt(NC) * (mass**2 / 2) ** 2)
        * sp.sqrt(6 / (4 * sp.pi * mass))
        * sp.sqrt(phi_prime2)
    )

    # For q1T=-q2T, the Minkowski product is +Q_t^2, fixing the amplitude sign
    # The coupling above uses the leading-power Q_t^2/M^2 -> 0 limit
    mbar_forward = sp.sqrt(sp.Rational(3, 2)) * c_chi * mass * qt2
    vertex = prove_chic0_superchic2_eq30_kernel()["forward"]
    require_equal("chi_c0 forward vertex sign", vertex.subs(sp.Symbol("c_chi"), c_chi), mbar_forward)
    # Dgg2chicJ uses the generated shat in DurhamMbarToQV, not the pole mass
    shat = sp.symbols("shat", positive=True)
    qv = shat * mbar_forward / 2
    require_equal("off-shell qV to Mbar", 2 * qv / shat, mbar_forward)
    qv_forward = qv.subs(shat, mass**2)
    loop_mbar_forward = sp.simplify(2 * qv_forward / mass**2)
    helicity_amp = sp.simplify(qv_forward / qt2)
    expected_helicity = sp.simplify(
        24 * sp.sqrt(sp.pi) * alpha_s * sp.sqrt(phi_prime2) / (sp.sqrt(NC) * mass ** sp.Rational(3, 2))
    )

    amp2 = sp.simplify(helicity_amp * sp.conjugate(helicity_amp))
    implied_width = sp.simplify(ncolor * amp2 / (16 * sp.pi * mass))
    expected_width = 96 * alpha_s**2 * phi_prime2 / mass**4

    require_equal("chi_c qV to Mbar roundtrip", loop_mbar_forward, mbar_forward)
    require_equal("chi_c Eq 30 qV gives the LO helicity amplitude", helicity_amp, expected_helicity)
    require_equal("chi_c Eq 30 common qV loop gives Eq 35 gg width", implied_width, expected_width)
    return {
        "LoopKernel": "QV",
        "qV_to_Mbar": loop_mbar_forward,
        "helicity_from_qV": helicity_amp,
        "Gamma_LO": implied_width,
    }


# Prove the off-shell scalar product used in the SuperChic chi_c coupling
def prove_charmonium_gluon_dot_kinematics():
    mass, q1_2, q2_2 = sp.symbols("M q1perp2 q2perp2", positive=True)
    qdot = sp.symbols("q1_dot_q2", positive=True)

    # q_i^2 = -|q_{i,t}|^2 and (q1 + q2)^2 = M^2
    on_shell_relation = sp.Eq(mass**2, -q1_2 - q2_2 + 2 * qdot)
    solved = sp.solve(on_shell_relation, qdot)[0]
    code = sp.Rational(1, 2) * (mass**2 + q1_2 + q2_2)

    require_equal("DurhamGluonDot", solved, code)
    return code


# Prove the SuperChic 2 Eq 33 charmonium coupling implemented by MDurham
def prove_charmonium_coupling_matches_superchic2_eq33():
    alpha_s, phi_prime2, mass, NC, qdot = sp.symbols("alpha_s phi_prime2 M N_C q1q2", positive=True)
    bw = sp.symbols("BW")

    code = (
        16 * sp.pi * alpha_s / (2 * sp.sqrt(NC) * qdot**2) * sp.sqrt(6 / (4 * sp.pi * mass)) * sp.sqrt(phi_prime2) * bw
    )
    # DurhamCharmoniumAlphaPhiPrime infers alpha_s |phi'_P| from chi_c0
    # [REFERENCE: Harland-Lang, Khoze and Ryskin, arXiv:1508.02718, Eqs. (33)-(35)]
    pole_mass, gamma_gg = sp.symbols("M0 Gamma_gg", positive=True)
    width = 96 * alpha_s**2 * phi_prime2 / pole_mass**4
    inferred_phi2 = sp.solve(sp.Eq(width, gamma_gg), phi_prime2)[0]
    code_from_width = 16 * sp.pi / (2 * sp.sqrt(NC) * qdot**2) * sp.sqrt(6 / (4 * sp.pi * mass))
    code_from_width *= pole_mass**2 * sp.sqrt(gamma_gg / 96) * bw
    require_equal("chi_c common pole width coupling", code.subs(phi_prime2, inferred_phi2), code_from_width)
    return code_from_width


# Prove the SuperChic 2 Eq 30 chi_c0 kernel implemented by MDurham
def prove_chic0_superchic2_eq30_kernel():
    q1x, q1y, q2x, q2y = sp.symbols("q1x q1y q2x q2y", real=True)
    mass = sp.symbols("M", positive=True)
    c_chi = sp.symbols("c_chi")

    q1 = transverse_four_vector(q1x, q1y)
    q2 = transverse_four_vector(q2x, q2y)
    q1t2 = lorentz_dot(q1, q1)
    q2t2 = lorentz_dot(q2, q2)
    q1q2t = lorentz_dot(q1, q2)

    code = sp.sqrt(sp.Rational(1, 6)) * c_chi / mass * (3 * mass**2 * q1q2t - q1q2t * (q1t2 + q2t2) - 2 * q1t2 * q2t2)
    q1e2 = q1x**2 + q1y**2
    q2e2 = q2x**2 + q2y**2
    q1e_dot_q2e = q1x * q2x + q1y * q2y
    superchic2 = (
        sp.sqrt(sp.Rational(1, 6))
        * c_chi
        / mass
        * (-3 * mass**2 * q1e_dot_q2e - q1e_dot_q2e * (q1e2 + q2e2) - 2 * q1e2 * q2e2)
    )

    require_equal("chi_c0 Eq 30 Minkowski form", code, superchic2)
    qt2 = sp.symbols("Q_t2", positive=True)
    forward = sp.simplify(code.subs({q1x: sp.sqrt(qt2), q1y: 0, q2x: -sp.sqrt(qt2), q2y: 0}))
    require_equal("chi_c0 forward Minkowski sign", forward, sp.sqrt(sp.Rational(3, 2)) * c_chi * mass * qt2)
    return {"minkowski": code, "euclidean": superchic2, "forward": forward}


# Prove the SuperChic 2 Eq 31 chi_c1 kernel implemented by MDurham
def prove_chic1_superchic2_eq31_kernel():
    # Arbitrary q_i here probe the projector algebra, physical rest events require q1T+q2T=0
    # Boosted physical spin states are tested in tests/cpp/test_durham.cc
    # [REFERENCE: Harland-Lang, Khoze and Ryskin, arXiv:1508.02718, Eq. (31)]
    q1x, q1y, q2x, q2y = sp.symbols("q1x q1y q2x q2y", real=True)
    shat_pp = sp.symbols("s", positive=True)
    c_chi = sp.symbols("c_chi")

    ebeam = sp.sqrt(shat_pp) / 2
    p1 = sp.Matrix([ebeam, 0, 0, ebeam])
    p2 = sp.Matrix([ebeam, 0, 0, -ebeam])
    q1 = transverse_four_vector(q1x, q1y)
    q2 = transverse_four_vector(q2x, q2y)
    q1_lo = lower_lorentz_components(q1)
    q2_lo = lower_lorentz_components(q2)
    qterm = q2_lo * lorentz_dot(q1, q1) - q1_lo * lorentz_dot(q2, q2)
    rotation = sp.symbols("rotation", real=True)
    rot = sp.rot_axis3(-rotation)[:2, :2]
    rotated_q = dict(zip((q1x, q1y, q2x, q2y), list(rot * q1[1:3, :]) + list(rot * q2[1:3, :]), strict=True))
    exchanged_qterm = qterm.subs({q1x: q2x, q1y: q2y, q2x: q1x, q2y: q1y}, simultaneous=True)

    out = {}
    for m in (-1, 0, 1):
        eps = spin1_rest_polarization(m)
        projector = axial_beam_projector(p1, p2, eps)
        contraction = sp.simplify(sum(qterm[mu] * projector[mu] for mu in range(4)))

        code = sp.simplify(-2 * sp.I * c_chi * contraction / shat_pp)
        require_equal_trig(
            f"chi_c1 rotation m={m}", code.subs(rotated_q, simultaneous=True), sp.exp(-sp.I * m * rotation) * code
        )
        exchanged = axial_beam_projector(p2, p1, eps)
        require_equal(
            f"chi_c1 beam exchange m={m}", sum(exchanged_qterm[mu] * exchanged[mu] for mu in range(4)), contraction
        )
        q1e2 = q1x**2 + q1y**2
        q2e2 = q2x**2 + q2y**2
        common_x = q1x * q2e2 - q2x * q1e2
        common_y = q1y * q2e2 - q2y * q1e2
        expected = {
            -1: -sp.sqrt(2) * sp.I * c_chi / 2 * (-common_y + sp.I * common_x),
            0: sp.Integer(0),
            1: -sp.sqrt(2) * sp.I * c_chi / 2 * (common_y + sp.I * common_x),
        }
        superchic2 = sp.simplify(expected[m])
        require_equal(f"chi_c1 Eq 31 m={m}", code, superchic2)

        forward = sp.simplify(code.subs({q2x: -q1x, q2y: -q1y}))
        mirrored = sp.simplify(forward.subs({q1x: -q1x, q1y: -q1y}))
        require_equal(f"chi_c1 forward oddness m={m}", forward + mirrored, 0)
        # The off-shell integrand is nonzero, forward suppression follows after the loop
        q, phi = sp.symbols("Q phi", real=True)
        forward_phi = sp.expand(forward.subs({q1x: q * sp.cos(phi), q1y: q * sp.sin(phi)}))
        require_zero(f"chi_c1 forward azimuthal integral m={m}", sp.integrate(forward_phi, (phi, 0, 2 * sp.pi)))
        if m != 0:
            require_nonzero(f"chi_c1 forward integrand m={m}", forward)
        out[m] = code

    require_equal("chi_c1 collinear longitudinal projector", out[0], 0)
    return out


# Prove the spin-2 basis used for chi_c2 is a physical massive tensor basis
def prove_spin2_rest_basis_properties():
    metric = sp.diag(1, -1, -1, -1)
    rest_p = sp.Matrix([sp.symbols("M", positive=True), 0, 0, 0])

    for m in (-2, -1, 0, 1, 2):
        eps = spin2_rest_polarization(m)
        require_matrix_equal(f"spin2 symmetry m={m}", eps, eps.T)
        trace = sp.simplify(sum(metric[mu, nu] * eps[mu, nu] for mu in range(4) for nu in range(4)))
        require_equal(f"spin2 tracelessness m={m}", trace, 0)
        for nu in range(4):
            trans = sp.simplify(sum(rest_p[mu] * metric[mu, mu] * eps[mu, nu] for mu in range(4)))
            require_equal(f"spin2 transversality m={m} nu={nu}", trans, 0)

    for m in (-2, -1, 0, 1, 2):
        eps_m = spin2_rest_polarization(m)
        eps_m_lo = lower_tensor2(eps_m)
        for n in (-2, -1, 0, 1, 2):
            eps_n = spin2_rest_polarization(n)
            overlap = 0
            for mu in range(4):
                for nu in range(4):
                    overlap += sp.conjugate(eps_m_lo[mu, nu]) * eps_n[mu, nu]
            expected = 1 if m == n else 0
            require_equal(f"spin2 basis overlap m={m} n={n}", sp.simplify(overlap), expected)
    return "spin-2 rest basis is symmetric, traceless, transverse, normalized"


# Prove the SuperChic 2 Eq 32 chi_c2 kernel implemented by MDurham
def prove_chic2_superchic2_eq32_kernel():
    # The m=0 zero is specific to this limit, not a general nonforward selection rule
    q1x, q1y, q2x, q2y = sp.symbols("q1x q1y q2x q2y", real=True)
    q, phi = sp.symbols("q phi", real=True)
    mass = sp.symbols("M", positive=True)
    shat_pp = sp.symbols("s", positive=True)
    c_chi = sp.symbols("c_chi")

    ebeam = sp.sqrt(shat_pp) / 2
    p1 = sp.Matrix([ebeam, 0, 0, ebeam])
    p2 = sp.Matrix([ebeam, 0, 0, -ebeam])
    q1 = transverse_four_vector(q1x, q1y)
    q2 = transverse_four_vector(q2x, q2y)
    q1_lo = lower_lorentz_components(q1)
    q2_lo = lower_lorentz_components(q2)
    q1q2t = lorentz_dot(q1, q2)
    rotation = sp.symbols("rotation", real=True)
    rot = sp.rot_axis3(-rotation)[:2, :2]
    rotated_q = dict(zip((q1x, q1y, q2x, q2y), list(rot * q1[1:3, :]) + list(rot * q2[1:3, :]), strict=True))

    out = {}
    for m in (-2, -1, 0, 1, 2):
        spin, beam = tensor_beam_projector(p1, p2, spin2_rest_polarization(m))
        q_contraction = 0
        for mu in range(4):
            for alpha in range(4):
                q_contraction += q1_lo[mu] * q2_lo[alpha] * spin[mu, alpha]

        code = sp.simplify(sp.sqrt(2) * c_chi * mass * (shat_pp * q_contraction + 2 * q1q2t * beam) / shat_pp)
        require_equal_trig(
            f"chi_c2 rotation m={m}", code.subs(rotated_q, simultaneous=True), sp.exp(-sp.I * m * rotation) * code
        )
        require_equal(
            f"chi_c2 beam exchange m={m}", code.subs({q1x: q2x, q1y: q2y, q2x: q1x, q2y: q1y}, simultaneous=True), code
        )
        expected = {
            -2: sp.sqrt(2) * mass * c_chi / 2 * (q1x * q2x + sp.I * q1x * q2y + sp.I * q1y * q2x - q1y * q2y),
            -1: sp.Integer(0),
            0: sp.Integer(0),
            1: sp.Integer(0),
            2: sp.sqrt(2) * mass * c_chi / 2 * (q1x * q2x - sp.I * q1x * q2y - sp.I * q1y * q2x - q1y * q2y),
        }
        superchic2 = sp.simplify(expected[m])
        require_equal(f"chi_c2 Eq 32 m={m}", code, superchic2)
        forward = sp.simplify(code.subs({q2x: -q1x, q2y: -q1y}))
        forward_phi = sp.simplify(forward.subs({q1x: q * sp.cos(phi), q1y: q * sp.sin(phi)}))
        forward_average = sp.integrate(forward_phi, (phi, 0, 2 * sp.pi))
        require_equal(f"chi_c2 forward azimuthal cancellation m={m}", forward_average, 0)
        out[m] = code

    return out


# Prove the meson wave functions cancel the hard kernel boundary factors
def prove_meson_open_grid_boundary_cancellation():
    # The octet kernel still has an integrable 1/r corner at x -> 0, y -> 1
    # Dgg2MMbar uses an open grid and hard angular cuts
    x, y = sp.symbols("x y", real=True)
    phi_x = x * (1 - x) * (1 - 2 * x) ** 2
    phi_y = y * (1 - y) * (1 - 2 * y) ** 2
    hard_boundary_factor = 1 / (x * y * (1 - x) * (1 - y))
    convolution_integrand = sp.cancel(phi_x * phi_y * hard_boundary_factor)
    expected = (1 - 2 * x) ** 2 * (1 - 2 * y) ** 2

    require_equal("meson boundary cancellation", convolution_integrand, expected)
    require_equal("meson x=0 boundary limit", sp.limit(convolution_integrand, x, 0), (1 - 2 * y) ** 2)
    require_equal("meson x=1 boundary limit", sp.limit(convolution_integrand, x, 1), (1 - 2 * y) ** 2)
    require_equal("meson y=0 boundary limit", sp.limit(convolution_integrand, y, 0), (1 - 2 * x) ** 2)
    require_equal("meson y=1 boundary limit", sp.limit(convolution_integrand, y, 1), (1 - 2 * x) ** 2)
    return convolution_integrand


# Prove the definite AP integrals used in MSudakov::AP_gg and MSudakov::AP_qg
def prove_sudakov_integrated_splitting_kernels():
    z, delta, ca, tr, nf = sp.symbols("z delta C_A T_R N_F", positive=True)
    upper = 1 - delta

    z_pgg = 2 * ca * (z**2 / (1 - z) + (1 - z) + z**2 * (1 - z))
    p_qg = tr * (z**2 + (1 - z) ** 2)

    integrated_gg = sp.integrate(z_pgg, (z, 0, upper))
    integrated_qg = nf * sp.integrate(p_qg, (z, 0, upper))

    expected_gg = 2 * ca * (-sp.log(delta) - (1 - delta) ** 2 * (3 * delta**2 - 2 * delta + 11) / 12)
    expected_qg = nf * tr * (-2 * delta**3 / 3 + delta**2 - delta + sp.Rational(2, 3))

    require_equal("Sudakov AP_gg integral", integrated_gg, expected_gg)
    require_equal("Sudakov AP_qg integral", integrated_qg, expected_qg)
    return {"AP_gg": expected_gg, "AP_qg": expected_qg}


# Prove Durham incoming color average with normalized final color projectors
def prove_colored_final_state_projection_tensors():
    incoming_gg = sp.Rational(1, 8)
    final_gg = sp.Rational(1, 1) / sp.sqrt(8)
    final_qqbar = sp.Rational(1, 1) / sp.sqrt(3)
    final_f = sp.Rational(1, 1) / sp.sqrt(24)
    final_d = sp.Rational(1, 1) / sp.sqrt(sp.Rational(40, 3))

    require_equal("gg Durham tensor normalization", incoming_gg * final_gg, sp.Rational(1, 8) / sp.sqrt(8))
    require_equal("qqbar singlet tensor normalization", incoming_gg * final_qqbar, sp.Rational(1, 8) / sp.sqrt(3))
    require_equal("ggg f singlet tensor normalization", incoming_gg * final_f, sp.Rational(1, 8) / sp.sqrt(24))
    require_equal(
        "ggg d singlet tensor normalization", incoming_gg * final_d, sp.Rational(1, 8) / sp.sqrt(sp.Rational(40, 3))
    )

    f_norm = 0
    d_norm = 0
    fd_overlap = 0
    for a in range(8):
        for b in range(8):
            for c in range(8):
                f = fabc_cached(a, b, c)
                d = dabc_cached(a, b, c)
                f_norm += f * f
                d_norm += d * d
                fd_overlap += f * d

    require_equal("fabc singlet norm", f_norm, 24)
    require_equal("dabc singlet norm", d_norm, sp.Rational(40, 3))
    require_equal("fabc dabc orthogonality", fd_overlap, 0)
    return {
        "gg": incoming_gg * final_gg,
        "qqbar": incoming_gg * final_qqbar,
        "ggg_f": incoming_gg * final_f,
        "ggg_d": incoming_gg * final_d,
    }


# Prove int_0^1 phi_CZ(x) dx = f_M/(2 sqrt(3))
def prove_cz_wave_function_normalization():
    x, fM = sp.symbols("x fM", positive=True)
    phi_cz = 5 * sp.sqrt(3) * fM * x * (1 - x) * (1 - 2 * x) ** 2
    integral = sp.integrate(phi_cz, (x, 0, 1))
    expected = fM / (2 * sp.sqrt(3))
    require_equal("CZ wave-function normalization", integral, expected)
    return expected


# Prove the scalar flavour-octet meson-pair hard kernels implemented by MDurham
def prove_meson_pair_scalar_octet_kernel():
    x, y, c2, NC = sp.symbols("x y c2 N_C", positive=True)
    CF = (NC**2 - 1) / (2 * NC)
    a = (1 - x) * (1 - y) + x * y
    b = (1 - x) * (1 - y) - x * y

    code_pp = 0
    code_mm = 0
    code_pm = (
        1
        / (x * y * (1 - x) * (1 - y))
        * (x * (1 - x) + y * (1 - y))
        / (a**2 - b**2 * c2)
        * (NC / 2)
        * (c2 - 2 * CF * a / NC)
    )
    paper_prefactor = 1 / (x * y * (1 - x) * (1 - y))
    paper_angular = (x * (1 - x) + y * (1 - y)) / (a**2 - b**2 * c2)
    paper_color = NC * c2 / 2 - CF * a
    expected_pm = paper_prefactor * paper_angular * paper_color

    radiation_zero = sp.simplify(code_pm.subs(c2, 2 * CF * a / NC))
    require_equal("scalar octet T+-", code_pm, expected_pm)
    require_equal("scalar octet radiation zero", radiation_zero, 0)
    require_equal("octet daughter exchange", code_pm.subs({x: y, y: x}, simultaneous=True), code_pm)
    require_equal("octet charge conjugation", code_pm.subs({x: 1 - x, y: 1 - y}, simultaneous=True), code_pm)
    return {"T++": code_pp, "T--": code_mm, "T+-": code_pm, "T-+": code_pm}


# Prove the scalar flavour-singlet eta0 eta0 ladder kernels implemented by MDurham
def prove_meson_pair_scalar_singlet_kernel():
    x, y, c2 = sp.symbols("x y c2", positive=True)
    denominator = x * y * (1 - x) * (1 - y)
    # [REFERENCE: Harland-Lang et al., arXiv:1105.1626, Eqs. (4.1)-(4.6)]
    # Dgg2MMbar: each eta0 carries (uubar+ddbar+ssbar)/sqrt(3)
    singlet = sp.ones(3, 1) / sp.sqrt(3)
    flavour_multiplicity = (singlet.T * sp.ones(3) * singlet)[0]

    one_flavour_pp = 1 / denominator * (1 + c2) / (1 - c2) ** 2
    one_flavour_pm = 1 / denominator * (1 + 3 * c2) / (2 * (1 - c2) ** 2)
    code_pp = flavour_multiplicity * one_flavour_pp
    code_pm = flavour_multiplicity * one_flavour_pm

    require_equal("scalar singlet flavour sum", flavour_multiplicity, 3)
    require_equal("scalar singlet angular ratio", code_pm / code_pp, (1 + 3 * c2) / (2 * (1 + c2)))
    for amplitude in (code_pp, code_pm):
        require_equal("singlet daughter exchange", amplitude.subs({x: y, y: x}, simultaneous=True), amplitude)
    return {"T++": code_pp, "T--": code_pp, "T+-": code_pm, "T-+": code_pm}


# Prove the eta and eta-prime octet-singlet decomposition used by Dgg2MMbar
def prove_eta_eta_prime_kernel_decomposition():
    theta8, theta1 = sp.symbols("theta8 theta1", real=True)
    phi8a, phi8b, phi0a, phi0b, t_base, t_ladder = sp.symbols("phi8a phi8b phi0a phi0b T_base T_ladder")

    mix = {"eta": {"8": sp.cos(theta8), "0": -sp.sin(theta1)}, "etap": {"8": sp.sin(theta8), "0": sp.cos(theta1)}}
    expected = {
        ("eta", "eta"): (
            sp.cos(theta8) ** 2 * phi8a * phi8b * t_base + sp.sin(theta1) ** 2 * phi0a * phi0b * (t_base + t_ladder)
        ),
        ("eta", "etap"): (
            sp.cos(theta8) * sp.sin(theta8) * phi8a * phi8b * t_base
            - sp.sin(theta1) * sp.cos(theta1) * phi0a * phi0b * (t_base + t_ladder)
        ),
        ("etap", "eta"): (
            sp.sin(theta8) * sp.cos(theta8) * phi8a * phi8b * t_base
            - sp.cos(theta1) * sp.sin(theta1) * phi0a * phi0b * (t_base + t_ladder)
        ),
        ("etap", "etap"): (
            sp.sin(theta8) ** 2 * phi8a * phi8b * t_base + sp.cos(theta1) ** 2 * phi0a * phi0b * (t_base + t_ladder)
        ),
    }

    out = {}
    for left in ("eta", "etap"):
        for right in ("eta", "etap"):
            octet_octet = mix[left]["8"] * mix[right]["8"] * phi8a * phi8b
            singlet_singlet = mix[left]["0"] * mix[right]["0"] * phi0a * phi0b
            code = (octet_octet + singlet_singlet) * t_base + singlet_singlet * t_ladder
            require_equal_trig(f"{left} {right} eta mixing decomposition", code, expected[(left, right)])

            omitted_base = sp.simplify(code - (octet_octet * t_base + singlet_singlet * t_ladder))
            require_equal_trig(f"{left} {right} ordinary singlet term", omitted_base, singlet_singlet * t_base)
            out[(left, right)] = sp.simplify(code)

    collapsed_product = (mix["eta"]["8"] * phi8a + mix["eta"]["0"] * phi0a) * (
        mix["eta"]["8"] * phi8b + mix["eta"]["0"] * phi0b
    )
    orthogonal_kernel_product = out[("eta", "eta")].subs({t_base: 1, t_ladder: 0})
    omitted_cross = sp.expand(collapsed_product - orthogonal_kernel_product)
    require_nonzero("eta collapsed wave-function cross terms", omitted_cross)
    return out


# Prove the meson-pair initial helicity content entering DHelProj
def prove_meson_pair_dhelproj_selection_rules():
    P0, M0, P2, M2, A8, A0, A2 = sp.symbols("P0 M0 P2 M2 A8 A0 A2")

    pp, mm, mp, pm = sp.symbols("A_pp A_mm A_mp A_pm")
    projected = P0 * (pp + mm) + M0 * (pp - mm) + P2 * mp + M2 * pm

    scalar_octet = projected.subs({pp: 0, mm: 0, mp: A8, pm: A8})
    scalar_singlet = projected.subs({pp: A0, mm: A0, mp: A2, pm: A2})

    require_equal("scalar octet has no Jz=0 term", sp.diff(scalar_octet, P0), 0)
    require_equal("scalar octet has no Jz=0- term", sp.diff(scalar_octet, M0), 0)
    require_equal("scalar singlet Jz=0 coefficient", sp.diff(scalar_singlet, P0), 2 * A0)
    require_equal("scalar singlet has no Jz=0- term", sp.diff(scalar_singlet, M0), 0)
    return {"octet": scalar_octet, "singlet": scalar_singlet}


# Compute HELAS vxxxxx spatial components for a massless vector
def helas_massless_vector(nsv, helicity, cos_theta, sin_theta, phi):
    h = sp.Integer(helicity)
    s = sp.Integer(nsv)
    sqh = 1 / sp.sqrt(2)
    return sp.Matrix(
        [
            (-h * sp.cos(phi) * cos_theta - sp.I * s * sp.sin(phi)) * sqh,
            (-h * sp.sin(phi) * cos_theta + sp.I * s * sp.cos(phi)) * sqh,
            h * sin_theta * sqh,
        ]
    )


# Compute CEP fixed helicity vector for an incoming gluon along +z
def cep_plus_beam_vector(helicity):
    h = sp.Integer(helicity)
    return sp.Matrix([-h, -sp.I, 0]) / sp.sqrt(2)


# Compute CEP fixed helicity vector for an incoming gluon along -z
def cep_minus_beam_vector(helicity):
    h = sp.Integer(helicity)
    return sp.Matrix([h, -sp.I, 0]) / sp.sqrt(2)


# Compute the transverse CEP vector for an incoming gluon along +z
def cep_plus_transverse_vector(helicity):
    h = sp.Integer(helicity)
    return sp.Matrix([-h, -sp.I]) / sp.sqrt(2)


# Compute the transverse CEP vector for an incoming gluon along -z
def cep_minus_transverse_vector(helicity):
    h = sp.Integer(helicity)
    return sp.Matrix([h, -sp.I]) / sp.sqrt(2)


# Require two symbolic vectors to be exactly equal
def require_vector_equal(name, left, right):
    for i, expr in enumerate(left - right):
        require_zero(f"{name}[{i}]", sp.expand_trig(sp.simplify(expr)).rewrite(sp.exp))


# Prove A_CEP = exp(i(lambda1 phi1 - lambda2 phi2)) A_HELAS for incoming gluons
def prove_incoming_helas_phase():
    phi1, phi2 = sp.symbols("phi1 phi2", real=True)
    helicities = (-1, 1)

    for h in helicities:
        plus_helas = helas_massless_vector(-1, h, 1, 0, phi1)
        plus_cep = cep_plus_beam_vector(h)
        require_vector_equal(f"HELAS +z phase h={h}", plus_helas, sp.exp(-sp.I * h * phi1) * plus_cep)

        minus_helas = helas_massless_vector(-1, h, -1, 0, phi2)
        minus_cep = cep_minus_beam_vector(h)
        require_vector_equal(f"HELAS -z phase h={h}", minus_helas, sp.exp(sp.I * h * phi2) * minus_cep)

    phase = {}
    for h1 in helicities:
        for h2 in helicities:
            phase[(h1, h2)] = sp.exp(sp.I * (h1 * phi1 - h2 * phi2))

    require_equal("MM incoming phase", phase[(-1, -1)], sp.exp(sp.I * (-phi1 + phi2)))
    require_equal("MP incoming phase", phase[(-1, 1)], sp.exp(-sp.I * (phi1 + phi2)))
    require_equal("PM incoming phase", phase[(1, -1)], sp.exp(sp.I * (phi1 + phi2)))
    require_equal("PP incoming phase", phase[(1, 1)], sp.exp(sp.I * (phi1 - phi2)))
    return phase


# Prove the MG5 transverse-source kernel and its collinear Durham limit
def prove_transverse_helicity_projection():
    q1x, q1y, q2x, q2y = sp.symbols("q1x q1y q2x q2y", real=True)
    q1 = sp.Matrix([q1x, q1y])
    q2 = sp.Matrix([q2x, q2y])

    upper = {
        helicity: sp.simplify((q1.T * cep_plus_transverse_vector(helicity).conjugate())[0]) for helicity in (-1, 1)
    }
    lower = {
        helicity: sp.simplify((q2.T * cep_minus_transverse_vector(helicity).conjugate())[0]) for helicity in (-1, 1)
    }

    amm, amp, apm, app = sp.symbols("A_mm A_mp A_pm A_pp")
    transverse_collinear = sp.expand(
        upper[-1] * lower[-1] * amm
        + upper[-1] * lower[1] * amp
        + upper[1] * lower[-1] * apm
        + upper[1] * lower[1] * app
    )
    projectors = durham_spin_projectors((q1x, q1y), (q2x, q2y))
    dhelproj = sp.expand(
        projectors["0+"] * (app + amm)
        + projectors["0-"] * (app - amm)
        + projectors["+2+"] * amp
        + projectors["-2+"] * apm
    )
    require_equal("transverse source projection collinear Durham limit", transverse_collinear, dhelproj)

    c1x, c1y, c2x, c2y = sp.symbols("c1x c1y c2x c2y")
    nodewise = sp.expand((q1x * c1x + q1y * c1y) * (q2x * c2x + q2y * c2y))
    tensor_form = q1x * q2x * c1x * c2x + q1x * q2y * c1x * c2y + q1y * q2x * c1y * c2x + q1y * q2y * c1y * c2y
    require_equal("integrated transverse tensor factorization", nodewise, tensor_form)
    return {
        "upper": upper,
        "lower": lower,
        "collinear_kernel": transverse_collinear,
        "tensor_factorization": tensor_form,
    }


# Prove the meson hard-frame azimuth is transported into the beam Jz basis
def prove_meson_beam_azimuth_transport():
    axis, phi, rotation = sp.symbols("axis phi rotation", real=True)
    p0, m0, p2, m2 = sp.symbols("P0 M0 P2 M2")
    app, amm, apm, amp = sp.symbols("A_pp A_mm A_pm A_mp")
    azimuth = phi + axis
    born = p0 * (app + amm) + m0 * (app - amm)
    born += p2 * amp * sp.exp(-2 * sp.I * azimuth)
    born += m2 * apm * sp.exp(2 * sp.I * azimuth)
    rotated = born.subs(
        {axis: axis + rotation, p2: p2 * sp.exp(2 * sp.I * rotation), m2: m2 * sp.exp(-2 * sp.I * rotation)},
        simultaneous=True,
    )
    require_equal("meson beam azimuth covariance", rotated, born)
    untransported = born.subs(axis, 0)
    rotated_untransported = untransported.subs(
        {p2: p2 * sp.exp(2 * sp.I * rotation), m2: m2 * sp.exp(-2 * sp.I * rotation)}, simultaneous=True
    )
    require_nonzero("missing meson beam azimuth", sp.simplify(rotated_untransported - untransported))
    return sp.simplify(rotated - born)


# Prove the covariant meson tensor reproduces the planar helicity amplitudes
def prove_meson_transverse_tensor():
    a0, a2, phi = sp.symbols("A0 A2 phi", real=True)
    energy = sp.symbols("E", positive=True)
    metric = sp.diag(1, -1, -1, -1)
    k1 = sp.Matrix([energy, 0, 0, energy])
    k2 = sp.Matrix([energy, 0, 0, -energy])
    direction = sp.Matrix([0, sp.cos(phi), sp.sin(phi), 0])
    kdot = (k1.T * metric * k2)[0]
    gperp = metric - metric * (k1 * k2.T + k2 * k1.T) * metric / kdot
    vertex = (a0 + a2) * gperp + 2 * a2 * metric * direction * direction.T * metric
    for component in vertex * k1:
        require_zero("meson Ward identity 1", component)
    for component in vertex * k2:
        require_zero("meson Ward identity 2", component)
    for h1 in (-1, 1):
        eps1 = sp.Matrix([0, -h1, -sp.I, 0]) / sp.sqrt(2)
        for h2 in (-1, 1):
            eps2 = sp.Matrix([0, h2, -sp.I, 0]) / sp.sqrt(2)
            value = (eps1.T * vertex * eps2)[0]
            expected = a0 if h1 == h2 else a2 * sp.exp(2 * sp.I * h1 * phi)
            require_equal_trig("meson tensor helicity", value, expected)

    x1, y1, x2, y2 = sp.symbols("q1x q1y q2x q2y", real=True)
    q1, q2 = sp.Matrix([0, x1, y1, 0]), sp.Matrix([0, x2, y2, 0])
    contracted = (q1.T * vertex * q2)[0]
    expected = -a0 * (x1 * x2 + y1 * y2) + a2 * (
        (x1 * x2 - y1 * y2) * sp.cos(2 * phi) + (x1 * y2 + y1 * x2) * sp.sin(2 * phi)
    )
    require_equal_trig("meson tensor Durham contraction", contracted, expected)
    return {"Ward": 0, "M_same": a0, "M_opposite": a2 * sp.exp(2 * sp.I * phi)}


# Compute the bilinear transverse tensor contraction left_i T_ij right_j
def tensor_contract(left, tensor, right):
    return sp.simplify((left.T * tensor * right)[0])


# Prove M_{+-}->exp(+2 i phi) M_{+-}, M_{-+}->exp(-2 i phi) M_{-+}
def prove_meson_pair_azimuthal_phases():
    phi = sp.symbols("phi", real=True)
    rotation = sp.Matrix([[sp.cos(phi), -sp.sin(phi)], [sp.sin(phi), sp.cos(phi)]])

    scalar_tensor_0 = sp.eye(2)
    helicity2_tensor_0 = sp.Matrix([[1, 0], [0, -1]])
    scalar_tensor_phi = rotation * scalar_tensor_0 * rotation.T
    helicity2_tensor_phi = rotation * helicity2_tensor_0 * rotation.T

    helicities = (-1, 1)
    phase = {}
    for h1 in helicities:
        for h2 in helicities:
            e1 = cep_plus_transverse_vector(h1)
            e2 = cep_minus_transverse_vector(h2)
            key = (h1, h2)

            scalar_plane = tensor_contract(e1, scalar_tensor_0, e2)
            scalar_azimuth = tensor_contract(e1, scalar_tensor_phi, e2)
            helicity2_plane = tensor_contract(e1, helicity2_tensor_0, e2)
            helicity2_azimuth = tensor_contract(e1, helicity2_tensor_phi, e2)

            if h1 == h2:
                require_nonzero(f"scalar meson-pair same-helicity base {key}", scalar_plane)
                require_equal_trig(f"scalar meson-pair phase {key}", scalar_azimuth, scalar_plane)
                require_equal_trig(f"helicity-2 same-helicity zero {key}", helicity2_azimuth, 0)
                phase[key] = sp.Integer(1)
            else:
                require_equal_trig(f"scalar meson-pair opposite-helicity zero {key}", scalar_azimuth, 0)
                require_nonzero(f"helicity-2 meson-pair base {key}", helicity2_plane)
                expected = helicity2_plane * sp.exp(sp.I * (h1 - h2) * phi)
                require_equal_trig(f"helicity-2 meson-pair phase {key}", helicity2_azimuth, expected)
                phase[key] = sp.exp(sp.I * (h1 - h2) * phi)

    require_equal("meson-pair PP phase", phase[(1, 1)], 1)
    require_equal("meson-pair MM phase", phase[(-1, -1)], 1)
    require_equal("meson-pair PM phase", phase[(1, -1)], sp.exp(2 * sp.I * phi))
    require_equal("meson-pair MP phase", phase[(-1, 1)], sp.exp(-2 * sp.I * phi))
    return phase


# Compute Durham spin projectors S_0+, S_0-, S_+2+, S_-2+
def durham_spin_projectors(q1, q2, signed_pseudoscalar=True):
    q1x, q1y = q1
    q2x, q2y = q2
    half = sp.Rational(1, 2)
    cross = q1x * q2y - q1y * q2x
    s0m = -sp.I * half * cross if signed_pseudoscalar else -sp.I * half * sp.Abs(cross)

    return {
        "0+": -half * (q1x * q2x + q1y * q2y),
        "0-": s0m,
        "+2+": half * ((q1x * q2x - q1y * q2y) + sp.I * (q1x * q2y + q1y * q2x)),
        "-2+": half * ((q1x * q2x - q1y * q2y) - sp.I * (q1x * q2y + q1y * q2x)),
    }


# Compute the helicity-amplitude coefficients inside each J_z^P sector
def durham_helicity_terms(projectors):
    return {
        "0+": {"++": projectors["0+"], "--": projectors["0+"]},
        "0-": {"++": projectors["0-"], "--": -projectors["0-"]},
        "+2+": {"-+": projectors["+2+"]},
        "-2+": {"+-": projectors["-2+"]},
    }


# Compute S_0+(A++ + A--) + S_0-(A++ - A--) + S_+2 A-+ + S_-2 A+-
def durham_projected_amplitude(projectors):
    app, amm, amp, apm = sp.symbols("A_pp A_mm A_mp A_pm")
    return sp.simplify(
        projectors["0+"] * (app + amm)
        + projectors["0-"] * (app - amm)
        + projectors["+2+"] * amp
        + projectors["-2+"] * apm
    )


# Prove the implemented q2 = -Q - p2 convention and isolate the M0 abs effect
def prove_loop_momentum_new_convention_and_m0_abs():
    qx, qy, p1x, p1y, p2x, p2y = sp.symbols("qx qy p1x p1y p2x p2y", real=True)
    q1 = (qx - p1x, qy - p1y)
    q2_kmr = (-qx - p2x, -qy - p2y)

    require_equal("KMR q2 denominator", q2_kmr[0] ** 2 + q2_kmr[1] ** 2, (qx + p2x) ** 2 + (qy + p2y) ** 2)

    kmr_signed = durham_spin_projectors(q1, q2_kmr, signed_pseudoscalar=True)
    kmr_abs = durham_spin_projectors(q1, q2_kmr, signed_pseudoscalar=False)
    for qn in ("0+", "+2+", "-2+"):
        require_equal(f"signed and abs projector agree for {qn}", kmr_signed[qn], kmr_abs[qn])

    q, phi = sp.symbols("q phi", real=True, positive=True)
    q1_forward = (q * sp.cos(phi), q * sp.sin(phi))
    q2_forward_kmr = (-q * sp.cos(phi), -q * sp.sin(phi))
    kmr_forward = durham_spin_projectors(q1_forward, q2_forward_kmr, signed_pseudoscalar=True)
    expected_forward = {
        "0+": q**2 / 2,
        "0-": 0,
        "+2+": -(q**2) * sp.exp(2 * sp.I * phi) / 2,
        "-2+": -(q**2) * sp.exp(-2 * sp.I * phi) / 2,
    }
    for qn, expected in expected_forward.items():
        require_equal_trig(f"KMR forward limit {qn}", kmr_forward[qn], expected)

    forward_terms = durham_helicity_terms(kmr_forward)
    expected_terms = {
        "0+:++": q**2 / 2,
        "0+:--": q**2 / 2,
        "0-:++": 0,
        "0-:--": 0,
        "+2+:-+": -(q**2) * sp.exp(2 * sp.I * phi) / 2,
        "-2+:+-": -(q**2) * sp.exp(-2 * sp.I * phi) / 2,
    }
    for qn, helicities in forward_terms.items():
        for helicity, value in helicities.items():
            require_equal_trig(f"KMR forward helicity term {qn} {helicity}", value, expected_terms[f"{qn}:{helicity}"])

    a, b, aodd = sp.symbols("a b A_odd", positive=True)
    m0_signed_pos = durham_spin_projectors((a, 0), (0, b), signed_pseudoscalar=True)["0-"]
    m0_signed_neg = durham_spin_projectors((a, 0), (0, -b), signed_pseudoscalar=True)["0-"]
    m0_abs_pos = durham_spin_projectors((a, 0), (0, b), signed_pseudoscalar=False)["0-"]
    m0_abs_neg = durham_spin_projectors((a, 0), (0, -b), signed_pseudoscalar=False)["0-"]

    require_equal("M0 abs equals signed for positive orientation", m0_abs_pos, m0_signed_pos)
    require_equal("M0 abs flips signed negative orientation", m0_abs_neg, -m0_signed_neg)
    require_equal("signed M0 mirrored orientations cancel", m0_signed_pos + m0_signed_neg, 0)
    require_equal("abs M0 mirrored orientations are equal", m0_abs_pos - m0_abs_neg, 0)
    require_equal("signed 0- odd-amplitude mirrored loop cancellation", aodd * (m0_signed_pos + m0_signed_neg), 0)
    require_equal("abs 0- odd-amplitude mirrored loop equality", aodd * (m0_abs_pos - m0_abs_neg), 0)

    return {
        "q1": "Q - p1",
        "q2": "-Q - p2",
        "forward_kmr": expected_forward,
        "forward_helicity_terms": expected_terms,
        "M0_abs": "same for mirrored orientations, signed M0 is opposite",
    }


# Prove gg -> gg Durham-color-averaged projection coefficients in the MG trace basis
def prove_gg_singlet_coefficients():
    kmr_incoming = 1 / sp.sqrt(8)
    perms = ((0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3), (0, 2, 3, 1), (0, 3, 1, 2), (0, 3, 2, 1))
    expected = (
        kmr_incoming * sp.Rational(2, 3),
        kmr_incoming * sp.Rational(2, 3),
        -kmr_incoming * sp.Rational(1, 12),
        kmr_incoming * sp.Rational(2, 3),
        -kmr_incoming * sp.Rational(1, 12),
        kmr_incoming * sp.Rational(2, 3),
    )

    coeffs = []
    for perm in perms:
        coeff = 0
        for a in range(8):
            for c in range(8):
                colors = (a, a, c, c)
                coeff += tr_word_cached(tuple(colors[i] for i in perm)) / (8 * sp.sqrt(8))
        coeffs.append(sp.simplify(coeff))

    for idx, (value, target) in enumerate(zip(coeffs, expected, strict=True)):
        require_equal(f"gg singlet coefficient {idx}", value, target)
    return coeffs


# Prove gg -> qqbar Durham-color-averaged projection coefficients in the MG color basis
def prove_qqbar_singlet_coefficients():
    coeff = sum(tr_word_cached((a, a)) for a in range(8)) / (8 * sp.sqrt(3))
    expected = 2 / (sp.sqrt(6) * sp.sqrt(8))
    require_equal("qqbar singlet coefficient", coeff, expected)
    return (sp.simplify(coeff), sp.simplify(coeff))


# Prove six normalized-incoming gg -> qqbar-g coefficients in native MG5 order
def prove_qqbarg_singlet_coefficients():
    # [REFERENCE: Del Duca, Dixon and Maltoni, arXiv:hep-ph/9910563]
    incoming_norm = sp.sqrt(8)
    final_norm = sp.sqrt(sp.Rational(1, 2) * 8)
    expected = (
        4 / (3 * sp.sqrt(2)),
        -1 / (6 * sp.sqrt(2)),
        4 / (3 * sp.sqrt(2)),
        -1 / (6 * sp.sqrt(2)),
        4 / (3 * sp.sqrt(2)),
        4 / (3 * sp.sqrt(2)),
    )

    # MG5 sorts the six permutations lexicographically before delta_ab contraction
    orderings = ("aac", "aca", "aac", "aca", "caa", "caa")
    coeffs = []
    for ordering in orderings:
        coeff = 0
        for a in range(8):
            for c in range(8):
                ordered = {"aac": (a, a, c), "caa": (c, a, a), "aca": (a, c, a)}[ordering]
                coeff += tr_word_cached((c,) + ordered) / (incoming_norm * final_norm)
        coeffs.append(sp.simplify(coeff))

    for idx, (value, target) in enumerate(zip(coeffs, expected, strict=True)):
        require_equal(f"qqbar-g singlet coefficient {idx}", value, target)

    final_tensor_norm = sum(sp.trace(T[c] * T[c]) for c in range(8)) / final_norm**2
    require_equal("qqbar-g final singlet norm", final_tensor_norm, 1)
    return coeffs


# Prove gg -> ggg Durham-color-averaged f_abc and d_abc coefficients in the MG trace basis
def prove_ggg_singlet_coefficients():
    kmr_incoming = 1 / sp.sqrt(8)
    inv_sqrt3 = kmr_incoming / sp.sqrt(3)
    inv8_sqrt3 = inv_sqrt3 / 8
    d_leading = kmr_incoming * sp.sqrt(sp.Rational(5, 27))
    d_sublead = d_leading / 8

    expected_f = (
        inv_sqrt3,
        -inv_sqrt3,
        -inv_sqrt3,
        inv_sqrt3,
        inv_sqrt3,
        -inv_sqrt3,
        -inv8_sqrt3,
        inv8_sqrt3,
        -inv8_sqrt3,
        inv_sqrt3,
        inv8_sqrt3,
        -inv_sqrt3,
        inv8_sqrt3,
        -inv8_sqrt3,
        inv8_sqrt3,
        -inv_sqrt3,
        -inv8_sqrt3,
        inv_sqrt3,
        -inv8_sqrt3,
        inv8_sqrt3,
        -inv8_sqrt3,
        inv_sqrt3,
        inv8_sqrt3,
        -inv_sqrt3,
    )
    expected_d = (
        d_leading,
        d_leading,
        d_leading,
        d_leading,
        d_leading,
        d_leading,
        -d_sublead,
        -d_sublead,
        -d_sublead,
        d_leading,
        -d_sublead,
        d_leading,
        -d_sublead,
        -d_sublead,
        -d_sublead,
        d_leading,
        -d_sublead,
        d_leading,
        -d_sublead,
        -d_sublead,
        -d_sublead,
        d_leading,
        -d_sublead,
        d_leading,
    )

    coeff_f = []
    coeff_d = []
    for tail in permutations((1, 2, 3, 4)):
        perm = (0,) + tail
        cf = 0
        cd = 0
        for a in range(8):
            for b in range(8):
                for c in range(8):
                    for d in range(8):
                        color = (a, a, b, c, d)
                        word = tuple(color[i] for i in perm)
                        trace = tr_word_cached(word)
                        cf += trace * fabc_cached(b, c, d) / (8 * sp.sqrt(24))
                        cd += trace * dabc_cached(b, c, d) / (8 * sp.sqrt(sp.Rational(40, 3)))
        coeff_f.append(sp.simplify(-sp.I * cf))
        coeff_d.append(sp.simplify(cd))

    for idx, (value, target) in enumerate(zip(coeff_f, expected_f, strict=True)):
        require_equal(f"ggg f-singlet coefficient {idx}", value, target)
    for idx, (value, target) in enumerate(zip(coeff_d, expected_d, strict=True)):
        require_equal(f"ggg d-singlet coefficient {idx}", value, target)
    return coeff_f, coeff_d


# Prove the generated finite-Nc MG5 Gram projector against exact SU(3) tensors
def prove_mg5_su3_gram_projection(qqbarg_coefficients, ggg_coefficients):
    # [REFERENCE: Alwall et al., JHEP 07 (2014) 079, arXiv:1405.0301]
    complex_gram = sp.Matrix([[1, sp.I], [-sp.I, 1]])
    check_color_factor(complex_gram, sp.Matrix([[1, sp.I]]))
    identity = sp.eye(3)
    casimir = sum((generator * generator for generator in T), sp.zeros(3))
    require_matrix_equal("SU(3) fundamental Casimir", casimir, sp.Rational(4, 3) * identity)

    # Prove the finite-Nc Fierz map on a basis spanning all complex 3 x 3 matrices
    for row in range(3):
        for column in range(3):
            unit = sp.zeros(3)
            unit[row, column] = 1
            sandwich = sum((generator * unit * generator for generator in T), sp.zeros(3))
            expected = sp.trace(unit) * identity / 2 - unit / 6
            require_matrix_equal(f"SU(3) Fierz map E{row}{column}", sandwich, expected)

    # Fix the trace, f_abc and d_abc phases used by the MG5 trace basis
    f_norm = 0
    d_norm = 0
    fd_overlap = 0
    for a in range(8):
        for b in range(8):
            for c in range(8):
                f = fabc_cached(a, b, c)
                d = dabc_cached(a, b, c)
                require_equal(f"SU(3) three-generator trace {a}{b}{c}", tr_word_cached((a, b, c)), (d + sp.I * f) / 4)
                f_norm += f * f
                d_norm += d * d
                fd_overlap += f * d
    require_equal("normalized ggg f tensor", f_norm, 24)
    require_equal("normalized ggg d tensor", d_norm, sp.Rational(40, 3))
    require_equal("orthogonal ggg f and d tensors", fd_overlap, 0)

    # The corrected six-flow q qbar g order collapses onto one normalized T^c/2
    qqbarg_orderings = ("aac", "aca", "aac", "aca", "caa", "caa")
    qqbarg_projector = sp.Matrix([qqbarg_coefficients])
    for index, ordering in enumerate(qqbarg_orderings):
        fierz_weight = sp.Rational(-1, 6) if ordering == "aca" else sp.Rational(4, 3)
        for color in range(8):
            direct = fierz_weight * T[color] / sp.sqrt(8)
            reconstructed = qqbarg_projector[0, index] * T[color] / 2
            require_matrix_equal(f"qqbar-g restricted color tensor {index},{color}", direct, reconstructed)

    qqbarg_gram = qqbarg_projector.conjugate().T * qqbarg_projector
    qqbarg_factor = check_color_factor(qqbarg_gram, qqbarg_projector)

    # Remove the later KMR 1/sqrt(8) average to recover delta_ab/sqrt(8)
    averaged_f, averaged_d = ggg_coefficients
    physical_projectors = sp.Matrix(
        [[sp.sqrt(8) * value for value in averaged_f], [sp.sqrt(8) * value for value in averaged_d]]
    )

    # Fierz reduction proves every restricted five-generator trace lies in f,d
    trace_permutations = [(0,) + tail for tail in permutations((1, 2, 3, 4))]
    for index, permutation in enumerate(trace_permutations):
        incoming_position = permutation.index(1)
        fierz_weight = sp.Rational(4, 3) if incoming_position in (1, 4) else sp.Rational(-1, 6)
        outgoing = [entry - 2 for entry in permutation if entry >= 2]
        inversions = sum(1 for left in range(3) for right in range(left + 1, 3) if outgoing[left] > outgoing[right])
        antisymmetric_sign = -1 if inversions % 2 else 1
        expected_f = antisymmetric_sign * fierz_weight * sp.sqrt(3) / 4
        expected_d = fierz_weight * sp.sqrt(sp.Rational(5, 3)) / 4
        require_equal(f"ggg restricted f coefficient {index}", physical_projectors[0, index], expected_f)
        require_equal(f"ggg restricted d coefficient {index}", physical_projectors[1, index], expected_d)

    # Since i f/sqrt(24) and d/sqrt(40/3) are orthonormal, this is the direct Gram
    ggg_gram = physical_projectors.conjugate().T * physical_projectors
    generated_factor = check_color_factor(ggg_gram)
    if generated_factor.shape != (2, 24):
        raise AssertionError(f"ggg incoming-singlet Gram rank is not two: {generated_factor.shape}")

    root42 = sp.sqrt(42)
    root210 = sp.sqrt(210)
    generated_cpp = sp.Matrix(
        [
            [
                root42 / 9,
                -2 * root42 / 63,
                -2 * root42 / 63,
                root42 / 9,
                root42 / 9,
                -2 * root42 / 63,
                -root42 / 72,
                root42 / 252,
                -root42 / 72,
                root42 / 9,
                root42 / 252,
                -2 * root42 / 63,
                root42 / 252,
                -root42 / 72,
                root42 / 252,
                -2 * root42 / 63,
                -root42 / 72,
                root42 / 9,
                -root42 / 72,
                root42 / 252,
                -root42 / 72,
                root42 / 9,
                root42 / 252,
                -2 * root42 / 63,
            ],
            [
                0,
                root210 / 21,
                root210 / 21,
                0,
                0,
                root210 / 21,
                0,
                -root210 / 168,
                0,
                0,
                -root210 / 168,
                root210 / 21,
                -root210 / 168,
                0,
                -root210 / 168,
                root210 / 21,
                0,
                0,
                0,
                -root210 / 168,
                0,
                0,
                -root210 / 168,
                root210 / 21,
            ],
        ]
    )
    check_color_factor(ggg_gram, generated_cpp)

    return {
        "qqbarg_rank": qqbarg_factor.shape[0],
        "ggg_rank": generated_factor.shape[0],
        "norm": "J^dagger G J = ||P J||^2",
    }


# Prove cyclic ggg shower color-flow tags are not an exact finite-N_C basis
def prove_ggg_color_flow_is_not_finite_nc_orthogonal():
    overlap = 0
    norm = 0
    for b in range(8):
        for c in range(8):
            for d in range(8):
                left = tr_word_cached((b, c, d))
                right = tr_word_cached((b, d, c))
                overlap += left * sp.conjugate(right)
                norm += left * sp.conjugate(left)

    overlap = sp.simplify(overlap)
    norm = sp.simplify(norm)
    if overlap == 0:
        raise AssertionError("ggg cyclic color flows unexpectedly orthogonal")
    require_equal("ggg trace basis norm", norm, sp.Rational(7, 3))
    return overlap


# Run all exact symbolic Durham checks
def main():
    checks = {
        "loop_measure": prove_loop_measure(),
        "gl_ir_loop_map": prove_ir_resolved_loop_quadrature_map(),
        "meson_color": prove_meson_color_factor(),
        "meson_pair_norm": prove_meson_pair_hard_kernel_normalization(),
        "charmonium_root_decay_normalization": prove_charmonium_root_decay_normalization(),
        "lonnblad_cross_section_bridge": prove_durham_bridge_from_lonnblad_cross_section(),
        "charmonium_kernel_convention": prove_charmonium_kernel_loop_convention(),
        "charmonium_gluon_dot": prove_charmonium_gluon_dot_kinematics(),
        "charmonium_coupling_eq33": prove_charmonium_coupling_matches_superchic2_eq33(),
        "chic0_eq30": prove_chic0_superchic2_eq30_kernel(),
        "chic1_eq31": prove_chic1_superchic2_eq31_kernel(),
        "spin2_basis": prove_spin2_rest_basis_properties(),
        "chic2_eq32": prove_chic2_superchic2_eq32_kernel(),
        "meson_open_grid_boundaries": prove_meson_open_grid_boundary_cancellation(),
        "sudakov_ap_kernels": prove_sudakov_integrated_splitting_kernels(),
        "colored_projectors": prove_colored_final_state_projection_tensors(),
        "cz_norm": prove_cz_wave_function_normalization(),
        "meson_octet_kernel": prove_meson_pair_scalar_octet_kernel(),
        "meson_singlet_kernel": prove_meson_pair_scalar_singlet_kernel(),
        "eta_kernel_decomposition": prove_eta_eta_prime_kernel_decomposition(),
        "meson_dhelproj_selection": prove_meson_pair_dhelproj_selection_rules(),
        "incoming_helas_phase": prove_incoming_helas_phase(),
        "transverse_projection": prove_transverse_helicity_projection(),
        "meson_pair_phase": prove_meson_pair_azimuthal_phases(),
        "meson_beam_azimuth": prove_meson_beam_azimuth_transport(),
        "meson_transverse_tensor": prove_meson_transverse_tensor(),
        "loop_q2_convention": prove_loop_momentum_new_convention_and_m0_abs(),
        "gg_color": prove_gg_singlet_coefficients(),
        "qqbar_color": prove_qqbar_singlet_coefficients(),
    }
    qqbarg_color = prove_qqbarg_singlet_coefficients()
    ggg_color = prove_ggg_singlet_coefficients()
    checks.update(
        {
            "qqbarg_color": qqbarg_color,
            "ggg_color": ggg_color,
            "mg5_su3_gram_projection": prove_mg5_su3_gram_projection(qqbarg_color, ggg_color),
            "ggg_flow_overlap": prove_ggg_color_flow_is_not_finite_nc_orthogonal(),
        }
    )
    for name, value in checks.items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
