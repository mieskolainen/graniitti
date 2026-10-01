# Symbolic Durham to GP glueball and meson helicity checks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import sympy as sp
from sympy.physics.wigner import clebsch_gordan


# Require an exact symbolic expression to vanish
def require_zero(name, expr):
    value = sp.simplify(sp.expand_trig(expr).rewrite(sp.exp))
    if value != 0:
        raise AssertionError(f"{name}: expected 0, got {value}")


# Require two exact symbolic expressions to be equal
def require_equal(name, expr, expected):
    require_zero(name, sp.simplify(expr - expected))


# Require equality after exact trigonometric simplification
def require_equal_trig(name, expr, expected):
    value = sp.trigsimp(sp.expand_trig(expr - expected))
    require_zero(name, value)


# Require an exact symbolic expression to be non-zero
def require_nonzero(name, expr):
    value = sp.simplify(expr)
    if value == 0:
        raise AssertionError(f"{name}: expected non-zero, got 0")


# Require a helicity dictionary to have unit norm
def require_unit_norm(name, amps):
    norm2 = helicity_norm2(amps)
    require_equal(name, norm2, 1)


# Compute the exact forward Durham spin projectors
def durham_projectors(q, phi):
    q1 = (q * sp.cos(phi), q * sp.sin(phi))
    q2 = (-q * sp.cos(phi), -q * sp.sin(phi))
    q1x, q1y = q1
    q2x, q2y = q2
    half = sp.Rational(1, 2)
    return {
        "0+": -half * (q1x * q2x + q1y * q2y),
        "0-": -sp.I * half * (q1x * q2y - q1y * q2x),
        "+2+": half * ((q1x * q2x - q1y * q2y) + sp.I * (q1x * q2y + q1y * q2x)),
        "-2+": half * ((q1x * q2x - q1y * q2y) - sp.I * (q1x * q2y + q1y * q2x)),
    }


# Contract hard gg helicity amplitudes with the Durham projectors
def durham_projected_amplitude(hard, q, phi):
    projectors = durham_projectors(q, phi)
    return sp.simplify(
        projectors["0+"] * (hard.get("++", 0) + hard.get("--", 0))
        + projectors["0-"] * (hard.get("++", 0) - hard.get("--", 0))
        + projectors["+2+"] * hard.get("-+", 0)
        + projectors["-2+"] * hard.get("+-", 0)
    )


# Compute the analytic Regge-helicity nonsense zero used by MReggeGP
def regge_helicity_nonsense_zero(alpha_t, m):
    factor = sp.Integer(1)
    for k in range(abs(m)):
        factor *= alpha_t - k
    return sp.simplify(factor)


# Compute the combined no-flip analytic Regge source vertex
def regge_residue(m, q, phi, alpha_t, s0, second_exchange_daughter=False, use_exchange_helicity_barrier=True):
    sign = -1 if second_exchange_daughter else 1
    barrier = (q / sp.sqrt(s0)) ** abs(m) if use_exchange_helicity_barrier else 1
    return sp.simplify(barrier * regge_helicity_nonsense_zero(alpha_t, m) * sp.exp(sp.I * sign * m * phi))


# Map a perturbative transverse gluon helicity into the GP m label
def gp_label_from_gluon_helicity(helicity):
    return -helicity


# Compute the forward GP no-flip vertex product for one gg helicity pair
def regge_forward_kernel(lambda1, lambda2, q, phi, alpha1, alpha2, s0):
    # src/Regge/MReggeGP.cc: ReggeFactors, NonsenseZero and gpom::Residue
    # This only compares azimuthal tensors, not the soft and perturbative dynamics
    m1 = gp_label_from_gluon_helicity(lambda1)
    m2 = gp_label_from_gluon_helicity(lambda2)
    phi1 = phi
    phi2 = phi + sp.pi
    return sp.simplify(regge_residue(m1, q, phi1, alpha1, s0, False) * regge_residue(m2, q, phi2, alpha2, s0, True))


# Prove the transverse GP source vertices reproduce Durham projectors
def prove_regge_residue_matches_durham_projectors():
    q, alpha1, alpha2, s0, phi = sp.symbols("q alpha1 alpha2 s0 phi", positive=True, real=True)
    projectors = durham_projectors(q, phi)
    expected = {
        "--": -2 * alpha1 * alpha2 * projectors["0+"] / s0,
        "++": -2 * alpha1 * alpha2 * projectors["0+"] / s0,
        "-+": 2 * alpha1 * alpha2 * projectors["+2+"] / s0,
        "+-": 2 * alpha1 * alpha2 * projectors["-2+"] / s0,
    }
    helicities = {"--": (-1, -1), "++": (1, 1), "-+": (-1, 1), "+-": (1, -1)}

    for label, (lambda1, lambda2) in helicities.items():
        kernel = regge_forward_kernel(lambda1, lambda2, q, phi, alpha1, alpha2, s0)
        require_equal(f"Regge residue bridge {label}", kernel, expected[label])

    alpha = sp.symbols("alpha", real=True)
    require_equal("Regge m=0 nonsense zero", regge_helicity_nonsense_zero(alpha, 0), 1)
    require_equal("Regge |m|=1 nonsense zero", regge_helicity_nonsense_zero(alpha, 1), alpha)
    require_equal("Regge |m|=2 nonsense zero", regge_helicity_nonsense_zero(alpha, 2), alpha * (alpha - 1))
    return expected


# Compute the raw GP scalar LS operator for a spin-two Pomeron pole
def scalar_regge_ls00_full(mmax=2, alpha=1):
    # MReggeGP::InitResonanceLS includes RawPoleLSNormalization
    # The rank-two scalar contraction has magnitude one, versus CG = 1/sqrt(5)
    pole_norm = 1 / clebsch_gordan(2, 2, 0, 2, -2, 0)
    rows = {(m, m): pole_norm * clebsch_gordan(alpha, alpha, 0, m, -m, 0) for m in range(-mmax, mmax + 1)}
    # Only integer trajectory values are tested here, not the Gamma continuation off pole
    # MWigner::W3jRegge vanishes at Gamma poles for |m| > 1 at alpha = 1
    return {key: value for key, value in rows.items() if value != 0}


# Compute a parity-even same-helicity gg ansatz, compatible with J=0
def scalar_durham_transverse_full():
    amp = sp.sqrt(sp.Rational(1, 2))
    return {(-1, -1): amp, (1, 1): amp}


# Compute a parity-even opposite-helicity gg ansatz, requiring J>=2
def tensor_durham_helicity2_full():
    amp = sp.sqrt(sp.Rational(1, 2))
    return {(-1, 1): amp, (1, -1): amp}


# Compute a possible J=2 helicity-zero ansatz, its strength is not predicted
def tensor_durham_helicity0_full():
    amp = sp.sqrt(sp.Rational(1, 2))
    return {(-1, -1): amp, (1, 1): amp}


# Compute the reduced color-singlet MHV gg -> gg hard helicity matrix
def gg_pair_spinor_helicity_matrix():
    # [REFERENCE: Harland-Lang, Khoze and Ryskin, arXiv:1508.02718, Eqs. (38)-(39)]
    shat, that, uhat, cgg = sp.symbols("s_hat t_hat u_hat C_gg", nonzero=True)
    phi = sp.symbols("phi", real=True)
    aux01 = shat**2 / (that * uhat)
    aux23 = uhat / that
    # In the Jacob-Wick basis the azimuthal phase is exp(i (lambda1-lambda2) phi)
    posphase = sp.exp(2 * sp.I * phi)
    negphase = sp.exp(-2 * sp.I * phi)
    return {
        "--": {"--": cgg * aux01},
        "++": {"++": cgg * aux01},
        "-+": {"-+": cgg * aux23 * negphase, "+-": cgg * posphase / aux23},
        "+-": {"-+": cgg * negphase / aux23, "+-": cgg * aux23 * posphase},
    }


# Compute the planar massless color-singlet gg -> q qbar matrix
def qqbar_pair_spinor_helicity_matrix():
    # [REFERENCE: Harland-Lang, Khoze and Ryskin, arXiv:1508.02718, Eq. (37)]
    theta = sp.symbols("theta_pair", positive=True, real=True)
    kq = sp.symbols("K_q", nonzero=True)
    x = sp.cot(theta / 2)
    y = sp.tan(theta / 2)
    return {"-+": {"-+": -kq * y, "+-": kq * x}, "+-": {"-+": kq * x, "+-": -kq * y}}


# Project a final two-parton wavefunction onto initial gg helicities
def project_initial_hard_matrix(matrix, final_wavefunction):
    out = {}
    for final_label, final_weight in final_wavefunction.items():
        for initial_label, value in matrix.get(final_label, {}).items():
            out[initial_label] = sp.simplify(out.get(initial_label, 0) + sp.conjugate(final_weight) * value)
    return {key: sp.simplify(value) for key, value in out.items() if sp.simplify(value) != 0}


# Compute fixed-angle hard coefficients for chosen final helicity superpositions without a bound-state projection
def partonic_pair_kernel_coefficients():
    theta = sp.symbols("theta_pair", positive=True, real=True)
    cgg, kq = sp.symbols("C_gg K_q", nonzero=True)
    s_hat, t_hat, u_hat = sp.symbols("s_hat t_hat u_hat", nonzero=True)
    gg_matrix = gg_pair_spinor_helicity_matrix()
    qqbar_matrix = qqbar_pair_spinor_helicity_matrix()
    phi = sp.symbols("phi", real=True)

    gg_scalar = project_initial_hard_matrix(
        gg_matrix, {"--": sp.sqrt(sp.Rational(1, 2)), "++": sp.sqrt(sp.Rational(1, 2))}
    )["--"]
    gg_tensor = project_initial_hard_matrix(
        gg_matrix, {"-+": sp.sqrt(sp.Rational(1, 2)), "+-": sp.sqrt(sp.Rational(1, 2))}
    )["-+"].subs(phi, 0)
    qqbar_tensor = project_initial_hard_matrix(
        qqbar_matrix, {"-+": sp.sqrt(sp.Rational(1, 2)), "+-": sp.sqrt(sp.Rational(1, 2))}
    )["-+"]

    massless_subs = {t_hat: -s_hat * (1 - sp.cos(theta)) / 2, u_hat: -s_hat * (1 + sp.cos(theta)) / 2}
    return {
        "gg_same": sp.trigsimp(gg_scalar.subs(massless_subs)),
        "gg_opposite": sp.trigsimp(gg_tensor.subs(massless_subs)),
        "qqbar_opposite": sp.trigsimp(qqbar_tensor),
    }


# Compute the exact helicity-vector norm squared
def helicity_norm2(amps):
    return sp.simplify(sum(value * sp.conjugate(value) for value in amps.values()))


# Compute the exact helicity-vector inner product
def helicity_inner(left, right):
    keys = set(left) | set(right)
    return sp.simplify(sum(sp.conjugate(left.get(key, 0)) * right.get(key, 0) for key in keys))


# Compute a coherent sum of scaled helicity dictionaries
def merge_amplitudes(*terms):
    out = {}
    for scale, amps in terms:
        for key, value in amps.items():
            out[key] = sp.simplify(out.get(key, 0) + scale * value)
    return {key: sp.simplify(value) for key, value in out.items() if sp.simplify(value) != 0}


# Compute a unit-normalized helicity dictionary
def normalize_amplitudes(amps):
    norm = sp.sqrt(helicity_norm2(amps))
    if norm == 0:
        raise AssertionError("normalize_amplitudes: zero norm")
    return {key: sp.simplify(value / norm) for key, value in amps.items()}


# Check the normalized Durham and absolute analytic Regge templates
def prove_scalar_templates():
    transverse = scalar_durham_transverse_full()
    soft = scalar_regge_ls00_full()

    require_unit_norm("0++ Durham transverse norm", transverse)
    require_equal("0++ soft Regge LS(0,0) absolute norm", helicity_norm2(soft), 5)
    require_equal("0++ Regge m=0 weight", soft[(0, 0)], -sp.sqrt(sp.Rational(5, 3)))
    require_equal("0++ Regge m=1 weight", soft[(1, 1)], sp.sqrt(sp.Rational(5, 3)))
    require_equal("0++ Regge m=2 pole zero", soft.get((2, 2), 0), 0)
    require_equal("0++ transverse-soft overlap", helicity_inner(transverse, soft), sp.sqrt(sp.Rational(10, 3)))
    require_equal("scalar pole adjacent-m phase", soft[(1, 1)] / soft[(0, 0)], -1)
    for (m, _), value in scalar_regge_ls00_full(alpha=2).items():
        require_equal("spin-two scalar pole contraction", value, sp.Integer(-1) ** m)

    q, phi = sp.symbols("q phi", positive=True, real=True)
    hard = {"++": sp.sqrt(sp.Rational(1, 2)), "--": sp.sqrt(sp.Rational(1, 2))}
    projected = durham_projected_amplitude(hard, q, phi)
    require_equal("0++ Durham hard amplitude", projected, sp.sqrt(2) * durham_projectors(q, phi)["0+"])

    return {"transverse": transverse, "soft_ls00": soft}


# Prove the tensor Durham helicity-zero and helicity-two templates
def prove_tensor_templates():
    helicity2 = tensor_durham_helicity2_full()
    helicity0 = tensor_durham_helicity0_full()

    require_unit_norm("2++ Durham helicity-two norm", helicity2)
    require_unit_norm("2++ Durham helicity-zero norm", helicity0)
    require_equal("2++ H2/H0 orthogonality", helicity_inner(helicity2, helicity0), 0)

    q, phi = sp.symbols("q phi", positive=True, real=True)
    hard_h2 = {"-+": sp.sqrt(sp.Rational(1, 2)), "+-": sp.sqrt(sp.Rational(1, 2))}
    projected_h2 = durham_projected_amplitude(hard_h2, q, phi)
    require_equal(
        "2++ Durham helicity-two hard amplitude",
        projected_h2,
        (durham_projectors(q, phi)["+2+"] + durham_projectors(q, phi)["-2+"]) / sp.sqrt(2),
    )

    hard_h0 = {"++": sp.sqrt(sp.Rational(1, 2)), "--": sp.sqrt(sp.Rational(1, 2))}
    projected_h0 = durham_projected_amplitude(hard_h0, q, phi)
    require_equal("2++ Durham helicity-zero hard amplitude", projected_h0, sp.sqrt(2) * durham_projectors(q, phi)["0+"])

    return {"helicity2": helicity2, "helicity0": helicity0}


# Check complex parton helicities and their fixed-angle projections
def prove_partonic_pair_matrix_elements():
    # A J projection requires angular integration with Wigner D and an orbital wavefunction
    # The forward poles also require physical cuts before angular integration
    coefficients = partonic_pair_kernel_coefficients()
    theta = sp.symbols("theta_pair", positive=True, real=True)
    cgg, kq = sp.symbols("C_gg K_q", nonzero=True)

    # Compare the half-angle expressions with the massless limit of Eq. (37)
    # [REFERENCE: Harland-Lang, Khoze and Ryskin, arXiv:1508.02718, Eq. (37)]
    for final, amplitudes in qqbar_pair_spinor_helicity_matrix().items():
        h = 1 if final[0] == "+" else -1
        for initial, value in amplitudes.items():
            lam = 1 if initial[0] == "+" else -1
            published = -lam * h * kq * (1 - lam * h * sp.cos(theta)) / sp.sin(theta)
            require_equal_trig(f"massless qqbar Eq. (37) {initial}->{final}", value, published)

    # Complex final-state coefficients enter as a bra, also check its overall phase
    phi, rotation = sp.symbols("phi rotation", real=True)
    matrix = gg_pair_spinor_helicity_matrix()
    wavefunction = {"-+": 1 / sp.sqrt(2), "+-": sp.I / sp.sqrt(2)}
    projected = project_initial_hard_matrix(matrix, wavefunction)
    phased = project_initial_hard_matrix(
        matrix, {key: sp.exp(sp.I * rotation) * value for key, value in wavefunction.items()}
    )
    for label, value in projected.items():
        require_equal("final-state bra phase", phased[label], sp.exp(-sp.I * rotation) * value)
        jz = sum((1 if h == "+" else -1) * sign for h, sign in zip(label, (1, -1), strict=True))
        require_equal(
            "hard gg azimuth covariance", value.subs(phi, phi + rotation), sp.exp(sp.I * jz * rotation) * value
        )

    a, b = sp.symbols("a b")
    qqbar_projected = project_initial_hard_matrix(qqbar_pair_spinor_helicity_matrix(), {"-+": a, "+-": b})
    require_equal("massless qqbar pair has no Jz=0 -- initial", qqbar_projected.get("--", 0), 0)
    require_equal("massless qqbar pair has no Jz=0 ++ initial", qqbar_projected.get("++", 0), 0)
    require_nonzero("massless qqbar pair has Jz=2 support", qqbar_projected["-+"])
    require_equal_trig(
        "gg same-helicity projection", coefficients["gg_same"], 2 * sp.sqrt(2) * cgg / sp.sin(theta) ** 2
    )
    require_equal_trig(
        "gg opposite-helicity projection",
        coefficients["gg_opposite"],
        sp.sqrt(2) * cgg * (1 + sp.cos(theta) ** 2) / sp.sin(theta) ** 2,
    )
    require_equal_trig(
        "qqbar opposite-helicity projection",
        coefficients["qqbar_opposite"],
        sp.sqrt(2) * kq * sp.cos(theta) / sp.sin(theta),
    )
    require_equal("symmetric qqbar projection at theta=pi/2", coefficients["qqbar_opposite"].subs(theta, sp.pi / 2), 0)
    require_nonzero(
        "gg and qqbar opposite-helicity kernels differ",
        coefficients["gg_opposite"].subs(cgg, 1) - coefficients["qqbar_opposite"].subs(kq, 1),
    )

    return coefficients


# Prove the Durham meson-pair kernels select transverse helicity sectors
def prove_meson_pair_selection_rules():
    q = sp.symbols("q", positive=True)
    phi, phi_pair = sp.symbols("phi phi_pair", real=True)
    a_sfo, a_sfs_pp, a_sfs_pm = sp.symbols("A_sfo A_sfs_pp A_sfs_pm")

    scalar_flavor_octet = {"-+": a_sfo * sp.exp(-2 * sp.I * phi_pair), "+-": a_sfo * sp.exp(2 * sp.I * phi_pair)}
    scalar_flavor_singlet = {
        "++": a_sfs_pp,
        "--": a_sfs_pp,
        "-+": a_sfs_pm * sp.exp(-2 * sp.I * phi_pair),
        "+-": a_sfs_pm * sp.exp(2 * sp.I * phi_pair),
    }

    sfo_projected = durham_projected_amplitude(scalar_flavor_octet, q, phi)
    sfs_projected = durham_projected_amplitude(scalar_flavor_singlet, q, phi)
    # The hard meson azimuth is fixed during the screening loop
    # [REFERENCE: Harland-Lang, Khoze and Ryskin, arXiv:1508.02718, Eqs. (5)-(6)]
    require_zero("SFO forward loop suppression", sp.integrate(sp.expand_trig(sfo_projected), (phi, 0, 2 * sp.pi)))
    require_equal(
        "SFS forward loop selects Jz=0",
        sp.integrate(sp.expand_trig(sfs_projected), (phi, 0, 2 * sp.pi)),
        2 * sp.pi * q**2 * a_sfs_pp,
    )
    require_nonzero("SFO meson pair has helicity-two term", sp.diff(sfo_projected, a_sfo))
    require_equal(
        "SFS meson pair has scalar same-helicity coefficient",
        sp.diff(sfs_projected, a_sfs_pp),
        2 * durham_projectors(q, phi)["0+"],
    )
    require_nonzero("SFS meson pair has helicity-two coefficient", sp.diff(sfs_projected, a_sfs_pm))

    return {"SFO": sfo_projected, "SFS": sfs_projected}


# Compute one analysis-normalized mixture, not a GP card normalization
def scalar_mixed_template(theta, delta):
    transverse = scalar_durham_transverse_full()
    soft = scalar_regge_ls00_full()
    mixed = merge_amplitudes((sp.cos(theta), transverse), (sp.exp(sp.I * delta) * sp.sin(theta), soft))
    return normalize_amplitudes(mixed)


# Compute 2++ mixing of helicity-two and helicity-zero Durham templates
def tensor_mixed_template(theta, delta):
    helicity2 = tensor_durham_helicity2_full()
    helicity0 = tensor_durham_helicity0_full()
    return merge_amplitudes((sp.cos(theta), helicity2), (sp.exp(sp.I * delta) * sp.sin(theta), helicity0))


# Prove the symbolic 0++ and 2++ mixing templates are normalized
def prove_mixing_templates():
    theta, delta = sp.symbols("theta delta", real=True)
    scalar_raw = merge_amplitudes(
        (sp.cos(theta), scalar_durham_transverse_full()),
        (sp.exp(sp.I * delta) * sp.sin(theta), scalar_regge_ls00_full()),
    )
    tensor_mixed = tensor_mixed_template(theta, delta)

    require_equal(
        "0++ raw mixed template norm",
        helicity_norm2(scalar_raw),
        sp.cos(theta) ** 2 + 5 * sp.sin(theta) ** 2 + sp.sqrt(sp.Rational(10, 3)) * sp.sin(2 * theta) * sp.cos(delta),
    )
    require_unit_norm("2++ mixed template norm", tensor_mixed)
    require_equal("0++ pure Durham limit", scalar_mixed_template(0, delta)[(-1, -1)], sp.sqrt(sp.Rational(1, 2)))
    require_equal("2++ pure helicity-two limit", tensor_mixed_template(0, delta)[(-1, 1)], sp.sqrt(sp.Rational(1, 2)))
    require_equal(
        "2++ pure helicity-zero limit", tensor_mixed_template(sp.pi / 2, 0)[(-1, -1)], sp.sqrt(sp.Rational(1, 2))
    )

    return {"scalar_raw": scalar_raw, "tensor": tensor_mixed}


# Run all symbolic glueball and meson helicity checks
def main():
    # These selected helicity wavefunctions do not determine a bound-state total J
    # Templates are hypotheses, not fitted couplings or derived steering predictions
    checks = {
        "regge_durham_bridge": prove_regge_residue_matches_durham_projectors(),
        "scalar_templates": prove_scalar_templates(),
        "tensor_templates": prove_tensor_templates(),
        "partonic_pair_matrix_elements": prove_partonic_pair_matrix_elements(),
        "meson_pair_selection": prove_meson_pair_selection_rules(),
        "mixing_templates": prove_mixing_templates(),
    }
    for name, value in checks.items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
