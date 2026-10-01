# Symbolic normalization checks for MGamma amplitudes and EPA factors
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import sympy as sp


# Require an exact symbolic expression to vanish
def require_zero(name, expr):
    value = sp.simplify(expr)
    if value != 0:
        raise AssertionError(f"{name}: expected 0, got {value}")


# Require two exact symbolic expressions to be equal
def require_equal(name, expr, expected):
    require_zero(name, sp.simplify(expr - expected))


# Derive the CP-even scalar helicity content from the gauge-invariant photon tensor
def prove_scalar_photon_tensor():
    # MGamma::yyMonopolium and yyHiggs use this parity-even on-shell tensor
    mass, coupling = sp.symbols("M g", positive=True)
    metric = sp.diag(1, -1, -1, -1)
    k1 = sp.Matrix([mass / 2, 0, 0, mass / 2])
    k2 = sp.Matrix([mass / 2, 0, 0, -mass / 2])
    vertex = coupling * ((k1.T * metric * k2)[0] * metric - (metric * k2) * (metric * k1).T)
    for value in list(k1.T * vertex) + list(vertex * k2):
        require_zero("scalar photon Ward identity", value)
    helicities = {}
    for h1 in (-1, 1):
        for h2 in (-1, 1):
            eps1 = sp.Matrix([0, -h1, -sp.I, 0]) / sp.sqrt(2)
            eps2 = sp.Matrix([0, h2, -sp.I, 0]) / sp.sqrt(2)
            amplitude = (eps1.T * vertex * eps2)[0]
            require_equal("scalar photon helicity", amplitude, coupling * mass**2 / 2 if h1 == h2 else 0)
            helicities[h1, h2] = amplitude
    return helicities


# Prove the scalar gamma-gamma Breit-Wigner coefficient from the decay width
def prove_scalar_resonance_normalization():
    shat, mass, gamma_tot, gamma_yy = sp.symbols("shat M Gamma Gamma_yy", positive=True)
    # [REFERENCE: PDG 2025, Kinematics review, two-body decay phase space]
    # Width matching assumes a constant pole coefficient and a fixed total width
    production_coupling2 = sp.symbols("g_yy^2", positive=True)

    denom = (shat - mass**2) ** 2 + (mass * gamma_tot) ** 2

    two_body_lips = 1 / (8 * sp.pi)
    identical_photons = sp.Rational(1, 2)
    same_helicity_states = 2
    width = sp.simplify(
        (1 / (2 * mass)) * identical_photons * two_body_lips * same_helicity_states * production_coupling2
    )
    solved_coupling2 = sp.solve(sp.Eq(width, gamma_yy), production_coupling2)[0]

    delta_bw = mass * gamma_tot / (sp.pi * denom)
    one_body_lips_bw = 2 * sp.pi * delta_bw
    same_helicity_amp2 = sp.simplify(solved_coupling2 * one_body_lips_bw)
    # src/Photon/MGamma.cc: ScalarYYResonanceAmplitude and ScalarDecayTensor
    # src/Particle/MResonance.cc: FixedWidthLineShape has no extra factor i
    amplitude = sp.sqrt(32 * sp.pi * mass**2 * gamma_yy * gamma_tot) / (shat - mass**2 + sp.I * mass * gamma_tot)
    require_equal("scalar complex pole residue", amplitude * sp.conjugate(amplitude), same_helicity_amp2)
    require_equal(
        "scalar pole phase", amplitude.subs(shat, mass**2), -sp.I * sp.sqrt(32 * sp.pi * gamma_yy / gamma_tot)
    )
    # Interference with a real continuum must change sign across the pole
    offset = sp.symbols("offset", real=True)
    residue = sp.sqrt(32 * sp.pi * mass**2 * gamma_yy * gamma_tot)
    require_equal(
        "scalar dispersive interference",
        sp.re(amplitude.subs(shat, mass**2 + offset)),
        residue * offset / (offset**2 + mass**2 * gamma_tot**2),
    )
    spin_avg_amp2 = sp.simplify((same_helicity_amp2 + same_helicity_amp2) / 4)
    required_norm = sp.simplify(spin_avg_amp2 / (gamma_yy * gamma_tot / denom))

    avg_xsec = sp.simplify(spin_avg_amp2 / (2 * shat))
    same_helicity_xsec = sp.simplify(same_helicity_amp2 / (2 * shat))
    # A constant pole coefficient fixes the on-shell width, not its off-shell continuation
    # The often quoted constant-numerator cross section agrees only at shat = M^2
    standard_avg_xsec = 8 * sp.pi * gamma_yy * gamma_tot / denom
    standard_same_helicity_xsec = 16 * sp.pi * gamma_yy * gamma_tot / denom
    off_shell_residual = sp.simplify(avg_xsec / standard_avg_xsec)

    require_equal("scalar yy decay coupling", solved_coupling2, 16 * sp.pi * mass * gamma_yy)
    require_equal(
        "scalar yy same-helicity amp2", same_helicity_amp2, 32 * sp.pi * mass**2 * gamma_yy * gamma_tot / denom
    )
    require_equal("scalar yy spin-averaged amp2", spin_avg_amp2, 16 * sp.pi * mass**2 * gamma_yy * gamma_tot / denom)
    require_equal("scalar yy normalization coefficient", required_norm, 16 * sp.pi * mass**2)
    require_equal("scalar yy off-shell residual", off_shell_residual, mass**2 / shat)
    require_equal(
        "scalar yy same-helicity off-shell residual", same_helicity_xsec / standard_same_helicity_xsec, mass**2 / shat
    )
    require_equal(
        "scalar yy on-shell averaged cross section",
        sp.simplify(avg_xsec.subs(shat, mass**2)),
        sp.simplify(standard_avg_xsec.subs(shat, mass**2)),
    )

    return {"required_norm": required_norm, "off_shell_shat_residual": off_shell_residual}


# Prove the qqbar charge-color scaling used after the MG5 unit-charge amplitude
def prove_quark_charge_color_factor():
    # MGamma::yy_ffbar multiplies the unit-charge amplitude by sqrt(3) Q^2
    # Two electromagnetic vertices supply Q^2 and delta_ij fixes the color state
    q = sp.symbols("Q", real=True)
    amplitude = sp.symbols("A")
    color_amplitude = q**2 * amplitude * sp.eye(3)
    explicit_sum = sp.trace(color_amplitude * color_amplitude.H)
    stored = sp.sqrt(3) * q**2 * amplitude
    require_equal("quark color trace", stored * sp.conjugate(stored), explicit_sum)
    require_equal("quark charge reversal", stored.subs(q, -q), stored)
    return {"amplitude_factor": sp.sqrt(3) * q**2, "squared_factor": 3 * q**4}


# Run all exact symbolic checks
def main():
    checks = {
        "scalar_tensor": prove_scalar_photon_tensor(),
        "scalar_resonance": prove_scalar_resonance_normalization(),
        "quark_charge_color": prove_quark_charge_color_factor(),
    }
    for name, value in checks.items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
