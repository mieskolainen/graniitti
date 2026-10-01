# Symbolic multichannel eikonal scattering amplitudes for arbitrary N
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import sympy as sp

# MMatrix::MixingReal, MEikonalMatrix::ProjectAmplitude and GoodWalkerSpace
# T[i,j] below is only the diagonal eigenchannel limit of the full pair matrix


# Construct the two-channel rotation in the C++ row convention
def make_U_for_N2():
    """
    Return a 2x2 orthonormal U matrix for the proton and one excitation

    The coefficients are alpha=cos(theta), beta=sin(theta), so the
    normalization constraint is exact rather than an external assumption
    """

    theta = sp.symbols("theta", real=True)

    U = sp.Matrix([[sp.cos(theta), sp.sin(theta)], [-sp.sin(theta), sp.cos(theta)]])

    return U, (theta,)


# Construct the three-angle CKM/PMNS-type matrix used by MixingReal for N = 3
def make_U_for_N3():
    """
    Construct the three-angle CKM/PMNS-type matrix used by MixingReal for N = 3

    Setting theta_13=theta_23=0 reduces to the conventional two-channel
    theta_12 rotation plus a decoupled third state
    """

    theta_12, theta_13, theta_23 = sp.symbols("theta_12 theta_13 theta_23", real=True)
    c12, s12 = sp.cos(theta_12), sp.sin(theta_12)
    c13, s13 = sp.cos(theta_13), sp.sin(theta_13)
    c23, s23 = sp.cos(theta_23), sp.sin(theta_23)

    # MMatrix::MixingReal multiplies R23 R13 R12 from the left
    r12 = sp.Matrix([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]])
    r13 = sp.Matrix([[c13, 0, s13], [0, 1, 0], [-s13, 0, c13]])
    r23 = sp.Matrix([[1, 0, 0], [0, c23, s23], [0, -s23, c23]])
    U = r23 * r13 * r12

    return U, (theta_12, theta_13, theta_23)


# Project a diagonal pair amplitude for real orthogonal proton mixing
def build_amplitude_matrix(U, T):
    weights = sp.diag(*U.row(0))
    return (U * weights * sp.Matrix(T) * weights * U.T).applyfunc(sp.simplify)


# Project a general pair operator exactly as MEikonalMatrix::ProjectAmplitude
# Rows specify physical-state coefficients in the pair basis
def pair_amplitudes(U, operator):
    rotation = sp.kronecker_product(U, U)
    proton = rotation.row(0).T
    return (rotation.conjugate() * operator * proton).reshape(U.rows, U.rows)


# Compute a squared complex amplitude as a polynomial in amplitudes and conjugates
def norm2(value):
    return sp.expand(value * sp.conjugate(value))


# Sum probabilities for orthogonal physical final states
def channel_weights(amp):
    return {
        "elastic": norm2(amp[0, 0]),
        "sd_up": sum(norm2(amp[m, 0]) for m in range(1, amp.rows)),
        "sd_down": sum(norm2(amp[0, n]) for n in range(1, amp.cols)),
        "dd": sum(norm2(amp[m, n]) for m in range(1, amp.rows) for n in range(1, amp.cols)),
    }


# Require exact scalar or matrix equality with complex amplitudes
def require_equal(name, left, right):
    residual = left - right
    values = list(residual) if isinstance(residual, sp.MatrixBase) else [residual]
    if any(sp.trigsimp(sp.expand(value)) != 0 for value in values):
        raise AssertionError(f"{name}: {residual}")


# Derive elastic, single and double dissociation from Good-Walker completeness
def check_multichannel_derivations():
    # [REFERENCE: Khoze, Martin and Ryskin, arXiv:1306.2149, Appendix A]
    common = sp.symbols("T_common")
    results = {}
    for channels, builder in ((2, make_U_for_N2), (3, make_U_for_N3)):
        U, angles = builder()
        require_equal(f"N={channels} orthogonality", U * U.T, sp.eye(channels))
        expected = sp.zeros(channels)
        expected[0, 0] = common
        require_equal("common eigenamplitude", build_amplitude_matrix(U, sp.ones(channels) * common), expected)

        rotation = U.subs(dict(zip(angles, (sp.pi / 6, sp.pi / 4, sp.pi / 3), strict=False)))
        eigen = sp.Matrix(channels, channels, sp.symbols(f"t0:{channels**2}"))
        exclusive = build_amplitude_matrix(rotation, eigen)
        require_equal("diagonal pair limit", pair_amplitudes(rotation, sp.diag(*eigen)), exclusive)
        weights = [rotation[0, i] ** 2 for i in range(channels)]
        mean = sum(weights[i] * weights[j] * eigen[i, j] for i in range(channels) for j in range(channels))
        mean2 = sum(weights[i] * weights[j] * norm2(eigen[i, j]) for i in range(channels) for j in range(channels))
        up_mean2 = sum(
            weights[i] * norm2(sum(weights[j] * eigen[i, j] for j in range(channels))) for i in range(channels)
        )
        down_mean2 = sum(
            weights[j] * norm2(sum(weights[i] * eigen[i, j] for i in range(channels))) for j in range(channels)
        )
        norms = channel_weights(exclusive)
        require_equal("elastic amplitude", exclusive[0, 0], mean)
        require_equal("upper SD variance", norms["sd_up"], up_mean2 - norm2(mean))
        require_equal("lower SD variance", norms["sd_down"], down_mean2 - norm2(mean))
        require_equal("DD covariance", norms["dd"], mean2 - up_mean2 - down_mean2 + norm2(mean))
        require_equal("inclusive completeness", sum(norms.values()), mean2)
        require_equal("beam exchange", build_amplitude_matrix(rotation, eigen.T), exclusive.T)
        results[f"N={channels}"] = "elastic, SD, DD and completeness"

        # A complex unitary change among excited states cannot change inclusive rates
        if channels == 3:
            phase = sp.symbols("phase", real=True)
            excited = sp.diag(1, 1, sp.exp(sp.I * phase))
            excited[1:3, 1:3] = sp.Matrix([[1, sp.I], [sp.I, 1]]) / sp.sqrt(2) * excited[1:3, 1:3]
            rotated = channel_weights(excited.conjugate() * exclusive * excited.H)
            for name, norm in norms.items():
                require_equal(f"excited basis invariant {name}", rotated[name], norm)

    U3, (t12, t13, t23) = make_U_for_N3()
    U2, (t2,) = make_U_for_N2()
    require_equal("three-to-two channel limit", U3.subs({t13: 0, t23: 0}), sp.diag(U2.subs(t2, t12), 1))
    return results


# Check the full C++ pair projection without assuming a diagonal interaction
def check_pair_operator():
    U, (theta,) = make_U_for_N2()
    U = U.subs(theta, sp.pi / 6)
    operator = sp.Matrix(4, 4, sp.symbols("A0:16"))
    amp = pair_amplitudes(U, operator)
    proton = sp.kronecker_product(U.row(0).T, U.row(0).T)
    require_equal(
        "general pair completeness", sum(channel_weights(amp).values()), (proton.H * operator.H * operator * proton)[0]
    )
    # A common diagonal amplitude is elastic, a common full matrix is not
    a = sp.symbols("a")
    require_equal("identity interaction", pair_amplitudes(U, a * sp.eye(4)), sp.diag(a, 0))
    return "general complex pair operator"


# Run physical channel checks without printing coherent sums of distinct states
def main():
    print(check_multichannel_derivations())
    print(check_pair_operator())


if __name__ == "__main__":
    main()
