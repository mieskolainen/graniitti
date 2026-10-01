# Generate and verify photon-Pomeron and Pomeron-Pomeron coupling structures
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

"""Generate complete photon-Pomeron and Pomeron-Pomeron meson pole vertex bases

Requires SymPy and the Python standard library, with no external input files

Conventions and scope
    g = diag(1,-1,-1,-1), eps_0123 = +1
    eps(v1,v2,v3,v4) is the determinant of their contravariant components
    The meson has p^2 = M^2 > 0 and a symmetric, p-transverse, traceless spin-J
    polarization tensor R. Each Pomeron P_munu is symmetric and Lorentz-traceless,
    with no momentum-transversality constraint unless explicitly selected
    m0 > 0 is a fixed reference mass. Couplings absorb changes of its value
    Each basis element has an independent scalar coefficient function of the
    central invariants. Completeness is at generic kinematics after projecting
    the meson onto its physical spin, modulo traces and scalar invariant factors
    P_X = +/-1 is physical parity. Natural parity means P_X = (-1)^J
    Charge conjugation requires C_X = -1 for gamma-P and C_X = +1 for P-P

Photon-Pomeron basis
    p = q+k, where q is the photon momentum, k the Pomeron momentum and e the
    photon polarization. Complex auxiliary vectors a,b encode the meson and
    Pomeron tensors through degrees J,2, with a.p = a^2 = b^2 = 0
    Differentiate a basis polynomial J times in a, twice in b and once in e,
    divide by J! 2!, then apply the meson and Pomeron tensor projectors
    H = k.q, E = H*e - (k.e)*q, so E is Ward invariant and k.E = 0
    H > 0 for spacelike incoming virtualities and timelike p
    z = a.q/m0, u = b.p/m0, v = b.q/m0, w = a.b
    x = E.a, y = E.b, ell = E.q/m0
    A = eps(E,p,q,a)/m0^2, B = eps(E,p,q,b)/m0^2, D = eps(E,q,a,b)/m0
    M(r,d) = {w^n z^(r-n) u^i v^(d-n-i), 0 <= n <= min(r,d), 0 <= i <= d-n}
    M(r,d) is empty for a negative index. The two parity bases are
        natural: x*M(J-1,2), y*M(J,1), ell*M(J,2)
        opposite: A*M(J-1,2), B*M(J,1), D*M(J-1,1)
    These span the Ward kernel after four-dimensional Schouten reduction
    Real photons obey q^2 = e.q = 0, ell = 0 and
        (p.q)*D = m0^2*(v*A-z*B)
    Thus --real retains only the x,y or A,B families. Additional Pomeron
    transversality gives b.k = 0, or u = v, implemented by --transverse
    All T_i have mass dimension two. No electromagnetic prefactor is included

Pomeron-Pomeron basis
    p = k1+k2, n = p/M, h_munu = g_munu-n_mu*n_nu
    Decompose each exchange into rotation spins 0,1,2 in the meson rest frame
        sigma = n^mu*n^nu*P_munu
        v_mu = h_mu^alpha*P_alphabeta*n^beta
        t_munu = (h_mu^alpha*h_nu^beta-h_munu*h^alphabeta/3)*P_alphabeta
        P_munu = t_munu+n_mu*v_nu+n_nu*v_mu+sigma*(n_mu*n_nu-h_munu/3)
    In any orthonormal spatial triad e_r, use v_r = e_r^mu*v_mu and
    t_rs = e_r^mu*e_s^nu*t_munu. Spherical components use Condon-Shortley phases
        U_0.v = v_z, U_(+1).v = -(v_x+i*v_y)/sqrt(2)
        U_(-1).v = (v_x-i*v_y)/sqrt(2)
        Pi_0_0 = sigma_i, Pi_1_m = sum_r U_mr*v_ir
        Pi_2_m = sum_(m1,m2,r,s) CG(1,m1,1,m2|2,m)*U_m1r*U_m2s*t_irs
    Pi means P1 or P2 in the printed symbols. Rb_J_M is the complex conjugate
    meson polarization in the same orthonormal spherical STF convention
    Q^mu = h^mu_nu*(k1-k2)^nu/(2*m0), with Cartesian Q_r = -e_r.Q
    Q_L_m = |Q|^L*sqrt(4*pi/(2*L+1))*Y_Lm(Qhat), Q_0_0 = 1
    These regular solid harmonics are expanded in Qx,Qy,Qz by --cartesian
    Ordered scalars are the sum over all magnetic projections of
        CG(s1,m1,s2,m2|S,mS)*CG(L,mL,S,mS|J,M)
        *P1_s1_m1*P2_s2_m2*Q_L_mL*Rb_J_M
    Allowed indices satisfy s1,s2 in {0,1,2}, |s1-s2| <= S <= s1+s2,
    |J-S| <= L <= J+S and (-1)^(s1+s2+L) = P_X
    The scalar, vector and tensor exchange components have parities +,-,+
    Bose exchange maps V_s1s2LS to eta*V_s2s1LS, eta = (-1)^(s1+s2-S+L)
    Printed V(s1,s2,L,S,xi) has exchange eigenvalue xi. For s1 < s2 it is
        (V_s1s2LS + xi*eta*V_s2s1LS)/sqrt(2)
    For s1 = s2 the ordered scalar is used directly and xi = eta
    Its coefficient obeys g_xi(t1,t2) = xi*g_xi(t2,t1), with ti = ki^2
    Analytic exchange-odd coefficients contain (t1-t2)/m0^2 times a symmetric
    function. --exchange even selects symmetric tensors. --spin2 retains only
    s1 = s2 = 2 in this rest-frame decomposition, not the constraint ki.Pi = 0

Checks and output
    --check verifies Ward and Schouten identities, exact photon ranks, tensor
    degrees, LS completeness, parity, Bose exchange and rotational covariance
    Synthetic rational momenta in these checks are algebraic test points
    --spin selects one spin, --max-spin all spins from zero to that value
    --format latex changes expression syntax in the labelled output listing
    --output writes that listing to the specified path
"""

import argparse
import sys
from itertools import product
from pathlib import Path
from random import Random

import sympy as sp
from sympy.physics.wigner import clebsch_gordan as cg


# Require an exact scalar or matrix identity
def require_zero(name, expr):
    values = list(expr) if isinstance(expr, sp.MatrixBase) else [expr]
    if any(sp.expand(value) != 0 for value in values):
        raise AssertionError(name)


# Enumerate scalar monomials with fixed meson and exchange degrees
def monomials(j, d):
    z, u, v, w = sp.symbols("z u v w")
    return [w**n * z ** (j - n) * u**i * v ** (d - n - i) for n in range(min(j, d) + 1) for i in range(d - n + 1)]


# Construct a Ward-invariant basis for physical parity +1 or -1
def photon_basis(j, parity, real=False, transverse=False):
    if j < 0 or parity not in (-1, 1):
        raise ValueError("Require integer spin >= 0 and parity +1 or -1")
    x, y, ell, a, b, d, u, v = sp.symbols("x y ell A B D u v")
    natural = parity == (-1) ** j
    families = [(x if natural else a, j - 1, 2), (y if natural else b, j, 1)]
    if not real:
        families.append((ell, j, 2) if natural else (d, j - 1, 1))
    terms = [prefactor * term for prefactor, degree, rank in families for term in monomials(degree, rank)]
    return list(dict.fromkeys(term.subs(u, v) if transverse else term for term in terms))


# Enumerate the exchange eigenbasis as (s1, s2, L, S, exchange sign)
def pp_indices(j, parity, exchange="all", spin2=False):
    if j < 0 or parity not in (-1, 1) or exchange not in ("all", "even", "odd"):
        raise ValueError("Invalid spin, parity or exchange symmetry")
    for a in [2] if spin2 else range(3):
        for b in range(a, 3):
            for s in range(abs(a - b), a + b + 1):
                for ell in range(abs(j - s), j + s + 1):
                    if (-1) ** (a + b + ell) != parity:
                        continue
                    signs = [(-1) ** (a + b - s + ell)] if a == b else [1, -1]
                    for sign in signs:
                        if exchange == "all" or sign == (1 if exchange == "even" else -1):
                            yield a, b, ell, s, sign


# Expand regular solid harmonics Q_lm in real Cartesian momentum components
def solid(ell, m):
    x, y, z = sp.symbols("Qx Qy Qz", real=True)
    if m < 0:
        return (-1) ** (-m) * sp.conjugate(solid(ell, -m))
    t = sp.Symbol("t")
    poly = sp.Poly(sp.diff(sp.legendre(ell, t), t, m), t)
    radial = sum(c * z ** n[0] * (x * x + y * y + z * z) ** ((ell - m - n[0]) // 2) for n, c in poly.terms())
    return sp.expand((-1) ** m * sp.sqrt(sp.factorial(ell - m) / sp.factorial(ell + m)) * (x + sp.I * y) ** m * radial)


# Construct one ordered LS scalar using Condon-Shortley coefficients
def pp_ordered(j, a, b, ell, s, cartesian=False):
    harmonics = {m: solid(ell, m) if cartesian else sp.Symbol(f"Q_{ell}_{m}") for m in range(-ell, ell + 1)}
    terms = []
    for m1, m2 in product(range(-a, a + 1), range(-b, b + 1)):
        ms = m1 + m2
        c1 = cg(a, b, s, m1, m2, ms)
        if c1 == 0:
            continue
        for ml in range(max(-ell, -j - ms), min(ell, j - ms) + 1):
            c2 = cg(ell, s, j, ml, ms, ml + ms)
            if c2 != 0:
                terms.append(
                    c1
                    * c2
                    * sp.Symbol(f"P1_{a}_{m1}")
                    * sp.Symbol(f"P2_{b}_{m2}")
                    * harmonics[ml]
                    * sp.Symbol(f"Rb_{j}_{ml + ms}")
                )
    return sp.Add(*terms)


# Combine the two leg orderings into an exchange eigenstate
def pp_vertex(j, index, cartesian=False):
    a, b, ell, s, sign = index
    direct = pp_ordered(j, a, b, ell, s, cartesian)
    if a == b:
        return direct
    reverse = pp_ordered(j, b, a, ell, s, cartesian)
    return sp.expand((direct + sign * (-1) ** (a + b - s + ell) * reverse) / sp.sqrt(2))


# Independently count photon amplitudes from magnetic projections and parity
def photon_count(j, parity, real, transverse):
    spins = [2] if transverse else range(3)
    photon = [-1, 1] if real else [-1, 0, 1]
    total = sum(abs(m + lam) <= j for s in spins for m in range(-s, s + 1) for lam in photon)
    fixed = 0 if real else len(spins)
    return (total + parity * (-1) ** j * fixed) // 2


# Prove the photon Ward map and the four-dimensional five-vector identity
def check_photon():
    metric = sp.diag(1, -1, -1, -1)
    e, p, q, a, b = [sp.Matrix(sp.symbols(f"{name}0:4")) for name in ("e", "p", "q", "a", "b")]
    k = p - q
    h = (k.T * metric * q)[0]
    electric = h * e - (k.T * metric * e)[0] * q
    require_zero("photon Ward identity", electric.xreplace(dict(zip(e, q, strict=True))))
    require_zero("E.k = 0", k.T * metric * electric)
    vectors = [e, p, q, a, b]
    schouten = sp.zeros(4, 1)
    for i, vector in enumerate(vectors):
        schouten += (-1) ** i * vector * sp.Matrix.hstack(*[v for n, v in enumerate(vectors) if n != i]).det()
    require_zero("five-vector Schouten identity", schouten)


# Evaluate the photon invariants at exact synthetic momenta with m0 = 1
def photon_sample(rng, real, transverse):
    t, r, v = [rng.randint(-11, 11) for _ in range(3)]
    metric = sp.diag(1, -1, -1, -1)
    p, q = sp.Matrix([8 if real else 5, 0, 0, 0]), sp.Matrix([5 if real else 2, 0, 0, 5])
    a = sp.Matrix([0, 1 - t * t, sp.I * (1 + t * t), 2 * t])
    b = sp.Matrix(
        [-5 * (1 + r * r), 4 * (1 - r * r), 8 * r, 3 * (1 + r * r)]
        if transverse
        else [1 + r * r + v * v, 2 * r, 2 * v, 1 - r * r - v * v]
    )
    e = sp.Matrix([rng.randint(-3, 3) for _ in range(4)])
    if real:
        e[0] = e[3] = 0
    k = p - q
    electric = (k.T * metric * q)[0] * e - (k.T * metric * e)[0] * q
    values = [
        (a.T * metric * q)[0],
        (b.T * metric * p)[0],
        (b.T * metric * q)[0],
        (a.T * metric * b)[0],
        (electric.T * metric * a)[0],
        (electric.T * metric * b)[0],
        (electric.T * metric * q)[0],
        sp.Matrix.hstack(electric, p, q, a).det(),
        sp.Matrix.hstack(electric, p, q, b).det(),
        sp.Matrix.hstack(electric, q, a, b).det(),
    ]
    require_zero("auxiliary traces", sp.Matrix([(a.T * metric * a)[0], (a.T * metric * p)[0], (b.T * metric * b)[0]]))
    if transverse:
        require_zero("exchange transversality", b.T * metric * k)
    return dict(zip(sp.symbols("z u v w x y ell A B D"), values, strict=True))


# Prove photon-basis independence by an exact nonzero minor at physical virtualities
def check_photon_rank(j, real, transverse):
    bases = [photon_basis(j, parity, real, transverse) for parity in (1, -1)]
    rng = Random(7919)
    samples = [photon_sample(rng, real, transverse) for _ in range(2 * max(map(len, bases)) + 5)]
    for basis in bases:
        matrix = sp.Matrix([[expr.xreplace(sample) for expr in basis] for sample in samples])
        if matrix.applyfunc(sp.expand).to_DM().convert_to(sp.QQ_I).rank() != len(basis):
            raise AssertionError(f"Photon rank J={j}, real={real}, transverse={transverse}")


# Check the LS-to-magnetic-projection map with exact rational square roots
def check_ls(j, a, b):
    indices = [(ell, s) for s in range(abs(a - b), a + b + 1) for ell in range(abs(j - s), j + s + 1)]
    projections = [(m1, m2) for m1 in range(-a, a + 1) for m2 in range(-b, b + 1) if abs(m1 + m2) <= j]
    matrix = sp.Matrix(
        [
            [cg(a, b, s, m1, m2, m1 + m2) * cg(ell, s, j, 0, m1 + m2, m1 + m2) for ell, s in indices]
            for m1, m2 in projections
        ]
    )
    norm = sp.diag(*[sp.Rational(2 * j + 1, 2 * ell + 1) for ell, _ in indices])
    require_zero(f"LS orthogonality J={j}, s1={a}, s2={b}", matrix.T * matrix - norm)
    if matrix.rows != matrix.cols:
        raise AssertionError("Incomplete LS basis")


# Construct parity, exchange and rotation transformations of spherical components
def transformations(expr):
    parity, exchange, raising = {}, {}, {}
    for symbol in expr.free_symbols:
        name, spin, m = str(symbol).split("_")
        spin, m = int(spin), int(m)
        parity[symbol] = (-1) ** spin * symbol if name != "Rb" else symbol
        if name in ("P1", "P2"):
            exchange[symbol] = sp.Symbol(f"{'P2' if name == 'P1' else 'P1'}_{spin}_{m}")
        elif name == "Q":
            exchange[symbol] = (-1) ** spin * symbol
        if name == "Rb":
            raising[symbol] = -sp.sqrt((spin - m) * (spin + m + 1)) * sp.Symbol(f"Rb_{spin}_{m + 1}")
        else:
            raising[symbol] = sp.sqrt((spin + m) * (spin - m + 1)) * sp.Symbol(f"{name}_{spin}_{m - 1}")
    return parity, exchange, raising


# Verify the Cartesian solid-harmonic phases and normalization under rotations
def check_harmonics(max_rank):
    x, y, z = sp.symbols("Qx Qy Qz", real=True)
    for ell in range(max_rank + 1):
        harmonics = {m: solid(ell, m) for m in range(-ell, ell + 1)}
        for m, expr in harmonics.items():
            require_zero("solid harmonic Laplacian", sum(sp.diff(expr, q, 2) for q in (x, y, z)))
            require_zero("solid harmonic normalization", expr.subs({x: 0, y: 0, z: 1}) - int(m == 0))
            require_zero("solid harmonic Lz", -sp.I * (x * sp.diff(expr, y) - y * sp.diff(expr, x)) - m * expr)
            raised = z * (sp.diff(expr, x) + sp.I * sp.diff(expr, y)) - (x + sp.I * y) * sp.diff(expr, z)
            require_zero(
                "solid harmonic raising", raised - sp.sqrt((ell - m) * (ell + m + 1)) * harmonics.get(m + 1, 0)
            )
    print(f"Exact Cartesian harmonic identities passed through L={max_rank}", flush=True)


# Check the actual coupling polynomials at amplitude level for all requested spins
def verify(max_spin):
    check_photon()
    print("Exact photon Ward and Schouten identities passed", flush=True)
    z, u, v, w, x, y, ell, aa, bb, dd, scale_a, scale_b = sp.symbols("z u v w x y ell A B D sa sb")
    degrees = {
        z: scale_a * z,
        u: scale_b * u,
        v: scale_b * v,
        w: scale_a * scale_b * w,
        x: scale_a * x,
        y: scale_b * y,
        aa: scale_a * aa,
        bb: scale_b * bb,
        dd: scale_a * scale_b * dd,
    }
    for j in range(max_spin + 1):
        for real, transverse in product((False, True), repeat=2):
            check_photon_rank(j, real, transverse)
        for parity, real, transverse in product((1, -1), (False, True), (False, True)):
            basis = photon_basis(j, parity, real, transverse)
            if len(basis) != photon_count(j, parity, real, transverse):
                raise AssertionError("Photon basis dimension")
            for expr in basis:
                require_zero("photon tensor degrees", expr.xreplace(degrees) - scale_a**j * scale_b**2 * expr)
        for a, b in product(range(3), repeat=2):
            check_ls(j, a, b)
        for parity in (1, -1):
            indices = list(pp_indices(j, parity))
            total = sum(
                abs(m + n) <= j for a in range(3) for b in range(3) for m in range(-a, a + 1) for n in range(-b, b + 1)
            )
            if len(indices) != (total + parity * (-1) ** j * 9) // 2:
                raise AssertionError("Pomeron-Pomeron basis dimension")
            for index in indices:
                expr = pp_vertex(j, index)
                inversion, exchange, raising = transformations(expr)
                require_zero("PP parity", expr.xreplace(inversion) - parity * expr)
                require_zero("PP Bose exchange", expr.xreplace(exchange) - index[-1] * expr)
                require_zero(
                    "PP rotational covariance", sum(sp.diff(expr, var) * value for var, value in raising.items())
                )
        print(
            f"J={j}: exact photon ranks, tensor degrees, dimensions, LS orthogonality, parity, exchange and rotations passed",
            flush=True,
        )

    check_harmonics(max_spin + 4)


# Print bases with explicit definitions of all symbolic inputs
def write_bases(args, stream):
    for line in __doc__.splitlines():
        print(f"# {line}".rstrip(), file=stream)
    spins = [args.spin] if args.spin is not None else range(args.max_spin + 1)
    parities = (1, -1) if args.parity == "both" else (1 if args.parity == "+" else -1,)
    printer = sp.latex if args.format == "latex" else sp.sstr
    for j, parity in product(spins, parities):
        if args.channel in ("gamma", "both"):
            basis = photon_basis(j, parity, args.real, args.transverse)
            print(
                f"\n# gamma-P: J={j}, P={parity:+d}, C=-1, real={args.real}, k-transverse={args.transverse}, N={len(basis)}",
                file=stream,
            )
            if not args.counts:
                for n, expr in enumerate(basis):
                    print(f"T{n} = {printer(expr)}", file=stream)
        if args.channel in ("pp", "both"):
            indices = list(pp_indices(j, parity, args.exchange, args.spin2))
            print(
                f"\n# P-P: J={j}, P={parity:+d}, C=+1, exchange={args.exchange}, spin2={args.spin2}, N={len(indices)}",
                file=stream,
            )
            if not args.counts:
                for index in indices:
                    print(f"V{index} = {printer(pp_vertex(j, index, args.cartesian))}", file=stream)


# Parse the requested spins, physical parities and optional physical restrictions
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    spin = parser.add_mutually_exclusive_group()
    spin.add_argument("--spin", type=int, help="one nonnegative integer spin")
    spin.add_argument("--max-spin", type=int, default=4, help="print all spins from zero to this value (default 4)")
    parser.add_argument("--parity", choices=("+", "-", "both"), default="both", help="physical meson parity")
    parser.add_argument("--channel", choices=("gamma", "pp", "both"), default="both")
    parser.add_argument("--real", action="store_true", help="real photon, q^2=0 and e.q=0")
    parser.add_argument("--transverse", action="store_true", help="add k.P=0 in the photon-Pomeron channel")
    parser.add_argument("--spin2", action="store_true", help="retain only s1=s2=2 in the PP rest-frame decomposition")
    parser.add_argument("--exchange", choices=("all", "even", "odd"), default="all")
    parser.add_argument("--cartesian", action="store_true", help="expand PP solid harmonics in Qx,Qy,Qz")
    parser.add_argument(
        "--format", choices=("text", "latex"), default="text", help="syntax of expressions in the labelled listing"
    )
    parser.add_argument("--counts", action="store_true", help="print dimensions only")
    parser.add_argument("--output", type=Path, help="write structures to this file instead of stdout")
    parser.add_argument(
        "--check", action="store_true", help="run exact identities through the requested spin, then exit"
    )
    args = parser.parse_args()
    last = args.spin if args.spin is not None else args.max_spin
    if last < 0:
        parser.error("Spin must be nonnegative")
    if args.check:
        verify(last)
    elif args.output:
        with args.output.open("w") as stream:
            write_bases(args, stream)
    else:
        write_bases(args, sys.stdout)


if __name__ == "__main__":
    main()
