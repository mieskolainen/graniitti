# Tests for tuning parameter topology utilities
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.tune.parameters import tools


# Check LS reflection acts on m while helicity parity also reverses both helicities
def test_gp_coupling_completion_weights():
    labels = ["2,0,0", "2,0,-1", "2,0,2"]
    assert tools.continuum_completion_weights("CON_GP|990:[211,211]/opposite:g_ls", labels) == [1.0, 2.0, 2.0]
    labels = ["0,0,0", "-1/2,-1/2,0", "0,0,-1"]
    assert tools.continuum_completion_weights("CON_GP|990:[2212,2212]/opposite:helicity", labels) == [1.0, 2.0, 2.0]


# Check compact helicity orbits use the same parity and leg exchange rules for all models
@pytest.mark.parametrize("model", ["MP", "XP", "GP"])
def test_continuum_completion_weights(model):
    base = f"CON_{model}|22:"
    assert tools.continuum_completion_weights(base + "[2212,2212]/opposite:helicity",
                                              ["0,0", "1/2,1/2", "1/2,-1/2"]) == [1.0, 2.0, 2.0]
    assert tools.continuum_completion_weights(base + "[113,113]/self:helicity",
                                              ["0,0", "1,1", "1,-1", "1,0"]) == [1.0, 2.0, 2.0, 4.0]
    assert tools.continuum_completion_weights(base + "[113,223]/self:helicity", ["1,0"]) == [2.0]
    assert tools.continuum_completion_weights(base + "[2212,2212]/opposite:g_ls", ["1,1", "3,1"]) == [1.0, 1.0]
    assert tools.continuum_completion_weights(f"RES|f0_980:{model}:g_ls", ["0,0", "2,2"]) is None


# Check the independent analytic projection also participates in the compact helicity orbit
def test_regge_helicity_completion_weights():
    assert tools.continuum_completion_weights("CON_GP|990:[113,113]/self:helicity",
                                              ["0,0,0", "1,0,0", "1,-1,0", "1,-1,1"]) == [1.0, 4.0, 2.0, 4.0]


@pytest.mark.parametrize(
    "phi",
    [
        0.0,
        1.0e-9,
        -1.0e-9,
        0.7,
        -1.3,
        math.pi - 1.0e-9,
        -math.pi + 1.0e-9,
    ],
)
def test_tools_phase_roundtrip_canonical_phase(phi):
    phase_u, phase_v = tools.encode_phase(phi)
    decoded = tools.decode_phase_components(phase_u, phase_v)

    assert decoded == pytest.approx(tools.canonicalize_phase(phi))
    assert -math.pi <= decoded < math.pi
    assert -1.0 <= phase_u <= 1.0
    assert -1.0 <= phase_v <= 1.0


def test_cayley_zero_anchor():
    assert tools.decode_phase_components(0.0, 0.0) == pytest.approx(0.0)
    assert tools.decode_phase_components(1.0e-20, -1.0e-20) == pytest.approx(0.0)


# Check signed real amplitudes use the canonical phase boundary
def test_tools_real_phase_is_canonical():
    assert tools.real_phase(1.0) == 0.0
    assert tools.real_phase(0.0) == 0.0
    assert tools.real_phase(-0.0) == 0.0
    assert tools.real_phase(-1.0) == -math.pi
    with pytest.raises(ValueError, match="finite"):
        tools.real_phase(math.inf)


@pytest.mark.parametrize(
    "u,v",
    [
        (0.0, 0.0),
        (0.35, 0.35),
        (-0.25, 0.6),
        (1.0, 0.0),
        (-1.0, 0.0),
        (1.0, 1.0),
        (-1.0, -1.0),
    ],
)
def test_cayley_product_angle(u, v):
    decoded = tools.decode_phase_components(u, v)
    expected = tools.canonicalize_phase(2.0 * math.atan(u) + 2.0 * math.atan(v))

    assert decoded == pytest.approx(expected)


# Check physical angular projective labels retain card row content and order
def test_projective_topology_angles():
    names = [
        "DECAY|9080225:[113,113]:alpha_ls(2,0)@PROJECTIVE",
        "DECAY|9080225:[113,113]:alpha_ls(2,2)@PROJECTIVE",
        "DECAY|9080225:[113,113]:alpha_ls(4,2)@PROJECTIVE",
    ]

    base, label = tools.parse_projective_angle_key(names[1])
    topology = tools.build_parameter_topology(names)

    assert base == "DECAY|9080225:[113,113]:alpha_ls"
    assert label == "2,2"
    assert topology["groups"][0]["parameters"] == names


# Check physical magnitude names have an explicit reversible marker
def test_physical_angular_magnitude_marker_roundtrip():
    base = "RES|f2_1270:GP:g_ls(2,0)"
    key = tools.magnitude_key(base)

    assert key == "RES|f2_1270:GP:g_ls(2,0)@MAG"
    assert tools.is_magnitude_key(key)
    assert tools.magnitude_base_key(key) == base


# Check a projective production norm retains its physical reference LS row
def test_projective_norm_marker_round_trip():
    key = tools.projective_norm_key("RES|f2_1270:GP:g_ls(0,2)")

    assert key == "RES|f2_1270:GP:g_ls(0,2)@NORM"
    assert tools.is_projective_norm_key(key)
    assert tools.parse_projective_norm_key(key) == ("RES|f2_1270:GP:g_ls", "0,2")


# Check the common Cartesian coefficient is not labeled as a complex norm
def test_projective_coeff_marker_roundtrip():
    reference = "RES|f2_1270:GP:g_ls(0,2)"
    re_key, im_key = tools.projective_coefficient_keys(reference)

    assert re_key == f"{reference}@COEFF_RE"
    assert im_key == f"{reference}@COEFF_IM"
    assert tools.is_projective_coefficient_key(re_key)
    assert tools.is_projective_coefficient_key(im_key)
    assert tools.parse_projective_coefficient_key(re_key) == (
        "RES|f2_1270:GP:g_ls",
        "0,2",
        "re",
    )
    assert tools.parse_projective_coefficient_key(im_key) == (
        "RES|f2_1270:GP:g_ls",
        "0,2",
        "im",
    )


# Check the projector embedding identifies opposite normalized directions
def test_projective_embedding_antipodes():
    vector = [0.2, -0.3, 0.7]
    opposite = [-value for value in vector]

    assert tools.projective_embedding_from_vector(vector) == pytest.approx(
        tools.projective_embedding_from_vector(opposite)
    )


# Check projector inner products equal squared physical direction overlaps
def test_projective_squared_overlap():
    left = tools.canonicalize_projective_vector([0.2, -0.3, 0.7])
    right = tools.canonicalize_projective_vector([-0.4, 0.8, 0.1])
    left_features = tools.projective_embedding_from_vector(left)
    right_features = tools.projective_embedding_from_vector(right)

    feature_overlap = sum(a * b for a, b in zip(left_features, right_features, strict=False))
    vector_overlap = sum(a * b for a, b in zip(left, right, strict=False))

    assert feature_overlap == pytest.approx(vector_overlap**2)


# Check spherical coordinates retain rather than identify the global sign
def test_spherical_angles_antipodal_vectors():
    vector = [0.2, -0.3, 0.7]
    opposite = [-value for value in vector]
    left = tools.spherical_vector_from_angles(tools.spherical_angles_from_vector(vector))
    right = tools.spherical_vector_from_angles(tools.spherical_angles_from_vector(opposite))

    assert left == pytest.approx(tools.normalize_spherical_vector(vector))
    assert right == pytest.approx(tools.normalize_spherical_vector(opposite))
    assert left != pytest.approx(right)


# Check physical spherical labels build one oriented sphere topology group
def test_spherical_topology_angles():
    names = [
        "RES|f0_980:GP:g_ls(2,2)@SPHERICAL",
        "RES|f0_980:GP:g_ls(4,4)@SPHERICAL",
    ]

    topology = tools.build_parameter_topology(names)

    assert topology["groups"] == [
        {
            "kind": "sphere",
            "base": "RES|f0_980:GP:g_ls",
            "parameters": names,
        }
    ]


# Check a resonance norm, phase and spherical direction form one physical group
def test_gp_res_polar_spherical_topology_joint():
    base = "RES|f0_980:GP:g_ls"
    norm = tools.projective_norm_key(f"{base}(0,0)")
    phase = tools.raw_phase_key("RES|f0_980:GP:phi")
    angles = [f"{base}({L},0){tools.SPHERICAL_VECTOR_SUFFIX}" for L in (2, 4)]

    topology = tools.build_parameter_topology([angles[0], phase, norm, angles[1]])

    assert topology["groups"] == [
        {
            "kind": "polar_sphere",
            "base": base,
            "parameters": [norm, phase, *angles],
            "period": 2.0 * math.pi,
        }
    ]


# Check direct GP resonance rows share their common physical phase
def test_gp_joint_magnitudes_phase():
    res_phase = tools.raw_phase_key("RES|f0_980:GP:phi")
    res_rows = [tools.magnitude_key(f"RES|f0_980:GP:g_ls({L},0)") for L in (0, 2)]
    con_rows = [tools.magnitude_key(f"CON_GP|990:[211,211]/opposite:helicity(0,0,{m})") for m in (-1, 0)]
    con_complex = [
        key
        for m in (-1, 0)
        for key in (
            tools.complex_re_key(f"CON_GP|990:[211,211]/opposite:helicity(0,0,{m})"),
            tools.complex_im_key(f"CON_GP|990:[211,211]/opposite:helicity(0,0,{m})"),
        )
    ]

    topology = tools.build_parameter_topology([*res_rows, res_phase, *con_rows, *con_complex])
    groups = {group["kind"]: group for group in topology["groups"]}

    assert groups["polar_components"]["parameters"] == [res_phase, *res_rows]
    assert all(name not in {parameter for group in topology["groups"] for parameter in group["parameters"]} for name in [*con_rows, *con_complex])
