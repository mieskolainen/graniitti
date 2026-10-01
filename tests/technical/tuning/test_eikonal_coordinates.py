# Test canonical eikonal optimizer coordinates
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest
from core.tune.drivers.graniitti.tunesetup import eikonal


# Check descending coupling coordinates are unique and reversible
def test_descending_coupling_round_trip():
    couplings = [9.0, 6.0, 1.5]
    encoded = eikonal.encode_descending_couplings(couplings)

    assert encoded == pytest.approx([9.0, 2.0 / 3.0, 0.25])
    assert eikonal.decode_descending_couplings(encoded) == pytest.approx(couplings)


# Check noncanonical coupling order is rejected by the inverse transform
def test_descending_coupling_order():
    with pytest.raises(ValueError, match="not descending"):
        eikonal.encode_descending_couplings([7.5, 9.0])


# Check ordered DPOW coordinates cover the bounded scale triangle
@pytest.mark.parametrize("scales", [(0.01, 0.01), (0.8, 1.2), (1.5, 1.75), (5.0, 5.0)])
def test_ordered_dpow_round_trip(scales):
    encoded = eikonal.encode_ordered_dpow_scales(*scales, limit=5.0)
    assert eikonal.decode_ordered_dpow_scales(*encoded, limit=5.0) == pytest.approx(scales)


# Check marked optimizer keys retain their flat steering-card targets
def test_canonical_coordinate_key_round_trip():
    coupling = "SOFT|MODEL.triple:EXCHANGE.P.g[1,1]"
    transition = "SOFT|MODEL.triple:EXCHANGE.O.g[0,2]"
    dpow = "SOFT|MODEL.triple:FF.P.param[2,0]"

    assert (
        eikonal.descending_coupling_target(eikonal.descending_coupling_key(coupling))
        == coupling
    )
    assert (
        eikonal.symmetric_coupling_target(eikonal.symmetric_coupling_key(transition))
        == transition
    )
    assert eikonal.ordered_dpow_target(eikonal.ordered_dpow_key(dpow, upper=5.0)) == dpow


# Decode saved independent DPOW groups using each group's own fit limit
def test_saved_dpow_limits_are_independent():
    from core.tune.drivers.graniitti.driver import GraniittiDriver

    values = {}
    for bank, limit in (("P", 7.0), ("O", 11.0)):
        for column, value in enumerate((1.0, 0.5)):
            values[eikonal.ordered_dpow_key(f"SOFT|MODEL.single:FF.{bank}.param[0,{column}]", limit)] = value
    driver = GraniittiDriver()
    decoded, groups = driver._prepare_eikonal_canonical_params(values)
    assert decoded["SOFT|MODEL.single:FF.P.param[0,1]"] == pytest.approx(4.0)
    assert decoded["SOFT|MODEL.single:FF.O.param[0,1]"] == pytest.approx(6.0)
    for group in groups:
        assert group.encode([decoded[key] for key in group.targets]) == pytest.approx((1.0, 0.5))
