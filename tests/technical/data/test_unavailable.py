# Check unavailable measurements cannot invoke publication parsing
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import pytest
from core.io.hepdata_reader import load_dataset_reader

from develop.tools.lib.hera import fit, inputs
from icepack.PHOTOPROD.LHCb_1373746.lhcb_upsilon import energy_fractions


# Reject unsupported rho measurements through the production reader API
@pytest.mark.parametrize("pack", ["CMS_1718344/rho", "STAR_1515028/rho_coherent", "LHCb_2931431/rho_coherent"])
def test_rho_unavailable(pack):
    root = Path(__file__).resolve().parents[3]
    reader = load_dataset_reader(f"icepack/UPC/PHOTOPROD/{pack}/dataset.json", cdir=root)
    with pytest.raises(ValueError, match="HEPData JSON"):
        reader.read(root / "HEPData/missing.pdf")


# Stop unsupported coupling fits before extracting any non-JSON measurements
@pytest.mark.parametrize("read", [inputs.psi2s_ratio, inputs.bottomonium, fit.h1_spectra,
                                 fit.h1_wt_spectra, fit.jpsi_energy, fit.jpsi_dissociation])
def test_fit_unavailable(read):
    with pytest.raises(ValueError, match="HEPData JSON"):
        read()


# Reject a combined beam sample without HEPData effective luminosities
def test_lhcb_fractions_unavailable():
    with pytest.raises(ValueError, match="HEPData JSON"):
        energy_fractions([])
