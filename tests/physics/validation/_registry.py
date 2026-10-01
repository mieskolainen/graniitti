# Physics-validation evidence for every registered GRANIITTI process
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from dataclasses import dataclass


@dataclass(frozen=True)
class ValidationClaim:
    """Describe comparison provenance and configured checks, without claiming a passing result"""

    path: str | None
    scope: str
    note: str


@dataclass(frozen=True)
class ValidationEvidence:
    """Describe integrated, differential and iceplot coverage for one process"""

    integrated: ValidationClaim
    differential: ValidationClaim
    iceplot: ValidationClaim
    references: tuple[str, ...]
    limitation: str | None


PHOTON_BENCHMARK = "icepack/GAMMA/continuum/dataset.json"
PHOTON_HIGGS_BENCHMARK = "icepack/GAMMA/resonances/higgs/dataset.json"
PHOTON_MONOPOLIUM_BENCHMARK = "icepack/GAMMA/resonances/monopolium/dataset.json"
HARDPOMERON_BENCHMARK = "icepack/HARDPOM/z/dataset.json"
REGGE_BENCHMARK = "icepack/SOFTCEP/processes/CON_all_pipi/dataset.json"
DURHAM_FLUX_BENCHMARK = "icepack/DURHAM/flux/dataset.json"
MINBIAS_BENCHMARK = "icepack/MINBIAS/nd/dataset.json"
DIFFRACTION_BENCHMARK = "icepack/MINBIAS/CMS_1356998/dataset.json"
UPSILON_BENCHMARK = "icepack/PHOTOPROD/LHCb_1373746/dataset.json"
CHARMONIUM_BENCHMARK = "icepack/PHOTOPROD/LHCb_1277076/jpsi/dataset.json"
ELASTIC_BENCHMARK = "icepack/ELASTIC/STAR_1791591/dataset.json"
CHIC_BENCHMARK = "icepack/DURHAM/chic/table8/dataset.json"
CHIC_DATA_BENCHMARK = "icepack/DURHAM/chic_oxford_thesis/dataset.json"
CHIC_DIFFERENTIAL_BENCHMARK = "icepack/DURHAM/chic/figure8a/dataset.json"
DURHAM_PARTON_BENCHMARK = "icepack/DURHAM/partons/dataset.json"
DURHAM_MESON_PAIR_BENCHMARK = "icepack/DURHAM/mmbar/dataset.json"
PHOTON_SUPERCHIC_BENCHMARK = "icepack/GAMMA/superchic2/table6/bare/dataset.json"
PHOTON_WW_BENCHMARK = "icepack/GAMMA/integrated/ATLAS_2021/dataset.json"
UPC_LIGHT_BY_LIGHT_BENCHMARK = "icepack/UPC/GAMMA/ATLAS_1811464/gammagamma/dataset.json"
TP_PHOTO_H1_BENCHMARK = "icepack/PHOTOPROD/H1_1798511/dataset.json"
VALIDATION_PREFIX = "tests/physics/validation"
VALIDATION_SCOPES = frozenset({"direct_data", "literature_model", "execution", "overlay", "unavailable"})


# Build one immutable evidence record
def evidence(
    path,
    *references,
    integrated_path=None,
    differential_path=None,
    iceplot_path=None,
    integrated_scope="execution",
    differential_scope="execution",
    integrated_note="Finite VEGAS rate and weighted-sample normalization closure",
    differential_note="Generated distribution support and internal shape closure",
    limitation=None,
):
    if path is None:
        claim = ValidationClaim(None, "unavailable", "No active icepack generates this process")
        note = "No active icepack validation is available"
        if limitation:
            note += f". {limitation}"
        return ValidationEvidence(claim, claim, claim, tuple(references), note)
    integrated_path = path if integrated_path is None else integrated_path
    differential_path = path if differential_path is None else differential_path
    iceplot_path = differential_path if iceplot_path is None else iceplot_path
    return ValidationEvidence(
        integrated=ValidationClaim(integrated_path, integrated_scope, integrated_note),
        differential=ValidationClaim(differential_path, differential_scope, differential_note),
        iceplot=ValidationClaim(
            iceplot_path,
            "execution",
            "iceplot reads the generated HepMC sample and writes the expected plots",
        ),
        references=tuple(references),
        limitation=limitation,
    )


VALIDATION_COVERAGE = {
    "yy[Higgs]": evidence(
        PHOTON_HIGGS_BENCHMARK,
        "d'Enterria, Martins and Rebello Teles, Phys. Rev. D 101 (2020) 033009, arXiv:1904.11936",
    ),
    "yy[monopolium(0)]": evidence(
        PHOTON_MONOPOLIUM_BENCHMARK,
        "Epele et al., Eur. Phys. J. C 62 (2009) 587, arXiv:0809.0272",
    ),
    "yy[EPA]": evidence(
        PHOTON_SUPERCHIC_BENCHMARK,
        "Harland-Lang, Khoze and Ryskin, Eur. Phys. J. C 76 (2016) 9, arXiv:1508.02718",
        integrated_scope="literature_model",
        integrated_note="Fiducial rate compared with the published SuperChic 2 table",
        differential_note="iceplot distribution only, since no matching data shape is used",
    ),
    "yy[QED]": evidence(
        PHOTON_BENCHMARK,
        "Budnev et al., Phys. Rept. 15 (1975) 181",
        integrated_note="Finite exact QED rate, without a matching external assertion",
        differential_note="iceplot distribution only, since no matching data shape is used",
    ),
    "yy[yy]": evidence(
        UPC_LIGHT_BY_LIGHT_BENCHMARK,
        "ATLAS, JHEP 03 (2021) 243, arXiv:2008.05355",
        "Hirschi et al., JHEP 05 (2011) 044, arXiv:1103.0621",
        integrated_scope="direct_data",
        differential_scope="direct_data",
        integrated_note="Fiducial PbPb rate compared with ATLAS Table 11",
        differential_note=(
            "Scattering angle, mean photon transverse momentum and rapidity "
            "compared bin by bin with ATLAS Tables 3, 5 and 7"
        ),
    ),
    "yy[Zjj]": evidence(
        HARDPOMERON_BENCHMARK,
        "Bowser-Chao and Cheung, Phys. Rev. D 48 (1993) 89, arXiv:hep-ph/9301252",
    ),
    "yy[jj]": evidence(
        None,
        limitation="No independent process-specific literature validation is available",
    ),
    "yy[WW]": evidence(
        PHOTON_WW_BENCHMARK,
        "ATLAS, Phys. Lett. B 816 (2021) 136190, arXiv:2010.04019",
        integrated_scope="direct_data",
        integrated_note="Fiducial e-mu-neutrino rate compared with the ATLAS 13 TeV measurement",
        limitation="The comparison uses the modeled elastic, single dissociative and double dissociative track veto acceptance",
    ),
    "yy[FLUX]": evidence(
        None,
        "Budnev et al., Phys. Rept. 15 (1975) 181",
    ),
    "yy_DZ[EPA]": evidence(
        PHOTON_BENCHMARK,
        "Drees and Zeppenfeld, Phys. Rev. D 39 (1989) 2536",
        "Budnev et al., Phys. Rept. 15 (1975) 181",
    ),
    "yy_DZ[FLUX]": evidence(
        None,
        "Drees and Zeppenfeld, Phys. Rev. D 39 (1989) 2536",
    ),
    "yy_DZ[Zjj]": evidence(
        HARDPOMERON_BENCHMARK,
        "Drees and Zeppenfeld, Phys. Rev. D 39 (1989) 2536",
        "Bowser-Chao and Cheung, Phys. Rev. D 48 (1993) 89, arXiv:hep-ph/9301252",
    ),
    "yy_DZ[jj]": evidence(
        None,
        "Drees and Zeppenfeld, Phys. Rev. D 39 (1989) 2536",
        limitation="No independent process-specific literature validation is available",
    ),
    "yy_DZ[WW]": evidence(
        None,
        "Drees and Zeppenfeld, Phys. Rev. D 39 (1989) 2536",
        limitation="No independent full decay process validation is available",
    ),
    "yy_LUX[EPA]": evidence(
        PHOTON_BENCHMARK,
        "Manohar et al., Phys. Rev. Lett. 117 (2016) 242002",
        "Budnev et al., Phys. Rept. 15 (1975) 181",
    ),
    "yy_LUX[Zjj]": evidence(
        HARDPOMERON_BENCHMARK,
        "Manohar et al., Phys. Rev. Lett. 117 (2016) 242002",
        "Bowser-Chao and Cheung, Phys. Rev. D 48 (1993) 89, arXiv:hep-ph/9301252",
    ),
    "yy_LUX[jj]": evidence(
        None,
        "Manohar et al., Phys. Rev. Lett. 117 (2016) 242002",
        limitation="No independent process-specific literature validation is available",
    ),
    "yy_LUX[WW]": evidence(
        None,
        "Manohar et al., Phys. Rev. Lett. 117 (2016) 242002",
        limitation="No independent full decay process validation is available",
    ),
    "ygg[Z]": evidence(
        HARDPOMERON_BENCHMARK,
        "Cisek, Schafer and Szczurek, Phys. Rev. D 80 (2009) 074013",
    ),
    "ygg[jpsi]": evidence(
        CHARMONIUM_BENCHMARK,
        "LHCb, J. Phys. G 41 (2014) 055002, arXiv:1401.3288",
        integrated_scope="direct_data",
        differential_scope="direct_data",
        integrated_note="Fiducial cross section compared with LHCb Table 1",
        differential_note="Rapidity bins compared with LHCb Table 2",
    ),
    "ygg[psi(2S)]": evidence(
        "icepack/PHOTOPROD/LHCb_1277076/psi2S/dataset.json",
        "LHCb, J. Phys. G 41 (2014) 055002, arXiv:1401.3288",
        integrated_scope="direct_data",
        differential_scope="direct_data",
        integrated_note="Fiducial cross section compared with LHCb Table 1",
        differential_note="Rapidity bins compared with LHCb Table 2",
    ),
    "ygg[Upsilon(1S)]": evidence(
        UPSILON_BENCHMARK,
        "LHCb, JHEP 09 (2015) 084, arXiv:1505.08139",
        integrated_scope="direct_data",
        differential_scope="direct_data",
        integrated_note="Fiducial cross section compared with LHCb Table 4",
        differential_note="Rapidity spectrum compared with LHCb Table 5",
    ),
    "ygg[Upsilon(2S)]": evidence(
        None,
        "LHCb, JHEP 09 (2015) 084, arXiv:1505.08139",
        integrated_scope="direct_data",
        integrated_note="Fiducial cross section compared with LHCb Table 4",
        differential_note="Generated spectrum only, since Table 5 is for Upsilon(1S)",
    ),
    "ygg[Upsilon(3S)]": evidence(
        None,
        "LHCb, JHEP 09 (2015) 084, arXiv:1505.08139",
        integrated_note="Generated rate only, since the published value is an upper limit",
        differential_note="Generated spectrum only, since Table 5 is for Upsilon(1S)",
        limitation=("No direct Upsilon(3S) normalization or differential assertion is implemented"),
    ),
    "MP[RES]": evidence(
        "icepack/SOFTCEP/processes/RES_all_f0_980/dataset.json",
        limitation="No independent process-specific literature validation is available",
    ),
    "MP[CON]": evidence(
        REGGE_BENCHMARK,
        limitation="No independent process-specific literature validation is available",
    ),
    "MP[RES+CON]": evidence(
        "icepack/SOFTCEP/processes/RES+CON_all_pipi/dataset.json",
        limitation="No independent process-specific literature validation is available",
    ),
    "XP[RES]": evidence(
        "icepack/SOFTCEP/processes/RES_all_f0_980/dataset.json",
        "Close and Schuler, Phys. Lett. B 458 (1999) 127",
    ),
    "XP[CON]": evidence(
        REGGE_BENCHMARK,
        "Close and Schuler, Phys. Lett. B 458 (1999) 127",
    ),
    "XP[RES+CON]": evidence(
        "icepack/SOFTCEP/processes/RES+CON_all_pipi/dataset.json",
        "Close and Schuler, Phys. Lett. B 458 (1999) 127",
    ),
    "GP[RES]": evidence(
        "icepack/SOFTCEP/processes/RES_all_f0_980/dataset.json",
        limitation=("No independent validation of the analytic GP resonance model is available"),
    ),
    "GP[CON]": evidence(
        REGGE_BENCHMARK,
        "Harland-Lang, Khoze and Ryskin, Eur. Phys. J. C 74 (2014) 2848",
    ),
    "GP[RES+CON]": evidence(
        "icepack/SOFTCEP/processes/RES+CON_all_pipi/dataset.json",
        "Harland-Lang, Khoze and Ryskin, Eur. Phys. J. C 74 (2014) 2848",
        limitation=("The reference covers the continuum, not the analytic GP resonance model"),
    ),
    "TP[RES]": evidence(
        "icepack/SOFTCEP/processes/RES_all_f0_980/dataset.json",
        "Lebiedowicz, Nachtmann and Szczurek, Phys. Rev. D 93 (2016) 054015, arXiv:1601.04537",
    ),
    "TP[CON]": evidence(
        REGGE_BENCHMARK,
        "Lebiedowicz, Nachtmann and Szczurek, Phys. Rev. D 93 (2016) 054015, arXiv:1601.04537",
    ),
    "TP[RES+CON]": evidence(
        "icepack/SOFTCEP/processes/RES+CON_all_pipi/dataset.json",
        "Lebiedowicz, Nachtmann and Szczurek, Phys. Rev. D 93 (2016) 054015, arXiv:1601.04537",
    ),
    "TP[PHOTO]": evidence(
        TP_PHOTO_H1_BENCHMARK,
        "Lebiedowicz, Nachtmann and Szczurek, arXiv:2508.06334",
        "Lebiedowicz, Nachtmann and Szczurek, arXiv:1804.04706",
        "H1, Eur. Phys. J. C 80 (2020) 1189, arXiv:2005.14471",
        integrated_scope="direct_data",
        differential_scope="direct_data",
        integrated_note="Fiducial dipion rate compared with the integral of the H1 mass spectrum",
        differential_note="Elastic dipion mass spectrum compared bin by bin with H1",
        limitation="This benchmark has no external pp, kaon or nuclear rate comparison",
    ),
    "X[EL]": evidence(
        ELASTIC_BENCHMARK,
        "STAR, Phys. Lett. B 808 (2020) 135663, arXiv:2003.12136",
        "D0, Phys. Rev. D 86 (2012) 012009, arXiv:1206.0687",
        "TOTEM, EPL 95 (2011) 41001",
        "TOTEM, EPL 101 (2013) 21002",
        "TOTEM, Eur. Phys. J. C 79 (2019) 861, arXiv:1812.08283",
        differential_scope="direct_data",
        integrated_note="Generator execution only, since the HEPData cases are differential",
        differential_note="Elastic distributions compared bin by bin with HEPData",
    ),
    "X[ND]": evidence(
        MINBIAS_BENCHMARK,
        "TOTEM, Eur. Phys. J. C 79 (2019) 103",
        "ALICE, Phys. Lett. B 753 (2016) 319",
        integrated_scope="overlay",
        differential_scope="overlay",
        integrated_note=(
            "Global eikonal total and inelastic rates are checked, not the X[ND] rate"
        ),
        differential_note=(
            "ALICE all-INEL primary-charged density is recorded only as a mismatched overlay"
        ),
        limitation="No process-matched X[ND] literature measurement is available",
    ),
    "X[SD]": evidence(
        DIFFRACTION_BENCHMARK,
        "CMS, Phys. Rev. D 92 (2015) 012003, arXiv:1503.08689",
        "Gribov, Sov. Phys. JETP 26 (1968) 414",
        integrated_scope="execution",
        differential_scope="direct_data",
        integrated_note="Generated rate only, no CMS Table 4 integrated measurement is selected",
        differential_note="Inclusive SD + DD + ND sum is tested against published CMS bins",
    ),
    "X[DD]": evidence(
        DIFFRACTION_BENCHMARK,
        "CMS, Phys. Rev. D 92 (2015) 012003, arXiv:1503.08689",
        "Khoze, Martin and Ryskin, Eur. Phys. J. C 73 (2013) 2503",
        integrated_scope="execution",
        differential_scope="direct_data",
        integrated_note="Generated rate only, no CMS Table 4 integrated measurement is selected",
        differential_note="Inclusive SD + DD + ND sum is tested against published CMS bins",
    ),
    "gg[chic(0)]": evidence(
        CHIC_BENCHMARK,
        "Harland-Lang et al., Eur. Phys. J. C 76 (2016) 9",
        "Harland-Lang et al., Eur. Phys. J. C 65 (2010) 433, arXiv:0909.4748",
        integrated_scope="overlay",
        differential_scope="overlay",
        differential_path=CHIC_DIFFERENTIAL_BENCHMARK,
        integrated_note="Fiducial spin rates are compared with the SuperChic 2 table",
        differential_note="Helicity acceptance is compared with SuperChic Figure 8a",
    ),
    "gg[chic(1)]": evidence(
        CHIC_DATA_BENCHMARK,
        "B. R. Gruberg Cazon, University of Oxford PhD thesis, CERN-THESIS-2021-222",
        "Harland-Lang et al., Eur. Phys. J. C 76 (2016) 9",
        "Harland-Lang et al., Eur. Phys. J. C 65 (2010) 433, arXiv:0909.4748",
        integrated_scope="direct_data",
        differential_scope="overlay",
        differential_path=CHIC_DIFFERENTIAL_BENCHMARK,
        integrated_note="Fiducial cross section compared with the Oxford LHCb thesis",
        differential_note="Helicity acceptance is compared with SuperChic Figure 8a",
    ),
    "gg[chic(2)]": evidence(
        CHIC_DATA_BENCHMARK,
        "B. R. Gruberg Cazon, University of Oxford PhD thesis, CERN-THESIS-2021-222",
        "Harland-Lang et al., Eur. Phys. J. C 76 (2016) 9",
        "Harland-Lang et al., Eur. Phys. J. C 65 (2010) 433, arXiv:0909.4748",
        integrated_scope="direct_data",
        differential_scope="overlay",
        differential_path=CHIC_DIFFERENTIAL_BENCHMARK,
        integrated_note="Fiducial cross section compared with the Oxford LHCb thesis",
        differential_note="Helicity acceptance is compared with SuperChic Figure 8a",
    ),
    "gg[QCD]": evidence(
        DURHAM_PARTON_BENCHMARK,
        "Harland-Lang et al., Eur. Phys. J. C 76 (2016) 9",
        "Harland-Lang et al., arXiv:1105.1626",
        integrated_scope="execution",
        differential_scope="overlay",
        integrated_note="Finite generated parton rate without an external table assertion",
        differential_note="Mass spectra are compared with SuperChic Figure 3",
        limitation="The removed SuperChic 2 Table 1 comparison is not claimed",
    ),
    "gg[MM]": evidence(
        DURHAM_MESON_PAIR_BENCHMARK,
        "Harland-Lang et al., arXiv:1105.1626",
        limitation="Meson spectra use internal normalization and shape closure",
    ),
    "gg[yy]": evidence(
        None,
        limitation="The analytic quark loop has amplitude tests but no active diphoton icepack",
    ),
    "gg[FLUX]": evidence(
        DURHAM_FLUX_BENCHMARK,
        "Khoze, Martin and Ryskin, Eur. Phys. J. C 23 (2002) 311",
    ),
    "IPp[Z]": evidence(
        HARDPOMERON_BENCHMARK,
        "Ceccopieri, Eur. Phys. J. C 77 (2017) 56, arXiv:1606.06134",
    ),
    "IPIP[Z]": evidence(
        HARDPOMERON_BENCHMARK,
        "CDF, Phys. Rev. D 82 (2010) 112004, arXiv:1007.5048",
    ),
    "IPp[Zj]": evidence(
        HARDPOMERON_BENCHMARK,
        limitation="No independent process-specific literature validation is available",
    ),
    "IPIP[Zj]": evidence(
        None,
        limitation="No independent process-specific literature validation is available",
    ),
    "IPp[jj]": evidence(
        None,
        limitation="No independent process-specific literature validation is available",
    ),
    "IPIP[jj]": evidence(
        "icepack/HARDPOM/jj/dataset.json",
        limitation="No independent process-specific literature validation is available",
    ),
    "IPp[W]": evidence(
        None,
        limitation="No independent process-specific literature validation is available",
    ),
    "IPIP[W]": evidence(
        None,
        limitation="No independent process-specific literature validation is available",
    ),
}

ALLOWED_EXCLUSIONS = {}


# Compute all claim paths with one external-comparison scope
def paths_with_scope(scope):
    return frozenset(
        claim.path
        for item in VALIDATION_COVERAGE.values()
        for claim in (item.integrated, item.differential)
        if claim.path is not None and claim.scope == scope
    )


PUBLICATION_DATASETS = frozenset(
    {
        "icepack/SOFTCEP/integrated/STAR_2020/dataset.json",
        "icepack/SOFTCEP/integrated/CMS_2017/dataset.json",
        "icepack/SOFTCEP/integrated/CMS_2019/dataset.json",
        "icepack/SOFTCEP/integrated/ATLAS_2017/dataset.json",
        "icepack/GAMMA/integrated/CDF_2007/dataset.json",
        "icepack/GAMMA/integrated/CDF_2012/dataset.json",
        "icepack/GAMMA/integrated/CMS_2012/dataset.json",
        "icepack/GAMMA/integrated/ATLAS_2015/dataset.json",
        "icepack/GAMMA/integrated/ATLAS_2018/dataset.json",
    }
)
DIRECT_DATA_VALIDATIONS = paths_with_scope("direct_data") | PUBLICATION_DATASETS
LITERATURE_MODEL_VALIDATIONS = paths_with_scope("literature_model") | frozenset(
    {DURHAM_MESON_PAIR_BENCHMARK}
)
LITERATURE_OVERLAY_VALIDATIONS = paths_with_scope("overlay")
LITERATURE_VALIDATIONS = (
    DIRECT_DATA_VALIDATIONS | LITERATURE_MODEL_VALIDATIONS | LITERATURE_OVERLAY_VALIDATIONS
)
CHECKED_IN_REFERENCE_VALIDATIONS = DIRECT_DATA_VALIDATIONS | LITERATURE_MODEL_VALIDATIONS
