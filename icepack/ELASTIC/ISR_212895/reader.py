# Read ISR elastic pp and ppbar differential cross sections in the diffraction dip
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.hist import bins2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, residual_uncertainty, uncertainty_source

from icepack._common.hepdata import error_pairs, print_rebin, read_table, rebin_common, table_bins, table_values


# Read the published bin edges and correlated pp and ppbar normalization errors
@cache
def read(filename: str, rebin_factor: int | None = None) -> dict:
    path = Path(filename)
    tables = {name: read_table(path.with_name(f"{name}.json")) for name in ("Table1", "Table2")}
    table = tables[path.stem]
    values = table_values(table)
    x, edges = table_bins(table)
    bins = doublet2linear(edges)
    y = np.asarray([value["value"] for value in values], dtype=float)
    errors = error_pairs(values, "error")
    dataset = {"x": x, "bins": bins, "binwidth": bins2binwidth(bins), "xlim": bins[[0, -1]], "y": y}

    # [REFERENCE: Breakstone et al., Phys. Rev. Lett. 54 (1985) 2180, Table I and normalization discussion on p. 2181]
    # Table 2 sys_1 is the relative ppbar/pp normalization, sys_2 is the absolute ppbar normalization
    pp_values, ppbar_values = (table_values(tables[name]) for name in ("Table1", "Table2"))
    pp = error_pairs(pp_values, "sys")[0, 0] / float(pp_values[0]["value"])
    ppbar = error_pairs(ppbar_values, "sys_2")[0, 0] / float(ppbar_values[0]["value"])
    relative = error_pairs(ppbar_values, "sys_1")[0, 0] / float(ppbar_values[0]["value"])
    # Var(delta_ppbar - delta_pp) fixes the covariance between the two fractional scales
    common = pp if path.stem == "Table1" else (pp**2 + ppbar**2 - relative**2) / (2 * pp)
    sources = [
        uncertainty_source(
            "published",
            errors[:, 0],
            down=errors[:, 1],
            category="combined",
            correlation="uncorrelated",
            provenance="Table I errors include statistics and point-to-point correction uncertainties",
        ),
        uncertainty_source(
            "normalization",
            common * y,
            category="systematic",
            correlation="collective",
            scope="ISR_212895:normalization",
            effect="multiplicative",
            provenance="Shared pp and ppbar scale fixed by published absolute and relative normalization errors",
        ),
    ]
    if path.stem == "Table2":
        sources.append(
            uncertainty_source(
                "ppbar_normalization",
                residual_uncertainty(ppbar * y, [common * y]),
                category="systematic",
                correlation="collective",
                scope="ISR_212895:ppbar_normalization",
                effect="multiplicative",
                provenance="Remaining ppbar normalization variance after the shared pp scale",
            )
        )
    finalize_uncertainties(dataset, sources)
    if rebin_factor is not None:
        print_rebin(__name__, rebin_factor)
        dataset = rebin_common(ds=dataset, rebin_factor=rebin_factor)
    return dataset
