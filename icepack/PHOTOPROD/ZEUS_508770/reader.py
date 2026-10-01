# ZEUS vector meson transfer spectra and correlated uncertainties
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.stats.hist import bins2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, number, read_table, resolve_group, table_bins, table_values


# Read original transfer points with common centre to edge bin conversion
# [REFERENCE: arXiv:hep-ex/9910038, Tables 6-7 and Figure 13]
def read(filename, file_filter=None, rebin_factor=None, **_kwargs):
    if rebin_factor is not None:
        raise ValueError("ZEUS transfer spectra cannot be rebinned")
    table = read_table(filename, dimensions=1)
    group = resolve_group(table, file_filter)
    values = table_values(table, group)
    x, edges = table_bins(table)
    bins = doublet2linear(edges)
    data = {"x": x, "binedges": edges, "bins": bins, "binwidth": bins2binwidth(bins), "xlim": bins[[0, -1]],
            "y": np.asarray([number(value["value"]) for value in values]), "group": group}
    sources = []
    for index, error in enumerate(values[0]["errors"]):
        label = error["label"]
        norm = "normalization" in label
        correlated = norm or label == "sys_2"
        up, down = error_pairs(values, index).T
        sources.append(uncertainty_source(
            "normalization" if norm else label, up, down=down,
            category="statistical" if label == "stat" else "systematic",
            correlation="collective" if correlated else "uncorrelated",
            effect="multiplicative" if norm else "additive",
            scope=f"{Path(filename).parent.name if norm else table['doi']}:{label}",
            provenance=("Published common normalization" if norm else
                        "Same-sign correlated proton dissociation model variation" if correlated else
                        "Published marginal errors, diagonal covariance without supplied bin correlations"),
        ))
    return finalize_uncertainties(data, sources)
