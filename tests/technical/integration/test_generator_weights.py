# End-to-end tests for iceplot overlays
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#
# Run with: pytest tests/technical/integration/test_generator_weights.py --run-integration -s -rP

import os

from tests.technical.support.iceplot import (
    N_EVENTS,
    PROCESS,
    assert_histograms_match,
    inputcard,
    read_cross_section_histogram,
    run_graniitti,
    run_iceplot,
)
from tests.technical.support.output import PhysicsOutput

OUTPUT = PhysicsOutput("generator_weights")


def test_weighted_unweighted_iceplot_overlay_matches():
    """
    Generate weighted and unweighted events, overlay them with iceplot, and
    verify that the dPhi_pp cross-section plot content agrees
    """
    suffix = f"{os.getpid()}"
    grid_tag = f"test_iceplot_grid_{suffix}"
    unweighted_tag = f"test_iceplot_unweighted_{suffix}"
    weighted_tag = f"test_iceplot_weighted_{suffix}"
    plot_tag = f"test_iceplot_weighted_unweighted_{suffix}"

    common_args = [
        "-i",
        inputcard,
        "-p",
        PROCESS,
        "-l",
        "false",
        "-h",
        "0",
        "-f",
        "hepmc3",
    ]

    # Initialize the rejection envelope needed by the unweighted sample
    run_graniitti(common_args + ["-w", "false", "-n", "0", "-o", grid_tag])
    run_graniitti(
        common_args
        + [
            "-w",
            "false",
            "-n",
            str(N_EVENTS),
            "-d",
            f"vgrid/{grid_tag}.vgrid",
            "-o",
            unweighted_tag,
            "--set",
            "GENERIC.RNDSEED=27183",
        ]
    )
    run_graniitti(
        common_args
        + [
            "-w", "true", "-n", str(N_EVENTS), "-d", f"vgrid/{grid_tag}.vgrid", "-o", weighted_tag,
            "--set", "GENERIC.RNDSEED=81937",
        ]
    )

    run_iceplot(
        [unweighted_tag, weighted_tag],
        plot_tag,
        ["unweighted", "weighted"],
        output=OUTPUT,
    )

    plot_path = OUTPUT.plot_dir(plot_tag) / "hplot__dPhi_pp.pdf"
    assert plot_path.is_file()

    assert_histograms_match(
        read_cross_section_histogram(unweighted_tag),
        read_cross_section_histogram(weighted_tag),
    )
