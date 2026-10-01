# Tests for pytest output paths and physics output helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
import pathlib

import pytest

from tests.technical.support.output import (
    PhysicsOutput,
    get_test_output_root,
    write_literature_table,
    write_physics_summary,
)


# Resolve relative output roots from the repository root
def test_output_root_relative_override(monkeypatch):
    monkeypatch.setenv("GRANIITTI_TEST_OUTPUT_DIR", "tmp/custom-results")
    project_root = pathlib.Path(__file__).resolve().parents[3]
    assert get_test_output_root() == project_root / "tmp" / "custom-results"


# Preserve absolute output roots selected by the caller
def test_output_root_absolute_override(monkeypatch, tmp_path):
    output_root = tmp_path / "results"
    monkeypatch.setenv("GRANIITTI_TEST_OUTPUT_DIR", str(output_root))
    assert get_test_output_root() == output_root


# Write summaries and literature tables with one common versioned schema
def test_standard_physics_output_writers(monkeypatch, tmp_path):
    monkeypatch.setenv("GRANIITTI_TEST_OUTPUT_DIR", str(tmp_path))
    output = PhysicsOutput("Example")
    summary_path = write_physics_summary(output, "Closure", {"pull": 0.5})
    json_path, text_path, latex_path = write_literature_table(
        output,
        "Rates",
        columns=[
            {"key": "process", "label": "process"},
            {"key": "cross_section", "label": "sigma", "unit": "pb"},
        ],
        rows=[{"process": r"$\gamma\gamma$", "cross_section": 1.25}],
        metadata={"reference": "example"},
    )

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    table = json.loads(json_path.read_text(encoding="utf-8"))
    assert summary["kind"] == "physics_validation_summary"
    assert summary["summary"] == {"pull": 0.5}
    assert table["kind"] == "literature_comparison_table"
    assert table["rows"] == [{"process": r"$\gamma\gamma$", "cross_section": 1.25}]
    assert text_path.read_text(encoding="utf-8").splitlines()[0] == ("process\tsigma [pb]")
    latex = latex_path.read_text(encoding="utf-8")
    assert r"\begin{tabular}{ll}" in latex
    assert r"$\gamma\gamma$ & 1.25 \\" in latex


# Reject non-standard non-finite JSON values from validation outputs
def test_output_nonfinite_values(
    monkeypatch,
    tmp_path,
):
    monkeypatch.setenv("GRANIITTI_TEST_OUTPUT_DIR", str(tmp_path))
    output = PhysicsOutput("Example")
    with pytest.raises(ValueError):
        write_physics_summary(output, "invalid", {"pull": math.nan})
