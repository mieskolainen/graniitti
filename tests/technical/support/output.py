# Shared output-path configuration for pytest-driven test results
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import os
import pathlib
import re
from dataclasses import dataclass

INTEGRATED_CROSS_SECTION_COLUMNS = [
    {"key": "state", "label": "state"},
    {"key": "phase_space", "label": "phase space"},
    {"key": "source", "label": "source"},
    {"key": "cross_section", "label": "sigma"},
    {"key": "unit", "label": "unit"},
    {"key": "statistical_error", "label": "stat"},
    {"key": "systematic_error", "label": "syst"},
    {"key": "luminosity_error", "label": "lumi"},
    {"key": "integration_error", "label": "integration"},
    {"key": "ratio_to_data", "label": "ratio to data"},
]


# Compute the configured test-output root with the repository default as fallback
def get_test_output_root():
    project_root = pathlib.Path(__file__).resolve().parents[3]
    configured_root = os.environ.get("GRANIITTI_TEST_OUTPUT_DIR")
    if configured_root:
        path = pathlib.Path(configured_root).expanduser()
        if not path.is_absolute():
            path = project_root / path
        return path.resolve()
    return project_root / "tmp" / "test-results"


# Convert one physics-suite or plot name into a stable directory component
def safe_output_name(value):
    name = re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")
    if not name:
        raise ValueError(f"Output name has no filesystem-safe characters: {value!r}")
    return name


@dataclass(frozen=True)
class PhysicsOutput:
    """Provide the standard plots, tables and summaries tree for one suite"""

    suite: str

    # Compute the suite root under the configured physics output directory
    @property
    def root(self):
        return get_test_output_root() / "physics" / safe_output_name(self.suite)

    # Compute the iceplot plot root for this suite
    @property
    def plots(self):
        return self.root / "plots"

    # Compute the literature and closure table root for this suite
    @property
    def tables(self):
        return self.root / "tables"

    # Compute the machine-readable summary root for this suite
    @property
    def summaries(self):
        return self.root / "summaries"

    # Compute the standard JSON report path for one iceplot invocation
    def report(self, plot_tag):
        return self.tables / f"{safe_output_name(plot_tag)}.json"

    # Compute one iceplot scale directory with optional HEPData set components
    def plot_dir(self, plot_tag, *dataset_parts, scale="linear"):
        return self.plots / plot_tag / pathlib.Path(*dataset_parts) / scale


# Convert one scalar table value to a stable text representation
def format_table_value(value):
    """Return one portable tab-separated table cell"""
    if value is None:
        return ""
    if isinstance(value, dict) and set(value) == {"down", "up"}:
        return f"+{value['up']:.4g}/-{value['down']:.4g}"
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


# Escape one formatted table cell for LaTeX tabular output
def format_latex_value(value):
    """Return one escaped LaTeX table cell"""
    formatted = format_table_value(value)
    if formatted.count("$") >= 2 and formatted.count("$") % 2 == 0:
        return formatted
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(character, character) for character in formatted)


# Write one standardized machine-readable physics-validation summary
def write_physics_summary(output, name, summary):
    """Write one versioned summary JSON under the canonical summary directory"""
    output.summaries.mkdir(parents=True, exist_ok=True)
    safe_name = safe_output_name(name)
    path = output.summaries / f"{safe_name}.json"
    payload = {
        "schema_version": 1,
        "kind": "physics_validation_summary",
        "suite": safe_output_name(output.suite),
        "name": safe_name,
        "summary": summary,
    }
    path.write_text(
        json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


# Write matching JSON, TSV and LaTeX representations of one literature table
def write_literature_table(output, name, *, columns, rows, metadata=None):
    """Write one versioned literature-comparison table in three forms"""
    keys = [str(column["key"]) for column in columns]
    if len(keys) != len(set(keys)):
        raise ValueError("Literature table column keys must be unique")
    for row in rows:
        missing = [key for key in keys if key not in row]
        if missing:
            raise ValueError(f"Literature table row is missing columns {missing}")

    output.tables.mkdir(parents=True, exist_ok=True)
    safe_name = safe_output_name(name)
    json_path = output.tables / f"{safe_name}.json"
    text_path = output.tables / f"{safe_name}.txt"
    latex_path = output.tables / f"{safe_name}.tex"
    payload = {
        "schema_version": 1,
        "kind": "literature_comparison_table",
        "suite": safe_output_name(output.suite),
        "name": safe_name,
        "metadata": {} if metadata is None else metadata,
        "columns": columns,
        "rows": [{key: row[key] for key in keys} for row in rows],
    }
    json_path.write_text(
        json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    labels = [
        (
            f"{column.get('label', column['key'])} [{column['unit']}]"
            if column.get("unit")
            else str(column.get("label", column["key"]))
        )
        for column in columns
    ]
    lines = ["\t".join(labels)]
    lines.extend("\t".join(format_table_value(row[key]) for key in keys) for row in rows)
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    alignment = "l" * len(keys)
    latex_lines = [
        rf"\begin{{tabular}}{{{alignment}}}",
        r"\hline",
        " & ".join(format_latex_value(label) for label in labels) + r" \\",
        r"\hline",
    ]
    latex_lines.extend(
        " & ".join(format_latex_value(row[key]) for key in keys) + r" \\"
        for row in rows
    )
    latex_lines.extend((r"\hline", r"\end{tabular}"))
    latex_path.write_text("\n".join(latex_lines) + "\n", encoding="utf-8")
    return json_path, text_path, latex_path
