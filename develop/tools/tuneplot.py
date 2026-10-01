#!/usr/bin/env python3
# Combine icetune linear and logarithmic plot PDFs into one titled PDF
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from core.io.files import ensure_dir


# Parse command-line input and output paths
def parse_args():
    parser = argparse.ArgumentParser(description='Combine icetune plot PDFs, pairing linear and log scales')
    parser.add_argument("input", type=Path, help="icetune plot directory")
    parser.add_argument("-o", "--output", type=Path, help="output PDF path")
    return parser.parse_args()


# Compute a natural-sort key for folder and observable names
def natural_key(value):
    return [
        int(token) if token.isdigit() else token.lower() for token in re.split(r"(\d+)", str(value))
    ]


# Collect plot paths by dataset folder, observable, and scale
def collect(input_dir, output_file):
    groups = {}

    for path in input_dir.rglob("hplot__*.pdf"):
        relative = path.relative_to(input_dir)
        if "default" in relative.parts or path.resolve() == output_file.resolve():
            continue

        scale = path.parent.name
        if scale not in {"linear", "log"}:
            continue

        folder = path.parent.parent.relative_to(input_dir)
        observable = path.stem.removeprefix("hplot__")
        groups.setdefault((folder, observable), {})[scale] = path

    return groups


# Sort observables in a compact physics-oriented order
def plot_order(item):
    folder, observable = item
    preferred = {
        "M": 0,
        "Rap": 1,
        "Abs_t1t2": 2,
        "dPhi_pp": 3,
        "costheta_CS": 4,
        "phi_CS": 5,
    }
    return natural_key(folder), preferred.get(observable, 100), natural_key(observable)


# Escape text for LaTeX titles
def latex_escape(value):
    replacements = {
        "\\": r"\textbackslash{}",
        "{": r"\{",
        "}": r"\}",
        "$": r"\$",
        "&": r"\&",
        "%": r"\%",
        "#": r"\#",
        "_": r"\_",
        "^": r"\textasciicircum{}",
        "~": r"\textasciitilde{}",
        "<": r"\textless{}",
        ">": r"\textgreater{}",
    }
    return "".join(replacements.get(character, character) for character in str(value))


# Link one source plot to a LaTeX-safe temporary filename
def link_plot(source, work_dir, page_index, scale):
    target = work_dir / f"plot_{page_index:04d}_{scale}.pdf"
    if source is not None:
        target.symlink_to(source.resolve())
    return target.name


# Compute LaTeX panel with a vector PDF or a missing marker
def latex_panel(label, source, work_dir, page_index):
    heading = rf"{{\sffamily\fontsize{{11}}{{13}}\selectfont {label}\par}}"
    if source is None:
        body = r"\vfill{\sffamily missing}\vfill"
    else:
        filename = link_plot(source, work_dir, page_index, label)
        body = (
            r"\includegraphics[width=\linewidth,height=4.9in,"
            rf"keepaspectratio]{{{filename}}}"
        )
    return "\n".join(
        [
            r"\begin{minipage}[t]{0.492\textwidth}",
            r"\centering",
            heading,
            r"\vspace{0.03in}",
            body,
            r"\end{minipage}",
        ]
    )


# Compute titled LaTeX page containing the linear-log pair
def latex_page(folder, observable, plots, work_dir, page_index):
    linear = latex_panel("linear", plots.get("linear"), work_dir, page_index)
    logarithmic = latex_panel("log", plots.get("log"), work_dir, page_index)
    return "\n".join(
        [
            r"\begin{center}",
            (
                r"{\sffamily\bfseries\fontsize{15}{17}\selectfont "
                f"{latex_escape(folder)}"
                r"\par}"
            ),
            r"\vspace{0.03in}",
            (
                r"{\sffamily\fontsize{12}{14}\selectfont "
                f"{latex_escape(observable)}"
                r"\par}"
            ),
            r"\vspace{0.06in}",
            linear,
            r"\hfill",
            logarithmic,
            r"\end{center}",
        ]
    )


# Write the complete vector PDF LaTeX document
def write_latex(tex_file, ordered, groups, work_dir):
    pages = []
    for page_index, key in enumerate(ordered):
        folder, observable = key
        pages.append(latex_page(folder, observable, groups[key], work_dir, page_index))

    preamble = "\n".join(
        [
            r"\documentclass{article}",
            r"\usepackage[T1]{fontenc}",
            r"\usepackage{graphicx}",
            r"\usepackage[paperwidth=16in,paperheight=6.4in,margin=0.18in]{geometry}",
            r"\pagestyle{empty}",
            r"\setlength{\parindent}{0pt}",
            r"\begin{document}",
        ]
    )
    document = preamble + "\n" + "\n\\newpage\n".join(pages) + "\n\\end{document}\n"
    tex_file.write_text(document, encoding="ascii")


# Compile the LaTeX book and expose useful diagnostics on failure
def compile_latex(tex_file, work_dir):
    command = [
        "pdflatex",
        "-interaction=nonstopmode",
        "-halt-on-error",
        "-output-directory",
        str(work_dir),
        str(tex_file),
    ]
    result = subprocess.run(
        command,
        cwd=work_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError("pdflatex failed:\n" + result.stdout[-4000:] + result.stderr[-4000:])
    return work_dir / f"{tex_file.stem}.pdf"


# Build the vector PDF and return the number of source plots and pages
def build_book(input_dir, output_file):
    groups = collect(input_dir, output_file)
    if not groups:
        raise RuntimeError(f"No non-default hplot PDFs found under {input_dir}")

    tmp_root = Path(__file__).resolve().parents[2] / "tmp"
    ensure_dir(tmp_root)
    ensure_dir(output_file.parent)

    source_count = sum(len(scales) for scales in groups.values())
    ordered = sorted(groups, key=plot_order)

    with tempfile.TemporaryDirectory(prefix="tuneplot_", dir=tmp_root) as tmp_name:
        work_dir = Path(tmp_name)
        tex_file = work_dir / "tuneplots.tex"
        write_latex(tex_file, ordered, groups, work_dir)
        compiled_pdf = compile_latex(tex_file, work_dir)
        shutil.copy2(compiled_pdf, output_file)

    return source_count, len(ordered)


# Validate input and run the PDF assembly
def main():
    args = parse_args()
    input_dir = args.input.resolve()
    if not input_dir.is_dir():
        raise NotADirectoryError(input_dir)

    output_file = args.output.resolve() if args.output else input_dir / f"{input_dir.name}.pdf"
    source_count, page_count = build_book(input_dir, output_file)
    print(
        f"Wrote {output_file} with {page_count} pages from {source_count} plots (default excluded)"
    )


if __name__ == "__main__":
    main()
