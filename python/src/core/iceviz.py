#!/usr/bin/env python
#
# iceviz: HepMC3 event-record visualization with Graphviz
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


from __future__ import annotations

import argparse
import math
import pathlib
import shutil
import subprocess
import sys
from collections.abc import Iterable
from typing import Any

from core.io.files import ensure_dir


# Load the first available HepMC3 Python binding
def import_hepmc() -> Any:
    try:
        import pyhepmc3  # type: ignore

        return pyhepmc3
    except ImportError:
        pass

    try:
        import pyhepmc  # type: ignore

        return pyhepmc
    except ImportError:
        pass

    try:
        from pyHepMC3 import HepMC3  # type: ignore

        return HepMC3
    except ImportError as exc:
        raise RuntimeError(
            "could not import pyhepmc3, pyhepmc, or pyHepMC3.HepMC3; run source install/setenv.sh first"
        ) from exc


# Read an attribute or method value from either HepMC3 binding
def value(obj: Any, name: str, default: Any = None) -> Any:
    attr = getattr(obj, name, default)
    if callable(attr):
        try:
            return attr()
        except TypeError:
            return default
    return attr


# Read a sequence from a HepMC3 object
def sequence(obj: Any, name: str) -> list[Any]:
    seq = value(obj, name, [])
    if seq is None:
        return []
    return list(seq)


# Read one event with the high-level pyhepmc interface
def read_event_pyhepmc(hepmc: Any, path: pathlib.Path, event_index: int) -> Any:
    reader = hepmc.open(str(path))
    try:
        for index, event in enumerate(reader):
            if index == event_index:
                return event
    finally:
        close = getattr(reader, "close", None)
        if callable(close):
            close()
    raise IndexError(f"event index {event_index} is outside {path}")


# Read one event with the lower-level pyHepMC3 interface
def read_event_pyhepmc3(hepmc: Any, path: pathlib.Path, event_index: int) -> Any:
    reader = hepmc.ReaderAscii(str(path))
    try:
        for _ in range(event_index + 1):
            event = hepmc.GenEvent()
            reader.read_event(event)
            if reader.failed():
                raise RuntimeError(f"failed to parse event {event_index} from {path}")
        return event
    finally:
        close = getattr(reader, "close", None)
        if callable(close):
            close()


# Read one event through the available HepMC3 Python API
def read_event(path: pathlib.Path, event_index: int) -> Any:
    hepmc = import_hepmc()
    if hasattr(hepmc, "open"):
        return read_event_pyhepmc(hepmc, path, event_index)
    if hasattr(hepmc, "ReaderAscii"):
        return read_event_pyhepmc3(hepmc, path, event_index)
    raise RuntimeError("HepMC3 Python binding has no supported ASCII reader")


# Compute a stable identifier for a HepMC3 object
def object_id(obj: Any, fallback: int) -> int:
    try:
        return int(value(obj, "id", fallback))
    except (TypeError, ValueError):
        return fallback


# Compute a compact string representation for one HepMC attribute value
def attribute_text(attr: Any, max_len: int = 80) -> str:
    text = str(attr).strip()
    if len(text) > max_len:
        return text[: max_len - 3] + "..."
    return text


# Collect event, particle and vertex attributes through the HepMC3 binding
def event_attributes(event: Any) -> dict[int, dict[str, str]]:
    owners = [(0, event)]
    owners.extend((object_id(obj, 0), obj) for name in ("particles", "vertices") for obj in sequence(event, name))
    table = {}
    for owner_id, owner in owners:
        if hasattr(event, "attribute_names"):
            attributes = {name: event.attribute_as_string(name, owner_id) for name in event.attribute_names(owner_id)}
        else:
            attributes = dict(value(owner, "attributes", {}).items())
        if attributes:
            table[owner_id] = {
                str(name): attribute_text(attr.astype(str) if hasattr(attr, "astype") else attr)
                for name, attr in attributes.items()
            }
    return table


# Look up HepMC3 attributes by owner id
def attributes_for(attr_table: dict[int, dict[str, str]], owner_id: int) -> dict[str, str]:
    return attr_table.get(owner_id, {})


# Compute the PDG display name used in graph labels
def pdg_name(pid: int) -> str:
    names = {
        1: "d",
        -1: "d~",
        2: "u",
        -2: "u~",
        3: "s",
        -3: "s~",
        4: "c",
        -4: "c~",
        5: "b",
        -5: "b~",
        6: "t",
        -6: "t~",
        11: "e-",
        -11: "e+",
        12: "nu_e",
        -12: "nu_e~",
        13: "mu-",
        -13: "mu+",
        14: "nu_mu",
        -14: "nu_mu~",
        15: "tau-",
        -15: "tau+",
        21: "g",
        22: "gamma",
        23: "Z0",
        24: "W+",
        -24: "W-",
        111: "pi0",
        211: "pi+",
        -211: "pi-",
        221: "eta",
        310: "K0S",
        311: "K0",
        -311: "K0~",
        321: "K+",
        -321: "K-",
        990: "Pomeron",
        1114: "Delta-",
        2114: "Delta0",
        2212: "p",
        -2212: "p~",
    }
    return names.get(pid, str(pid))


# Compute the particle category used for graph styling
def particle_category(pid: int, status: int) -> str:
    apid = abs(pid)
    if status == 4:
        return "beam"
    if pid == 990:
        return "pomeron"
    if apid == 2212:
        return "proton"
    if apid in {11, 13, 15}:
        return "lepton"
    if apid in {12, 14, 16}:
        return "neutrino"
    if apid in {1, 2, 3, 4, 5, 6, 21}:
        return "parton"
    if apid in {22, 23, 24}:
        return "boson"
    if status == 1:
        return "stable"
    return "other"


# Identify particles in the Pythia hardest subprocess
def is_hard_particle(particle: Any) -> bool:
    status = abs(int(value(particle, "status", 0)))
    return 21 <= status <= 29


# Identify the central hard-scattering vertex
def is_hard_scattering_vertex(vertex: Any) -> bool:
    incoming = sequence(vertex, "particles_in")
    outgoing = sequence(vertex, "particles_out")
    if not incoming or not outgoing:
        return False
    has_hard_in = any(abs(int(value(p, "status", 0))) == 21 for p in incoming)
    has_hard_out = any(22 <= abs(int(value(p, "status", 0))) <= 29 for p in outgoing)
    return has_hard_in and has_hard_out


# Identify vertices attached to the hard subprocess
def touches_hard_branch(vertex: Any) -> bool:
    particles = sequence(vertex, "particles_in") + sequence(vertex, "particles_out")
    return any(is_hard_particle(p) for p in particles)


# Compute the Graphviz node style for one particle category
def particle_style(
    particle: Any, category: str, status: int, highlight_hard: bool
) -> dict[str, str]:
    styles = {
        "beam": {"fillcolor": "#d9e8ff", "shape": "box"},
        "pomeron": {"fillcolor": "#c8f1ff", "shape": "diamond"},
        "proton": {"fillcolor": "#d4edda", "shape": "box"},
        "lepton": {"fillcolor": "#fff3bf", "shape": "box"},
        "neutrino": {"fillcolor": "#f1f3f5", "shape": "box"},
        "parton": {"fillcolor": "#ffd8a8", "shape": "box"},
        "boson": {"fillcolor": "#e5dbff", "shape": "box"},
        "stable": {"fillcolor": "#f8f9fa", "shape": "box"},
        "other": {"fillcolor": "#ffffff", "shape": "box"},
    }
    style = dict(styles.get(category, styles["other"]))
    if highlight_hard and is_hard_particle(particle):
        style["color"] = "#9c36b5"
        style["fillcolor"] = "#f3d9fa"
        style["penwidth"] = "2.4"
        style["style"] = "rounded,filled,bold"
        return style
    if status != 1 and category not in {"beam", "pomeron"}:
        style["style"] = "rounded,filled,dashed"
    else:
        style["style"] = "rounded,filled"
    return style


# Compute a four-momentum sum for one particle sequence
def four_momentum_sum(particles: Iterable[Any]) -> tuple[float, float, float, float]:
    px = py = pz = energy = 0.0
    for particle in particles:
        momentum = value(particle, "momentum")
        px += momentum_value(momentum, "px")
        py += momentum_value(momentum, "py")
        pz += momentum_value(momentum, "pz")
        energy += momentum_value(momentum, "e")
    return px, py, pz, energy


# Compute the four-momentum residual for one vertex
def vertex_balance(vertex: Any) -> tuple[tuple[float, float, float, float], float, float, bool]:
    in_particles = sequence(vertex, "particles_in")
    out_particles = sequence(vertex, "particles_out")
    incoming = four_momentum_sum(in_particles)
    outgoing = four_momentum_sum(out_particles)
    residual = tuple(outgoing[i] - incoming[i] for i in range(4))
    max_abs = max(abs(x) for x in residual)
    scale = max(1.0, abs(incoming[3]), abs(outgoing[3]))
    tolerance = 1.0e-2 + 2.0e-4 * scale
    if len(in_particles) == 0 or len(out_particles) == 0:
        return residual, max_abs, tolerance, True
    return residual, max_abs, tolerance, max_abs <= tolerance


# Identify one-to-one nonlocal Pythia copy vertices
def is_nonlocal_copy_vertex(vertex: Any) -> bool:
    incoming = sequence(vertex, "particles_in")
    outgoing = sequence(vertex, "particles_out")
    if len(incoming) != 1 or len(outgoing) != 1:
        return False
    if int(value(incoming[0], "pid", 0)) != int(value(outgoing[0], "pid", 0)):
        return False
    _, _, _, balanced = vertex_balance(vertex)
    return not balanced


# Compute the compact vertex-balance diagnostic text
def vertex_balance_label(vertex: Any, force: bool) -> list[str]:
    incoming = sequence(vertex, "particles_in")
    outgoing = sequence(vertex, "particles_out")
    if len(incoming) == 0 or len(outgoing) == 0:
        if not force:
            return []
        return ["event source" if len(incoming) == 0 else "event sink"]

    residual, max_abs, tolerance, balanced = vertex_balance(vertex)
    copy_vertex = is_nonlocal_copy_vertex(vertex)
    if not force and not copy_vertex:
        return []
    if balanced and not force:
        return []
    prefix = "copy/nonlocal" if copy_vertex else "nonlocal"
    if balanced:
        prefix = "balanced"
    return [
        prefix,
        f"max|dp|={max_abs:.3g}",
        f"tol={tolerance:.3g}",
        f"dpt={math.hypot(residual[0], residual[1]):.3g}",
    ]


# Compute the Graphviz node style for one vertex
def vertex_style(vertex: Any, show_vertex_balance: bool, highlight_hard: bool) -> dict[str, str]:
    _, _, _, balanced = vertex_balance(vertex)
    if highlight_hard and is_hard_scattering_vertex(vertex):
        return {
            "shape": "circle",
            "style": "filled,bold",
            "fillcolor": "#e599f7",
            "color": "#9c36b5",
            "penwidth": "2.8",
            "width": "0.42",
            "height": "0.42",
        }
    if highlight_hard and touches_hard_branch(vertex):
        return {
            "shape": "circle",
            "style": "filled",
            "fillcolor": "#f8f0fc",
            "color": "#9c36b5",
            "penwidth": "1.6",
            "width": "0.35",
            "height": "0.35",
        }
    if balanced or (not show_vertex_balance and not is_nonlocal_copy_vertex(vertex)):
        return {
            "shape": "circle",
            "style": "filled",
            "fillcolor": "#e9ecef",
            "color": "#495057",
            "width": "0.35",
            "height": "0.35",
        }
    return {
        "shape": "circle",
        "style": "filled,bold",
        "fillcolor": "#ffe3e3",
        "color": "#c92a2a",
        "width": "0.35",
        "height": "0.35",
    }


# Read a finite HepMC four-vector component
def momentum_value(momentum: Any, name: str) -> float:
    try:
        return float(value(momentum, name, 0.0))
    except (TypeError, ValueError):
        return 0.0


# Compute a compact kinematic label for one momentum
def momentum_label(momentum: Any) -> str:
    px = momentum_value(momentum, "px")
    py = momentum_value(momentum, "py")
    pz = momentum_value(momentum, "pz")
    e = momentum_value(momentum, "e")
    pt = math.hypot(px, py)
    p = math.sqrt(max(0.0, px * px + py * py + pz * pz))
    if p > abs(pz):
        eta = 0.5 * math.log((p + pz) / max(p - pz, 1.0e-300))
        eta_text = f"{eta:.2f}"
    else:
        eta_text = "inf" if pz >= 0.0 else "-inf"
    return f"pT={pt:.3g} eta={eta_text} E={e:.3g}"


# Quote and escape a DOT string
def dot_quote(text: str) -> str:
    escaped = text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
    return f'"{escaped}"'


# Compute one DOT attribute list
def dot_attributes(attrs: dict[str, str]) -> str:
    return ", ".join(f"{key}={dot_quote(value)}" for key, value in attrs.items())


# Compute DOT edge attributes for hard and nonlocal links
def edge_attributes(particle: Any, vertex: Any, highlight_hard: bool) -> str:
    if is_nonlocal_copy_vertex(vertex):
        return ' [color="#c92a2a", style=dashed, penwidth=1.4]'
    if highlight_hard and (is_hard_particle(particle) or is_hard_scattering_vertex(vertex)):
        return ' [color="#9c36b5", penwidth=2.0]'
    return ""


# Compute the particle label shown in the graph
def particle_label(
    particle: Any,
    attr_table: dict[int, dict[str, str]],
    show_momentum: bool,
    show_attributes: bool,
    highlight_hard: bool,
) -> str:
    pid = int(value(particle, "pid", 0))
    status = int(value(particle, "status", 0))
    lines = [f"{pdg_name(pid)}  [{object_id(particle, 0)}]", f"PDG={pid} st={status}"]
    if highlight_hard and is_hard_particle(particle):
        lines.append("hard subprocess")
    if show_momentum:
        lines.append(momentum_label(value(particle, "momentum")))
    if show_attributes:
        for key, attr in sorted(attributes_for(attr_table, object_id(particle, 0)).items()):
            lines.append(f"{key}={attr}")
    return "\n".join(lines)


# Compute the vertex label shown in the graph
def vertex_label(
    vertex: Any,
    attr_table: dict[int, dict[str, str]],
    show_attributes: bool,
    show_vertex_balance: bool,
    highlight_hard: bool,
) -> str:
    lines = [f"V[{abs(object_id(vertex, 0))}]"]
    if highlight_hard and is_hard_scattering_vertex(vertex):
        lines.append("hard scatter")
    lines.extend(vertex_balance_label(vertex, show_vertex_balance))
    if show_attributes:
        for key, attr in sorted(attributes_for(attr_table, object_id(vertex, 0)).items()):
            lines.append(f"{key}={attr}")
    return "\n".join(lines)


# Compute the graph-level label for one event
def event_label(
    event: Any,
    input_path: pathlib.Path,
    event_index: int,
    attr_table: dict[int, dict[str, str]],
    highlight_hard: bool,
) -> str:
    event_number = value(event, "event_number", event_index)
    particles = sequence(event, "particles")
    vertices = sequence(event, "vertices")
    lines = [
        f"{input_path.name} event {event_index} (HepMC event {event_number})",
        f"particles={len(particles)} vertices={len(vertices)}",
    ]
    event_attrs = attributes_for(attr_table, 0)
    closure = event_attrs.get("graniitti_closure_max_abs")
    if closure is not None:
        tolerance = event_attrs.get("graniitti_closure_tolerance")
        passed = event_attrs.get("graniitti_closure_pass")
        status = "pass" if passed == "1" else "fail" if passed == "0" else "unknown"
        closure_line = f"energy-momentum closure: {status}, max|dP4|={closure} GeV"
        if tolerance is not None:
            closure_line += f", tolerance={tolerance} GeV"
        lines.append(closure_line)
    if highlight_hard:
        lines.append("purple=hard subprocess status 21-29; red=dashed nonlocal copy")
    return "\n".join(lines)


# Compute the selected particles after applying the graph-size limit
def selected_particles(event: Any, max_particles: int) -> list[Any]:
    particles = sequence(event, "particles")
    if max_particles <= 0 or len(particles) <= max_particles:
        return particles
    beams = [p for p in particles if int(value(p, "status", 0)) == 4]
    unstable = [p for p in particles if int(value(p, "status", 0)) != 1]
    stable = [p for p in particles if int(value(p, "status", 0)) == 1]
    ordered = []
    seen = set()
    for particle in beams + unstable + stable:
        key = object_id(particle, id(particle))
        if key in seen:
            continue
        ordered.append(particle)
        seen.add(key)
        if len(ordered) >= max_particles:
            break
    return ordered


# Compute the vertices connected to the selected particle set
def connected_vertices(particles: Iterable[Any]) -> list[Any]:
    vertices = []
    seen = set()
    for particle in particles:
        for vertex in (value(particle, "production_vertex"), value(particle, "end_vertex")):
            if vertex is None:
                continue
            key = object_id(vertex, id(vertex))
            if key in seen:
                continue
            vertices.append(vertex)
            seen.add(key)
    return vertices


# Build one DOT graph for a HepMC3 event
def build_dot(
    event: Any,
    input_path: pathlib.Path,
    event_index: int,
    attr_table: dict[int, dict[str, str]],
    max_particles: int,
    show_momentum: bool,
    show_attributes: bool,
    show_vertex_balance: bool,
    highlight_hard: bool,
) -> str:
    particles = selected_particles(event, max_particles)
    vertices = connected_vertices(particles)
    vertex_ids = {object_id(v, id(v)) for v in vertices}

    lines = [
        "digraph HepMC3Event {",
        '  graph [rankdir=LR, bgcolor="white", pad="0.2", nodesep="0.25", ranksep="0.55",',
        f'         label={dot_quote(event_label(event, input_path, event_index, attr_table, highlight_hard))}, labelloc=t, fontsize=18, fontname="Helvetica"];',
        '  node [fontname="Helvetica", fontsize=10, color="#495057", margin="0.06,0.04"];',
        '  edge [fontname="Helvetica", fontsize=8, color="#495057", arrowsize=0.7];',
    ]

    for vertex in vertices:
        node = f"v{abs(object_id(vertex, id(vertex)))}"
        attrs = vertex_style(vertex, show_vertex_balance, highlight_hard)
        attrs["label"] = vertex_label(
            vertex, attr_table, show_attributes, show_vertex_balance, highlight_hard
        )
        lines.append(f"  {node} [{dot_attributes(attrs)}];")

    for particle in particles:
        pid = int(value(particle, "pid", 0))
        status = int(value(particle, "status", 0))
        category = particle_category(pid, status)
        node = f"p{abs(object_id(particle, id(particle)))}"
        attrs = particle_style(particle, category, status, highlight_hard)
        attrs["label"] = particle_label(
            particle, attr_table, show_momentum, show_attributes, highlight_hard
        )
        lines.append(f"  {node} [{dot_attributes(attrs)}];")

    for particle in particles:
        pid = object_id(particle, id(particle))
        pnode = f"p{abs(pid)}"
        prod = value(particle, "production_vertex")
        end = value(particle, "end_vertex")
        if prod is not None and object_id(prod, id(prod)) in vertex_ids:
            attrs = edge_attributes(particle, prod, highlight_hard)
            lines.append(f"  v{abs(object_id(prod, id(prod)))} -> {pnode}{attrs};")
        elif prod is None:
            source = f"s{abs(pid)}"
            lines.append(f'  {source} [label="", shape=point, width=0.05, color="#adb5bd"];')
            lines.append(f'  {source} -> {pnode} [style=dotted, color="#adb5bd"];')
        if end is not None and object_id(end, id(end)) in vertex_ids:
            attrs = edge_attributes(particle, end, highlight_hard)
            lines.append(f"  {pnode} -> v{abs(object_id(end, id(end)))}{attrs};")

    hidden = len(sequence(event, "particles")) - len(particles)
    if hidden > 0:
        lines.append(
            f'  hidden [label={dot_quote(f"{hidden} particles hidden by --max-particles")}, shape=note, fillcolor="#fff9db", style=filled];'
        )

    lines.append("}")
    return "\n".join(lines) + "\n"


# Resolve the output path selected by CLI options
def output_path(args: argparse.Namespace) -> pathlib.Path:
    if args.output:
        return pathlib.Path(args.output)
    suffix = args.format if args.format else "svg"
    stem = pathlib.Path(args.input).stem
    return pathlib.Path("output") / "iceviz" / f"{stem}_event{args.event}.{suffix}"


# Write one text file with robust parent-directory creation
def write_text(path: pathlib.Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


# Render one DOT file with Graphviz
def render_dot(dot_path: pathlib.Path, rendered_path: pathlib.Path, fmt: str) -> None:
    dot = shutil.which("dot")
    if dot is None:
        raise RuntimeError("Graphviz dot executable was not found in PATH")
    ensure_dir(rendered_path.parent)
    subprocess.run([dot, f"-T{fmt}", str(dot_path), "-o", str(rendered_path)], check=True)


# Parse command-line arguments
def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize one HepMC3 event-record graph with Graphviz"
    )
    parser.add_argument("input", help="input HepMC3 ASCII file")
    parser.add_argument("-e", "--event", type=int, default=0, help="zero-based event index")
    parser.add_argument("-o", "--output", help="output .svg, .png, .pdf, or .dot file")
    parser.add_argument(
        "-f",
        "--format",
        choices=("svg", "png", "pdf", "dot"),
        help="rendering format when --output has no recognized suffix",
    )
    parser.add_argument(
        "--max-particles",
        type=int,
        default=250,
        help="maximum particles to draw; use 0 for all particles",
    )
    parser.add_argument(
        "--no-momentum",
        action="store_true",
        help="omit compact pT/eta/E labels",
    )
    parser.add_argument(
        "--show-attributes",
        action="store_true",
        help="include particle and vertex attributes in node labels",
    )
    parser.add_argument(
        "--show-vertex-balance",
        action="store_true",
        help="show momentum-balance diagnostics for every vertex",
    )
    parser.add_argument(
        "--no-hard-highlight",
        action="store_true",
        help="disable purple highlighting for Pythia hard-subprocess particles",
    )
    return parser.parse_args(argv)


# Run the command-line visualizer
def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    input_path = pathlib.Path(args.input)
    if args.event < 0:
        raise ValueError("--event must be non-negative")
    event = read_event(input_path, args.event)
    attr_table = event_attributes(event)
    dot = build_dot(
        event=event,
        input_path=input_path,
        event_index=args.event,
        attr_table=attr_table,
        max_particles=args.max_particles,
        show_momentum=not args.no_momentum,
        show_attributes=args.show_attributes,
        show_vertex_balance=args.show_vertex_balance,
        highlight_hard=not args.no_hard_highlight,
    )

    out = output_path(args)
    suffix = out.suffix.lower().lstrip(".")
    fmt = args.format or (suffix if suffix in {"svg", "png", "pdf", "dot"} else "svg")
    if fmt == "dot":
        dot_path = out
        write_text(dot_path, dot)
        print(f"DOT output: {dot_path}")
        return 0

    dot_path = out.with_suffix(".dot")
    write_text(dot_path, dot)
    render_dot(dot_path, out, fmt)
    print(f"Graphviz output: {out}")
    print(f"DOT output: {dot_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"iceviz.py: error: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
