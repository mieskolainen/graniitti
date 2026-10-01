#!/usr/bin/env python3
#
# Convert and build source-only MadLoop process runtimes
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import posixpath
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath
from typing import Any

from core.io.files import ensure_dir
from core.io.serialize import load_json_file

MG5_VERSION = "2.9.27"
SOURCE_DATA = "source.json"
SOURCE_LICENSE = "MADGRAPH_LICENSE"
SOURCE_ARCHIVES = {
    "generated": "madloop.tar.gz",
    "cuttools": "cuttools.tar.gz",
    "ninja": "ninja.tar.gz",
    "oneloop": "oneloop.tar.gz",
}
CARD_VALUES = {
    # Preserve complex helicity phases with quadruple precision loop reduction
    "CTModeRun": "4",
    "CheckCycle": "5",
    "DoubleCheckHelicityFilter": ".FALSE.",
    "HelicityFilterLevel": "0",
    "HelInitStartOver": ".TRUE.",
    "ImprovePSPoint": "-1",
    "LoopInitStartOver": ".TRUE.",
    "MLReductionLib": "6|1",
    "WriteOutFilters": ".FALSE.",
}
BINARY_SUFFIXES = {".a", ".mod", ".o", ".so"}


# Describe one generated loop runtime at the MG2GRA boundary
@dataclass(frozen=True)
class MadLoopMetadata:
    name: str
    process: str
    model: str
    mg5_version: str
    loop_diagrams: int
    counterterms: int
    helicities: int
    denominator: int
    helicity_average: int
    nexternal: int
    ncolor: int
    color_norm: float
    final_pdgs: list[int]
    final_color_representations: list[int]
    qcd_order: int
    qed_order: int
    process_directory: str
    prefix: str
    helicity_table: list[list[int]]
    default_alpha_s: float
    default_alpha_qed: float
    parameters: list[dict[str, Any]]
    particles: list[dict[str, Any]]
    charge: str | None
    alpha_zero: bool


# Compute the SHA256 digest of one file
def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            value.update(block)
    return value.hexdigest()


# Extract one required integer parameter from fixed-form Fortran
def integer_parameter(source: str, name: str) -> int:
    match = re.search(rf"PARAMETER\s*\([^)]*\b{name}\s*=\s*(\d+)", source)
    if match is None:
        raise RuntimeError(f"MadLoop output has no {name} parameter")
    return int(match.group(1))


# Extract one required initialized integer from fixed-form Fortran
def data_integer(source: str, name: str) -> int:
    match = re.search(rf"DATA\s+{name}\s*/\s*(\d+)\s*/", source)
    if match is None:
        raise RuntimeError(f"MadLoop output has no {name} value")
    return int(match.group(1))


# Compute the sole virtual loop subprocess directory
def process_directory(output: Path) -> Path:
    candidates = sorted((output / "SubProcesses").glob("PV*"))
    candidates = [path for path in candidates if (path / "loop_matrix.f").is_file()]
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one MadLoop subprocess in {output}")
    return candidates[0]


# Parse the generated physical helicity configurations
def helicity_table(output: Path, prefix: str, nexternal: int) -> list[list[int]]:
    resource = output / "SubProcesses" / "MadLoop5_resources"
    path = resource / f"{prefix}HelConfigs.dat"
    if not path.is_file():
        raise RuntimeError("MadLoop output has no helicity configurations")
    rows = [[int(value) for value in line.split()] for line in path.read_text().splitlines()]
    if not rows or any(len(row) != nexternal for row in rows):
        raise RuntimeError("MadLoop helicity configurations have inconsistent dimensions")
    return rows


# Validate Durham ordering of the two incoming helicities
def validate_durham_helicities(rows: list[list[int]]) -> None:
    if len(rows) % 4 != 0:
        raise RuntimeError("Durham MadLoop helicities do not factor into four incoming states")
    final_count = len(rows) // 4
    expected = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    finals = [tuple(row[2:]) for row in rows[:final_count]]
    for block, incoming in enumerate(expected):
        selected = rows[block * final_count : (block + 1) * final_count]
        if any(tuple(row[:2]) != incoming for row in selected):
            raise RuntimeError("MadLoop incoming helicity ordering is not lexicographic")
        if [tuple(row[2:]) for row in selected] != finals:
            raise RuntimeError("MadLoop final helicity ordering changes between incoming states")


# Read the scalar loop color-flow norm supported by the Durham adapter
def loop_color_norm(output: Path, prefix: str) -> float:
    path = output / "SubProcesses" / "MadLoop5_resources" / f"{prefix}LoopColorFlowMatrix.dat"
    values = [line.strip() for line in path.read_text().splitlines() if line.strip() != "EOF"]
    if len(values) != 2:
        raise RuntimeError("MG2GRA currently requires one MadLoop color flow")
    numerator = int(values[0])
    denominator = int(values[1])
    if denominator == 0:
        raise RuntimeError("MadLoop color-flow norm has a zero denominator")
    value = numerator / abs(denominator)
    if denominator < 0:
        raise RuntimeError("MG2GRA does not support an imaginary scalar color norm")
    if not math.isfinite(value) or value <= 0.0:
        raise RuntimeError("MadLoop color-flow norm is not positive")
    return value


# Evaluate the selected loop UFO and its actual amplitude orders
def loop_model(mg5_root: Path, output: Path, entry: dict, model: dict) -> dict:
    from . import mg5_color
    from .mg5_family import chdir, load_madgraph_family_api

    with chdir(output):
        command_class, _ = load_madgraph_family_api(mg5_root)
        command = command_class()
        command.no_notification()
        command.exec_cmd("set complex_mass_scheme " +
                         ("True --allow_qed" if model.get("complex_mass_scheme", False) else "False"),
                         printcmd=False)
        command.exec_cmd(f"import model {model['import']}", printcmd=False)
        command.exec_cmd(f"generate {entry['process']}", printcmd=False)
        ufo = command._curr_model
        orders = {tuple((name, diagram.get("orders").get(name, 0)) for name in ("QCD", "QED"))
                  for amplitude in command._curr_amps
                  for diagram in amplitude.get("loop_diagrams")}
        if len(orders) != 1:
            raise RuntimeError("MadLoop coupling rescaling requires one actual amplitude order")
        charge = model.get("alpha_qed", {}).get("charge")
        alpha_zero = charge is not None
        alpha_s, alpha_qed = mg5_color.model_couplings(
            ufo, (output / "Cards/param_card.dat").read_text(), charge)
        parameters = [{"name": p.name, "block": p.lhablock.lower(), "indices": list(p.lhacode)}
                      for p in ufo.get("parameters")[("external",)]]
        particles = [{"pdg": abs(p.get("pdg_code")), "mass": p.get("mass"), "width": p.get("width"),
                      "signed_mass": p.is_fermion() and p.get("self_antipart")}
                     for p in ufo.get("particles")]
        if charge is None:
            inverse = next((p["name"] for p in parameters
                            if p["block"] == "sminputs" and p["indices"] == [1]), None)
            charge = f"sqrt(4d0 * acos(-1d0) / {inverse})" if inverse else None
        return {"orders": dict(orders.pop()), "parameters": parameters, "particles": particles,
                "charge": charge, "alpha_zero": alpha_zero, "alpha_s": alpha_s, "alpha_qed": alpha_qed}


# Read the sole concrete Les Houches external state
def external_state(output: Path, nexternal: int) -> tuple[list[int], list[int]]:
    candidates = sorted((output / "SubProcesses").glob("P*/leshouche.inc"))
    candidates = [path for path in candidates if not path.parent.name.startswith("PV")]
    if len(candidates) != 1:
        raise RuntimeError("MG2GRA requires one concrete MadLoop external state")
    source = candidates[0].read_text()
    pdg_match = re.search(r"DATA\s*\(IDUP\([^/]+/([^/]+)/", source)
    color_matches = re.findall(r"DATA\s*\(ICOLUP\([^/]+/([^/]+)/", source)
    if pdg_match is None or len(color_matches) != 2:
        raise RuntimeError("MadLoop output has no complete Les Houches external state")
    pdgs = [int(value) for value in pdg_match.group(1).split(",")]
    colors = [[int(value) for value in row.split(",")] for row in color_matches]
    if len(pdgs) != nexternal or any(len(row) != nexternal for row in colors):
        raise RuntimeError("MadLoop Les Houches external state has inconsistent dimensions")
    if pdgs[:2] != [21, 21]:
        raise RuntimeError("Durham MadLoop conversion requires two incoming gluons")
    if any(colors[row][leg] != 0 for row in range(2) for leg in range(2, nexternal)):
        raise RuntimeError("MG2GRA currently supports colorless MadLoop final states")
    return pdgs[2:], [1] * (nexternal - 2)


# Inspect and validate one generated loop-induced Durham process
def inspect_output(output: Path, entry: dict[str, Any], model: dict) -> MadLoopMetadata:
    directory = process_directory(output)
    source = (directory / "loop_matrix.f").read_text()
    version_match = re.search(r"MadGraph5_aMC@NLO v\.\s*([0-9.]+)", source)
    if version_match is None or version_match.group(1) != MG5_VERSION:
        raise RuntimeError(f"MadLoop output must use MG5 {MG5_VERSION}")
    process_match = re.search(r"Process:\s*(.*?)\s*$", source, re.M)
    if process_match is None or "noborn" not in process_match.group(1).lower():
        raise RuntimeError("MadLoop output is not loop induced")
    prefix = (directory / "proc_prefix.txt").read_text().strip()
    if re.fullmatch(r"ML5_[0-9_]+_", prefix) is None:
        raise RuntimeError("MadLoop output has an unsafe subprocess prefix")
    nexternal = integer_parameter(source, "NEXTERNAL")
    rows = helicity_table(output, prefix, nexternal)
    validate_durham_helicities(rows)
    color_source = (directory / "compute_color_flows.f").read_text()
    if integer_parameter(color_source, "NLOOPAMPSO") != 1:
        raise RuntimeError("MG2GRA requires one selected MadLoop squared order")
    ncolor = integer_parameter(color_source, "NLOOPFLOWS")
    if ncolor != 1:
        raise RuntimeError("MG2GRA currently requires one MadLoop color flow")
    coefficients = (
        output / "SubProcesses" / "MadLoop5_resources" / f"{prefix}LoopColorFlowCoefs.dat"
    ).read_text()
    if "Coefficient for flow number 1 with expr. 1 Tr(1,2)" not in coefficients:
        raise RuntimeError("MadLoop color flow is not the incoming gluon trace")
    final_pdgs, final_representations = external_state(output, nexternal)
    return MadLoopMetadata(
        name=entry["name"],
        process=entry["process"],
        model=entry["model"],
        mg5_version=MG5_VERSION,
        loop_diagrams=integer_parameter(source, "NLOOPS"),
        counterterms=integer_parameter(source, "NCTAMPS"),
        helicities=len(rows),
        denominator=data_integer(source, "IDEN"),
        helicity_average=data_integer(source, "HELAVGFACTOR"),
        nexternal=nexternal,
        ncolor=ncolor,
        color_norm=loop_color_norm(output, prefix),
        final_pdgs=final_pdgs,
        final_color_representations=final_representations,
        qcd_order=model["orders"].get("QCD", 0),
        qed_order=model["orders"].get("QED", 0),
        process_directory=directory.name,
        prefix=prefix,
        helicity_table=rows,
        default_alpha_s=model["alpha_s"],
        default_alpha_qed=model["alpha_qed"],
        parameters=model["parameters"],
        particles=model["particles"],
        charge=model["charge"],
        alpha_zero=model["alpha_zero"],
    )


# Materialize all immediate MadLoop card values for immutable runtime use
def runtime_card(source: str) -> str:
    lines = source.splitlines()
    seen: set[str] = set()
    for index, line in enumerate(lines[:-1]):
        if not line.startswith("#"):
            continue
        name = line[1:].strip()
        value_index = index + 1
        value = lines[value_index]
        if value.startswith("!") and not value.startswith("! "):
            value = value[1:]
        if name in CARD_VALUES:
            value = CARD_VALUES[name]
            seen.add(name)
        if not value or value.startswith(("!", "#")):
            raise RuntimeError(f"MadLoop card has no immediate value for {name}")
        lines[value_index] = value
    missing = set(CARD_VALUES) - seen
    if missing:
        raise RuntimeError(f"MadLoop card has no required settings {sorted(missing)}")
    return "\n".join(lines) + "\n"


# Compute a fixed-helicity C ABI adapter for one generated subprocess
def adapter_source(metadata: MadLoopMetadata) -> str:
    fixed_norm = metadata.color_norm * metadata.helicity_average / metadata.denominator
    return f"""! Expose complex fixed-helicity MadLoop amplitudes through a C ABI
!
! (c) 2026 Mikael Mieskolainen
! Licensed under the MIT License <http://opensource.org/licenses/MIT>.

module graniitti_madloop
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  implicit none
contains

  ! Expose generated dimensions before passing model or amplitude buffers
  subroutine gra_madloop_sizes(sizes) bind(C, name='gra_madloop_sizes')
    integer(c_int), intent(out) :: sizes(5)
    sizes = [{metadata.nexternal}, {metadata.helicities}, {metadata.ncolor}, &
             {len(metadata.parameters)}, {len(metadata.particles)}]
  end subroutine gra_madloop_sizes

  ! Copy a bounded C character array into one Fortran path
  subroutine copy_path(source, count, destination, status)
    character(kind=c_char), intent(in) :: source(*)
    integer(c_int), value, intent(in) :: count
    character(len=512), intent(out) :: destination
    integer(c_int), intent(out) :: status
    integer :: i

    destination = ' '
    status = 0_c_int
    if (count < 1_c_int .or. count > len(destination)) then
      status = 1_c_int
      return
    end if
    do i = 1, count
      if (source(i) == c_null_char) exit
      destination(i:i) = source(i)
    end do
  end subroutine copy_path

  ! Initialize the model card and immutable reduction resources
  subroutine gra_madloop_init(param_c, param_len, resource_c, resource_len, status) &
      bind(C, name='gra_madloop_init')
    character(kind=c_char), intent(in) :: param_c(*)
    integer(c_int), value, intent(in) :: param_len
    character(kind=c_char), intent(in) :: resource_c(*)
    integer(c_int), value, intent(in) :: resource_len
    integer(c_int), intent(out) :: status
    character(len=512) :: param_path
    character(len=512) :: resource_path
    integer(c_int) :: path_status
    logical :: forbid_doublecheck

    call copy_path(param_c, param_len, param_path, status)
    if (status /= 0_c_int) return
    call copy_path(resource_c, resource_len, resource_path, path_status)
    if (path_status /= 0_c_int) then
      status = path_status
      return
    end if
    call setparamlog(.false.)
    call setpara(param_path)
    call setmadlooppath(resource_path)
    forbid_doublecheck = .true.
    call set_forbid_hel_doublecheck(forbid_doublecheck)
  end subroutine gra_madloop_init

  ! Evaluate all physical helicities and preserve each complex color flow
  subroutine gra_madloop_eval(momentum, amp_re, amp_im, pole_re, pole_im, &
      amp2, ret_code, status) bind(C, name='gra_madloop_eval')
    real(c_double), intent(in) :: momentum(4, {metadata.nexternal})
    real(c_double), intent(out) :: amp_re({metadata.ncolor}, {metadata.helicities})
    real(c_double), intent(out) :: amp_im({metadata.ncolor}, {metadata.helicities})
    real(c_double), intent(out) :: pole_re(2, {metadata.ncolor}, {metadata.helicities})
    real(c_double), intent(out) :: pole_im(2, {metadata.ncolor}, {metadata.helicities})
    real(c_double), intent(out) :: amp2({metadata.helicities})
    integer(c_int), intent(out) :: ret_code({metadata.helicities})
    integer(c_int), intent(out) :: status
    complex(c_double_complex) :: jampl(3, {metadata.ncolor}, 1)
    real(c_double) :: answer(0:3, 0:1)
    real(c_double) :: accuracy(0:1)
    real(c_double) :: norm_error
    real(c_double) :: norm_value
    integer :: color
    integer :: h
    integer :: hel
    integer :: t
    integer :: u
    common /{metadata.prefix}JAMPL/ jampl
    common /{metadata.prefix}ACC/ accuracy, h, t, u

    status = 0_c_int
    do hel = 1, {metadata.helicities}
      call {metadata.prefix}SLOOPMATRIXHEL(momentum, hel, answer)
      norm_value = 0.0_c_double
      do color = 1, {metadata.ncolor}
        amp_re(color, hel) = real(jampl(1, color, 1), kind=c_double)
        amp_im(color, hel) = aimag(jampl(1, color, 1))
        pole_re(1, color, hel) = real(jampl(2, color, 1), kind=c_double)
        pole_im(1, color, hel) = aimag(jampl(2, color, 1))
        pole_re(2, color, hel) = real(jampl(3, color, 1), kind=c_double)
        pole_im(2, color, hel) = aimag(jampl(3, color, 1))
        norm_value = norm_value + amp_re(color, hel)**2 + amp_im(color, hel)**2
      end do
      amp2(hel) = answer(1, 0)
      ret_code(hel) = 100_c_int * h + 10_c_int * t + u
      norm_error = abs(amp2(hel) - {fixed_norm:.17e}_c_double * norm_value)
      if (.not. ieee_is_finite(amp2(hel))) status = 2_c_int
      if (norm_error > 1.0e-8_c_double * max(1.0_c_double, abs(amp2(hel)))) &
          status = 3_c_int
    end do
  end subroutine gra_madloop_eval

end module graniitti_madloop
"""


# Expose evaluated UFO parameters without assumptions about their SLHA blocks
def model_source(metadata: MadLoopMetadata) -> str:
    assignments = "\n".join(f"        {p['name']} = values({i})"
                            for i, p in enumerate(metadata.parameters, 1))
    assignments += "\n" + "\n".join(
        f"        {p['width']}=sign(abs({p['width']}),{p['mass']})"
        for p in metadata.particles if p["signed_mass"] and p["width"].upper() != "ZERO")
    assignments += "\n" + "\n".join(f"        MP__{p['name']} = {p['name']}" for p in metadata.parameters)
    values = "\n".join(f"      values({i}) = {p['name']}"
                       for i, p in enumerate(metadata.parameters, 1))
    poles = "\n".join(f"      poles({2 * i + j + 1}) = {p[key]}"
                      for i, p in enumerate(metadata.particles)
                      for j, key in enumerate(("mass", "width")))
    strong = "AS" if any(p["name"].lower() == "as" for p in metadata.parameters) else "0d0"
    qed = f"abs({metadata.charge})**2 / (4d0 * acos(-1d0))" if metadata.charge else "0d0"
    return f"""C     Expose model parameters and physical poles through the C ABI
      subroutine gra_madloop_model(values,poles,alpha,update)
     $ bind(C,name='gra_madloop_model')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double) values({len(metadata.parameters)})
      real(c_double) poles({2 * len(metadata.particles)}),alpha(2)
      integer(c_int), value :: update
      double precision ZERO
      parameter (ZERO=0d0)
      include 'input.inc'
      include 'coupl.inc'
      include 'mp_input.inc'
      include 'mp_coupl.inc'
      if (update.ne.0) then
{assignments}
        call COUP()
        call {metadata.prefix}CLEAR_CACHES()
      endif
{values}
{poles}
      alpha(1) = {strong}
      alpha(2) = {qed}
      end
"""


# Compute link stubs for the unused IREGI implementation
def iregi_stub_source() -> str:
    return """C     IREGI is not selected by the portable MadLoop runtime
      SUBROUTINE INITIREGI()
      RETURN
      END

      SUBROUTINE IMLOOP()
      RETURN
      END

      SUBROUTINE IREGI_FREE_PS()
      RETURN
      END
"""


# Extract portable regular files while rejecting archive path traversal
def extract(archive: Path, output: Path) -> None:
    ensure_dir(output)
    links: list[tarfile.TarInfo] = []
    with tarfile.open(archive, "r:gz") as source:
        for member in source.getmembers():
            path = PurePosixPath(member.name)
            if path.is_absolute() or ".." in path.parts:
                raise RuntimeError(f"Unsafe source archive member {member.name}")
            if any(part.startswith("._") for part in path.parts):
                continue
            if member.issym() or member.islnk():
                links.append(member)
                continue
            source.extract(member, output, filter="data")
    for member in links:
        path = PurePosixPath(member.name)
        raw_target = (
            path.parent / member.linkname if member.issym() else PurePosixPath(member.linkname)
        )
        target = PurePosixPath(posixpath.normpath(str(raw_target)))
        if target.is_absolute() or ".." in target.parts:
            continue
        destination = output / path
        source_path = output / target
        if not source_path.exists():
            raise RuntimeError(f"Source archive link has no target {member.name}")
        ensure_dir(destination.parent)
        if member.issym():
            destination.symlink_to(member.linkname, target_is_directory=source_path.is_dir())
        else:
            os.link(source_path, destination)


# Execute one checked build command in its source directory
def run(command: list[str], directory: Path, environment: dict[str, str] | None = None) -> None:
    subprocess.run(command, cwd=directory, env=environment, check=True)


# Replace one checked source fragment in an extracted dependency
def replace_once(path: Path, old: str, new: str) -> None:
    source = path.read_text()
    if source.count(old) != 1:
        raise RuntimeError(f"Unexpected MadLoop dependency source {path}")
    path.write_text(source.replace(old, new))


# Silence third-party banners emitted once for every isolated runtime
def silence_runtime(generated: Path, dependencies: Path) -> None:
    loop = process_directory(generated) / "loop_matrix.f"
    replace_once(
        loop, "        CALL PRINT_MADLOOP_BANNER()", "C       MadLoop banner disabled by GRANIITTI"
    )
    replace_once(
        loop,
        "        CALL MADLOOPPARAMREADER(PARAMFN,.TRUE.)",
        "        CALL MADLOOPPARAMREADER(PARAMFN,.FALSE.)",
    )
    replace_once(
        dependencies / "OneLOop-3.6/src/avh_olo_version.f90",
        "logical ,save :: done=.false.",
        "logical ,save :: done=.true.",
    )
    replace_once(
        dependencies / "CutTools/src/avh/avh_olo.f90",
        "logical ,save :: done=.false.",
        "logical ,save :: done=.true.",
    )
    cuttools = dependencies / "CutTools/src/cts/cts_cuttools.f90"
    source = cuttools.read_text()
    source, count = re.subn(
        r"   write \(\*,\*\) ' '\n.*?   write \(\*,\*\) '  '\n", "", source, count=1, flags=re.S
    )
    if count != 1:
        raise RuntimeError(f"Unexpected MadLoop dependency source {cuttools}")
    cuttools.write_text(source)
    replace_once(
        dependencies / "ninja-1.2.0/src/ninja.cc",
        "bool Options::quiet = false;",
        "bool Options::quiet = true;",
    )


# Build the scalar OneLOop library used by Ninja
def build_oneloop(source: Path, prefix: Path, compiler: str) -> None:
    config = source / "Config"
    text = config.read_text()
    text = re.sub(r"^FC\s*=.*$", f"FC = {compiler}", text, flags=re.M)
    text = re.sub(r"^FFLAGS\s*=.*$", "FFLAGS = -O2 -fPIC", text, flags=re.M)
    text = re.sub(r"^#?QPKIND\s*=.*$", "QPKIND = 16", text, flags=re.M)
    config.write_text(text)
    run([sys.executable, "create.py"], source)
    ensure_dir(prefix / "lib")
    ensure_dir(prefix / "include")
    shutil.copy2(source / "libavh_olo.a", prefix / "lib" / "libavh_olo.a")
    for module in source.glob("*.mod"):
        shutil.copy2(module, prefix / "include" / module.name)


# Build CutTools with its quadruple precision fallback
def build_cuttools(source: Path, prefix: Path, compiler: str) -> None:
    arguments = [f"FC={compiler}", "FFLAGS=-O2 -fPIC -fno-automatic"]
    run(["make", "-j1", *arguments, "cpqp"], source)
    run(["make", "-j1", *arguments, "default"], source)
    include = source / "includects"
    shutil.copy2(include / "libcts.a", prefix / "lib" / "libcts.a")
    shutil.copy2(include / "mpmodule.mod", prefix / "include" / "mpmodule.mod")


# Build Ninja against the local OneLOop copy
def build_ninja(source: Path, prefix: Path, compiler: str, cxx: str) -> None:
    environment = os.environ.copy()
    environment.update({"FC": compiler, "CXX": cxx})
    run(
        [
            "./configure",
            f"--prefix={prefix}",
            "--enable-higher_rank",
            f"--with-avholo=-L{prefix / 'lib'} -lavh_olo",
            f"FCINCLUDE=-I{prefix / 'include'}",
            "CXXFLAGS=-O2 -fPIC -fcx-fortran-rules -fno-exceptions -fno-rtti",
            "CPPFLAGS=-DNINJA_NO_EXCEPTIONS -fPIC",
            "--enable-quadninja",
            "--disable-shared",
            "--enable-static",
            "LIBS=-lstdc++",
        ],
        source,
        environment,
    )
    run(["make", "-j1"], source, environment)
    run(["make", "install"], source, environment)


# Build generated model support and the complete MadLoop shared library
def build_generated(generated: Path, prefix: Path, compiler: str, cxx: str) -> None:
    ensure_dir(generated / "lib")
    definitions = generated / "SubProcesses" / "MadLoop_makefile_definitions"
    definitions.write_text(
        f"LINK_LOOP_LIBS = -L{prefix / 'lib'} -lninja -lavh_olo -lcts\n"
        f"LOOP_LIBS = {prefix / 'lib/libninja.a'} {prefix / 'lib/libavh_olo.a'} "
        f"{prefix / 'lib/libcts.a'}\n"
        f"LOOP_INCLUDE = -I {prefix / 'include'}\n"
        "LOOP_PREFIX = PV\nDOTO = %.o\nDOTF = %.f\n"
        "LINK_MADLOOP_LIB = -L$(LIBDIR) -lMadLoop\n"
        "MADLOOP_LIB = $(LIBDIR)libMadLoop.$(libext)\n\n"
        "$(MADLOOP_LIB):\n\tcd ..; make -f makefile_MadLoop OLP_static\n"
    )
    subprocess_dir = generated / "SubProcesses"
    makefile = subprocess_dir / "makefile_MadLoop"
    original = makefile.read_text()
    marker = "OLP_PROCESS= MadLoopParamReader.o MadLoopCommons.o " + "\\\n"
    if original.count(marker) != 1:
        raise RuntimeError("Unexpected generated MadLoop OLP object list")
    makefile.write_text(
        original.replace(
            marker, "OLP_PROCESS= MadLoopParamReader.o MadLoopCommons.o iregi_stub.o " + "\\\n"
        )
    )
    (subprocess_dir / "iregi_stub.f").write_text(iregi_stub_source())
    variables = [f"FC={compiler}", f"CXX={cxx}"]
    run(["make", "-j4", *variables, "libmodel", "libdhelas"], generated / "Source")
    run(["make", "-f", "makefile_MadLoop", "-j4", *variables, "OLP"], subprocess_dir)


# Load and validate one source bundle manifest
def load_source_data(source: Path) -> dict[str, Any]:
    data = load_json_file(source / SOURCE_DATA)
    required = {"version", "metadata", "archives", "license_sha256"}
    if not isinstance(data, dict) or set(data) != required:
        raise RuntimeError("Invalid MadLoop source data")
    if data["version"] != 1:
        raise RuntimeError("Unsupported MadLoop source data version")
    metadata = data["metadata"]
    archives = data["archives"]
    if not isinstance(metadata, dict) or not isinstance(archives, dict):
        raise RuntimeError("Invalid MadLoop source data")
    if set(archives) != set(SOURCE_ARCHIVES.values()):
        raise RuntimeError("Incomplete MadLoop source archives")
    for filename, expected in archives.items():
        if digest(source / filename) != expected:
            raise RuntimeError(f"MadLoop source checksum failed for {filename}")
    if digest(source / SOURCE_LICENSE) != data["license_sha256"]:
        raise RuntimeError("MadLoop source license checksum failed")
    return data


# Build the MadLoop shared library and its C ABI adapter
def build_output(output: Path, metadata: MadLoopMetadata, compiler: str) -> None:
    subprocess_dir = output / "SubProcesses"
    library = subprocess_dir / "libMadLoop.so"
    if not library.is_file():
        run(["make", "-f", "makefile_MadLoop", "-j4", "OLP"], subprocess_dir)
    adapter = process_directory(output) / "graniitti_madloop.f90"
    adapter.write_text(adapter_source(metadata))
    model = adapter.with_name("graniitti_model.f")
    model.write_text(model_source(metadata))
    run(
        [
            compiler,
            "-O2",
            "-fPIC",
            "-shared",
            str(adapter),
            str(model),
            "-ffixed-line-length-none",
            "-I" + str(output / "Source/MODEL"),
            "-L.",
            "-Wl,--no-as-needed",
            "-lMadLoop",
            "-lstdc++",
            "-Wl,-rpath,$ORIGIN",
            "-o",
            str(subprocess_dir / "libMadLoopAdapter.so"),
        ],
        subprocess_dir,
    )


# Preserve one previous runtime before replacing it
def preserve_output(output: Path) -> None:
    if not output.exists():
        return
    index = 0
    while True:
        suffix = "._old" if index == 0 else f"._old.{index}"
        previous = output.with_name(output.name + suffix)
        if not previous.exists():
            output.rename(previous)
            return
        index += 1


# Install runtime libraries, cards and immutable reduction resources
def install_output(generated: Path, output: Path, metadata: MadLoopMetadata, compiler: str) -> None:
    build_output(generated, metadata, compiler)
    subprocess_dir = generated / "SubProcesses"
    resource_source = subprocess_dir / "MadLoop5_resources"
    preserve_output(output)
    resource_output = output / "MadLoop5_resources"
    ensure_dir(resource_output, exist_ok=False)
    for source in sorted(resource_source.iterdir()):
        if source.is_file():
            shutil.copy2(source.resolve(), resource_output / source.name)
    card_source = generated / "Cards" / "param_card.dat"
    shutil.copy2(card_source, output / "param_card.dat")
    shutil.copy2(card_source, resource_output / "param_card.dat")
    ident_card = generated / "Cards" / "ident_card.dat"
    if ident_card.is_file():
        shutil.copy2(ident_card, resource_output / "ident_card.dat")
    shutil.copy2(subprocess_dir / "libMadLoop.so", output / "libMadLoop.so")
    shutil.copy2(subprocess_dir / "libMadLoopAdapter.so", output / "libMadLoopAdapter.so")
    resource_card = resource_output / "MadLoopParams.dat"
    if not resource_card.is_file():
        shutil.copy2(generated / "Cards" / "MadLoopParams.dat", resource_card)
    resource_card.write_text(runtime_card(resource_card.read_text()))
    runtime = asdict(metadata)
    runtime["lib_madloop_sha256"] = digest(output / "libMadLoop.so")
    runtime["adapter_sha256"] = digest(output / "libMadLoopAdapter.so")
    (output / "runtime.json").write_text(json.dumps(runtime, indent=2) + "\n")


# Build and install a checked source bundle without invoking MadGraph
def build_bundle(source: Path, work: Path, output: Path, compiler: str, cxx: str) -> None:
    source = source.resolve()
    work = work.resolve()
    output = output.resolve()
    data = load_source_data(source)
    metadata = MadLoopMetadata(**data["metadata"])
    preserve_output(work)
    generated = work / "generated"
    dependencies = work / "dependencies"
    prefix = work / "prefix"
    extract(source / SOURCE_ARCHIVES["generated"], generated)
    extract(source / SOURCE_ARCHIVES["cuttools"], dependencies)
    extract(source / SOURCE_ARCHIVES["ninja"], dependencies)
    extract(source / SOURCE_ARCHIVES["oneloop"], dependencies)
    silence_runtime(generated, dependencies)
    build_oneloop(dependencies / "OneLOop-3.6", prefix, compiler)
    build_cuttools(dependencies / "CutTools", prefix, compiler)
    build_ninja(dependencies / "ninja-1.2.0", prefix, compiler, cxx)
    build_generated(generated, prefix, compiler, cxx)
    install_output(generated, output, metadata, compiler)


# Copy one source tree while excluding compiled files
def copy_source_tree(source: Path, destination: Path) -> None:
    shutil.copytree(
        source,
        destination,
        symlinks=False,
        ignore=shutil.ignore_patterns("._*", "*.a", "*.mod", "*.o", "*.so", "compiler_version.log"),
    )


# Write one deterministic source-only gzip tar archive
def write_archive(source: Path, archive: Path) -> None:
    ensure_dir(archive.parent)
    with (
        archive.open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed,
        tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as target,
    ):
        for path in sorted(source.rglob("*"), key=lambda item: item.relative_to(source).as_posix()):
            relative = path.relative_to(source).as_posix()
            info = target.gettarinfo(str(path), relative)
            info.uid = 0
            info.gid = 0
            info.uname = ""
            info.gname = ""
            info.mtime = 0
            if info.isfile():
                if path.suffix in BINARY_SUFFIXES:
                    raise RuntimeError(f"Compiled file in MadLoop source bundle: {relative}")
                with path.open("rb") as stream:
                    target.addfile(info, stream)
            else:
                target.addfile(info)


# Stage the minimal generated sources needed by the runtime build
def stage_generated_source(output: Path, staging: Path) -> None:
    cards = staging / "Cards"
    ensure_dir(cards, exist_ok=False)
    for name in ("MadLoopParams.dat", "ident_card.dat", "param_card.dat", "proc_card_mg5.dat"):
        source = output / "Cards" / name
        if source.is_file():
            shutil.copy2(source, cards / name)
    source_out = staging / "Source"
    ensure_dir(source_out, exist_ok=False)
    for name in ("make_opts", "makefile", "param_card.inc"):
        source = output / "Source" / name
        if source.is_file():
            shutil.copy2(source, source_out / name)
    for name in ("DHELAS", "MODEL"):
        copy_source_tree(output / "Source" / name, source_out / name)
    subprocess_out = staging / "SubProcesses"
    ensure_dir(subprocess_out, exist_ok=False)
    for name in (
        "MGVersion.txt",
        "MadLoopCommons.f",
        "MadLoopParamReader.f",
        "MadLoopParams.dat",
        "MadLoopParams.inc",
        "MadLoop_makefile_definitions",
        "coupl.inc",
        "cts_mpc.h",
        "cts_mprec.h",
        "global_specs.inc",
        "makefile_MadLoop",
        "mp_coupl.inc",
        "mp_coupl_same_name.inc",
    ):
        source = output / "SubProcesses" / name
        if source.is_file():
            shutil.copy2(source.resolve(), subprocess_out / name)
    copy_source_tree(process_directory(output), subprocess_out / process_directory(output).name)
    copy_source_tree(
        output / "SubProcesses" / "MadLoop5_resources",
        subprocess_out / "MadLoop5_resources",
    )
    resource = subprocess_out / "MadLoop5_resources"
    for name in ("MadLoopParams.dat", "ident_card.dat", "param_card.dat"):
        source = cards / name
        if source.is_file() and not (resource / name).is_file():
            shutil.copy2(source, resource / name)


# Create a checked source-only runtime bundle from one MG5 export
def create_bundle(
    output: Path, destination: Path, mg5_root: Path, metadata: MadLoopMetadata
) -> dict[str, Path]:
    output = output.resolve()
    environment = os.environ.copy()
    environment["PYTHONWARNINGS"] = "ignore"
    run([str(output / "bin" / "madevent"), "treatcards", "param"], output, environment)
    ensure_dir(destination)
    with tempfile.TemporaryDirectory(prefix="madloop_bundle_", dir=destination.parent) as temporary:
        staging = Path(temporary)
        generated = staging / "generated"
        stage_generated_source(output, generated)
        generated_archive = destination / SOURCE_ARCHIVES["generated"]
        write_archive(generated, generated_archive)
        cuttools = staging / "cuttools" / "CutTools"
        copy_source_tree(mg5_root / "vendor" / "CutTools", cuttools)
        cuttools_archive = destination / SOURCE_ARCHIVES["cuttools"]
        write_archive(staging / "cuttools", cuttools_archive)
    shutil.copy2(mg5_root / "vendor" / "ninja.tar.gz", destination / SOURCE_ARCHIVES["ninja"])
    shutil.copy2(mg5_root / "vendor" / "oneloop.tar.gz", destination / SOURCE_ARCHIVES["oneloop"])
    license_path = destination / SOURCE_LICENSE
    shutil.copy2(mg5_root / "LICENSE", license_path)
    archives = {filename: digest(destination / filename) for filename in SOURCE_ARCHIVES.values()}
    source_data = {
        "version": 1,
        "metadata": asdict(metadata),
        "archives": archives,
        "license_sha256": digest(license_path),
    }
    (destination / SOURCE_DATA).write_text(json.dumps(source_data, indent=2) + "\n")
    return {name: destination / filename for name, filename in SOURCE_ARCHIVES.items()}


# Compute exact finite-Nc Durham data for one colorless scalar loop flow
def durham_color_data(entry: dict[str, Any], metadata: MadLoopMetadata) -> dict[str, Any]:
    if metadata.ncolor != 1 or metadata.color_norm <= 0.0:
        raise RuntimeError("Durham loop conversion requires one positive color flow")
    basis = f"{metadata.prefix}Tr(1,2)\n{metadata.color_norm:.17g}"
    return {
        "version": 1,
        "process": entry["process"],
        "incoming_pdgs": [21, 21],
        "final_pdgs": metadata.final_pdgs,
        "final_color_representations": metadata.final_color_representations,
        "ncolor": 1,
        "rank": 1,
        "basis_sha256": hashlib.sha256(basis.encode()).hexdigest(),
        "factorization_residual": 0.0,
        "projectors": [[[math.sqrt(metadata.color_norm), 0.0]]],
        "flow_candidates": [],
        "flow_weights": [],
    }


# Format one C++ floating-point literal
def cpp_float(value: float) -> str:
    return format(float(value), ".17g")


# Generate one concrete Durham MadLoop process header
def durham_header(entry: dict[str, Any], metadata: MadLoopMetadata) -> str:
    name = entry["name"]
    guard = f"AMP_MG5_{name.upper()}_H"
    final_pdgs = ", ".join(str(value) for value in metadata.final_pdgs)
    final_representations = ", ".join(str(value) for value in metadata.final_color_representations)
    return f"""// Generated Durham MadLoop process
// MadGraph5_aMC@NLO v. {metadata.mg5_version}, 2026-01-05
// @@@@ MadGraph to GRANIITTI conversion done @@@@
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef {guard}
#define {guard}

#include <atomic>
#include <complex>
#include <cstddef>
#include <mutex>
#include <string>
#include <vector>

#include "Graniitti/Amplitude/MG5/Durham/AMP_MG5_DurhamRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_MadLoop.h"

namespace gra {{

// Evaluate one loop-induced Durham hard subprocess
class AMP_MG5_{name} final : public DurhamMG5Process {{
 public:
  // Load the standard generated MadLoop runtime
  AMP_MG5_{name}();

  // Load one explicit generated MadLoop runtime
  explicit AMP_MG5_{name}(const std::string &directory);

  // Initialize the evaluated UFO parameters and invalidate the hard cache
  void InitParameters(SLHAReader card) override;

  // Compute the physical pole parameters used by phase space
  mg5::ParticleMap Particles() const override {{ return madloop_.Particles(); }}

  // Compute the stable MG2GRA process name
  const std::string &Name() const override;

  // Compute final-state PDG codes in generated order
  const std::vector<int> &FinalPDGs() const override;

  // Compute final-state SU(3) representations
  const std::vector<int> &FinalColorRepresentations() const override;

  // Compute the number of loop color-flow tensors
  std::size_t ColorCount() const override;

  // Compute the incoming-singlet color rank
  std::size_t ColorRank() const override;

  // Compute the complete generated helicity count
  std::size_t HelicityCount() const override;

  // Compute the one loop subprocess represented by this wrapper
  std::size_t SubprocessCount() const override;

  // Compute exact incoming-singlet projectors
  const std::vector<std::complex<double>> &ExactProjectors() const override;

  // Compute the empty colorless shower-flow candidates
  const std::vector<std::vector<MColorFlow>> &FlowCandidates() const override;

  // Evaluate exact projected fixed-helicity amplitudes
  DurhamMG5Evaluation Evaluate(LORENTZSCALAR &lts, double alpha_s,
                               M4Vec *hard_k1, M4Vec *hard_k2) override;

  // Compute whether the generated loop runtime is ready
  bool Ready() const noexcept;

  // Compute the generated loop runtime diagnostic
  std::string Error() const;

  // Compute the number of expensive loop evaluations
  std::size_t HardCallCount() const noexcept;

  // Compute the number of exact hard-point cache reuses
  std::size_t HardReuseCount() const noexcept;

 private:
  mg5::MadLoop madloop_;
  std::string name_ = "{name}";
  std::vector<int> final_pdgs_ = {{{final_pdgs}}};
  std::vector<int> final_color_representations_ = {{{final_representations}}};
  std::vector<std::complex<double>> exact_projectors_ = {{
      std::complex<double>({cpp_float(math.sqrt(metadata.color_norm))}, 0.0)}};
  std::vector<std::vector<MColorFlow>> flow_candidates_;
  std::vector<M4Vec> final_buffer_;
  std::vector<double> hard_point_;
  mg5::MadLoopResult hard_result_;
  bool hard_ready_ = false;
  std::mutex evaluation_mutex_;
  std::atomic<std::size_t> hard_calls_ = 0;
  std::atomic<std::size_t> hard_reuses_ = 0;
}};

}}  // namespace gra

#endif
"""


# Generate one concrete Durham MadLoop process source
def durham_source(entry: dict[str, Any], metadata: MadLoopMetadata) -> str:
    name = entry["name"]
    helicities = ",\n    ".join(
        "{" + ", ".join(str(value) for value in row) + "}" for row in metadata.helicity_table
    )
    env_name = f"GRANIITTI_MADLOOP_{name.upper()}"
    return f"""// Generated Durham MadLoop process
// MadGraph5_aMC@NLO v. {metadata.mg5_version}, 2026-01-05
// @@@@ MadGraph to GRANIITTI conversion done @@@@
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// [REFERENCE: Hirschi et al., JHEP 05 (2011) 044, arXiv:1103.0621]
// [REFERENCE: Alwall et al., JHEP 07 (2014) 079, arXiv:1405.0301]

#include "Graniitti/Amplitude/MG5/Durham/AMP_MG5_{name}.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Kinematics.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {{
namespace {{

constexpr std::size_t kExternal = {metadata.nexternal};
constexpr std::size_t kHelicity = {metadata.helicities};
constexpr std::size_t kColor = {metadata.ncolor};
constexpr int kHelicities[kHelicity][kExternal] = {{
    {helicities}}};

// Resolve an explicit runtime followed by the configured build path
std::string RuntimeDirectory() {{
  const char *configured = std::getenv("{env_name}");
  if (configured != nullptr && configured[0] != '\\0') {{ return configured; }}
#if defined(GRANIITTI_MADLOOP_ROOT)
  return std::string(GRANIITTI_MADLOOP_ROOT) + "/{name}/runtime";
#else
  return aux::ResolveProjectPath("MG5cards/Durham/{name}");
#endif
}}

// Store one four-vector in MadGraph E, px, py, pz order
void StoreMomentum(const M4Vec &momentum, std::vector<double> &storage) {{
  storage.push_back(momentum.E());
  storage.push_back(momentum.Px());
  storage.push_back(momentum.Py());
  storage.push_back(momentum.Pz());
}}

// Compare canonical hard points without an approximate cache key
bool SameHardPoint(const std::vector<double> &first,
                   const std::vector<double> &second) {{
  return first.size() == second.size() &&
         std::equal(first.cbegin(), first.cend(), second.cbegin(),
                    [](double lhs, double rhs) {{
                      return std::bit_cast<std::uint64_t>(lhs) ==
                             std::bit_cast<std::uint64_t>(rhs);
                    }});
}}

// Compute the event-to-card coupling rescaling for this loop order
double CouplingScale(double alpha_s, double alpha_s_card,
                     double alpha_qed_card, const std::string &qed_scheme) {{
  const double alpha_qed =
      qed_scheme == "ZERO" ? qed::alpha_QED() : alpha_qed_card;
  return std::pow({"alpha_s / alpha_s_card" if metadata.qcd_order else "1.0"}, {cpp_float(metadata.qcd_order / 2.0)}) *
         std::pow({"alpha_qed / alpha_qed_card" if metadata.qed_order and metadata.alpha_zero else "1.0"},
                  {cpp_float(metadata.qed_order / 2.0)});
}}

}}  // namespace

// Load the standard generated MadLoop runtime
AMP_MG5_{name}::AMP_MG5_{name}() : AMP_MG5_{name}(RuntimeDirectory()) {{}}

// Load one explicit generated MadLoop runtime
AMP_MG5_{name}::AMP_MG5_{name}(const std::string &directory)
    : DurhamMG5Process(amplitude::Processes("DURHAM", "{name}")),
      madloop_(directory, "{name}", kExternal, kHelicity, kColor) {{}}

// Initialize the evaluated UFO parameters and invalidate the hard cache
void AMP_MG5_{name}::InitParameters(SLHAReader card) {{
  std::lock_guard<std::mutex> guard(evaluation_mutex_);
  hard_ready_ = false;
  madloop_.InitParameters(card);
}}

// Compute the stable MG2GRA process name
const std::string &AMP_MG5_{name}::Name() const {{ return name_; }}

// Compute final-state PDG codes in generated order
const std::vector<int> &AMP_MG5_{name}::FinalPDGs() const {{ return final_pdgs_; }}

// Compute final-state SU(3) representations
const std::vector<int> &AMP_MG5_{name}::FinalColorRepresentations() const {{
  return final_color_representations_;
}}

// Compute the number of loop color-flow tensors
std::size_t AMP_MG5_{name}::ColorCount() const {{ return kColor; }}

// Compute the incoming-singlet color rank
std::size_t AMP_MG5_{name}::ColorRank() const {{ return 1; }}

// Compute the complete generated helicity count
std::size_t AMP_MG5_{name}::HelicityCount() const {{ return kHelicity; }}

// Compute the one loop subprocess represented by this wrapper
std::size_t AMP_MG5_{name}::SubprocessCount() const {{ return 1; }}

// Compute exact incoming-singlet projectors
const std::vector<std::complex<double>> &AMP_MG5_{name}::ExactProjectors() const {{
  return exact_projectors_;
}}

// Compute the empty colorless shower-flow candidates
const std::vector<std::vector<MColorFlow>> &AMP_MG5_{name}::FlowCandidates() const {{
  return flow_candidates_;
}}

// Evaluate exact projected fixed-helicity amplitudes
DurhamMG5Evaluation AMP_MG5_{name}::Evaluate(LORENTZSCALAR &lts, double alpha_s,
                                             M4Vec *hard_k1, M4Vec *hard_k2) {{
  std::lock_guard<std::mutex> guard(evaluation_mutex_);
  DurhamMG5Evaluation evaluation;
  lts.hamp.clear();
  if (hard_k1 != nullptr) {{ *hard_k1 = M4Vec(); }}
  if (hard_k2 != nullptr) {{ *hard_k2 = M4Vec(); }}
  if (!Ready() ||
      !std::isfinite(alpha_s) || alpha_s < 0.0) {{
    evaluation.status = mg5helas::EvaluationStatus::AmplitudeFailure;
    return evaluation;
  }}

  const auto leaves = mg5::StableDecayLeaves(lts.decaytree);
  final_buffer_ = mg5::StableLeafMomenta(leaves);
  if (final_buffer_.size() + 2 != kExternal) {{
    evaluation.status = mg5helas::EvaluationStatus::KinematicsFailure;
    return evaluation;
  }}
  M4Vec p1;
  M4Vec p2;
  if (!mg5helas::PrepareOnShellKinematics(lts, final_buffer_, p1, p2)) {{
    evaluation.status = mg5helas::EvaluationStatus::KinematicsFailure;
    return evaluation;
  }}

  std::vector<double> momentum;
  momentum.reserve(4 * kExternal);
  StoreMomentum(p1, momentum);
  StoreMomentum(p2, momentum);
  for (const auto &p4 : final_buffer_) {{ StoreMomentum(p4, momentum); }}
  if (hard_ready_ && SameHardPoint(momentum, hard_point_)) {{
    ++hard_reuses_;
  }} else {{
    ++hard_calls_;
    hard_ready_ = false;
    if (!madloop_.Evaluate(momentum, hard_result_)) {{
      evaluation.status = mg5helas::EvaluationStatus::AmplitudeFailure;
      return evaluation;
    }}
    hard_point_ = momentum;
    hard_ready_ = true;
  }}
  if (lts.model_cache == nullptr) {{
    evaluation.status = mg5helas::EvaluationStatus::AmplitudeFailure;
    return evaluation;
  }}

  const double coupling =
      CouplingScale(alpha_s, madloop_.AlphaS(), madloop_.AlphaQED(),
                    lts.model_cache->Tune().Structure().QED_alpha);
  evaluation.projected.assign(kHelicity, 0.0);
  for (std::size_t hel = 0; hel < kHelicity; ++hel) {{
    const std::complex<double> value =
        coupling * exact_projectors_[0] *
        hard_result_.amplitude[hel * kColor];
    evaluation.projected[hel] = value;
  }}
  if (hard_k1 != nullptr) {{ *hard_k1 = p1; }}
  if (hard_k2 != nullptr) {{ *hard_k2 = p2; }}
  return evaluation;
}}

// Compute whether the generated loop runtime is ready
bool AMP_MG5_{name}::Ready() const noexcept {{ return madloop_.Ready(); }}

// Compute the generated loop runtime diagnostic
std::string AMP_MG5_{name}::Error() const {{ return madloop_.Error(); }}

// Compute the number of expensive loop evaluations
std::size_t AMP_MG5_{name}::HardCallCount() const noexcept {{
  return hard_calls_.load(std::memory_order_relaxed);
}}

// Compute the number of exact hard-point cache reuses
std::size_t AMP_MG5_{name}::HardReuseCount() const noexcept {{
  return hard_reuses_.load(std::memory_order_relaxed);
}}

}}  // namespace gra
"""


# Parse and execute the standalone source-bundle builder
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    bundle = subparsers.add_parser("bundle")
    bundle.add_argument("--source", type=Path, required=True)
    bundle.add_argument("--work", type=Path, required=True)
    bundle.add_argument("--output", type=Path, required=True)
    bundle.add_argument("--compiler", required=True)
    bundle.add_argument("--cxx", required=True)
    arguments = parser.parse_args()
    build_bundle(
        arguments.source,
        arguments.work,
        arguments.output,
        arguments.compiler,
        arguments.cxx,
    )


if __name__ == "__main__":
    main()
