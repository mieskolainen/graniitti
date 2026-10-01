#!/usr/bin/env python3
#
# Generate exact finite-Nc color data from MadGraph color bases
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import ast
import hashlib
import io
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

from core.io.files import ensure_dir

NC = 3
INCOMING_GLUON_DIMENSION = NC * NC - 1


# Compute actual amplitude coupling orders from generated HELAS diagrams
def coupling_orders(matrix_element: Any) -> list[dict[str, int]]:
    rows = {tuple(sorted(diagram.calculate_orders().items())) for diagram in matrix_element.get("diagrams")}
    return [dict(row) for row in sorted(rows)]


# Evaluate reference couplings through the selected UFO model and parameter card
def model_couplings(model: Any, card: str, charge: str | None = None) -> tuple[float, float]:
    from models.check_param_card import ParamCard
    from models.model_reader import ModelReader

    reader = ModelReader(model)
    reader.set_parameters_and_couplings(ParamCard(io.StringIO(card)))
    parameters = reader.get("parameter_dict")
    alpha_s = float(complex(parameters.get("aS", 0.0)).real)
    if charge is not None:
        if charge not in parameters:
            raise RuntimeError(f"Configured electromagnetic charge {charge} is not a UFO parameter")
        alpha_qed = abs(parameters[charge]) ** 2 / (4.0 * math.pi)
    else:
        inverse = next((p for p in model.get("parameters").get(("external",), [])
                        if p.lhablock.lower() == "sminputs" and tuple(p.lhacode) == (1,)), None)
        alpha_qed = 1.0 / float(complex(parameters[inverse.name]).real) if inverse else 0.0
    if not all(math.isfinite(value) and value >= 0.0 for value in (alpha_s, alpha_qed)):
        raise RuntimeError("Invalid generated UFO reference coupling")
    return alpha_s, alpha_qed


# Compute a stable signature for one simplified MadGraph color factor
def color_factor_signature(factor: Any) -> tuple[Any, ...]:
    terms = []
    for term in factor:
        terms.append(
            (
                term.to_immutable(),
                term.coeff.numerator,
                term.coeff.denominator,
                bool(term.is_imaginary),
                int(term.Nc_power),
                int(term.loop_Nc_power),
            )
        )
    return tuple(sorted(terms))


# Contract the two incoming adjoint indices into a normalized SU(3) singlet
def project_incoming_singlet(basis_keys: list[Any], color_algebra: Any) -> list[Any]:
    projected = []
    for immutable in basis_keys:
        color_string = color_algebra.ColorString()
        color_string.from_immutable(immutable)
        color_string.replace_indices({2: 1})
        factor = color_algebra.ColorFactor([color_string]).full_simplify()
        for term in factor:
            remaining = [index for obj in term for index in obj]
            if 1 in remaining or 2 in remaining:
                raise RuntimeError("MadGraph did not fully contract the incoming gluon indices")
        projected.append(factor)
    return projected


# Rename summed indices in the conjugate factor before forming an inner product
def independent_conjugate(color_string: Any, occupied: set[int]) -> Any:
    conjugate = color_string.create_copy().complex_conjugate()
    negative = sorted({index for obj in conjugate for index in obj if index < 0})
    if not negative:
        return conjugate

    used = set(occupied)
    replacement = {}
    candidate = min(used | set(negative) | {-1}) - 1
    for index in negative:
        while candidate in used:
            candidate -= 1
        replacement[index] = candidate
        used.add(candidate)
        candidate -= 1
    conjugate.replace_indices(replacement)
    return conjugate


# Evaluate the exact SU(3) inner product of two color factors
def color_factor_inner(
    left: Any, right: Any, color_algebra: Any, normalization: float = 1.0
) -> complex:
    if not math.isfinite(normalization) or normalization <= 0.0:
        raise RuntimeError("Color-factor normalization must be positive")
    value = 0.0j
    for left_term in left:
        for right_term in right:
            occupied = {index for obj in right_term for index in obj}
            product = independent_conjugate(left_term, occupied)
            product.product(right_term.create_copy())
            reduced = color_algebra.ColorFactor([product]).full_simplify()
            real, imaginary = reduced.set_Nc(NC)
            value += complex(float(real), float(imaginary))
    return value / normalization


# Build one normalized finite-Nc color Gram matrix
def color_gram(
    factors: list[Any], color_algebra: Any, normalization: float = 1.0
) -> list[list[complex]]:
    size = len(factors)
    gram = [[0.0j for _ in range(size)] for _ in range(size)]
    cache: dict[tuple[Any, Any], complex] = {}

    for row in range(size):
        left_signature = color_factor_signature(factors[row])
        for column in range(row, size):
            right_signature = color_factor_signature(factors[column])
            key = (left_signature, right_signature)
            reverse = (right_signature, left_signature)
            if key in cache:
                value = cache[key]
            elif reverse in cache:
                value = cache[reverse].conjugate()
            else:
                value = color_factor_inner(
                    factors[row], factors[column], color_algebra, normalization
                )
                cache[key] = value
            gram[row][column] = value
            gram[column][row] = value.conjugate()
    return gram


# Build the normalized incoming-singlet restricted color Gram matrix
def restricted_gram(projected: list[Any], color_algebra: Any) -> list[list[complex]]:
    return color_gram(projected, color_algebra, INCOMING_GLUON_DIMENSION)


# Factor a Hermitian positive-semidefinite matrix as gram = P^dagger P
def pivoted_cholesky(
    gram: list[list[complex]], relative_tolerance: float = 1.0e-12
) -> list[list[complex]]:
    size = len(gram)
    if size == 0 or any(len(row) != size for row in gram):
        raise RuntimeError("Durham color Gram matrix must be non-empty and square")

    scale = max(abs(gram[index][index].real) for index in range(size))
    if not math.isfinite(scale) or scale <= 0.0:
        raise RuntimeError("Durham color Gram matrix has no positive support")
    tolerance = relative_tolerance * scale

    permutation = list(range(size))
    lower = [[0.0j for _ in range(size)] for _ in range(size)]
    rank = 0
    for column in range(size):
        pivot_position = column
        pivot_value = -math.inf
        for row in range(column, size):
            original = permutation[row]
            residual = gram[original][original].real
            residual -= sum(abs(lower[row][k]) ** 2 for k in range(column))
            if residual > pivot_value:
                pivot_value = residual
                pivot_position = row

        if pivot_value <= tolerance:
            break
        if pivot_position != column:
            permutation[column], permutation[pivot_position] = (
                permutation[pivot_position],
                permutation[column],
            )
            lower[column], lower[pivot_position] = (
                lower[pivot_position],
                lower[column],
            )

        lower[column][column] = math.sqrt(pivot_value)
        pivot = permutation[column]
        for row in range(column + 1, size):
            original = permutation[row]
            residual = gram[original][pivot]
            residual -= sum(lower[row][k] * lower[column][k].conjugate() for k in range(column))
            lower[row][column] = residual / lower[column][column]
        rank += 1

    projectors = [[0.0j for _ in range(size)] for _ in range(rank)]
    for component in range(rank):
        for row, original in enumerate(permutation):
            projectors[component][original] = lower[row][component].conjugate()

    residual = 0.0
    for row in range(size):
        for column in range(size):
            reconstructed = sum(
                projector[row].conjugate() * projector[column] for projector in projectors
            )
            residual = max(residual, abs(reconstructed - gram[row][column]))
    if residual > 100.0 * tolerance:
        raise RuntimeError(f"Unstable Durham Gram factorization: residual {residual:.3e}")
    return projectors


# Compute the leading-color flow carried by one projected color string
def leading_color_flow(
    term: Any,
    final_representations: dict[int, int],
    color_amp: Any,
) -> tuple[tuple[int, int], ...]:
    basis = color_amp.ColorBasis()
    basis[term.to_immutable()] = []
    decomposition = basis.color_flow_decomposition(final_representations, 0)
    if len(decomposition) != 1:
        raise RuntimeError("MadGraph returned an ambiguous leading-color flow")
    flow = decomposition[0]
    return tuple(tuple(flow[leg]) for leg in sorted(final_representations))


# Partition each exact projector into MadGraph leading-color shower candidates
def shower_flow_data(
    projected: list[Any],
    final_representations: dict[int, int],
    color_amp: Any,
) -> tuple[list[list[list[int]]], list[list[float]]]:
    colored = [
        leg for leg, representation in final_representations.items() if abs(representation) != 1
    ]
    if not colored:
        return [], []

    column_flows: list[list[tuple[tuple[tuple[int, int], ...], float]]] = []
    candidates: list[tuple[tuple[int, int], ...]] = []
    for factor in projected:
        nonzero = [term for term in factor if term.coeff]
        if not nonzero:
            column_flows.append([])
            continue
        leading_power = max(term.Nc_power for term in nonzero)
        weights: dict[tuple[tuple[int, int], ...], float] = {}
        for term in nonzero:
            if term.Nc_power != leading_power:
                continue
            flow = leading_color_flow(term, final_representations, color_amp)
            magnitude = abs(float(term.coeff)) * float(NC**term.Nc_power)
            weights[flow] = weights.get(flow, 0.0) + magnitude
            if flow not in candidates:
                candidates.append(flow)
        column_flows.append(list(weights.items()))

    flow_weights = [[0.0 for _ in range(len(projected))] for _ in range(len(candidates))]
    for column, contributions in enumerate(column_flows):
        normalization = sum(weight for _, weight in contributions)
        if normalization <= 0.0:
            continue
        for flow, weight in contributions:
            flow_weights[candidates.index(flow)][column] = weight / normalization

    serialized = [
        [[int(color), int(anticolor)] for color, anticolor in candidate] for candidate in candidates
    ]
    return serialized, flow_weights


# Serialize one complex matrix without relying on JSON implementation details
def serialize_complex_matrix(matrix: list[list[complex]]) -> list[list[list[float]]]:
    return [[[float(value.real), float(value.imag)] for value in row] for row in matrix]


# Compute the maximum reconstruction residual of one projector factorization
def factorization_residual(gram: list[list[complex]], projectors: list[list[complex]]) -> float:
    residual = 0.0
    for row in range(len(gram)):
        for column in range(len(gram)):
            reconstructed = sum(
                projector[row].conjugate() * projector[column] for projector in projectors
            )
            residual = max(residual, abs(reconstructed - gram[row][column]))
    return residual


# Load MadGraph lazily from the installation owning the requested executable
def load_madgraph(mg5_root: Path) -> tuple[Any, Any, Any]:
    root = str(mg5_root.resolve())
    if root not in sys.path:
        sys.path.insert(0, root)
    from madgraph.core import color_algebra, color_amp
    from madgraph.interface.master_interface import MasterCmd

    return MasterCmd, color_amp, color_algebra


# Generate complete Durham color data for one concrete MadGraph process
def generate_durham_data(
    mg5_root: Path,
    model_import: str,
    process: str,
    work_dir: Path | None = None,
) -> dict[str, Any]:
    previous_directory = Path.cwd()
    parser_directory = (work_dir or previous_directory).resolve()
    ensure_dir(parser_directory)
    try:
        os.chdir(parser_directory)
        MasterCmd, color_amp, color_algebra = load_madgraph(mg5_root)
        command = MasterCmd()
        command.no_notification()
        command.exec_cmd(f"import model {model_import}", printcmd=False, precmd=True, postcmd=True)
        command.exec_cmd(f"generate {process}", printcmd=False, precmd=True, postcmd=True)
        amplitudes = command._curr_amps
    finally:
        os.chdir(previous_directory)
    if len(amplitudes) != 1:
        raise RuntimeError(
            "--durham requires one concrete MadGraph subprocess, "
            f"but generation produced {len(amplitudes)}"
        )

    amplitude = amplitudes[0]
    process_object = amplitude.get("process")
    legs = list(process_object.get("legs"))
    incoming = [leg for leg in legs if not leg.get("state")]
    outgoing = [leg for leg in legs if leg.get("state")]
    incoming_pdgs = [int(leg.get("id")) for leg in incoming]
    if incoming_pdgs != [21, 21]:
        raise RuntimeError(
            f"--durham requires the concrete incoming state g g, not {incoming_pdgs}"
        )

    basis = color_amp.ColorBasis(amplitude)
    basis_keys = sorted(basis.keys())
    if not basis_keys:
        raise RuntimeError("MadGraph generated no color basis for the Durham process")

    projected = project_incoming_singlet(basis_keys, color_algebra)
    gram = restricted_gram(projected, color_algebra)
    projectors = pivoted_cholesky(gram)

    model = process_object.get("model")
    final_representations = {
        int(leg.get("number")): int(model.get_particle(leg.get("id")).get_color())
        for leg in outgoing
    }
    flow_candidates, flow_weights = shower_flow_data(projected, final_representations, color_amp)

    basis_text = [str(key) for key in basis_keys]
    basis_hash = hashlib.sha256("\n".join(basis_text).encode()).hexdigest()
    return {
        "version": 1,
        "process": process,
        "incoming_pdgs": incoming_pdgs,
        "final_pdgs": [int(leg.get("id")) for leg in outgoing],
        "final_color_representations": [
            final_representations[int(leg.get("number"))] for leg in outgoing
        ],
        "ncolor": len(basis_keys),
        "rank": len(projectors),
        "basis_sha256": basis_hash,
        "factorization_residual": factorization_residual(gram, projectors),
        "projectors": serialize_complex_matrix(projectors),
        "flow_candidates": flow_candidates,
        "flow_weights": flow_weights,
    }



# Parse one C++ brace initializer as a Python scalar array
def parse_initializer(text: str) -> Any:
    return ast.literal_eval(text.replace("{", "[").replace("}", "]"))


# Extract one complete standalone helicity table from generated C++
def extract_helicities(source: str) -> list[list[int]]:
    match = re.search(
        r"(?:static )?const int helicities\[ncomb\]\[nexternal\]\s*=\s*(\{\{.*?\}\});",
        source,
        flags=re.S,
    )
    if match is None:
        raise RuntimeError("Could not extract the generated helicity table")
    table = parse_initializer(match.group(1))
    if not table or any(len(row) != len(table[0]) for row in table):
        raise RuntimeError("Generated helicity table is empty or ragged")
    return [[int(value) for value in row] for row in table]


# Extract the sole generated spin, color and symmetry denominator
def extract_process_denominator(source: str) -> int:
    match = re.search(
        r"const int denominators\[nprocesses\]\s*=\s*(\{.*?\});",
        source,
        flags=re.S,
    )
    if match is None:
        raise RuntimeError("Could not extract the generated process denominator")
    values = parse_initializer(match.group(1))
    if len(values) != 1 or int(values[0]) <= 0:
        raise RuntimeError("Photon registry requires one positive process denominator")
    return int(values[0])


# Compute the factorial symmetry factor for stable final-state PDG ids
def final_state_symmetry_factor(final_pdgs: list[int]) -> int:
    multiplicities: dict[int, int] = {}
    for pdg in final_pdgs:
        multiplicities[pdg] = multiplicities.get(pdg, 0) + 1
    return math.prod(math.factorial(count) for count in multiplicities.values())


# Convert one immutable MadGraph basis tensor into a simplified color factor
def basis_factor(immutable: Any, color_algebra: Any) -> Any:
    color_string = color_algebra.ColorString()
    color_string.from_immutable(immutable)
    return color_algebra.ColorFactor([color_string]).full_simplify()


# Compute stable external PDGs in generated momentum order
def external_pdgs(matrix_element: Any) -> list[int]:
    external = sorted(
        matrix_element.get_external_wavefunctions(),
        key=lambda wavefunction: int(wavefunction.get("number_external")),
    )
    return [int(wavefunction.get("pdg_code")) for wavefunction in external]


# Compute the raw MadGraph color metric embedded in standalone C++
def generated_color_gram(color_structure: dict[str, str | int]) -> list[list[complex]]:
    denominators = [float(value) for value in parse_initializer(str(color_structure["denom"]))]
    factors = parse_initializer(str(color_structure["cf"]))
    size = int(color_structure["ncolor"])
    if len(denominators) != size or len(factors) != size:
        raise RuntimeError("Generated standalone color metric has inconsistent dimensions")
    gram: list[list[complex]] = []
    for row in range(size):
        if denominators[row] <= 0.0 or len(factors[row]) != size:
            raise RuntimeError("Generated standalone color metric is invalid")
        gram.append(
            [complex(float(factors[row][column]) / denominators[row]) for column in range(size)]
        )
    return gram


# Compute the maximum elementwise difference of two square matrices
def matrix_residual(left: list[list[complex]], right: list[list[complex]]) -> float:
    if len(left) != len(right) or any(
        len(left[row]) != len(right[row]) for row in range(len(left))
    ):
        return math.inf
    return max(
        abs(left[row][column] - right[row][column])
        for row in range(len(left))
        for column in range(len(left[row]))
    )


# Generate complete color and helicity data for one incoming-photon subprocess
def generate_photon_data(
    mg5_root: Path,
    model_import: str,
    process: str,
    raw_source: str,
    param_card: str,
    color_structure: dict[str, str | int],
    work_dir: Path | None = None,
    charge: str | None = None,
    complex_mass_scheme: bool = False,
) -> dict[str, Any]:
    previous_directory = Path.cwd()
    parser_directory = (work_dir or previous_directory).resolve()
    ensure_dir(parser_directory)
    try:
        os.chdir(parser_directory)
        MasterCmd, color_amp, color_algebra = load_madgraph(mg5_root)
        from madgraph.core import helas_objects

        command = MasterCmd()
        command.no_notification()
        command.exec_cmd("set complex_mass_scheme " + ("True --allow_qed" if complex_mass_scheme else "False"), printcmd=False)
        command.exec_cmd(f"import model {model_import}", printcmd=False, precmd=True, postcmd=True)
        command.exec_cmd(f"generate {process}", printcmd=False, precmd=True, postcmd=True)
        matrix_elements = helas_objects.HelasMultiProcess.generate_matrix_elements(
            command._curr_amps
        )
    finally:
        os.chdir(previous_directory)

    if len(matrix_elements) != 1:
        raise RuntimeError(
            "Photon automation requires one concrete MadGraph subprocess, "
            f"but generation produced {len(matrix_elements)}"
        )

    matrix_element = matrix_elements[0]
    pdgs = external_pdgs(matrix_element)
    if pdgs[:2] != [22, 22]:
        raise RuntimeError(
            f"Photon automation requires the concrete incoming state a a, not {pdgs[:2]}"
        )

    basis = matrix_element.get("color_basis")
    basis_keys = sorted(basis.keys())
    if not basis_keys:
        raise RuntimeError("MadGraph generated no color basis for the photon process")
    factors = [basis_factor(key, color_algebra) for key in basis_keys]
    gram = color_gram(factors, color_algebra)
    standalone_gram = generated_color_gram(color_structure)
    residual = matrix_residual(gram, standalone_gram)
    scale = max(abs(value) for row in standalone_gram for value in row)
    if residual > 1.0e-10 * max(1.0, scale):
        raise RuntimeError(
            "MadGraph API and standalone color metrics disagree for "
            f"{process}: residual {residual:.3e}"
        )

    projectors = pivoted_cholesky(gram)
    final_pdgs = pdgs[2:]
    model = matrix_element.get("processes")[0].get("model")
    orders = coupling_orders(matrix_element)
    alpha_s, alpha_qed = model_couplings(model, param_card, charge)
    final_representations = [int(model.get_particle(pdg).get_color()) for pdg in final_pdgs]
    representation_map = {
        index + 3: representation for index, representation in enumerate(final_representations)
    }
    flow_candidates, flow_weights = shower_flow_data(
        factors, representation_map, color_amp
    )

    helicities = extract_helicities(raw_source)
    if len(helicities) == 0 or len(helicities[0]) != len(pdgs):
        raise RuntimeError("Generated helicities do not match the photon external state")
    process_denominator = extract_process_denominator(raw_source)
    final_symmetry_factor = final_state_symmetry_factor(final_pdgs)
    if process_denominator % final_symmetry_factor != 0:
        raise RuntimeError(
            "Photon process denominator does not factor into initial-state "
            "average and final-state symmetry"
        )
    basis_hash = hashlib.sha256("\n".join(str(key) for key in basis_keys).encode()).hexdigest()
    return {
        "version": 1,
        "process": process,
        "incoming_pdgs": pdgs[:2],
        "final_pdgs": final_pdgs,
        "final_color_representations": final_representations,
        "has_decay_chain": "// *   Decay:" in raw_source,
        "ncolor": len(basis_keys),
        "rank": len(projectors),
        "process_denominator": process_denominator,
        "final_symmetry_factor": final_symmetry_factor,
        "initial_state_denominator": process_denominator // final_symmetry_factor,
        "default_alpha_s": alpha_s,
        "default_alpha_qed": alpha_qed,
        "coupling_orders": orders,
        "alpha_s_power": max(row.get("QCD", 0) for row in orders),
        "alpha_qed_power": max(row.get("QED", 0) for row in orders),
        "helicities": helicities,
        "basis_sha256": basis_hash,
        "color_metric_residual": residual,
        "factorization_residual": factorization_residual(gram, projectors),
        "projectors": serialize_complex_matrix(projectors),
        "flow_candidates": flow_candidates,
        "flow_weights": flow_weights,
    }
