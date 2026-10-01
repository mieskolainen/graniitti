#!/usr/bin/env python3
# Parse mapped SOFT exchanges and project their physical proton vertices
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from dataclasses import dataclass

REGGE_ROLES = frozenset({"pomeron", "reggeon", "odderon"})
TRAJECTORY_MODES = frozenset({"pion_loop", "linear"})
TRANSITION_FORM_FACTORS = frozenset({"geometric", "arithmetic", "diagonal"})
ETA_MODES = frozenset({"raw", "rotating_t0", "rotating"})


@dataclass(frozen=True)
class Exchange:
    """Regge alias and SOFT exchange mapping"""

    role: str
    pdg: tuple[int, ...]
    pole_spin: int
    soft_exchange: str


@dataclass(frozen=True)
class ProtonVertex:
    """Proton coupling and form factor"""

    beam_residue_per_gev: float
    amplitude_slope_per_gev2: float
    form_factor: str
    parameters: tuple[tuple[float, ...], ...]
    source: str

    # Compute twice the amplitude form factor slope
    @property
    def cross_section_slope_per_gev2(self) -> float:
        return 2.0 * self.amplitude_slope_per_gev2


@dataclass(frozen=True)
class Config:
    """Validated SOFT model and Regge mappings"""

    model_name: str
    channels: int
    forward_excitation_exchange: str
    regge_exchanges: tuple[Exchange, ...]
    definitions: dict[str, object]
    model: dict[str, object]


# Convert one JSON number to a finite float
def finite_number(value: object, context: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{context} must be numeric")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{context} must be finite")
    return out


# Convert one JSON number to an exact integer
def exact_integer(value: object, context: str) -> int:
    numeric = finite_number(value, context)
    if not numeric.is_integer():
        raise ValueError(f"{context} must be an exact integer")
    return int(numeric)


# Parse the typed central Regge exchange rows
def exchange_rows(regge: object) -> tuple[Exchange, ...]:
    if not isinstance(regge, dict):
        raise TypeError("PARAM_REGGE must be an object")
    entries = regge.get("EXCHANGES")
    if not isinstance(entries, list) or not entries:
        raise ValueError("PARAM_REGGE.EXCHANGES must be a nonempty array")

    rows = []
    seen_pdgs: set[int] = set()
    seen_soft: set[str] = set()
    for index, entry in enumerate(entries):
        context = f"PARAM_REGGE.EXCHANGES[{index}]"
        if not isinstance(entry, dict):
            raise TypeError(f"{context} must be an object")
        if set(entry) != {"role", "pdg", "pole_spin", "soft_exchange"}:
            raise ValueError(
                f"{context} must contain only role, pdg, pole_spin, and soft_exchange"
            )
        role = entry.get("role")
        if role not in REGGE_ROLES:
            raise ValueError(f"{context}.role must be pomeron, reggeon, or odderon")
        pdg_values = entry.get("pdg")
        if not isinstance(pdg_values, list) or not pdg_values:
            raise ValueError(f"{context}.pdg must be a nonempty array")
        pdgs = []
        for offset, value in enumerate(pdg_values):
            pdg = exact_integer(value, f"{context}.pdg[{offset}]")
            if pdg == 0:
                raise ValueError(f"{context}.pdg[{offset}] must be nonzero")
            if pdg in seen_pdgs:
                raise ValueError(f"PARAM_REGGE PDG alias {pdg} occurs in multiple rows")
            seen_pdgs.add(pdg)
            pdgs.append(pdg)
        pole_spin = exact_integer(entry.get("pole_spin"), f"{context}.pole_spin")
        if pole_spin <= 0:
            raise ValueError(f"{context}.pole_spin must be positive")
        soft_exchange = entry.get("soft_exchange")
        if not isinstance(soft_exchange, str) or not soft_exchange:
            raise ValueError(f"{context}.soft_exchange must be a nonempty name")
        if soft_exchange in seen_soft:
            raise ValueError(
                f"PARAM_REGGE soft exchange {soft_exchange} is mapped by multiple rows"
            )
        seen_soft.add(soft_exchange)
        rows.append(Exchange(str(role), tuple(pdgs), pole_spin, soft_exchange))
    return tuple(rows)


# Compute the physical proton row of the Good Walker rotation
def proton_state(theta: object, channels: int) -> list[float]:
    if channels < 1:
        raise ValueError("the SOFT channel count must be positive")
    if not isinstance(theta, list):
        raise TypeError("PARAM_SOFT.GW.theta must be an array")
    expected = channels * (channels - 1) // 2
    if len(theta) != expected:
        raise ValueError("PARAM_SOFT.GW.theta does not match the exchange coupling matrices")
    rotation = [
        [1.0 if row == column else 0.0 for column in range(channels)] for row in range(channels)
    ]
    angle_index = 0
    for row in range(channels - 1):
        for column in range(row + 1, channels):
            angle = finite_number(theta[angle_index], f"PARAM_SOFT.GW.theta[{angle_index}]")
            angle_index += 1
            cosine = math.cos(angle)
            sine = math.sin(angle)
            updated = [entry[:] for entry in rotation]
            for state in range(channels):
                updated[row][state] = cosine * rotation[row][state] + sine * rotation[column][state]
                updated[column][state] = (
                    -sine * rotation[row][state] + cosine * rotation[column][state]
                )
            rotation = updated
    return rotation[0]


# Validate and normalize the resolved Good Walker excitation direction
def direction(a_c: object, channels: int) -> list[float]:
    if channels < 1:
        raise ValueError("the SOFT channel count must be positive")
    if not isinstance(a_c, list):
        raise TypeError("PARAM_SOFT.GW.a_c must be an array")
    if len(a_c) != channels - 1:
        raise ValueError("PARAM_SOFT.GW.a_c must contain N-1 coefficients")
    values = [
        finite_number(value, f"PARAM_SOFT.GW.a_c[{index}]") for index, value in enumerate(a_c)
    ]
    if channels == 1:
        return values
    scale = max(abs(value) for value in values)
    if scale <= 0.0:
        raise ValueError("PARAM_SOFT.GW.a_c direction must be nonzero")
    values = [value / scale for value in values]
    norm = math.hypot(*values)
    return [value / norm for value in values]


# Compute the local amplitude slope of one SOFT form factor row
def form_slope(
    form_factor: str,
    parameters: tuple[float, ...],
) -> float:
    values = tuple(
        finite_number(value, f"{form_factor} form factor parameter") for value in parameters
    )
    if form_factor == "EXP":
        if len(values) != 1 or values[0] <= 0.0:
            raise ValueError("EXP proton form factor requires positive [B]")
        slope = 0.5 * values[0]
    elif form_factor in {"DPOW", "EXPOW"}:
        if len(values) != 3 or any(value <= 0.0 for value in values):
            raise ValueError(f"{form_factor} proton form factor requires positive [b,c,d]")
        b_value, c_value, power = values
        slope = (
            power * (1.0 / b_value + 1.0 / c_value)
            if form_factor == "DPOW"
            else power * b_value * (b_value * c_value) ** (power - 1.0)
        )
    elif form_factor == "MIXEXP":
        if len(values) != 4:
            raise ValueError("MIXEXP proton form factor requires [a,B1,B2,t_zero]")
        weight, slope1, slope2, zero = values
        if not 0.0 <= weight <= 1.0 or slope1 < 0.0 or slope2 < 0.0 or zero < 0.0:
            raise ValueError("MIXEXP parameters are outside their physical ranges")
        slope = 0.5 * ((1.0 - weight) * slope1 + weight * slope2)
        if zero > 0.0:
            slope += 1.0 / zero
    elif form_factor == "ODD3G_NODE":
        if len(values) != 3:
            raise ValueError("ODD3G_NODE proton form factor requires [B,t_zero,m2]")
        slope0, zero, mass2 = values
        if slope0 < 0.0 or zero <= 0.0 or mass2 <= 0.0:
            raise ValueError("ODD3G_NODE parameters are outside their physical ranges")
        slope = 0.5 * slope0 + 1.0 / zero + 3.0 / mass2
    elif form_factor == "3G":
        if len(values) != 1 or values[0] <= 0.0:
            raise ValueError("3G proton form factor requires positive [m2]")
        slope = 2.0 / values[0]
    elif form_factor == "GKERNEL":
        block_size = 4
        if len(values) > 1 and len(values) % block_size == 1:
            if values[0] < 0.0:
                raise ValueError("GKERNEL prefactor must be nonnegative")
            slope = -values[0]
            offset = 1
        elif values and len(values) % block_size == 0:
            slope = 0.0
            offset = 0
        else:
            raise ValueError("GKERNEL requires 4*N or 1+4*N parameters")
        for index in range(offset, len(values), block_size):
            scale, power, deformation, shift = values[index : index + block_size]
            if scale < 0.0 or power <= 0.0 or deformation < 0.0 or shift < 0.0:
                raise ValueError("GKERNEL parameters are outside their physical ranges")
            if shift == 0.0 and power < 1.0:
                raise ValueError("GKERNEL has a nonfinite forward slope")
            slope += (
                scale
                if shift == 0.0 and power == 1.0
                else (0.0 if shift == 0.0 else scale * power * shift ** (power - 1.0))
            )
    else:
        raise ValueError(f"unsupported proton form factor: {form_factor}")
    if not math.isfinite(slope):
        raise ValueError("proton form factor must have a finite forward slope")
    return slope


# Validate one symmetric coupling matrix
def _coupling_matrix(value: object, channels: int | None, context: str) -> list[list[float]]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{context} must be a nonempty matrix")
    size = len(value)
    if channels is not None and size != channels:
        raise ValueError(f"{context} must use the common N channel count")
    matrix = []
    for row_index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != size:
            raise ValueError(f"{context} must be an N x N matrix")
        matrix.append(
            [
                finite_number(entry, f"{context}[{row_index}][{column}]")
                for column, entry in enumerate(row)
            ]
        )
    for row in range(size):
        for column in range(row + 1, size):
            scale = max(1.0, abs(matrix[row][column]), abs(matrix[column][row]))
            if abs(matrix[row][column] - matrix[column][row]) > 1.0e-12 * scale:
                raise ValueError(f"{context} must be symmetric")
    return matrix


# Validate one helicity transition scalar or symmetric channel matrix
def _helicity_transition(
    value: object,
    channels: int,
    context: str,
    require_nonnegative: bool,
) -> None:
    if isinstance(value, list):
        values = [entry for row in _coupling_matrix(value, channels, context) for entry in row]
    else:
        values = [finite_number(value, context)]
    if require_nonnegative and any(entry < 0.0 for entry in values):
        raise ValueError(f"{context} must be nonnegative")


# Validate one SOFT exchange definition
def _exchange_definition(value: object, name: str) -> dict[str, object]:
    context = f"PARAM_SOFT.EXCHANGE_DEF.{name}"
    if not isinstance(value, dict):
        raise TypeError(f"{context} must be an object")
    expected = {"role", "trajectory_mode", "tau", "crossing"}
    if set(value) != expected:
        raise ValueError(f"{context} must contain only role, trajectory_mode, tau, and crossing")
    role = value.get("role")
    if role not in REGGE_ROLES:
        raise ValueError(f"{context}.role must be pomeron, reggeon, or odderon")
    trajectory_mode = value.get("trajectory_mode")
    if trajectory_mode not in TRAJECTORY_MODES:
        raise ValueError(f"{context}.trajectory_mode must be pion_loop or linear")
    if trajectory_mode == "pion_loop" and role != "pomeron":
        raise ValueError(f"{context} uses pion_loop for a nonpomeron exchange")
    tau = exact_integer(value.get("tau"), f"{context}.tau")
    crossing = exact_integer(value.get("crossing"), f"{context}.crossing")
    if tau not in {-1, 1} or crossing not in {-1, 1}:
        raise ValueError(f"{context} has invalid signature quantum numbers")
    if tau != crossing:
        raise ValueError(f"{context}.tau and crossing must agree")
    if role == "pomeron" and tau != 1:
        raise ValueError(f"{context} pomeron requires positive signature")
    if role == "odderon" and (tau != -1 or trajectory_mode != "linear"):
        raise ValueError(f"{context} odderon requires negative signature and linear mode")
    if role == "reggeon" and trajectory_mode != "linear":
        raise ValueError(f"{context} reggeon requires linear mode")
    return value


# Validate one SOFT exchange parameter and form factor block
def _validate_exchange_model(
    model: dict[str, object],
    name: str,
    channels: int,
    model_name: str,
    definition: dict[str, object],
    helicity_enabled: bool,
) -> None:
    exchanges = model.get("EXCHANGE")
    if not isinstance(exchanges, dict) or name not in exchanges:
        raise ValueError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE must contain {name}")
    exchange = exchanges[name]
    context = f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE.{name}"
    if not isinstance(exchange, dict):
        raise TypeError(f"{context} must be an object")
    alpha = exchange.get("alpha")
    if not isinstance(alpha, list) or len(alpha) != 2:
        raise ValueError(f"{context}.alpha must contain intercept and slope")
    finite_number(alpha[0], f"{context}.alpha[0]")
    finite_number(alpha[1], f"{context}.alpha[1]")
    if not isinstance(exchange.get("on"), bool):
        raise TypeError(f"{context}.on must be boolean")
    sign = exact_integer(exchange.get("sign"), f"{context}.sign")
    if sign not in {-1, 1}:
        raise ValueError(f"{context}.sign must be -1 or 1")
    eta_mode = exchange.get("eta_mode")
    if eta_mode not in ETA_MODES:
        raise ValueError(f"{context}.eta_mode has an unknown mode")
    _coupling_matrix(exchange.get("g"), channels, f"{context}.g")
    transition = exchange.get("transition_ff")
    if transition not in TRANSITION_FORM_FACTORS:
        raise ValueError(f"{context}.transition_ff has an unknown mode")
    bank_name = exchange.get("ff")
    banks = model.get("FF")
    if not isinstance(bank_name, str) or not isinstance(banks, dict) or bank_name not in banks:
        raise ValueError(f"{context}.ff selects an unknown form factor bank")
    bank = banks[bank_name]
    if definition["role"] == "pomeron" and bank["type"] in {"3G", "ODD3G_NODE"}:
        raise ValueError(f"{context} pomeron cannot use a three gluon form factor")
    if helicity_enabled:
        helicity = exchange.get("helicity")
        if not isinstance(helicity, dict):
            raise TypeError(f"{context}.helicity must be an object")
        _helicity_transition(
            helicity.get("kappa"),
            channels,
            f"{context}.helicity.kappa",
            False,
        )
        _helicity_transition(
            helicity.get("B_kappa"),
            channels,
            f"{context}.helicity.B_kappa",
            True,
        )


# Validate every named SOFT form factor bank
def _validate_form_factor_banks(model: dict[str, object], channels: int, model_name: str) -> None:
    banks = model.get("FF")
    if not isinstance(banks, dict) or not banks:
        raise ValueError(f"PARAM_SOFT.MODEL.{model_name}.FF must be a nonempty object")
    for bank_name, bank in banks.items():
        context = f"PARAM_SOFT.MODEL.{model_name}.FF.{bank_name}"
        if not isinstance(bank, dict):
            raise TypeError(f"{context} must be an object")
        form_factor = bank.get("type")
        if not isinstance(form_factor, str):
            raise TypeError(f"{context}.type must be a form factor name")
        parameters = bank.get("param")
        if not isinstance(parameters, list) or len(parameters) != channels:
            raise ValueError(f"{context}.param must have N rows")
        for row_index, row in enumerate(parameters):
            if not isinstance(row, list):
                raise TypeError(f"{context}.param[{row_index}] must be an array")
            form_slope(form_factor, tuple(row))


# Validate the eikonal proton helicity switch
def _validate_helicity_controls(model: dict[str, object], model_name: str) -> bool:
    eikonal = model.get("EIKONAL")
    helicity = eikonal.get("helicity") if isinstance(eikonal, dict) else None
    context = f"PARAM_SOFT.MODEL.{model_name}.EIKONAL.helicity"
    if not isinstance(helicity, bool):
        raise TypeError(f"{context} must be a boolean")
    return helicity


# Validate the complete matrix eikonal control block
def _validate_eikonal_controls(model: dict[str, object], model_name: str) -> None:
    eikonal = model.get("EIKONAL")
    context = f"PARAM_SOFT.MODEL.{model_name}.EIKONAL"
    if not isinstance(eikonal, dict):
        raise TypeError(f"{context} must be an object")
    unitarization = eikonal.get("unitarization")
    if unitarization not in {"exp", "q_exp"}:
        raise ValueError(f"{context}.unitarization must be exp or q_exp")
    q_value = finite_number(eikonal.get("q"), f"{context}.q")
    if not 0.0 < q_value <= 1.0:
        raise ValueError(f"{context}.q must be in (0,1]")
    if unitarization == "exp" and abs(q_value - 1.0) > 1.0e-12:
        raise ValueError(f"{context} exp unitarization requires q equal to one")


# Expand one explicit exchange selection or its wildcard
def _expand_exchange_selection(
    selection: object,
    exchanges: dict[str, object],
    definitions: dict[str, object],
    context: str,
    pomeron_only: bool,
) -> tuple[str, ...]:
    if not isinstance(selection, list) or not selection:
        raise ValueError(f"{context} must be a nonempty array")
    if selection == ["*"]:
        selected = tuple(
            name
            for name, exchange in exchanges.items()
            if isinstance(exchange, dict)
            and exchange.get("on") is True
            and (
                not pomeron_only
                or isinstance(definitions.get(name), dict)
                and definitions[name].get("role") == "pomeron"
            )
        )
        if not selected:
            raise ValueError(f"{context} wildcard selects no exchanges")
        return selected
    if "*" in selection:
        raise ValueError(f"{context} wildcard must be used alone")
    if any(not isinstance(name, str) or not name for name in selection):
        raise TypeError(f"{context} must contain exchange names")
    if len(set(selection)) != len(selection):
        raise ValueError(f"{context} must not contain duplicates")
    for name in selection:
        if name not in exchanges or name not in definitions:
            raise ValueError(f"{context} selects unknown exchange {name}")
        if exchanges[name].get("on") is not True:
            raise ValueError(f"{context} selects disabled exchange {name}")
        if pomeron_only and definitions[name].get("role") != "pomeron":
            raise ValueError(f"{context} must select pomeron exchanges")
    return tuple(selection)


# Validate named exchange selections owned by the active SOFT model
def _validate_selections(
    model: dict[str, object],
    definitions: dict[str, object],
    model_name: str,
) -> str:
    exchanges = model.get("EXCHANGE")
    if not isinstance(exchanges, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE must be an object")

    eikonal = model["EIKONAL"]
    screening_context = f"PARAM_SOFT.MODEL.{model_name}.EIKONAL.screening_exchanges"
    _expand_exchange_selection(
        eikonal.get("screening_exchanges"),
        exchanges,
        definitions,
        screening_context,
        False,
    )

    selected = eikonal.get("excitation_exchanges")
    context = f"PARAM_SOFT.MODEL.{model_name}.EIKONAL.excitation_exchanges"
    excitation = _expand_exchange_selection(
        selected, exchanges, definitions, context, True
    )
    if len(excitation) != 1:
        raise ValueError(f"{context} must resolve to exactly one pomeron exchange")
    return excitation[0]


# Load and validate one active mapped SOFT configuration
def load(general: object) -> Config:
    if not isinstance(general, dict):
        raise TypeError("GENERAL.json root must be an object")
    rows = exchange_rows(general.get("PARAM_REGGE"))
    soft = general.get("PARAM_SOFT")
    if not isinstance(soft, dict):
        raise TypeError("PARAM_SOFT must be an object")
    model_name = soft.get("active_model")
    models = soft.get("MODEL")
    if not isinstance(model_name, str) or not isinstance(models, dict) or model_name not in models:
        raise ValueError("PARAM_SOFT.active_model must select a defined model")
    model = models[model_name]
    if not isinstance(model, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name} must be an object")
    definitions = soft.get("EXCHANGE_DEF")
    if not isinstance(definitions, dict) or not definitions:
        raise ValueError("PARAM_SOFT.EXCHANGE_DEF must be a nonempty object")
    for name, definition in definitions.items():
        _exchange_definition(definition, str(name))

    pomeron_rows = [row for row in rows if row.role == "pomeron"]
    if len(pomeron_rows) != 1:
        raise ValueError("PARAM_REGGE.EXCHANGES must contain exactly one pomeron row")
    pomeron_name = pomeron_rows[0].soft_exchange
    exchanges = model.get("EXCHANGE")
    if not isinstance(exchanges, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE must be an object")
    definition_names = {str(name) for name in definitions}
    exchange_names = {str(name) for name in exchanges}
    missing_exchanges = definition_names - exchange_names
    unknown_exchanges = exchange_names - definition_names
    if missing_exchanges:
        missing = sorted(missing_exchanges)[0]
        raise ValueError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE is missing {missing}")
    if unknown_exchanges:
        unknown = sorted(unknown_exchanges)[0]
        raise ValueError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE has unknown exchange {unknown}")
    if pomeron_name not in exchanges:
        raise ValueError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE must contain {pomeron_name}")
    pomeron = exchanges[pomeron_name]
    if not isinstance(pomeron, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE.{pomeron_name} must be an object")
    channels = len(
        _coupling_matrix(
            pomeron.get("g"), None, f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE.{pomeron_name}.g"
        )
    )
    gw = model.get("GW")
    if not isinstance(gw, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.GW must be an object")
    proton_state(gw.get("theta"), channels)
    direction(gw.get("a_c"), channels)

    _validate_eikonal_controls(model, model_name)
    helicity_enabled = _validate_helicity_controls(model, model_name)
    _validate_form_factor_banks(model, channels, model_name)
    for name, definition in definitions.items():
        _validate_exchange_model(
            model,
            str(name),
            channels,
            model_name,
            definition,
            helicity_enabled,
        )
    forward_excitation_exchange = _validate_selections(model, definitions, model_name)
    for row in rows:
        if row.soft_exchange not in definitions:
            raise ValueError(f"PARAM_REGGE maps unknown PARAM_SOFT exchange {row.soft_exchange}")
        definition = definitions[row.soft_exchange]
        if definition["role"] != row.role:
            raise ValueError(f"PARAM_REGGE role for {row.soft_exchange} disagrees with PARAM_SOFT")
        exchange = exchanges[row.soft_exchange]
        if exchange.get("on") is not True:
            raise ValueError(f"PARAM_REGGE maps disabled PARAM_SOFT exchange {row.soft_exchange}")
    return Config(
        model_name=model_name,
        channels=channels,
        forward_excitation_exchange=forward_excitation_exchange,
        regge_exchanges=rows,
        definitions=definitions,
        model=model,
    )


# Select one unique mapped exchange by its physics role
def exchange(
    configuration: Config,
    role: str,
) -> Exchange:
    matches = [row for row in configuration.regge_exchanges if row.role == role]
    if len(matches) != 1:
        raise ValueError(f"PARAM_REGGE must map exactly one {role} exchange")
    return matches[0]


# Project one complete SOFT exchange vertex onto the physical proton
def proton_vertex(
    model: dict[str, object],
    exchange_name: str,
    model_name: str = "active",
) -> ProtonVertex:
    exchanges = model.get("EXCHANGE")
    if not isinstance(exchanges, dict) or exchange_name not in exchanges:
        raise ValueError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE must contain {exchange_name}")
    exchange = exchanges[exchange_name]
    if not isinstance(exchange, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE.{exchange_name} must be an object")
    context = f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE.{exchange_name}"
    matrix = _coupling_matrix(exchange.get("g"), None, f"{context}.g")
    channels = len(matrix)
    gw = model.get("GW")
    if not isinstance(gw, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.GW must be an object")
    proton = proton_state(gw.get("theta"), channels)
    transition = exchange.get("transition_ff")
    if transition not in TRANSITION_FORM_FACTORS:
        raise ValueError(f"{context}.transition_ff has an unknown mode")
    bank_name = exchange.get("ff")
    banks = model.get("FF")
    if not isinstance(bank_name, str) or not isinstance(banks, dict) or bank_name not in banks:
        raise ValueError(f"{context}.ff selects an unknown form factor bank")
    bank = banks[bank_name]
    if not isinstance(bank, dict):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.FF.{bank_name} must be an object")
    form_factor = bank.get("type")
    parameters = bank.get("param")
    if not isinstance(form_factor, str) or not isinstance(parameters, list):
        raise TypeError(f"PARAM_SOFT.MODEL.{model_name}.FF.{bank_name} has invalid fields")
    if len(parameters) != channels or any(not isinstance(row, list) for row in parameters):
        raise ValueError(f"PARAM_SOFT.MODEL.{model_name}.FF.{bank_name}.param must have N rows")
    slopes = [form_slope(form_factor, tuple(row)) for row in parameters]

    coupling = 0.0
    derivative = 0.0
    for row in range(channels):
        for column in range(channels):
            if transition == "diagonal" and row != column:
                continue
            projected = proton[row] * matrix[row][column] * proton[column]
            coupling += projected
            derivative += projected * (
                slopes[row] if row == column else 0.5 * (slopes[row] + slopes[column])
            )
    if not math.isfinite(coupling) or coupling <= 0.0:
        raise ValueError(f"PARAM_SOFT exchange {exchange_name} has no positive proton residue")
    slope = derivative / coupling
    if not math.isfinite(slope):
        raise ValueError(f"PARAM_SOFT exchange {exchange_name} has a nonfinite proton slope")
    return ProtonVertex(
        beam_residue_per_gev=coupling,
        amplitude_slope_per_gev2=slope,
        form_factor=form_factor,
        parameters=tuple(tuple(float(value) for value in row) for row in parameters),
        source=f"PARAM_SOFT.MODEL.{model_name}.EXCHANGE.{exchange_name}",
    )


# Compute only the projected physical beam coupling of one SOFT exchange
def coupling(model: dict[str, object], exchange_name: str) -> float:
    return proton_vertex(model, exchange_name).beam_residue_per_gev
