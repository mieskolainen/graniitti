# Regge continuum model parameters
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
import re

import pyjson5
from jsonpointer import JsonPointerException
from ray import tune

from core.io import jsonref, steering
from core.tune.drivers.graniitti.tunesetup.domains import (
    add_magnitude_phase_params,
    add_residue_params,
    angular_row_name,
    load_model_json,
    model_dir,
    settings,
    source_root,
    tp_con_names,
)
from core.tune.drivers.graniitti.tunesetup.domains import load_general_data as _load_general_data
from core.tune.parameters import tools


# Load one model-specific central continuum card
def _load_production_data(model: str) -> dict:
    model = str(model)
    if model not in {"MP", "XP", "GP", "TP"}:
        raise ValueError(f"Unknown continuum model {model}")
    return load_model_json(model_dir() / f"CON_{model}.json")


def _normalize_models(continuum_models) -> tuple[str, ...]:
    if isinstance(continuum_models, str):
        continuum_models = (continuum_models,)

    out = tuple(str(model) for model in continuum_models)
    if not out:
        raise ValueError("At least one continuum model must be selected")
    unknown = [model for model in out if model not in {"MP", "GP", "XP", "TP"}]
    if unknown:
        raise ValueError(f"Unknown continuum model(s): {unknown}")
    return out


# Compute true for elastic pp datasets which are handled by the eikonal setup
def _is_elastic_proton_pair(pair: tuple[int, int]) -> bool:
    return pair == (2212, 2212) or pair == (-2212, -2212)


def discover_final_states_from_datacards(datacards: list[dict]) -> list[tuple[int, int]]:
    """Discover two-body final-state PID pairs from the selected icetune datacards."""

    final_states = []
    seen = set()
    for datacard in datacards:
        data, _ = steering.load_dataset(datacard["datacard"], cdir=source_root())
        for dataset in data.get("sets", []):
            pid = dataset.get("pid", [])
            if len(pid) != 2:
                continue
            pair = (int(pid[0]), int(pid[1]))
            if _is_elastic_proton_pair(pair):
                continue
            key = tuple(sorted(pair))
            if key in seen:
                continue
            seen.add(key)
            final_states.append(pair)
    if not final_states:
        raise ValueError("No two-body PID final states discovered from datacards")
    return final_states


def _normalize_final_state_selection(
    final_state_pdgs: list[tuple[int | None, int | None] | list[int | None] | None] | None,
) -> list[tuple[int, int]] | None:
    if not final_state_pdgs:
        raise ValueError("Explicit continuum final-state selection must not be empty")

    final_states = []
    for pdg_pair in final_state_pdgs:
        if pdg_pair is None or tuple(pdg_pair) == (None, None):
            raise ValueError("Continuum final states require an explicit PDG pair")
        if len(pdg_pair) != 2:
            raise ValueError(f"Expected a two-body final-state PID pair, got {pdg_pair}")
        pair = (pdg_pair[0], pdg_pair[1])
        if pair[0] is None or pair[1] is None:
            raise ValueError("Continuum final states require two integer PDGs")
        final_states.append(pair)
    return final_states


# Compute magnitude bounds centered on one default coupling value
def _magnitude_bounds(mag: float, *, relative_range: float) -> tuple[float, float]:
    relative_range = float(relative_range)
    if not math.isfinite(mag) or mag < 0.0:
        raise ValueError(f"Continuum coupling magnitude must be non-negative, got {mag}")
    if not math.isfinite(relative_range) or relative_range < 0.0:
        raise ValueError(f"Continuum coupling relative range must be non-negative, got {relative_range}")
    return max(0.0, (1.0 - relative_range) * mag), (1.0 + relative_range) * mag


# Compute Cartesian bounds centered on the default continuum coupling
def _channel_cartesian_bounds(
    channel: list, *, relative_range: float
) -> tuple[tuple[float, float], tuple[float, float]]:
    re_val, im_val = tools.encode_polar(float(channel[2]), float(channel[3]))
    span = relative_range * float(channel[2])
    return ((re_val - span, re_val + span), (im_val - span, im_val + span))


# Add one normalized continuum coupling shape and its common coupling
def _add_projective_params(
    p: dict, *, model: str, exchange: str, pair_key: str, sector: str, rows: list, field: str, mode: str, relative_range: float
) -> None:
    if mode not in {"magnitude", "mag_phase_cartesian"}:
        raise ValueError("Projective continuum tuning requires magnitude or mag_phase_cartesian mode")
    if any(abs(math.sin(tools.canonicalize_phase(float(row[-1])))) > 1.0e-8 for row in rows):
        raise ValueError("Projective continuum rows must share a real signed phase convention")
    row_names = [angular_row_name(field, row) for row in rows]
    base = f"CON_{model}|{exchange}:{pair_key}/{sector}"
    weights = tools.continuum_completion_weights(
        f"{base}:{field}", [name.partition("(")[2][:-1] for name in row_names]
    ) or [1.0] * len(rows)
    norm = math.sqrt(sum(weight * float(row[-2])**2 for weight, row in zip(weights, rows, strict=True)))
    if norm <= tools.COMPLEX_EPS:
        raise ValueError("Projective continuum table must have nonzero completed norm")
    add_residue_params(
        p,
        rows=[f"{base}:{row}" for row in row_names],
        bounds=_magnitude_bounds(norm, relative_range=relative_range),
        mode=mode,
    )


# Compute an order-independent key for one explicit continuum PDG pair
def _canonical_continuum_pair(first: int, second: int, *, context: str) -> tuple[int, int]:
    if not isinstance(first, int) or isinstance(first, bool) or not isinstance(second, int) or isinstance(second, bool):
        raise ValueError(f"{context} must contain two integer PDGs")
    return min(first, second), max(first, second)


# Compute validated continuum entries and reject order-equivalent pair keys
def _continuum_entries(
    model_block: dict, *, context: str = "PARAM_CON", allow_default: bool = False
) -> list[tuple[str, tuple[int | None, int | None], list]]:
    if not isinstance(model_block, dict) or not model_block:
        raise ValueError(f"{context} must contain explicit entries")

    entries = []
    seen = {}
    for entry, pairs in model_block.items():
        if entry == "[*]":
            if not allow_default:
                raise ValueError(f"{context} does not support a [*] entry")
            entries.append((entry, (None, None), pairs))
            continue
        try:
            row = json.loads(entry)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{context} has invalid PDG pair key {entry}") from exc
        if not isinstance(row, list) or len(row) != 2:
            raise ValueError(f"{context} PDG pair key must contain two integers: {entry}")
        pair = (row[0], row[1])
        key = _canonical_continuum_pair(*pair, context=f"{context}.{entry}")
        if key in seen:
            raise ValueError(f"{context}.{entry} duplicates the order-equivalent pair key {seen[key]}")
        seen[key] = entry
        if not isinstance(pairs, list) or not pairs:
            raise ValueError(f"{context}.{entry} must contain a nonempty exchange pair list")
        if any(
            not isinstance(exchange_pair, list)
            or len(exchange_pair) != 2
            or any(not isinstance(value, int) or isinstance(value, bool) for value in exchange_pair)
            for exchange_pair in pairs
        ):
            raise ValueError(f"{context}.{entry} exchange pairs must be [exchange_top, exchange_bottom]")
        entries.append((entry, pair, pairs))
    return entries


def _continuum_entry_for_pair(
    model_block: dict, pdg_pair: tuple[int | None, int | None] | None, *, allow_default: bool = False
) -> tuple[str, tuple[int | None, int | None], list]:
    default = None
    for entry, entry_pdg_pair, pairs in _continuum_entries(model_block, allow_default=allow_default):
        if entry == "[*]":
            default = (entry, entry_pdg_pair, pairs)
            continue
        if pdg_pair is not None and (entry_pdg_pair == pdg_pair or (entry_pdg_pair[1], entry_pdg_pair[0]) == pdg_pair):
            return entry, entry_pdg_pair, pairs

    if default is not None:
        return default
    raise ValueError(f"No explicit PARAM_CON entry for final-state pair {pdg_pair}")


# Compute exchange ids used by one continuum model and final-state pair
def _continuum_exchange_ids(
    model_block: dict,
    pdg_pair: tuple[int | None, int | None] | None,
    *,
    tune_all_channels: bool,
    allow_default: bool = False,
) -> list[str]:
    _, pdg_pair, pairs = _continuum_entry_for_pair(model_block, pdg_pair, allow_default=allow_default)
    selected_pairs = pairs if tune_all_channels else pairs[:1]
    return list(dict.fromkeys(str(exchange) for pair in selected_pairs for exchange in pair))


# Iterate configured exchanges for each selected explicit continuum final state
def _model_channels(model: str, final_states: list, tune_all_channels: bool, *, allow_default: bool = False):
    general = _load_general_data()
    section = "PARAM_TENSORPOM" if model == "TP" else "PARAM_REGGE"
    model_block = general[section]["PARAM_CON"][model]
    for pair in final_states:
        if pair is not None:
            for exchange in _continuum_exchange_ids(
                model_block, pair, tune_all_channels=tune_all_channels, allow_default=allow_default
            ):
                yield exchange, pair


# Resolve one continuum exchange alias to its unique PARAM_REGGE trajectory row
def _regge_index_for_exchange(regge: dict, exchange: str) -> int:
    matches = [
        index for index, row in enumerate(regge["EXCHANGES"]) if int(exchange) in {int(alias) for alias in row["pdg"]}
    ]
    if len(matches) != 1:
        raise ValueError(f"Continuum exchange {exchange} must match exactly one PARAM_REGGE exchange row")
    return matches[0]


# Compute trajectory rows used by the selected continuum final states
def _continuum_regge_indices(
    *,
    continuum_models: tuple[str, ...],
    final_states: list[tuple[int, int]],
    tune_all_channels: bool,
) -> tuple[int, ...]:
    general = _load_general_data()
    regge_continuum = general["PARAM_REGGE"]["PARAM_CON"]
    tensor_continuum = general["PARAM_TENSORPOM"]["PARAM_CON"]
    regge = general["PARAM_REGGE"]
    indices = []
    for model in continuum_models:
        model_block = tensor_continuum["TP"] if model == "TP" else regge_continuum[model]
        for pdg_pair in final_states:
            for exchange in _continuum_exchange_ids(model_block, pdg_pair, tune_all_channels=tune_all_channels):
                indices.append(_regge_index_for_exchange(regge, exchange))
    return tuple(dict.fromkeys(indices))


# Add mapped SOFT trajectory parameters while keeping secondary exchanges fixed
def _add_regge_params(p: dict, indices: tuple[int, ...]) -> None:
    general = _load_general_data()
    regge = general["PARAM_REGGE"]
    model = str(general["PARAM_SOFT"]["active_model"])
    for index in indices:
        row = regge["EXCHANGES"][index]
        role = str(row["role"])
        if role not in {"pomeron", "odderon"}:
            continue
        exchange = str(row["soft_exchange"])
        for column, bounds in enumerate(settings()["bounds"]["eikonal"]["EXCHANGE"][role]["alpha"]):
            p[f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.alpha[{column}]"] = tune.uniform(*bounds)


# Compute the central continuum pair-key block matching one final state
def _production_entry_for_pair(
    production_block: dict, pdg_pair: tuple[int | None, int | None], *, card_name: str
) -> tuple[str, dict]:
    for key, value in production_block.items():
        if not isinstance(key, str) or not isinstance(value, dict):
            raise ValueError(f"Invalid {card_name} pair entry {key}")
        parts = key[1:-1].split(",") if key.startswith("[") and key.endswith("]") else []
        try:
            stored = tuple(int(part) for part in parts)
        except ValueError as exc:
            raise ValueError(f"Invalid {card_name} pair key {key}") from exc
        canonical = (
            len(stored) == 2
            and all(value > 0 for value in stored)
            and stored == tuple(sorted(stored))
            and key == f"[{stored[0]},{stored[1]}]"
        )
        if not canonical:
            raise ValueError(f"Invalid {card_name} pair key {key}")

    target = tuple(sorted(abs(int(value)) for value in pdg_pair))
    pair_key = f"[{target[0]},{target[1]}]"
    block = production_block.get(pair_key)
    if isinstance(block, dict):
        return pair_key, block
    raise ValueError(f"No {card_name} entry for final-state pair {pdg_pair}")


# Compute the sector selected by the physical signed continuum production legs
def _production_sector_for_pair(pair_block: dict, pdg_pair: tuple[int | None, int | None], *, card_name: str) -> str:
    if "self" in pair_block:
        if "same" in pair_block or "opposite" in pair_block:
            raise ValueError(f"{card_name} self sector cannot coexist with same or opposite sectors")
        return "self"

    first_positive = int(pdg_pair[0]) > 0
    second_positive = int(pdg_pair[1]) > 0
    sector = "same" if first_positive == second_positive else "opposite"
    if sector not in pair_block:
        raise ValueError(f"No {sector} sign {card_name} sector for physical final-state pair {pdg_pair}")
    return sector


# Compute one model-specific central continuum block and charge sector
def _production_block_for_pair(
    production: dict, exchange: str, pdg_pair: tuple[int | None, int | None], *, card_name: str
) -> tuple[str, str, dict]:
    if exchange not in production:
        raise KeyError(f"{card_name} is missing active continuum exchange {exchange}")
    pair_key, pair_block = _production_entry_for_pair(production[exchange], pdg_pair, card_name=card_name)
    sector = _production_sector_for_pair(pair_block, pdg_pair, card_name=card_name)
    return pair_key, sector, pair_block[sector]


# Add one active central continuum amplitude table
def _add_continuum_production_block_params(
    p: dict,
    *,
    continuum_model: str,
    exchange: str,
    pair_key: str,
    sector: str,
    block: dict,
    mode: str,
    relative_range: float,
    equal_gp_m: bool,
    geometry: str,
    active_only: bool,
    projections: list[int] | None,
) -> None:
    basis = block.get("basis")
    allowed_bases = {
        "MP": {
            "crossed_auto_min_L",
            "crossed_auto_min_S",
            "crossed_auto_equal_ls",
            "crossed_auto_equal_helicity",
            "crossed_ls",
            "crossed_helicity",
        },
        "XP": {
            "crossed_auto_min_L",
            "crossed_auto_min_S",
            "crossed_auto_equal_ls",
            "crossed_auto_equal_helicity",
            "crossed_ls",
            "crossed_helicity",
        },
        "GP": {"crossed_ls", "crossed_helicity"},
    }
    if basis not in allowed_bases[continuum_model]:
        allowed = sorted(allowed_bases[continuum_model])
        raise ValueError(f"{continuum_model} production block {exchange}:{pair_key}/{sector} must use basis {allowed}")
    if basis in {"crossed_auto_min_L", "crossed_auto_min_S", "crossed_auto_equal_ls", "crossed_auto_equal_helicity"}:
        g = block.get("g")
        if not isinstance(g, list) or len(g) != 2:
            raise ValueError(f"{continuum_model} production block {exchange}:{pair_key}/{sector} has no g")
        magnitude_bounds = _magnitude_bounds(float(g[0]), relative_range=relative_range)
        re_bounds, im_bounds = _channel_cartesian_bounds([None, None, *g], relative_range=relative_range)
        base_key = f"CON_{continuum_model}|{exchange}:{pair_key}/{sector}:g"
        add_magnitude_phase_params(
            p,
            base_key=base_key,
            magnitude_key=f"{base_key}[0]",
            phase_key=f"{base_key}[1]",
            mode=mode,
            magnitude_bounds=magnitude_bounds,
            cartesian_bounds=(re_bounds, im_bounds),
            mode_name="continuum_production_mode",
        )
        return
    field = "g_ls" if basis == "crossed_ls" else "helicity"
    rows = block.get(field)
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"{continuum_model} production block {exchange}:{pair_key}/{sector} has no {field} rows")
    if block.get("CP") != [True, True]:
        if continuum_model == "GP":
            raise ValueError(
                f"GP production block {exchange}:{pair_key}/{sector} must declare "
                "CP = [true, true]. C parity is applied by exchange antiparticle signs "
                "and the continuum t/u projector"
            )
        raise ValueError(
            f"{continuum_model} production block {exchange}:{pair_key}/{sector} must enforce C and parity symmetry"
        )

    nonphoton_gp = continuum_model == "GP" and int(exchange) != 22
    row_size = 5 if nonphoton_gp else 4
    if any(not isinstance(row, list) or len(row) != row_size for row in rows):
        raise ValueError(
            f"{continuum_model} production block {exchange}:{pair_key}/{sector} {field} rows must contain {row_size} columns"
        )
    if projections is not None and nonphoton_gp:
        rows = [row for row in rows if abs(row[2]) in projections]
    mode = tools.validate_tune_mode(mode, name="continuum_production_mode")
    if mode == "none":
        return
    magnitude_column = row_size - 2
    phase_column = magnitude_column + 1
    magnitudes = [float(row[magnitude_column]) for row in rows]
    if any(value < 0.0 or not math.isfinite(value) for value in magnitudes):
        raise ValueError("Continuum coupling magnitudes must be finite and non-negative")
    if geometry == "projective":
        if equal_gp_m:
            raise ValueError("Projective GP continuum geometry is incompatible with equal_gp_m")
        _add_projective_params(
            p,
            model=continuum_model,
            exchange=exchange,
            pair_key=pair_key,
            sector=sector,
            rows=[row for row in rows if not active_only or float(row[magnitude_column]) > 0.0],
            field=field,
            mode=mode,
            relative_range=relative_range,
        )
        return
    if equal_gp_m and nonphoton_gp:
        if mode != "magnitude":
            raise ValueError("Equal GP m coupling tuning requires magnitude mode")
        groups = {}
        for row in rows:
            if float(row[magnitude_column]) <= 0.0:
                continue
            groups.setdefault((row[0], row[1]), []).append(row)
        for group in groups.values():
            if sorted(int(row[2]) for row in group) != [-2, -1, 0]:
                raise ValueError("Equal GP m coupling tuning requires m = -2,-1,0 rows")
            if not all(
                math.isclose(
                    float(row[magnitude_column]), float(group[0][magnitude_column]), rel_tol=0.0, abs_tol=1.0e-12
                )
                for row in group[1:]
            ):
                raise ValueError("Equal GP m coupling tuning requires equal card magnitudes")
            magnitude_bounds = _magnitude_bounds(float(group[0][magnitude_column]), relative_range=relative_range)
            first = angular_row_name(field, group[0]).rsplit(",", maxsplit=1)[0]
            row_key = f"CON_{continuum_model}|{exchange}:{pair_key}/{sector}:{first},m=-2,-1,0)"
            p[tools.magnitude_key(row_key)] = tune.uniform(*magnitude_bounds)
        return
    for row_index, row in enumerate(rows):
        if not magnitudes[row_index] > 0.0:
            continue
        magnitude_bounds = _magnitude_bounds(float(row[magnitude_column]), relative_range=relative_range)
        re_bounds, im_bounds = _channel_cartesian_bounds(
            [None, None, row[magnitude_column], row[phase_column]], relative_range=relative_range
        )
        row_key = f"CON_{continuum_model}|{exchange}:{pair_key}/{sector}:{angular_row_name(field, row)}"
        add_magnitude_phase_params(
            p,
            base_key=row_key,
            magnitude_key=tools.magnitude_key(row_key),
            phase_key=row_key,
            mode=mode,
            magnitude_bounds=magnitude_bounds,
            cartesian_bounds=(re_bounds, im_bounds),
            mode_name="continuum_production_mode",
        )


# Compute a relative interval for one real tensor continuum coupling
def _tensor_coupling_bounds(value: float, *, relative_range: float) -> tuple[float, float]:
    value = float(value)
    relative_range = float(relative_range)
    if not math.isfinite(value):
        raise ValueError(f"Tensor continuum coupling must be finite, got {value}")
    if not math.isfinite(relative_range) or relative_range < 0.0:
        raise ValueError(f"Continuum coupling relative range must be finite and non-negative, got {relative_range}")
    bounds = ((1.0 - relative_range) * value, (1.0 + relative_range) * value)
    return min(bounds), max(bounds)


# Add active real tensor continuum couplings from CON_TP.json
def _add_tensor_continuum_params(p: dict, *, block: dict, prefix: str, relative_range: float) -> None:
    couplings = block.get("g_tensor")
    if not isinstance(couplings, list) or not couplings:
        raise ValueError(f"TP production block {prefix} has no g_tensor couplings")
    for name, coupling in zip(tp_con_names(couplings), couplings, strict=True):
        lower, upper = _tensor_coupling_bounds(coupling, relative_range=relative_range)
        if lower < upper:
            p[f"{prefix}:{name}"] = tune.uniform(lower, upper)


# Add active model-specific central continuum amplitude parameters
def _add_continuum_production_params(
    p: dict,
    *,
    continuum_models: tuple[str, ...],
    final_states: list[tuple[int, int]],
    tune_all_channels: bool,
    mode: str,
    relative_range: float,
    equal_gp_m: bool,
    geometry: str,
    active_only: bool,
    projections: list[int] | None,
) -> None:
    mode = tools.validate_tune_mode(mode, name="continuum_production_mode")
    seen = set()
    for model in continuum_models:
        production = _load_production_data(model)
        card_name = f"CON_{model}.json"
        for exchange, pair in _model_channels(model, final_states, tune_all_channels, allow_default=model == "TP"):
            if exchange not in production:
                raise KeyError(f"{card_name} is missing active continuum exchange {exchange}")
            try:
                pair_key, block = _production_entry_for_pair(production[exchange], pair, card_name=card_name)
            except ValueError as exc:
                if model == "TP" or not str(exc).startswith(f"No {card_name} entry"):
                    raise
                raise ValueError(f"{card_name} requires an explicit subvertex for {pair}") from exc
            sector = None if model == "TP" else _production_sector_for_pair(block, pair, card_name=card_name)
            key = (model, exchange, pair_key, sector)
            if key in seen:
                continue
            seen.add(key)
            if model == "TP":
                _add_tensor_continuum_params(
                    p, block=block, prefix=f"CON_TP|{exchange}:{pair_key}", relative_range=relative_range
                )
            else:
                _add_continuum_production_block_params(
                    p,
                    continuum_model=model,
                    exchange=exchange,
                    pair_key=pair_key,
                    sector=sector,
                    block=block[sector],
                    mode=mode,
                    relative_range=relative_range,
                    equal_gp_m=equal_gp_m,
                    geometry=geometry,
                    active_only=active_only,
                    projections=projections,
                )


# Add one propagator-specific continuum off-shell form factor
def form_factor(
    *, model: str, exchange: str, pdg_pair: tuple[int | None, int | None], ff_type: str | None = None, kernels: int = 1
):
    p = {}
    production = _load_production_data(model)
    pair_key, pair_block = _production_entry_for_pair(
        production[str(exchange)], pdg_pair, card_name=f"CON_{model}.json"
    )
    tag = f"CON_{model}|{exchange}:{pair_key}:"
    if ff_type is None:
        ff_type = str(pair_block["FF_offshell"]["type"])
    bounds = settings()["bounds"]["continuum"]["FF_offshell"]
    if ff_type in bounds and ff_type != "gkernel":
        for field, limits in bounds[ff_type].items():
            p[f"{tag}FF_offshell.{field}"] = tune.uniform(*limits)
    elif ff_type == "gkernel":
        if kernels < 1:
            raise Exception(__name__ + ".form_factor: gkernel requires kernels >= 1")
        for i in range(kernels):
            term_bounds = bounds["gkernel"]
            for field, limits in term_bounds.items():
                p[f"{tag}FF_offshell.terms[{i}].{field}"] = tune.uniform(*limits)
    else:
        raise Exception(__name__ + f".form_factor: Unknown ff_type = {ff_type}")

    aux = {f"{tag}FF_offshell.type": ff_type, f"{tag}FF_offshell.norm": "pole"}

    return p, aux


# Add Poisson zero-secondary veto parameters for one exchange and hadron pair
def poisson_veto(*, model: str, exchange: str, pdg_pair: tuple[int, int], active=False):
    if not active:
        return {}, {}
    production = _load_production_data(model)
    pair, block = _production_entry_for_pair(production[str(exchange)], pdg_pair, card_name=f"CON_{model}.json")
    default_m0 = float(block["pveto"]["M0"])
    tag = f"CON_{model}|{exchange}:{pair}:pveto"
    p = {f"{tag}.M0": tune.uniform(*(scale * default_m0 for scale in settings()["icetune"]["continuum"]["pveto"]["M0"]["relative_bounds"])), f"{tag}.c": tune.uniform(*settings()["bounds"]["continuum"]["pveto"]["c"])}
    return p, {f"{tag}.active": True}


# Add inverse transfer scales for active continuum exchange and hadron vertices
def transfer_form_factor(
    *,
    continuum_models: tuple[str, ...],
    final_states: list[tuple[int, int]],
    tune_all_channels: bool,
) -> dict:
    p = {}
    for model in continuum_models:
        production = _load_production_data(model)
        card_name = f"CON_{model}.json"
        selected = [pair for pair in final_states
                    if model != "TP" or sorted(abs(int(value)) for value in pair) != [2212, 2212]]
        for exchange, pair in _model_channels(model, selected, tune_all_channels, allow_default=model == "TP"):
            pair_key, block = _production_entry_for_pair(production[exchange], pair, card_name=card_name)
            form = block["FF_transfer"]
            if form["type"] in {"none", "dirac"}:
                continue
            if form["type"] != "power":
                raise ValueError(f"{card_name} {exchange}:{pair_key} transfer tuning requires a power form")
            p[f"CON_{model}|{exchange}:{pair_key}:FF_transfer.LambdaInv2"] = tune.uniform(*settings()["icetune"]["FF_transfer"]["LambdaInv2"])
    return p


# Keep one optimizer coordinate for each referenced source parameter
def _shared_parameters(parameters):
    reader = jsonref.JsonReader(pyjson5.load)
    output = {}
    seen = set()
    for key, value in parameters.items():
        card, _, suffix = key.partition("|")
        if not card.startswith("CON_"):
            output[key] = value
            continue
        exchange, pair, field = suffix.split(":", 2)
        if field.split(".", 1)[0] not in {"FF_transfer", "FF_offshell", "pveto", "reggeize"}:
            output[key] = value
            continue
        tokens = [exchange, pair, *re.findall(r"[^.\[\]]+", field)]
        if tokens[-1] == "LambdaInv2":
            tokens[-1] = "Lambda2"
        remaining = []
        while True:
            try:
                target, parts = reader.origin(model_dir() / f"{card}.json", tokens)
                break
            except JsonPointerException:
                # A selected form family can introduce fields absent from the initial card
                if not tokens:
                    raise
                remaining.insert(0, tokens.pop())
        source = target, tuple(parts + remaining)
        if source not in seen:
            output[key] = value
            seen.add(source)
    return output


# Construct continuum trajectory, coupling and form-factor tuning parameters
def setup(
    *, final_state_pdgs,
    REGGE: bool = False,
    FF_offshell: bool = True,
    pveto: bool = False,
    FF_transfer: bool = False,
    ff_type: str | None = None,
    continuum_models="MP",
    tune_couplings: bool = True,
    tune_all_channels: bool = True,
    production_mode: str = "magnitude",
    production_geometry: str = "direct",
    production_active_only: bool = True,
    production_m: list[int] | None = None,
    production_relative_range: float | None = None,
    equal_gp_m: bool = False,
    tune_reggeize: bool = False,
    freeze_scale2_range: tuple[float, float] | None = None,
):
    p = {}
    p_aux = {}
    production_relative_range = (settings()["icetune"]["continuum"]["production"]["relative_range"]
                                 if production_relative_range is None else production_relative_range)
    continuum_models = _normalize_models(continuum_models)
    if production_m is not None and (continuum_models != ("GP",) or not isinstance(production_m, list)
            or not production_m or any(type(m) is not int or m < 0 for m in production_m)
            or len(set(production_m)) != len(production_m)):
        raise ValueError("production_m requires distinct nonnegative integer |m| sectors of GP")
    if pveto and all(model == "TP" for model in continuum_models):
        raise ValueError("Continuum Poisson veto tuning is not supported by TP")
    production_geometry = str(production_geometry)
    if production_geometry not in {"direct", "projective"}:
        raise ValueError(f'Unknown continuum production geometry "{production_geometry}"')
    if not production_active_only and production_geometry != "projective":
        raise ValueError("Tuning zero continuum rows requires projective production geometry")
    final_states = _normalize_final_state_selection(final_state_pdgs)

    # Tune active Pomeron or odderon trajectories while secondary rows stay fixed
    if REGGE:
        indices = _continuum_regge_indices(
            continuum_models=continuum_models, final_states=final_states, tune_all_channels=tune_all_channels
        )
        _add_regge_params(p, indices)

    if tune_couplings:
        _add_continuum_production_params(
            p,
            continuum_models=continuum_models,
            final_states=final_states,
            tune_all_channels=tune_all_channels,
            mode=production_mode,
            relative_range=production_relative_range,
            equal_gp_m=equal_gp_m,
            geometry=production_geometry,
            active_only=production_active_only,
            projections=production_m,
        )

    if FF_transfer:
        p.update(
            transfer_form_factor(
                continuum_models=continuum_models, final_states=final_states, tune_all_channels=tune_all_channels
            )
        )

    if FF_offshell:
        form_entries_seen = set()
        for model in continuum_models:
            for exchange, pair in _model_channels(model, final_states, tune_all_channels, allow_default=model == "TP"):
                key = (model, exchange, pair)
                if key in form_entries_seen:
                    continue
                t, t_aux = form_factor(model=model, exchange=exchange, pdg_pair=pair, ff_type=ff_type)
                p.update(t)
                p_aux.update(t_aux)
                form_entries_seen.add(key)

    if freeze_scale2_range is not None:
        low, high = freeze_scale2_range
        if not (0 < low < high and math.isfinite(high)):
            raise ValueError("freeze_scale2_range requires finite 0 < low < high in GeV^2")
    if tune_reggeize or pveto or freeze_scale2_range is not None:
        for model in continuum_models:
            if model == "TP":
                continue
            production = _load_production_data(model)
            for exchange, pdg_pair in _model_channels(model, final_states, tune_all_channels):
                pair, _ = _production_entry_for_pair(production[exchange], pdg_pair, card_name=f"CON_{model}.json")
                if tune_reggeize:
                    p[f"CON_{model}|{exchange}:{pair}:reggeize.active"] = tune.choice([False, True])
                if freeze_scale2_range is not None:
                    p[f"CON_{model}|{exchange}:{pair}:reggeize.freeze_scale2"] = tune.uniform(*freeze_scale2_range)
                t, t_aux = poisson_veto(model=model, exchange=exchange, pdg_pair=pdg_pair, active=pveto)
                p.update(t)
                p_aux.update(t_aux)

    return _shared_parameters(p), _shared_parameters(p_aux)
