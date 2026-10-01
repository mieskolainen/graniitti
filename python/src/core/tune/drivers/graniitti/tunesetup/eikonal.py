# Eikonal Pomeron (elastic / screening) parameters
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from ray import tune

from core.numerics import array
from core.tune.drivers.graniitti.tunesetup.domains import load_general_data as _load_general_data
from core.tune.drivers.graniitti.tunesetup.domains import settings

# Canonical eikonal optimizer coordinate transforms
DESCENDING_COUPLING_SUFFIX = "@DESCENDING"
ORDERED_DPOW_SUFFIX = "@ORDERED_DPOW="
SYMMETRIC_COUPLING_SUFFIX = "@SYMMETRIC"

# Mark one physical coupling target as a descending optimizer coordinate
def descending_coupling_key(raw_key: str) -> str:
    return f"{raw_key}{DESCENDING_COUPLING_SUFFIX}"


# Compute true when a key is a descending coupling coordinate
def is_descending_coupling_key(key: str) -> bool:
    return str(key).endswith(DESCENDING_COUPLING_SUFFIX)


# Compute the physical coupling target of a descending coordinate
def descending_coupling_target(key: str) -> str:
    key = str(key)
    if not is_descending_coupling_key(key):
        raise ValueError(f'descending_coupling_target: invalid key "{key}"')
    return key[: -len(DESCENDING_COUPLING_SUFFIX)]


# Mark one upper-triangle coupling target as a symmetric matrix coordinate
def symmetric_coupling_key(raw_key: str) -> str:
    return f"{raw_key}{SYMMETRIC_COUPLING_SUFFIX}"


# Compute true when a key is a symmetric matrix coupling coordinate
def is_symmetric_coupling_key(key: str) -> bool:
    return str(key).endswith(SYMMETRIC_COUPLING_SUFFIX)


# Compute the upper-triangle target of a symmetric coupling coordinate
def symmetric_coupling_target(key: str) -> str:
    key = str(key)
    if not is_symmetric_coupling_key(key):
        raise ValueError(f'symmetric_coupling_target: invalid key "{key}"')
    return key[: -len(SYMMETRIC_COUPLING_SUFFIX)]


# Decode one leading coupling and successive ratios into descending couplings
def decode_descending_couplings(coordinates: list[float]) -> list[float]:
    values = [array.asarray(value, dtype=float) for value in coordinates]
    if not values or not all(array.namespace(value).isfinite(value) for value in values):
        raise ValueError("decode_descending_couplings: expected finite coordinates")
    if values[0] < 0.0 or any(value < 0.0 or value > 1.0 for value in values[1:]):
        raise ValueError("decode_descending_couplings: invalid leading value or ratio")

    couplings = [values[0]]
    for ratio in values[1:]:
        couplings.append(couplings[-1] * ratio)
    return couplings


# Encode descending nonnegative couplings as one leading value and ratios
def encode_descending_couplings(couplings: list[float]) -> list[float]:
    values = [float(value) for value in couplings]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("encode_descending_couplings: expected finite couplings")
    if values[-1] < 0.0 or any(left < right for left, right in zip(values[:-1], values[1:], strict=True)):
        raise ValueError("encode_descending_couplings: couplings are not descending")

    coordinates = [values[0]]
    for previous, value in zip(values[:-1], values[1:], strict=True):
        coordinates.append(0.0 if previous == 0.0 else value / previous)
    return coordinates


# Mark one physical DPOW scale target as an ordered optimizer coordinate
def ordered_dpow_key(raw_key: str, upper: float) -> str:
    return f"{raw_key}{ORDERED_DPOW_SUFFIX}{upper:.17g}"


# Compute true when a key is an ordered DPOW scale coordinate
def is_ordered_dpow_key(key: str) -> bool:
    return ORDERED_DPOW_SUFFIX in str(key)


# Compute the physical DPOW target of an ordered coordinate
def ordered_dpow_target(key: str) -> str:
    key = str(key)
    if not is_ordered_dpow_key(key):
        raise ValueError(f'ordered_dpow_target: invalid key "{key}"')
    return key.rsplit(ORDERED_DPOW_SUFFIX, 1)[0]


# Decode a lower scale and unit-interval fraction into an ordered scale pair
def decode_ordered_dpow_scales(lower: float, fraction: float, *, limit: float) -> tuple[float, float]:
    lower = array.asarray(lower, dtype=float)
    fraction = array.asarray(fraction, dtype=float)
    if not array.namespace(lower).isfinite(lower) or not math.isfinite(limit) or lower < 0.0 or lower > limit or limit <= 0.0:
        raise ValueError("decode_ordered_dpow_scales: lower scale outside bounds")
    if not array.namespace(fraction).isfinite(fraction) or fraction < 0.0 or fraction > 1.0:
        raise ValueError("decode_ordered_dpow_scales: fraction outside [0, 1]")
    upper = lower + fraction * (limit - lower)
    return lower, upper


# Encode one bounded ascending DPOW scale pair into canonical coordinates
def encode_ordered_dpow_scales(lower: float, upper: float, *, limit: float) -> tuple[float, float]:
    lower = float(lower)
    upper = float(upper)
    if not all(math.isfinite(value) for value in (lower, upper)):
        raise ValueError("encode_ordered_dpow_scales: expected finite scales")
    if not math.isfinite(limit) or limit <= 0.0 or lower < 0.0 or upper > limit or lower > upper:
        raise ValueError("encode_ordered_dpow_scales: scales are not bounded and ascending")
    width = limit - lower
    fraction = 0.0 if width == 0.0 else (upper - lower) / width
    return lower, fraction


# Encode a symmetric pair through its single independent coupling
def encode_symmetric_coupling(values: list[float]) -> tuple[float]:
    value, mirror = values
    if not math.isclose(value, mirror, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("Asymmetric pair in eikonal coupling defaults")
    return (value,)


# Compute the requested or source-card eikonal model name
def source_model(model: str | None = None) -> str:
    if model is not None:
        return model
    return str(_load_general_data()["PARAM_SOFT"]["active_model"])


# Compute the source-card block for one eikonal model
def source_model_block(model: str) -> dict:
    soft = _load_general_data()["PARAM_SOFT"]
    models = soft["MODEL"]
    if model not in models:
        raise Exception(__name__ + f".source_model_block: Unknown model = {model}")
    return models[model]


# Compute true when the source card enables the Odderon
def source_odderon_active(model: str) -> bool:
    exchange = source_model_block(model)["EXCHANGE"]
    if "O" not in exchange:
        return False
    enabled = exchange["O"]["on"]
    if not isinstance(enabled, bool):
        raise TypeError(__name__ + ".source_odderon_active: EXCHANGE.O.on must be boolean")
    return enabled


# Compute true when the source card enables hard three-gluon exchange
def source_o3g_active(model: str) -> bool:
    exchange = source_model_block(model)["EXCHANGE"]
    if "O3g" not in exchange:
        return False
    enabled = exchange["O3g"]["on"]
    if not isinstance(enabled, bool):
        raise TypeError(__name__ + ".source_o3g_active: EXCHANGE.O3g.on must be boolean")
    return enabled


# Compute enabled secondary Reggeon exchanges from the source card
def source_reggeon_exchanges(model: str) -> tuple[str, ...]:
    soft = _load_general_data()["PARAM_SOFT"]
    definitions = soft["EXCHANGE_DEF"]
    exchanges = source_model_block(model)["EXCHANGE"]
    active = []
    for name, entry in exchanges.items():
        if str(definitions[name]["role"]) != "reggeon":
            continue
        enabled = entry["on"]
        if not isinstance(enabled, bool):
            raise TypeError(__name__ + f".source_reggeon_exchanges: EXCHANGE.{name}.on must be boolean")
        if enabled:
            active.append(name)
    return tuple(active)


# Compute true when the source card selects q-exponential unitarization
def source_q_exp_enabled(model: str) -> bool:
    eikonal = source_model_block(model)["EIKONAL"]
    return str(eikonal["unitarization"]) == "q_exp"


# Compute the source-card form-factor bank for one requested soft exchange
def source_form_factor_bank(model: str, exchange: str, value: str | None = None) -> str:
    block = source_model_block(model)
    selected = str(block["EXCHANGE"][exchange]["ff"])
    if value is None or str(block["FF"][selected]["type"]) == value:
        return selected
    matches = [name for name, entry in block["FF"].items() if str(entry["type"]) == value]
    if not matches:
        raise Exception(__name__ + f".source_form_factor_bank: Unknown FF = {value}")
    return matches[0]


# Compute the single form factor bank shared by selected exchanges
def source_common_form_factor_bank(model: str, exchanges: tuple[str, ...]) -> str:
    if not exchanges:
        raise ValueError(__name__ + ".source_common_form_factor_bank: no exchanges selected")
    block = source_model_block(model)
    banks = {str(block["EXCHANGE"][exchange]["ff"]) for exchange in exchanges}
    if len(banks) != 1:
        raise ValueError(__name__ + ".source_common_form_factor_bank: selected exchanges do not share one form factor")
    return banks.pop()


# Compute the explicit forward-excitation Pomeron selected by one model
def source_forward_exchange(model: str) -> str:
    general = _load_general_data()
    block = source_model_block(model)
    selected_exchanges = block["EIKONAL"]["excitation_exchanges"]
    if selected_exchanges == ["*"]:
        definitions = general["PARAM_SOFT"]["EXCHANGE_DEF"]
        selected_exchanges = [
            name
            for name, exchange in block["EXCHANGE"].items()
            if exchange["on"] and definitions[name]["role"] == "pomeron"
        ]
    if not isinstance(selected_exchanges, list) or len(selected_exchanges) != 1:
        raise ValueError(__name__ + ".source_forward_exchange: excitation_exchanges must resolve to one Pomeron")
    selected = str(selected_exchanges[0])
    if selected not in block["EXCHANGE"]:
        raise Exception(__name__ + f".source_forward_exchange: Unknown exchange = {selected}")
    definition = general["PARAM_SOFT"]["EXCHANGE_DEF"][selected]
    if str(definition["role"]) != "pomeron":
        raise Exception(__name__ + f".source_forward_exchange: Exchange {selected} is not a Pomeron")
    enabled = block["EXCHANGE"][selected]["on"]
    if not isinstance(enabled, bool):
        raise TypeError(__name__ + f".source_forward_exchange: EXCHANGE.{selected}.on must be boolean")
    if not enabled:
        raise Exception(__name__ + f".source_forward_exchange: Exchange {selected} is disabled")
    return selected


# Compute the Good-Walker channel count from the selected Pomeron matrix
def source_channel_count(model: str) -> int:
    selected = source_forward_exchange(model)
    matrix = source_model_block(model)["EXCHANGE"][selected]["g"]
    if not isinstance(matrix, list) or not matrix:
        raise TypeError(__name__ + f".source_channel_count: EXCHANGE.{selected}.g must be a matrix")
    channels = len(matrix)
    if any(not isinstance(row, list) or len(row) != channels for row in matrix):
        raise ValueError(__name__ + f".source_channel_count: EXCHANGE.{selected}.g must be square")
    return channels


# Compute a positive diagonal coupling interval around the source tune value
def source_diagonal_coupling_bounds(model: str, exchange: str, row: int) -> tuple[float, float]:
    matrix = source_model_block(model)["EXCHANGE"][exchange]["g"]
    try:
        source_value = float(matrix[row][row])
    except (IndexError, TypeError, ValueError) as exc:
        raise ValueError(
            __name__ + f".source_diagonal_coupling_bounds: invalid EXCHANGE.{exchange}.g[{row},{row}]"
        ) from exc
    if not math.isfinite(source_value) or source_value <= 0.0:
        raise ValueError(
            __name__ + f".source_diagonal_coupling_bounds: EXCHANGE.{exchange}.g[{row},{row}] must be positive"
        )
    lower_scale, upper_scale = settings()["icetune"]["eikonal"]["EXCHANGE"]["g"]["diagonal_relative_bounds"]
    return lower_scale * source_value, upper_scale * source_value


# Dispatch one configured eikonal form-factor family
def _form_factor_parameters(name: str, *, model: str, bank: str, allowed: set[str], kernels: int = 1) -> dict:
    if name not in allowed:
        raise Exception(__name__ + f".setup: Unknown FF = {name}")
    if kernels < 1:
        raise Exception(__name__ + f".{name}: kernels must be >= 1")
    bounds = settings()["bounds"]["eikonal"]["FF"][name]["param"]
    kernel_count = kernels if name == "GKERNEL" else 1
    parameters = {}
    for row in range(source_channel_count(model)):
        for kernel in range(kernel_count):
            for column, limits in enumerate(bounds):
                raw_key = f"SOFT|MODEL.{model}:FF.{bank}.param[{row},{len(bounds) * kernel + column}]"
                key = ordered_dpow_key(raw_key, bounds[1][1]) if name == "DPOW" and column < 2 else raw_key
                domain = (0.0, 1.0) if name == "DPOW" and column == 1 else limits
                parameters[key] = tune.uniform(*domain)
    return parameters


# Add secondary Reggeon trajectory and Good Walker coupling tune parameters
def _add_reggeon_parameters(parameters: dict, *, model: str, exchanges: tuple[str, ...], rows: int) -> None:
    for exchange in exchanges:
        parameters.update(
            {
                f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.alpha[0]": tune.uniform(*settings()["bounds"]["eikonal"]["EXCHANGE"]["reggeon"]["alpha"][0]),
                f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.alpha[1]": tune.uniform(*settings()["bounds"]["eikonal"]["EXCHANGE"]["reggeon"]["alpha"][1]),
            }
        )
        for row in range(rows):
            limits = source_diagonal_coupling_bounds(model, exchange, row)
            parameters[f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.g[{row},{row}]"] = tune.uniform(*limits)
            for column in range(row + 1, rows):
                raw_key = f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.g[{row},{column}]"
                parameters[symmetric_coupling_key(raw_key)] = tune.uniform(*settings()["bounds"]["eikonal"]["EXCHANGE"]["g[i,j]"])


# Add trajectory and Good Walker coupling tune parameters
def general(
    model: str = "single", tune_odderon: bool | None = None, tune_o3g: bool | None = None, tune_reggeons: bool = False
):
    exchange = source_model_block(model)["EXCHANGE"]
    pomeron = source_forward_exchange(model)
    o_index = "O" if "O" in exchange else None
    o3g_index = "O3g" if "O3g" in exchange else None
    if tune_odderon is None:
        tune_odderon = source_odderon_active(model)
    if tune_o3g is None:
        tune_o3g = source_o3g_active(model)
    reggeons = source_reggeon_exchanges(model) if tune_reggeons else ()
    if tune_reggeons and not reggeons:
        raise Exception(__name__ + ".general: Reggeon tuning requires an enabled exchange")

    p = {
        f"SOFT|MODEL.{model}:EXCHANGE.{pomeron}.alpha[0]": tune.uniform(*settings()["bounds"]["eikonal"]["EXCHANGE"]["pomeron"]["alpha"][0]),
        f"SOFT|MODEL.{model}:EXCHANGE.{pomeron}.alpha[1]": tune.uniform(*settings()["bounds"]["eikonal"]["EXCHANGE"]["pomeron"]["alpha"][1]),
        f"SOFT|MODEL.{model}:pion_loop_scale2": tune.uniform(*settings()["bounds"]["eikonal"]["pion_loop_scale2"]),
    }
    if source_q_exp_enabled(model):
        p[f"SOFT|MODEL.{model}:EIKONAL.q"] = tune.uniform(*settings()["bounds"]["eikonal"]["EIKONAL"]["q"])

    secondary = []
    for name, enabled, index in (("O", tune_odderon, o_index), ("O3g", tune_o3g, o3g_index)):
        if not enabled:
            continue
        if index is None:
            label = "Odderon" if name == "O" else name
            raise Exception(__name__ + f".general: {label} tuning requires exchange {name}")
        secondary.append(name)
        for column, bounds in enumerate(settings()["bounds"]["eikonal"]["EXCHANGE"]["odderon"]["alpha"]):
            p[f"SOFT|MODEL.{model}:EXCHANGE.{name}.alpha[{column}]"] = tune.uniform(*bounds)

    rows = source_channel_count(model)
    _add_reggeon_parameters(parameters=p, model=model, exchanges=reggeons, rows=rows)
    # Tune only rotations that control the physical proton row
    theta_count = rows - 1
    for index in range(theta_count):
        p[f"SOFT|MODEL.{model}:GW.theta[{index}]"] = tune.uniform(*settings()["bounds"]["eikonal"]["GW"]["theta"])
    for index in range(rows):
        raw_key = f"SOFT|MODEL.{model}:EXCHANGE.{pomeron}.g[{index},{index}]"
        key = descending_coupling_key(raw_key)
        p[key] = tune.uniform(*settings()["bounds"]["eikonal"]["EXCHANGE"]["pomeron"]["g[0,0]"]) if index == 0 else tune.uniform(*settings()["icetune"]["eikonal"]["EXCHANGE"]["g"]["pomeron_diagonal_ratio_bounds"])
        for name in secondary:
            limits = source_diagonal_coupling_bounds(model, name, index)
            p[f"SOFT|MODEL.{model}:EXCHANGE.{name}.g[{index},{index}]"] = tune.uniform(*limits)

    for row in range(rows):
        for column in range(row + 1, rows):
            for name in secondary:
                raw_key = f"SOFT|MODEL.{model}:EXCHANGE.{name}.g[{row},{column}]"
                p[symmetric_coupling_key(raw_key)] = tune.uniform(*settings()["bounds"]["eikonal"]["EXCHANGE"]["g[i,j]"])

    return p


# Fix every Pomeron transition coupling to zero in the Good Walker eigenbasis
def _diagonal_pomeron_constraints(model: str, pomeron: str) -> dict:
    constraints = {}
    rows = source_channel_count(model)
    for row in range(rows):
        for column in range(row + 1, rows):
            raw_key = f"SOFT|MODEL.{model}:EXCHANGE.{pomeron}.g[{row},{column}]"
            constraints[symmetric_coupling_key(raw_key)] = 0.0
    return constraints


# Construct the eikonal tune parameter and auxiliary steering dictionaries
def setup(
    model: str | None = None,
    pFF: str | None = None,
    oFF: str | None = None,
    o3gFF: str | None = None,
    tune_reggeons: bool = False,
):
    model = source_model(model)
    pomeron = source_forward_exchange(model)
    tune_odderon = source_odderon_active(model)
    tune_o3g = source_o3g_active(model)
    reggeons = source_reggeon_exchanges(model) if tune_reggeons else ()
    if tune_reggeons and not reggeons:
        raise Exception(__name__ + ".setup: Reggeon tuning requires an enabled exchange")
    soft_forms = {"EXP", "DPOW", "EXPOW", "GKERNEL"}
    requests = [((pomeron,), pFF, soft_forms)]
    if tune_odderon:
        requests.append((("O",), oFF, soft_forms | {"ODD3G_NODE"}))
    if tune_o3g:
        requests.append((("O3g",), o3gFF, {"3G"}))
    if reggeons:
        requests.append((reggeons, None, soft_forms))
    banks = []
    for exchanges, value, allowed in requests:
        bank = (
            source_common_form_factor_bank(model, exchanges)
            if exchanges == reggeons
            else source_form_factor_bank(model=model, exchange=exchanges[0], value=value)
        )
        banks.append((exchanges, bank, allowed))

    # Get general
    p = general(model=model, tune_odderon=tune_odderon, tune_o3g=tune_o3g, tune_reggeons=tune_reggeons)

    # Get form factor
    for _, bank, allowed in banks:
        name = str(source_model_block(model)["FF"][bank]["type"])
        p.update(_form_factor_parameters(name, model=model, bank=bank, allowed=allowed))

    # Constants
    p_aux = {"SOFT|active_model": model, f"SOFT|MODEL.{model}:EXCHANGE.{pomeron}.ff": banks[0][1]}
    p_aux.update(_diagonal_pomeron_constraints(model, pomeron))
    for exchanges, bank, _ in banks[1:]:
        p_aux.update({f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.ff": bank for exchange in exchanges})
    return p, p_aux
