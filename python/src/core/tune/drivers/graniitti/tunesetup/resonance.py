# Resonance tuning helpers for MP, XP, GP and TP models
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re

import numpy as np
from ray import tune

from core.io import steering
from core.tune.drivers.graniitti.tunesetup.domains import (
    add_direction_params,
    add_magnitude_phase_params,
    add_residue_params,
    angular_row_name,
    load_model_json,
    model_dir,
    mp_spin_sectors,
    settings,
    source_root,
    tp_res_names,
)
from core.tune.parameters import tools


# Load decay branching fractions from the source model card
def _load_branching_data() -> dict:
    return load_model_json(model_dir() / "DECAYS.json")


# Compute the parameter block from one resonance source card
def _load_resonance_param(resonance: str) -> dict:
    return load_model_json(model_dir() / "RES" / f"{resonance}.json")["PARAM_RES"]



def _has_model_channels(resonance: str, model: str) -> bool:
    model_block = _load_resonance_param(resonance)["MODELS"].get(model)
    return model_block is not None and any(re.fullmatch(r"\[-?\d+,-?\d+\]", key) for key in model_block)


# Compute true when one production model contains a nonzero physical coupling
def _has_model_coupling(resonance: str, model: str) -> bool:
    model_block = _load_resonance_param(resonance)["MODELS"].get(model)
    if model_block is None:
        return False
    for key, block in model_block.items():
        if not re.fullmatch(r"\[-?\d+,-?\d+\]", key):
            continue
        if model == "TP":
            if any(float(value) != 0.0 for value in block.get("g_tensor", [])):
                return True
            continue
        basis = block.get("basis")
        if basis in {"auto_min_L", "auto_min_S", "auto_equal_ls", "auto_equal_helicity"}:
            rows = [block.get("g", [])]
            magnitude_index = 0
        else:
            rows = block.get("g_ls" if basis == "g_ls" else "helicity", [])
            magnitude_index = 2
        if any(
            len(row) > magnitude_index and (row[magnitude_index] is None or float(row[magnitude_index]) != 0.0)
            for row in rows
        ):
            return True
    return False



def discover_resonances(*, model: str | None = None, production_required: bool = False) -> list[str]:
    """Discover source resonance cards usable by the requested production model."""

    resonances = []
    for path in sorted((model_dir() / "RES").glob("*.json")):
        resonance = path.stem
        if model is not None and not _has_model_channels(resonance, model):
            continue
        if production_required and not _has_model_channels(resonance, "XP"):
            continue
        resonances.append(resonance)

    if not resonances:
        detail = f" model={model}" if model is not None else ""
        raise ValueError(f"No source resonance cards discovered for{detail}")
    return resonances


def discover_resonances_from_datacards(
    datacards: list[dict], *, model: str | None = None, production_required: bool = False
) -> list[str]:
    """Discover resonance names from the selected datacards and their gencards."""

    resonances = []
    seen = set()
    for datacard in datacards:
        data, datacard_path = steering.load_dataset(datacard["datacard"], cdir=source_root())
        for sample in data["samples"]:
            gencard_path = steering.resolve_gencard_reference(
                sample["gencard"], dataset_path=datacard_path, cdir=source_root()
            )
            gencard = load_model_json(gencard_path)
            for resonance in gencard["SCATTERING"].get("RES", []):
                if resonance in seen:
                    continue
                if model is not None and not _has_model_channels(resonance, model):
                    continue
                if production_required and not _has_model_channels(resonance, "XP"):
                    continue
                seen.add(resonance)
                resonances.append(resonance)

    return resonances


# Compute default K+K- branching-ratio targets present in selected resonances
def default_kk_branching_ratio_targets(resonances: list[str]) -> list[tuple[int, list[int]]]:
    selected = set(resonances)
    return [target for resonance, target in settings()["icetune"]["resonance"]["decay"]["branching"].items() if resonance in selected]


# Compute the resonance spin from its doubled card value
def _resonance_spin(resonance: str) -> float:
    return 0.5 * float(_load_resonance_param(resonance)["spinX2"])


def _discover_tensor_resonances(resonances: list[str]) -> list[str]:
    return [resonance for resonance in resonances if abs(_resonance_spin(resonance) - 2.0) < 1e-12]


def _discover_coherent_ajzp_resonances(resonances: list[str]) -> list[str]:
    selected = []
    for resonance in resonances:
        spin = _resonance_spin(resonance)
        if spin >= 1.0 and abs(spin - round(spin)) < 1e-12:
            selected.append(resonance)
    return selected


# Select resonances with MP production channels for independent spin steering
def _mp_resonances(resonances: list[str]) -> list[str]:
    return [name for name in resonances if _has_model_channels(name, "MP")]


# Require all active MP blocks to enable C and parity symmetry
def _validate_mp_symmetry(resonances: list[str]) -> None:
    for resonance in resonances:
        for _, key, block in _active_inline_model_blocks(resonance, "MP"):
            if block.get("CP") != [True, True]:
                raise ValueError(
                    f'MP.{key} for resonance "{resonance}" must enforce C and parity symmetry with CP=[true,true]'
                )


def _normalize_channel_pair(values: list) -> tuple[int, int]:
    if len(values) < 2:
        raise ValueError(f"Expected a 2-body production channel, got {values}")
    return tuple(sorted((int(values[0]), int(values[1]))))


# Compute active eta-style inline central coupling blocks for a resonance model
def _active_inline_model_blocks(resonance: str, model: str) -> list[tuple[int, str, dict]]:
    model_block = _load_resonance_param(resonance)["MODELS"].get(model)
    if model_block is None:
        raise KeyError(f'Missing {model} model block for resonance "{resonance}"')

    blocks = []
    keys = sorted(key for key in model_block if re.fullmatch(r"\[-?\d+,-?\d+\]", key))
    for channel_index, key in enumerate(keys):
        cp = model_block[key].get("CP")
        if cp != [True, True]:
            raise ValueError(
                f'{model}.{key} for resonance "{resonance}" must enforce C and parity symmetry with CP=[true,true]'
            )
        blocks.append((channel_index, key, model_block[key]))
    if not blocks:
        raise ValueError(f'No {model} pair entries found for resonance "{resonance}"')
    return blocks


def _normalize_final_state(values: list[int]) -> tuple[int, ...]:
    return tuple(sorted(int(value) for value in values))


# Check whether one DECAYS.json key identifies a physical decay channel
def _is_decay_channel_key(channel_id: str) -> bool:
    return re.fullmatch(r"\[-?\d+(?:,-?\d+)*\]", channel_id) is not None


# Parse one compact DECAYS.json PDG-array key
def _decay_channel_pdgs(channel_id: str) -> list[int]:
    if not _is_decay_channel_key(channel_id):
        raise ValueError(f"Invalid DECAYS.json channel key {channel_id}")
    return [int(value) for value in channel_id[1:-1].split(",")]


# Discover the selected two-body decay-coupling basis and its row indices
def _discover_decay_coupling_specs(resonances: list[str]) -> list[tuple[str, str, str, list[str]]]:
    branching_data = _load_branching_data()
    specs = []
    seen = set()

    for resonance in resonances:
        pdg_id = str(_load_resonance_param(resonance)["PDG"])
        if pdg_id not in branching_data:
            raise KeyError(f"Missing DECAYS.json block for resonance PDG {pdg_id} ({resonance})")

        for channel_id, channel_block in branching_data[pdg_id].items():
            channel_id = str(channel_id)
            if not _is_decay_channel_key(channel_id):
                continue
            basis = channel_block.get("basis")
            if basis not in {"alpha_ls", "helicity"}:
                continue
            rows = channel_block.get(basis)
            if not isinstance(rows, list) or not rows:
                raise ValueError(f"DECAYS.json channel {pdg_id}:{channel_id} has no active {basis} rows")
            row_names = [angular_row_name(basis, row) for row in rows]
            if len(row_names) != len(set(row_names)):
                raise ValueError(f"DECAYS.json channel {pdg_id}:{channel_id} has duplicate {basis} rows")
            spec_key = (pdg_id, channel_id, basis, tuple(row_names))
            if row_names and spec_key not in seen:
                specs.append((pdg_id, channel_id, basis, row_names))
                seen.add(spec_key)

    return specs


def _normalize_branching_ratio_targets(
    targets: list[tuple[int, list[int] | tuple[int, int]]],
) -> list[tuple[str, tuple[int, int]]]:
    normalized = []
    seen = set()
    for resonance_pdg, final_state in targets:
        spec = (str(int(resonance_pdg)), _normalize_channel_pair(list(final_state)))
        if spec not in seen:
            normalized.append(spec)
            seen.add(spec)
    return normalized


# Compute whether one resonance keeps its decay phases and branching fractions fixed
def _decay_phase_and_branching_are_fixed(pdg_id: str) -> bool:
    return str(pdg_id) == "9000221"


# Discover explicitly selected free two-body branching fractions
def _discover_decay_branching_ratio_specs(
    *, resonances: list[str], branching_ratio_targets: list[tuple[int, list[int] | tuple[int, int]]]
) -> list[tuple[str, str]]:
    normalized_targets = [
        target
        for target in _normalize_branching_ratio_targets(branching_ratio_targets)
        if not _decay_phase_and_branching_are_fixed(target[0])
    ]
    selected_pdgs = {str(_load_resonance_param(resonance)["PDG"]) for resonance in resonances}
    missing_pdgs = sorted({pdg_id for pdg_id, _ in normalized_targets if pdg_id not in selected_pdgs})
    if missing_pdgs:
        raise ValueError(
            f"Requested branching-ratio tuning for resonance PDGs not present in selected resonances: {missing_pdgs}"
        )

    branching_data = _load_branching_data()
    specs = []
    for pdg_id, final_state in normalized_targets:
        if pdg_id not in branching_data:
            raise KeyError(f"Missing DECAYS.json block for resonance PDG {pdg_id}")

        matched_channels = [
            str(channel_id)
            for channel_id in branching_data[pdg_id]
            if _is_decay_channel_key(str(channel_id))
            and len(_decay_channel_pdgs(str(channel_id))) == 2
            and _normalize_channel_pair(_decay_channel_pdgs(str(channel_id))) == final_state
        ]
        if not matched_channels:
            raise ValueError(f"No 2-body DECAYS.json channel matched resonance PDG {pdg_id} and pair {final_state}")
        if len(matched_channels) > 1:
            raise ValueError(
                f"Ambiguous DECAYS.json 2-body channels matched resonance PDG {pdg_id} and pair {final_state}: {matched_channels}"
            )
        specs.append((pdg_id, matched_channels[0]))

    return specs


# Discover selected decay phases while retaining fixed resonance decays
def _discover_decay_zeta_specs_filtered(
    *, resonances: list[str], allowed_final_states: set[tuple[int, ...]] | None, allowed_resonances: set[str] | None
) -> list[tuple[str, str]]:
    branching_data = _load_branching_data()
    specs = []
    seen = set()

    for resonance in resonances:
        if allowed_resonances is not None and resonance not in allowed_resonances:
            continue
        pdg_id = str(_load_resonance_param(resonance)["PDG"])
        if pdg_id not in branching_data:
            raise KeyError(f"Missing DECAYS.json block for resonance PDG {pdg_id} ({resonance})")
        if _decay_phase_and_branching_are_fixed(pdg_id):
            continue

        for channel_id, channel_block in branching_data[pdg_id].items():
            if not _is_decay_channel_key(str(channel_id)):
                continue
            if "zeta" not in channel_block:
                continue
            if allowed_final_states is not None:
                channel_pdg = _decay_channel_pdgs(str(channel_id))
                if _normalize_final_state(channel_pdg) not in allowed_final_states:
                    continue
            spec_key = (pdg_id, str(channel_id))
            if spec_key not in seen:
                specs.append(spec_key)
                seen.add(spec_key)

    return specs


# Add one scalar phase in raw or Cayley coordinates
def _add_phase_param(p: dict, base_key: str, mode: str) -> None:
    mode = tools.validate_tune_mode(mode, name="phase_mode")
    if mode == "none":
        return
    if mode == "phase_raw":
        p[tools.raw_phase_key(base_key)] = tune.uniform(-np.pi, np.pi)
    elif mode == "phase_cayley":
        p[tools.phase_u_key(base_key)] = tune.uniform(-1.0, 1.0)
        p[tools.phase_v_key(base_key)] = tune.uniform(-1.0, 1.0)
    else:
        raise ValueError(f'Unknown phase mode "{mode}"; expected phase_raw or phase_cayley')


# Add one normalized signed real-projective coupling vector
def _add_projective_vector_params(p: dict, *, base_key: str, row_names: list[str], mode: str) -> None:
    mode = tools.validate_tune_mode(mode, name="projective_vector_mode")
    if mode == "none":
        return
    if mode != "magnitude":
        raise ValueError(
            "Normalized internal coupling vectors use the real-projective baseline; complex row phases require a separate complex-projective model"
        )
    prefix = base_key.rsplit(":", 1)[0]
    add_direction_params(p, [f"{prefix}:{row}{tools.PROJECTIVE_VECTOR_SUFFIX}" for row in row_names[1:]])


# Add observable decay phases and gauge-fixed active-basis coupling vectors
def _add_decay_coupling_params(
    p: dict,
    *,
    resonances: list[str],
    model: str,
    zeta_phase_mode: str,
    decay_coupling_mode: str,
    zeta_decay_pdgs: list[list[int]] | None,
    zeta_decay_resonances: list[str] | None,
) -> None:
    zeta_phase_mode = tools.validate_tune_mode(zeta_phase_mode, name="decay.zeta.mode")
    decay_coupling_mode = tools.validate_tune_mode(decay_coupling_mode, name="decay.mode")
    allowed_final_states = None
    if zeta_decay_pdgs is not None:
        allowed_final_states = {_normalize_final_state(final_state) for final_state in zeta_decay_pdgs}
    allowed_resonances = None
    if zeta_decay_resonances is not None:
        allowed_resonances = set(zeta_decay_resonances)

    for pdg_id, channel in _discover_decay_zeta_specs_filtered(
        resonances=resonances, allowed_final_states=allowed_final_states, allowed_resonances=allowed_resonances
    ):
        _add_phase_param(p, f"DECAY|{pdg_id}:{channel}:zeta.{model}", zeta_phase_mode)

    if decay_coupling_mode != "none":
        specs = _discover_decay_coupling_specs(resonances=resonances)
        for pdg_id, channel, basis, row_names in specs:
            _add_projective_vector_params(
                p, base_key=f"DECAY|{pdg_id}:{channel}:{basis}", row_names=row_names, mode=decay_coupling_mode
            )


# Add selected free decay branching fractions
def _add_branching_ratio_params(
    p: dict, *, resonances: list[str], branching_ratio_targets: list[tuple[int, list[int] | tuple[int, int]]]
) -> None:
    specs = _discover_decay_branching_ratio_specs(
        resonances=resonances, branching_ratio_targets=branching_ratio_targets
    )
    targets_by_pdg = {}
    for pdg_id, channel in specs:
        targets_by_pdg.setdefault(pdg_id, []).append(channel)
    for pdg_id, channels in targets_by_pdg.items():
        if len(channels) != 1:
            raise ValueError(
                f"Multiple free branching fractions for PDG {pdg_id} require an explicit simplex parametrization"
            )
        channel = channels[0]
        fixed_sum = sum(
            float(block.get("BR", 0.0))
            for name, block in _load_branching_data()[pdg_id].items()
            if str(name) != channel
        )
        upper = 1.0 - fixed_sum
        if upper <= 0.0:
            raise ValueError(f"No branching-fraction remainder is available for PDG {pdg_id}")
        p[f"DECAY|{pdg_id}:{channel}:BR"] = tune.uniform(0.0, upper)


# Compute the production tune mode after applying fixed-phase overrides
def _production_tune_mode(resonance: str, production_mode: str, magnitude_only_resonances: set[str]) -> str:
    production_mode = tools.validate_tune_mode(production_mode, name="production.mode")
    if production_mode != "none" and resonance in magnitude_only_resonances:
        return "magnitude"
    return production_mode


# Add one direct production coupling in the selected coordinate mode
def _add_direct_coupling_params(
    p: dict,
    base_key: str,
    max_abs: float,
    *,
    production_mode: str,
    derived_magnitude: bool = False,
    angular: bool = False,
    bounds: tuple[float, float] | None = None,
    initial: float | None = None,
    relative_bounds=None,
) -> None:
    production_mode = tools.validate_tune_mode(production_mode, name="production.mode")
    if derived_magnitude and tools.tune_mode_is_cartesian(production_mode):
        return
    if relative_bounds is not None and initial is not None:
        bounds = tuple(abs(initial) * scale for scale in relative_bounds)
    if angular:
        magnitude_key = tools.magnitude_key(base_key)
        phase_key = base_key
    elif base_key.endswith("]"):
        magnitude_key = f"{base_key[:-1]},2]"
        phase_key = f"{base_key[:-1]},3]"
    else:
        magnitude_key = f"{base_key}[0]"
        phase_key = f"{base_key}[1]"
    add_magnitude_phase_params(
        p,
        base_key=base_key,
        magnitude_key=magnitude_key,
        phase_key=phase_key,
        mode=production_mode,
        magnitude_bounds=(1.0, 1.0) if derived_magnitude else (bounds or (0.0, max_abs)),
    )


# Add normalized signed real MP spin amplitude direction coordinates
def _add_coherent_ajz_params(p: dict, resonances: list[str], geometry: str, frame: str = "CS") -> None:
    if geometry not in {"projective", "spherical"}:
        raise ValueError(f'Unknown coherent a_Jz geometry "{geometry}"')
    for res in resonances:
        J = int(round(_resonance_spin(res)))
        vector_size = len(mp_spin_sectors(J, frame))
        base = f"RES|{res}:MP:polarization.a_Jz"
        angle_key = tools.ajzp_angle_key if geometry == "projective" else tools.spherical_angle_key
        add_direction_params(p, [angle_key(base, k) for k in range(vector_size - 1)], geometry=geometry)


# Add an overall norm and signed direction for active production LS rows
def _add_production_ls_params(
    p: dict, *, res: str, model: str, rows: list, max_abs: float, mode: str, geometry: str, relative_bounds=None
) -> bool:
    allowed = {"magnitude", "mag_phase_cartesian"} if geometry == "projective" else {"magnitude"}
    if mode not in allowed:
        modes = "magnitude or mag_phase_cartesian" if geometry == "projective" else "magnitude"
        raise ValueError(f"{geometry.capitalize()} production LS geometry requires {modes} mode")
    if any(row[2] is None for row in rows):
        return False
    active = [row for row in rows if abs(float(row[2])) > 0.0]
    if not active:
        return True
    if any(abs(np.sin(float(row[3]))) > 1.0e-8 for row in active):
        raise ValueError(f'{geometry.capitalize()} production LS rows must be real for resonance "{res}"')
    norm = float(np.linalg.norm([row[2] for row in active]))
    bounds = tuple(norm * scale for scale in relative_bounds) if relative_bounds is not None else (0.0, max_abs)
    prefix = f"RES|{res}:{model}:"
    if len(active) == 1 and mode == "magnitude":
        _add_direct_coupling_params(
            p, f"{prefix}{angular_row_name('g_ls', active[0])}", max_abs, production_mode=mode, angular=True, bounds=bounds
        )
    else:
        add_residue_params(
            p,
            rows=[f"{prefix}{angular_row_name('g_ls', row)}" for row in active],
            bounds=bounds,
            mode=mode,
            geometry=geometry,
        )
    return True


# Add parity-symmetric incoherent spin-population simplex coordinates
def _add_incoherent_density_params(p: dict, resonances: list) -> None:
    for res in resonances:
        J = int(round(_resonance_spin(res)))
        base_key = f"RES|{res}:MP:polarization.rho_population"
        for k in range(J):
            p[tools.simplex_theta_key(base_key, k)] = tune.uniform(0.0, 0.5 * np.pi)


# Add direct resonance couplings for one named production model
def _add_resonance_model_couplings(
    p: dict,
    resonances: list,
    *,
    model: str,
    production_mode: str,
    production_magnitude_only_resonances: list[str] | tuple[str, ...] | None = None,
    production_max_abs: float | None = None,
    active_only: bool = False,
    geometry: str = "direct",
    relative_bounds=None,
) -> None:
    magnitude_only_resonances = set(production_magnitude_only_resonances or [])
    production_max_abs = float(settings()["icetune"]["resonance"]["production"]["max_abs"][model] if production_max_abs is None else production_max_abs)
    if not np.isfinite(production_max_abs) or production_max_abs <= 0.0:
        raise ValueError("Resonance production maximum magnitude must be positive and finite")
    if geometry not in {"direct", "projective", "spherical"}:
        raise ValueError(f'Unknown production LS geometry "{geometry}"')
    for res in resonances:
        scales = relative_bounds.get(res, relative_bounds["default"]) if isinstance(relative_bounds, dict) else relative_bounds
        model_block = _load_resonance_param(res)["MODELS"].get(model)
        if model_block is None:
            raise KeyError(f'Missing {model} model block for resonance "{res}"')
        direct_mode = _production_tune_mode(res, production_mode, magnitude_only_resonances)
        if direct_mode == "none":
            continue
        for _, _, block in _active_inline_model_blocks(res, model):
            basis = block.get("basis")
            if model == "GP" and basis not in {"g_ls", "helicity"}:
                raise ValueError(f'GP block basis for resonance "{res}" must be g_ls or helicity')
            if model != "GP" and basis not in {
                "auto_min_L",
                "auto_min_S",
                "auto_equal_ls",
                "auto_equal_helicity",
                "g_ls",
                "helicity",
            }:
                raise ValueError(f'{model} block basis for resonance "{res}" is not supported')
            if basis in {
                "auto_min_L",
                "auto_min_S",
                "auto_equal_ls",
                "auto_equal_helicity",
            }:
                coupling = block.get("g")
                if not isinstance(coupling, list) or len(coupling) != 2:
                    raise ValueError(f'Missing {model} block g for resonance "{res}"')
                if active_only and coupling[0] == 0.0:
                    continue
                _add_direct_coupling_params(
                    p,
                    f"RES|{res}:{model}:g",
                    production_max_abs,
                    production_mode=direct_mode,
                    derived_magnitude=coupling[0] is None,
                    initial=coupling[0],
                    relative_bounds=scales,
                )
                continue
            field = "g_ls" if basis == "g_ls" else "helicity"
            rows = block.get(field)
            if not isinstance(rows, list) or not rows:
                raise ValueError(f'Missing {model} block {field} for resonance "{res}"')
            if (
                geometry in {"projective", "spherical"}
                and field == "g_ls"
                and _add_production_ls_params(
                    p, res=res, model=model, rows=rows, max_abs=production_max_abs, mode=direct_mode, geometry=geometry,
                    relative_bounds=scales,
                )
            ):
                continue
            for row in rows:
                if active_only and row[2] == 0.0:
                    continue
                _add_direct_coupling_params(
                    p,
                    f"RES|{res}:{model}:{angular_row_name(field, row)}",
                    production_max_abs,
                    production_mode=direct_mode,
                    derived_magnitude=row[2] is None,
                    angular=True,
                    initial=row[2],
                    relative_bounds=scales,
                )


# Add standalone direct resonance production coupling parameters
def production_setup(
    *,
    resonances: list[str],
    model: str,
    mode: str = "magnitude",
    magnitude_only: list[str] | tuple[str, ...] | None = (),
    max_abs: float | None = None,
    active_only: bool = False,
    geometry: str = "direct",
):
    p = {}
    p_aux = {}
    selected_resonances = list(resonances)

    if model not in {"MP", "XP", "GP"}:
        raise ValueError(f'Unknown resonance production model "{model}"')
    if model == "MP" and geometry != "direct":
        raise ValueError("Projective production LS geometry applies only to XP and GP")
    _add_resonance_model_couplings(
        p,
        selected_resonances,
        model=model,
        production_mode=mode,
        production_magnitude_only_resonances=magnitude_only,
        production_max_abs=max_abs,
        active_only=active_only,
        geometry=geometry,
    )
    return p, p_aux


# Add one periodic model level coupling phase for each active resonance
def phase_setup(*, resonances: list[str], model: str, mode: str = "phase_raw"):
    if model not in {"MP", "XP", "GP", "TP"}:
        raise ValueError(f'Unknown resonance phase model "{model}"')
    if mode not in {"none", "phase_raw"}:
        raise ValueError('Resonance phi mode must be "none" or "phase_raw"')

    p = {}
    for resonance in resonances:
        if not _has_model_channels(resonance, model):
            raise KeyError(f'Missing {model} model block for resonance "{resonance}"')
        if mode == "phase_raw" and _has_model_coupling(resonance, model):
            key = f"RES|{resonance}:{model}:phi"
            p[tools.raw_phase_key(key)] = tune.uniform(-np.pi, np.pi)
    return p, {}


# Compute TP coefficient names and the canonical nonzero reference row
def tensor_names(resonance: str, channel: str) -> tuple[tuple[str, ...], int]:
    param = _load_resonance_param(resonance)
    couplings = param["MODELS"].get("TP", {}).get(channel, {}).get("g_tensor")
    names = tp_res_names(param)
    if not isinstance(couplings, list) or len(names) != len(couplings):
        raise ValueError(f'Invalid TP g_tensor size for resonance "{resonance}"')
    reference = next((i for i, value in enumerate(couplings) if abs(float(value)) > 1.0e-15), None)
    if reference is None:
        raise ValueError(f'TP production model is entirely zero for resonance "{resonance}"')
    return names, reference


# Add TP resonance tensor couplings in direct or normalized direction coordinates
def tensor_setup(
    *,
    resonances: list[str],
    channel: str = "[995,995]",
    geometry: str = "direct",
    max_abs: float | None = None,
) -> tuple[dict, dict]:
    if geometry not in {"direct", "projective", "spherical"}:
        raise ValueError(f'Unknown TP tensor geometry "{geometry}"')
    max_abs = settings()["icetune"]["resonance"]["production"]["max_abs"]["TP"] if max_abs is None else max_abs
    if not np.isfinite(max_abs) or max_abs <= 0.0:
        raise ValueError("TP tensor maximum magnitude must be positive and finite")
    p = {}
    for res in resonances:
        names, reference = tensor_names(res, channel)
        selected = names[reference:]
        prefix = f"RES|{res}:TP:{channel}:"
        if geometry == "direct" or len(selected) == 1:
            for name in selected:
                p[f"{prefix}{name}"] = tune.uniform(0.0, max_abs)
            continue
        add_residue_params(
            p, rows=[f"{prefix}{name}" for name in selected], bounds=(0.0, max_abs), mode="magnitude", geometry=geometry
        )
    return p, {}


# Compute fixed MP spin basis steering for selected resonances
def spin_basis(*, resonances: list[str], basis: str):
    if basis not in {"a_Jz", "rho"}:
        raise ValueError(f'Unknown MP spin basis "{basis}"')
    selected_resonances = _mp_resonances(list(resonances))
    _validate_mp_symmetry(selected_resonances)
    return {}, {f"RES|{resonance}:MP:polarization.mode": basis for resonance in selected_resonances}


# Add explicit resonance and decay form factor coordinates
def res_form_factors(model: str, *, resonances: list[str], fields: set[str]) -> dict:
    decays = _load_branching_data() if "FF_decay" in fields else {}
    p = {}
    for name in resonances:
        res = _load_resonance_param(name)
        forms = res["MODELS"][model]
        if "FF_prod" in fields and forms["FF_prod"]["type"] != "none":
            p[f"RES|{name}:{model}:FF_prod.Lambda2"] = tune.uniform(*settings()["bounds"]["resonance"]["MODELS"][model]["FF_prod"]["Lambda2"])
        if "FF_transfer" in fields and forms["FF_transfer"]["type"] != "none":
            transfer = forms["FF_transfer"]
            if transfer["type"] == "power":
                p[f"RES|{name}:{model}:FF_transfer.LambdaInv2"] = tune.uniform(*settings()["icetune"]["FF_transfer"]["LambdaInv2"])
            elif transfer["type"] == "exp":
                p[f"RES|{name}:{model}:FF_transfer.b"] = tune.uniform(*(scale * transfer["b"] for scale in settings()["icetune"]["FF_transfer"]["b"]["relative_bounds"]))
            else:
                raise ValueError(f"{name}:{model} transfer tuning requires a power or exp form")
        if "FF_decay" in fields:
            for channel, block in decays[str(res["PDG"])].items():
                if _is_decay_channel_key(channel) and block["FF_decay"][model]["type"] != "none":
                    p[f"DECAY|{res['PDG']}:{channel}:FF_decay.{model}.Lambda2"] = tune.uniform(
                        *settings()["bounds"]["decay"]["FF_decay"][model]["Lambda2"])
    return p


# Initialize direct resonance couplings for one process-level model
def _model_production_setup(
    model: str,
    resonances: list | None,
    production_mode: str,
    production_magnitude_only_resonances,
    production_max_abs: float,
    active_only: bool,
    geometry: str,
    relative_bounds=None,
) -> tuple[dict, list[str]]:
    selected = list(resonances) if resonances is not None else discover_resonances(model=model, production_required=model == "XP")
    parameters = {}
    _add_resonance_model_couplings(
        parameters,
        selected,
        model=model,
        production_mode=production_mode,
        production_magnitude_only_resonances=production_magnitude_only_resonances,
        production_max_abs=production_max_abs,
        active_only=active_only,
        geometry=geometry,
        relative_bounds=relative_bounds,
    )
    return parameters, selected


# Build one MP parameter block in the selected production basis
def _mp_model_setup(resonances, production_mode, magnitude_only, max_abs, active_only, discover, add, basis, geometry, frame,
                    relative_bounds=None):
    p, selected = _model_production_setup(
        "MP", resonances, production_mode, magnitude_only, max_abs, active_only, "direct", relative_bounds=relative_bounds
    )
    selected = _mp_resonances(selected)
    _validate_mp_symmetry(selected)
    if basis == "a_Jz":
        add(p, discover(selected), geometry, frame)
    else:
        add(p, discover(selected))
    return p, spin_basis(resonances=selected, basis=basis)[1]


# Merge and validate one grouped setup configuration
def _config(name: str, values: dict | None, defaults: dict) -> dict:
    if values is None:
        return dict(defaults)
    if not isinstance(values, dict):
        raise TypeError(f"{name} configuration must be a dictionary")
    unknown = sorted(set(values) - set(defaults))
    if unknown:
        raise ValueError(f"Unknown {name} configuration keys: {unknown}")
    return {**defaults, **values}


# Validate the explicitly varied resonance parameter families
def _vary_names(values) -> set[str]:
    if values is None:
        return set()
    if isinstance(values, str):
        raise TypeError("vary must be a sequence of parameter-family names")
    selected = {str(value) for value in values}
    unknown = sorted(selected - {"FF_prod", "FF_decay", "FF_transfer", "omega", "mass_width"})
    if unknown:
        raise ValueError(f"Unknown varied resonance parameter families: {unknown}")
    return selected


# Add the full resonance tuning block for one process-level model
def setup(
    *,
    model: str = "MP",
    resonances: list | None = None,
    production: dict | None = None,
    decay: dict | None = None,
    vary=(),
    frame: str | None = None,
):
    p = {}
    if model == "MP":
        p_aux = {"REGGE|MP_FRAME": "CS" if frame is None else frame}
    else:
        if frame is not None:
            raise ValueError('frame applies only to model "MP"')
        p_aux = {}
    selected_resonances = list(resonances) if resonances is not None else discover_resonances(model=model, production_required=model == "XP")
    production_values = dict(production or {})
    default_coherence = "incoherent"
    coherence_value = production_values.get("coherence", default_coherence)
    default_production_mode = (
        "magnitude" if model == "TP" or (model == "MP" and coherence_value == "incoherent") else "mag_phase_cartesian"
    )
    production_config = _config(
        "production",
        production,
        {
            "basis": None,
            "coherence": default_coherence,
            "geometry": "direct",
            "max_abs": settings()["icetune"]["resonance"]["production"]["max_abs"][model],
            "magnitude_only": (),
            "mode": default_production_mode,
            "active_only": False,
            "relative_bounds": None,
        },
    )
    relative_bounds = production_config["relative_bounds"]
    if relative_bounds is not None:
        if "max_abs" in production_values:
            raise ValueError("Choose production.relative_bounds or production.max_abs, not both")
        if isinstance(relative_bounds, dict) and (
                "default" not in relative_bounds or set(relative_bounds) - {"default", *selected_resonances}):
            raise ValueError("production.relative_bounds requires a default and selected resonance names")
        for bounds in relative_bounds.values() if isinstance(relative_bounds, dict) else [relative_bounds]:
            if (np.shape(bounds) != (2,) or not np.all(np.isfinite(bounds))
                    or not 0.0 <= bounds[0] < bounds[1]):
                raise ValueError("production.relative_bounds requires two finite, nonnegative, increasing scales")
        if model not in {"GP", "XP", "MP"} or production_config["mode"] != "magnitude":
            raise ValueError("Relative production bounds require GP, XP or MP magnitude tuning")
    decay_config = _config("decay", decay, {"branching": None, "mode": "none", "zeta": None})
    zeta_config = _config(
        "decay.zeta", decay_config["zeta"], {"final_states": None, "mode": "none", "resonances": None}
    )
    if zeta_config["mode"] != "none":
        p_aux["TENSORPOM|use_zeta" if model == "TP" else f"REGGE|use_zeta.{model}"] = True
    varied = _vary_names(vary)
    production_max_abs = float(production_config["max_abs"])
    if not np.isfinite(production_max_abs) or production_max_abs <= 0.0:
        raise ValueError('production["max_abs"] must be a positive finite number')

    if model == "MP" and production_config["basis"] not in {
        None,
        "auto_min_L",
        "auto_min_S",
        "auto_equal_ls",
        "auto_equal_helicity",
        "g_ls",
        "helicity",
    }:
        raise ValueError("MP production basis is not supported")
    if model != "MP" and "basis" in production_values:
        raise ValueError('production["basis"] applies only to model "MP"')
    if model == "MP" and production_config["basis"] is not None:
        basis = production_config["basis"]
        field = "g" if basis.startswith("auto_") else basis
        for name in selected_resonances:
            for _, _, block in _active_inline_model_blocks(name, "MP"):
                if field not in block:
                    raise ValueError(f"MP basis {basis} requires {field} in the source card for {name}")
    if model == "MP" and production_config["geometry"] not in {"direct", "projective", "spherical"}:
        raise ValueError('MP production["geometry"] must be direct, projective or spherical')
    if model != "MP":
        res_only = {"coherence"}.intersection(production_values)
        if res_only:
            raise ValueError(f"MP-only production configuration keys: {sorted(res_only)}")

    if "mass_width" in varied:
        p.update({f"RES|{name}:{model}:{field}": tune.uniform(*limits) for name in selected_resonances if name in settings()["icetune"]["resonance"]["mass_width"]
                  for field, limits in settings()["bounds"]["resonance"]["PARAM_RES"].get(name, {}).items()})
    coherence = production_config["coherence"] if model == "MP" else None
    if model in {"XP", "GP", "TP"} or coherence == "coherent":
        _add_decay_coupling_params(
            p,
            resonances=selected_resonances,
            model=model,
            zeta_phase_mode=zeta_config["mode"],
            decay_coupling_mode=decay_config["mode"],
            zeta_decay_pdgs=zeta_config["final_states"],
            zeta_decay_resonances=zeta_config["resonances"],
        )

    if model in {"XP", "GP"}:
        model_p, _ = _model_production_setup(
            model,
            selected_resonances,
            production_config["mode"],
            production_config["magnitude_only"],
            production_max_abs,
            production_config["active_only"],
            production_config["geometry"],
            relative_bounds=relative_bounds,
        )
        model_aux = {}
    elif model == "TP":
        if production_config["mode"] not in {"none", "magnitude"}:
            raise ValueError("TP production tuning requires real tensor couplings")
        model_p, model_aux = tensor_setup(
            resonances=selected_resonances if production_config["mode"] != "none" else [],
            geometry=production_config["geometry"],
            max_abs=production_max_abs,
        )
    elif model == "MP":
        if coherence == "incoherent":
            if zeta_config["mode"] != "none" or decay_config["mode"] != "none":
                raise ValueError("Incoherent MP tuning cannot contain decay amplitudes or phases")
            if production_config["mode"] not in {"none", "magnitude"}:
                raise ValueError("Incoherent MP tuning can vary resonance magnitudes but not phases")
        if coherence == "none":
            model_p, _ = _model_production_setup("MP", selected_resonances, production_config["mode"],
                                                production_config["magnitude_only"], production_max_abs,
                                                production_config["active_only"], production_config["geometry"],
                                                relative_bounds=relative_bounds)
            model_aux = {f"RES|{name}:MP:polarization.mode": "none" for name in selected_resonances}
        else:
            res_modes = {
                "coherent": (_discover_coherent_ajzp_resonances, _add_coherent_ajz_params, "a_Jz"),
                "incoherent": (_discover_tensor_resonances, _add_incoherent_density_params, "rho"),
            }
            if coherence not in res_modes:
                raise ValueError(f'Unknown MP production coherence "{coherence}"')
            model_p, model_aux = _mp_model_setup(
                selected_resonances,
                production_config["mode"],
                production_config["magnitude_only"],
                production_max_abs,
                production_config["active_only"],
                *res_modes[coherence],
                ("projective" if production_config["geometry"] == "direct" else production_config["geometry"]),
                p_aux["REGGE|MP_FRAME"],
                relative_bounds=relative_bounds,
            )
        if production_config["basis"] is not None:
            model_aux.update(
                {f"RES|{resonance}:MP:basis": production_config["basis"] for resonance in selected_resonances}
            )
    else:
        raise ValueError(f'Unknown resonance tuning model "{model}"')

    if decay_config["branching"] is not None:
        _add_branching_ratio_params(
            p, resonances=selected_resonances, branching_ratio_targets=decay_config["branching"]
        )
    if "omega" in varied:
        if model == "TP":
            raise ValueError("omega applies only to MP, XP and GP")
        p[f"REGGE|omega.{model}"] = tune.uniform(*settings()["bounds"]["resonance"]["PARAM_REGGE"]["omega"][model])
    if varied & {"FF_prod", "FF_decay", "FF_transfer"}:
        p.update(res_form_factors(model, resonances=selected_resonances, fields=varied))

    p.update(model_p)
    p_aux.update(model_aux)
    return p, p_aux
