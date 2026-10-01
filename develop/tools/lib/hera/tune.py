# Write an isolated tune with all derived photoproduction parameters
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import shutil
from contextlib import redirect_stdout
from pathlib import Path

import pyjson5
from core.io.serialize import load_json_file

from .. import common, push
from . import couplings, model


# Resolve source card references so that the candidate tune is independent of its source
def copy(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if destination == source or source in destination.parents:
        raise ValueError("An isolated tune must be outside the input tune")
    shutil.copytree(source, destination)
    for path in source.rglob("*.json"):
        (destination / path.relative_to(source)).write_text(common.dumps(load_json_file(path, loader=pyjson5.load)) + "\n")


# Regenerate elastic couplings, transfer profiles and dissociation rows together
def write(source, destination, selected):
    copy(source, destination)
    with (Path(destination) / "HERA_update.log").open("w") as log, redirect_stdout(log):
        push_cards(Path(destination) / "GENERAL.json", selected, confirm=lambda rows: True)


# Compute the unique PARAM_REGGE photoproduction row for one vector meson
def _row_index(regge: dict, field: str, pdg: int) -> int:
    rows = regge[field]
    columns = {"photoprod": 5, "photoprod_diss": 8, "photoprod_dlog": 3}[field]
    if not isinstance(rows, list):
        raise TypeError(f"PARAM_REGGE.{field} must be an array")
    matches = [
        index for index, row in enumerate(rows) if isinstance(row, list) and len(row) == columns and row[0] == pdg
    ]
    if len(matches) != 1:
        raise ValueError(f"PARAM_REGGE.{field} requires one row for PDG {pdg}")
    return matches[0]


# Compute direct production coupling block
def _helicity_block(
    theory: object, first_pdg: int, second_pdg: int
) -> str:
    if not isinstance(theory, dict):
        raise TypeError("PARAM_RES.MODELS production model must be an object")
    key = f"[{first_pdg},{second_pdg}]"
    block = theory.get(key)
    if not isinstance(block, dict) or block.get("basis") != "helicity":
        raise ValueError(f"PARAM_RES.MODELS requires helicity basis in block {key}")
    return key


# Find the unique resonance card implementing one HERA production channel
def _resonance_card(res_directory: Path, channel: model.HERAChannel) -> tuple[Path, dict]:
    matches = []
    for path in sorted(res_directory.glob("*.json")):
        with path.open(encoding="utf-8") as stream:
            card = load_json_file(stream.name, loader=pyjson5.load)
        param_res = card.get("PARAM_RES", {})
        if param_res.get("PDG") != channel.pdg:
            continue
        models = param_res.get("MODELS", {})
        valid = True
        for theory, first_pdg, second_pdg in channel.production_channels:
            key = f"[{first_pdg},{second_pdg}]"
            valid = valid and key in models.get(theory, {})
        if valid:
            matches.append((path, card))
    if len(matches) != 1:
        raise ValueError(
            f"RES directory requires one HERA resonance card for {channel.name}, found {len(matches)}"
        )
    return matches[0]


# Build scalar updates for selected HERA trajectory and resonance parameters
def card_updates(
    general_path: Path, selected: list[tuple[model.HERAChannel, model.DerivedChannel]]
) -> dict[Path, list[push.ScalarUpdate]]:
    with general_path.open(encoding="utf-8") as stream:
        general = load_json_file(stream.name, loader=pyjson5.load)
    general_updates = []
    updates_by_card: dict[Path, list[push.ScalarUpdate]] = {general_path: general_updates}
    fields = {
        "photoprod": ("W0", "B_vector", "alpha0", "alpha_prime"),
        "photoprod_diss": ("W0", "ratio", "delta", "b", "n", "epsilon", "M_max"),
    }
    for channel, derived in selected:
        cards = couplings.channel_output(channel, derived)["cards"]
        for field, names in fields.items():
            if field not in cards:
                continue
            row_index = _row_index(general["PARAM_REGGE"], field, channel.pdg)
            for column, name in enumerate(names, start=1):
                general_updates.append(
                    push.ScalarUpdate(
                        path=("PARAM_REGGE", field, row_index, column),
                        label=f"GENERAL.json:PARAM_REGGE.{field}[{channel.pdg}].{name}",
                        value=cards[field][column],
                    )
                )

        resonance_path, resonance = _resonance_card(general_path.parent / "RES", channel)
        resonance_updates = updates_by_card.setdefault(resonance_path, [])
        models = resonance["PARAM_RES"]["MODELS"]
        desired_models = couplings.production(channel, derived.gamma_pomeron_vector_coupling_per_gev)
        for theory, first_pdg, second_pdg in channel.production_channels:
            key = _helicity_block(models[theory], first_pdg, second_pdg)
            rows = models[theory][key]["helicity"]
            desired_rows = desired_models[theory][key]["helicity"]
            desired = {(row[0], row[1]): row[2] for row in desired_rows}
            if not isinstance(rows, list) or len(rows) != len(desired) or any(
                not isinstance(row, list) or len(row) != 4 for row in rows
            ):
                raise ValueError(f"{resonance_path}:{theory}.{key} requires one HERA helicity row")
            for index, row in enumerate(rows):
                coordinate = (row[0], row[1])
                if coordinate not in desired:
                    raise ValueError(f"missing HERA helicity coordinate {coordinate}")
                resonance_updates.append(
                    push.ScalarUpdate(
                        path=("PARAM_RES", "MODELS", theory, key, "helicity", index, 2),
                        label=(
                            f"RES/{resonance_path.name}:PARAM_RES.MODELS.{theory}.{key}.helicity[{coordinate}].magnitude"
                        ),
                        value=desired[coordinate],
                    )
                )
    return updates_by_card


# Preview and push selected HERA photoproduction parameters
def push_cards(
    general_path: Path,
    selected: list[tuple[model.HERAChannel, model.DerivedChannel]],
    *,
    confirm=None,
) -> bool:
    updates = card_updates(general_path, selected)
    general = load_json_file(general_path, loader=pyjson5.load)
    payloads = [couplings.channel_output(channel, derived) for channel, derived in selected]
    pdgs = {channel.pdg for channel, _ in selected}
    arrays = {}
    for field in ("photoprod_dlog", "photoprod_diss_dlog"):
        rows = [row for row in general["PARAM_REGGE"].get(field, []) if row[0] not in pdgs]
        for payload in payloads:
            value = payload["cards"].get(field)
            if value is not None:
                rows.extend(value if field == "photoprod_diss_dlog" else [value])
        if rows or field in general["PARAM_REGGE"]:
            arrays[("PARAM_REGGE", field)] = rows
    options = {} if confirm is None else {"confirm": confirm}
    return push.push_json5_updates(updates, array_updates={general_path: arrays}, **options)
