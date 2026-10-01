# Construct GRANIITTI parameter spaces from explicit JSON study selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

from core.tune.drivers.graniitti.tunesetup import continuum, eikonal, resonance
from core.tune.drivers.graniitti.tunesetup.card import card
from core.tune.parameters.space import normalize_param_space


# Combine explicitly selected resonance and continuum parameters
def central(*, model: str, datacards: list[dict], production: dict | None = None, continuum_options: dict | None = None,
    decay: dict | None = None, phi: dict | None = None, select=(), exclude=(), vary=(), res_fields=(),
    spin_basis: str | None = None, production_required: bool = False, ) -> tuple[list[dict], dict, dict]:
    discovered = resonance.discover_resonances_from_datacards(
        datacards, model=model, production_required=production_required)
    unknown = set(select) - set(discovered)
    if unknown:
        raise ValueError(f"Requested resonances are not active: {sorted(unknown)}")
    selected = [name for name in discovered if (not select or name in select) and name not in exclude]
    production = dict(production or {"mode": "none"})
    phi = dict(phi or {"mode": "none"})
    unknown = set(phi.get("resonances", ())) - set(discovered)
    if unknown:
        raise ValueError(f"Requested resonance phi values are not active: {sorted(unknown)}")
    phase_resonances = [name for name in discovered if phi["mode"] != "none"
        and (not phi.get("resonances") or name in phi["resonances"]) and name not in phi.get("exclude", ())]
    if model in {"XP", "GP"} and production.get("mode") == "mag_phase_cartesian":
        cartesian = set(phase_resonances) & set(selected) - set(production.get("magnitude_only", ()))
        phase_resonances = [name for name in phase_resonances if name not in cartesian]
        production["magnitude_only"] = sorted(set(selected) - cartesian)

    parameters, auxiliary = resonance.setup(
        model=model, resonances=selected, production=production, decay=decay, vary=vary)
    phase_parameters, phase_auxiliary = resonance.phase_setup(
        model=model, resonances=phase_resonances, mode=phi["mode"])
    parameters.update(phase_parameters)
    auxiliary.update(phase_auxiliary)
    if continuum_options:
        options = {"REGGE": False, "FF_offshell": False, "FF_transfer": False,
            "pveto": False, "tune_couplings": False, "tune_reggeize": False, **continuum_options, }
        con_parameters, con_auxiliary = continuum.setup(continuum_models=model,
            final_state_pdgs=continuum.discover_final_states_from_datacards(datacards), **options, )
        parameters.update(con_parameters)
        auxiliary.update(con_auxiliary)
    if spin_basis is not None:
        _, basis_auxiliary = resonance.spin_basis(resonances=discovered, basis=spin_basis)
        auxiliary.update(basis_auxiliary)
    unknown = set(res_fields) - {"FF_prod"}
    if unknown:
        raise ValueError(f"Unknown RES MODELS fields: {sorted(unknown)}")
    if "FF_prod" in res_fields:
        parameters.update(resonance.res_form_factors(fields={"FF_prod"}, model=model, resonances=selected))
    return datacards, parameters, auxiliary


# Build and combine model selections, requiring identical shared parameter definitions
def build(config):
    if not isinstance(config["models"], list) or not config["models"]:
        raise ValueError("GRANIITTI tuning setup requires a non-empty models array")
    datacards, parameters, fixed = [], {}, {}
    for model in config["models"]:
        options = copy.deepcopy(model)
        datasets = [card(item["datacard"], **{k: v for k, v in item.items() if k != "datacard"})
                    for item in options.pop("datasets")]
        form_factors = options.pop("form_factors", None)
        if options["model"] in {"single", "double", "triple"}:
            varied, auxiliary = eikonal.setup(**options)
        else:
            _, varied, auxiliary = central(datacards=datasets, **options)
        if form_factors:
            varied.update(resonance.res_form_factors(model=options["model"], **form_factors))
        previous, current = normalize_param_space(parameters), normalize_param_space(varied)
        conflicts = {key for key in previous.keys() & current.keys() if previous[key] != current[key]}
        conflicts |= {key for key in fixed.keys() & auxiliary.keys() if fixed[key] != auxiliary[key]}
        conflicts |= (parameters.keys() & auxiliary.keys()) | (varied.keys() & fixed.keys())
        if conflicts:
            raise ValueError(f"Conflicting shared fit parameters: {sorted(conflicts)}")
        datacards.extend(datasets)
        parameters.update(varied)
        fixed.update(auxiliary)
    return datacards, parameters, fixed
