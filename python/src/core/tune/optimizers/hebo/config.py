# Validated HEBO JSON configuration loading
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import math
from pathlib import Path

from core import resource
from core.io.serialize import load_json_file

DEFAULT_HEBO_CONFIG = resource("tune/settings/hebo.json")
HEBO_SECTIONS = {
    "output": {"transform"},
    "gp": {
        "lr",
        "num_epochs",
        "optimizer",
        "noise_lb",
        "noise_guess",
        "pred_likeli",
        "ard_kernel",
        "max_exact_trials",
        "verbose",
    },
    "sparse_gp": {"lr", "lr_vp", "num_epochs", "batch_size", "num_inducing", "learn_u"},
    "topology": {"enabled", "nu", "lengthscale_scale", "outputscale_prior_shape", "outputscale_prior_rate"},
    "acquisition": {"population", "generations", "eps", "kappa_scale", "delta"},
    "pending": {"lie", "penalty_scale"},
}


# Require a finite real HEBO setting above its lower limit
def _real(value, name: str, lower: float = 0.0, *, inclusive: bool = False) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"HEBO {name} must be finite and numeric")
    if value < lower or (not inclusive and value <= lower):
        raise ValueError(f"HEBO {name} must be {'>=' if inclusive else '>'} {lower}")


# Require an integer HEBO setting without accepting Boolean values
def _integer(value, name: str, minimum: int = 1) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"HEBO {name} must be an integer >= {minimum}")


# Validate complete HEBO settings before starting any trial
def validate_hebo_config(payload: object) -> dict:
    if (
        not isinstance(payload, dict)
        or type(payload.get("schema_version")) is not int
        or payload["schema_version"] != 1
    ):
        raise ValueError("Unsupported HEBO configuration schema")
    if set(payload) != {"schema_version", *HEBO_SECTIONS}:
        raise ValueError(f"HEBO configuration sections must be exactly {sorted(HEBO_SECTIONS)}")
    for section, keys in HEBO_SECTIONS.items():
        if not isinstance(payload[section], dict) or set(payload[section]) != keys:
            raise ValueError(f"HEBO {section} keys must be exactly {sorted(keys)}")
    for section, keys in {
        "gp": ("lr", "noise_lb", "noise_guess"),
        "sparse_gp": ("lr", "lr_vp"),
        "topology": ("lengthscale_scale", "outputscale_prior_shape", "outputscale_prior_rate"),
        "acquisition": ("kappa_scale", "delta"),
    }.items():
        for key in keys:
            _real(payload[section][key], f"{section}.{key}")
    for section, keys in {
        "gp": ("num_epochs",),
        "sparse_gp": ("num_epochs", "batch_size", "num_inducing"),
        "acquisition": ("generations",),
    }.items():
        for key in keys:
            _integer(payload[section][key], f"{section}.{key}")
    _integer(payload["gp"]["max_exact_trials"], "gp.max_exact_trials", 0)
    _integer(payload["acquisition"]["population"], "acquisition.population", 3)
    _real(payload["pending"]["penalty_scale"], "pending.penalty_scale", inclusive=True)
    _real(payload["acquisition"]["eps"], "acquisition.eps", inclusive=True)
    for section, keys in {
        "gp": ("pred_likeli", "ard_kernel", "verbose"),
        "sparse_gp": ("learn_u",),
        "topology": ("enabled",),
    }.items():
        for key in keys:
            if not isinstance(payload[section][key], bool):
                raise ValueError(f"HEBO {section}.{key} must be Boolean")
    for section, key, choices in (
        ("output", "transform", ("power", "standardize")),
        ("gp", "optimizer", ("psgld", "adam", "lbfgs")),
        ("pending", "lie", ("worst", "mean", "best")),
    ):
        if payload[section][key] not in choices:
            raise ValueError(f"HEBO {section}.{key} must be one of {choices}")
    _real(payload["topology"]["nu"], "topology.nu")
    if not any(math.isclose(payload["topology"]["nu"], nu, rel_tol=0.0, abs_tol=1e-12) for nu in (1.5, 2.5)):
        raise ValueError("HEBO topology.nu must be 1.5 or 2.5")
    if payload["acquisition"]["delta"] >= 1.0:
        raise ValueError("HEBO acquisition.delta must be < 1")
    if payload["gp"]["noise_guess"] <= payload["gp"]["noise_lb"]:
        raise ValueError("HEBO gp.noise_guess must exceed gp.noise_lb")
    return copy.deepcopy(payload)


# Load the shared HEBO card or an explicitly selected configuration
def load_hebo_config(path: str | Path = DEFAULT_HEBO_CONFIG) -> dict:
    return validate_hebo_config(load_json_file(path))
