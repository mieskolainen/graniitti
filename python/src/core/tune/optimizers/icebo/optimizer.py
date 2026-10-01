# Standalone exact Bayesian optimization
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import hashlib
import json
import math
import pathlib
from collections.abc import Sequence

import numpy as np
import torch

from core.inference.gp import ExactGaussianProcess
from core.tune.optimizers.icebo.config import validate_icebo_config
from core.tune.optimizers.icebo.output import StandardizeOutput
from core.tune.optimizers.icebo.qlognei import QLogNEIOptimizer, QLogNEISettings
from core.tune.optimizers.icebo.space import IceboSpace, resolve_icebo_device, resolve_icebo_dtype
from core.tune.optimizers.icebo.turbo import TurboState, turbo_bounds


class ICEBO:
    """Torch based optimizer for noisy mixed-domain black-box functions"""

    # Initialize the search policy and empty observation state
    def __init__(
        self,
        bounds: dict,
        *,
        settings: dict,
        seed: int,
        parameter_topology: dict | None = None,
        warmup: int | None = None,
    ) -> None:
        config = validate_icebo_config(copy.deepcopy(settings), source=pathlib.Path("<ICEBO settings>"))
        search = config["search"]
        acquisition = config["acquisition"]
        turbo = config["turbo"]
        self.device = resolve_icebo_device(search["device"])
        self.dtype = resolve_icebo_dtype(self.device, search["dtype"])
        self.space = IceboSpace(bounds, parameter_topology=parameter_topology, device=self.device, dtype=self.dtype)
        self.direction = search["direction"]
        self.duplicate_tolerance = search["duplicate_tolerance"]
        self.seed = int(seed)
        self.warmup = max(4, 2 * self.space.dimension + 1) if warmup is None else max(0, int(warmup))
        self.gp_settings = copy.deepcopy(config["gp"])
        # Scale the lengthscale prior with the kernel dimension [REFERENCE: arXiv:2402.02229]
        for name in ("initial_lengthscale", "lengthscale_prior_median"):
            self.gp_settings[name] *= math.sqrt(self.space.geometry.effective_dim)
        self.acquisition_settings = acquisition
        self.turbo_settings = turbo
        self.trust_region = turbo["enabled"]
        self.settings = config
        self._physical = torch.empty((0, self.space.dimension), device=self.device, dtype=self.dtype)
        self._values = torch.empty((0,), device=self.device, dtype=self.dtype)
        self._errors = torch.empty((0,), device=self.device, dtype=self.dtype)
        self._pending = torch.empty((0, self.space.dimension), device=self.device, dtype=self.dtype)
        self._excluded = self._pending.clone()
        self._surrogate = None
        self._transform = None
        self._fit_stats = {}
        self._fit_state = None
        self._dirty = True
        self._proposal_count = 0
        self._turbo = TurboState.from_settings(self.space.dimension, turbo)
        self._restart_remaining = 0
        self._turbo_active_start = 0
        self._turbo_center = None
        self._turbo_updates = []
        self._turbo_batches = set()

    # Compute the number of finite observations retained by the optimizer
    @property
    def observation_count(self) -> int:
        return len(self._values)

    # Observe one or more completed objective evaluations
    def observe(
        self,
        configs: dict | Sequence[dict],
        values: float | Sequence[float] | np.ndarray,
        errors: float | Sequence[float] | np.ndarray | None = None,
        *,
        update_turbo: bool = True,
    ) -> None:
        physical = self.space.configs_to_tensor(configs)
        objective = torch.as_tensor(values, device=self.device, dtype=self.dtype).reshape(-1)
        if len(physical) != len(objective):
            raise ValueError("ICEBO configurations and objective values must have equal length")
        uncertainty = self._observation_errors(errors, len(objective))
        self._remove_observed_pending(physical)
        finite = torch.isfinite(objective) & torch.isfinite(uncertainty) & (uncertainty >= 0.0)
        self._excluded = self._merge_physical(self._excluded, physical[~finite])
        if not bool(torch.any(finite)):
            return
        adaptive = self.observation_count - self._turbo_active_start >= self._design_size()
        physical = physical[finite]
        objective = objective[finite]
        uncertainty = uncertainty[finite]
        if self.direction == "maximize":
            objective = -objective
        self._physical = torch.cat((self._physical, physical), dim=0)
        self._values = torch.cat((self._values, objective), dim=0)
        self._errors = torch.cat((self._errors, uncertainty), dim=0)
        if update_turbo and adaptive:
            self._schedule_turbo_update(physical)
        self._dirty = True

    # Compute aligned non-negative observation errors
    def _observation_errors(self, errors: float | Sequence[float] | np.ndarray | None, count: int) -> torch.Tensor:
        if errors is None:
            return torch.zeros(count, device=self.device, dtype=self.dtype)
        values = torch.as_tensor(errors, device=self.device, dtype=self.dtype).reshape(-1)
        if len(values) == 1 and count > 1:
            values = values.expand(count)
        if len(values) != count:
            raise ValueError("ICEBO observation errors must be scalar or aligned with values")
        return values

    # Schedule one posterior-mean TuRBO update after the next GP fit
    def _schedule_turbo_update(self, physical: torch.Tensor, batch_id: str | None = None) -> None:
        if not self.trust_region:
            return
        self._turbo_updates.append((physical.detach().clone(), batch_id))

    # Update TuRBO from denoised posterior means after fitting all new observations
    def _apply_turbo_update(self) -> None:
        if not self.trust_region:
            return
        for physical, batch_id in self._turbo_updates:
            if self._turbo.restart_triggered:
                break
            if len(physical) == 0 or batch_id is not None and batch_id in self._turbo_batches:
                continue
            with torch.no_grad():
                points = physical if self._turbo_center is None else torch.cat((self._turbo_center, physical))
                unit = self.space.to_unit(points)
                transformed = self._surrogate.predict_mean(self.space.surrogate_inputs_from_unit(unit)).flatten()
                objective = self._transform.inverse(transformed)
                if self._turbo_center is not None:
                    self._turbo.best_value = float(objective[0])
                    objective = objective[1:]
                best = int(torch.argmin(objective))
                if float(objective[best]) < self._turbo.best_value:
                    self._turbo_center = physical[best : best + 1].detach().clone()
            self._turbo.batch_size = len(objective)
            self._turbo.update(objective)
        self._turbo_batches.update(batch_id for _, batch_id in self._turbo_updates if batch_id is not None)
        self._turbo_updates.clear()
        if self._turbo.restart_triggered and len(self._pending) == 0:
            self._turbo.restart()
            self._restart_remaining = self._design_size()
            self._turbo_active_start = self.observation_count
            self._turbo_center = None
            self._dirty = True

    # Compute the required number of completed initial observations in this region
    def _design_size(self) -> int:
        if not self.trust_region or self._turbo.restart_count == 0:
            return max(1, self.warmup)
        return max(self.turbo_settings["restart_min_points"],
                   self.turbo_settings["restart_points_per_dimension"] * self.space.dimension)

    # Update TuRBO once after an externally scheduled proposal batch completes
    def update_turbo_batch(self, configs: Sequence[dict], *, batch_id: str | None = None) -> None:
        physical = self.space.configs_to_tensor(configs)
        if len(physical) == 0:
            return
        self._schedule_turbo_update(physical, batch_id)

    # Remove completed configurations from internal pending reservations
    def _remove_observed_pending(self, observed: torch.Tensor) -> None:
        if len(self._pending) == 0:
            return
        pending_unit = self.space.to_unit(self._pending)
        observed_unit = self.space.to_unit(observed)
        distance = self.space.minimum_topology_distance(pending_unit, observed_unit)
        keep = distance > self.duplicate_tolerance
        self._pending = self._pending[keep]

    # Replace internal pending reservations from an external scheduler
    def set_pending(self, configs: Sequence[dict] | None) -> None:
        self._pending = self.space.configs_to_tensor([] if configs is None else configs)

    # Exclude failed configurations without fantasizing pending objective values
    def set_excluded(self, configs: Sequence[dict]) -> None:
        self._excluded = self.space.configs_to_tensor(configs)

    # Add scheduler reservations without discarding existing pending points
    def reserve(self, configs: dict | Sequence[dict]) -> None:
        physical = self.space.configs_to_tensor(configs)
        self._pending = self._merge_physical(self._pending, physical)

    # Suggest one or more pending-aware configurations
    def suggest(
        self, n_suggestions: int = 1, *, pending: Sequence[dict] | None = None, return_diagnostics: bool = False
    ):
        count = max(1, int(n_suggestions))
        external = self.space.configs_to_tensor([] if pending is None else pending)
        pending_physical = self._merge_pending(external)
        self._pending = pending_physical
        if self._turbo.restart_triggered:
            self._apply_turbo_update()
            if self._turbo.restart_triggered:
                return ([], []) if return_diagnostics else []
        if self.observation_count - self._turbo_active_start >= self._design_size():
            self._fit_if_needed()
        if self._turbo.restart_triggered:
            return ([], []) if return_diagnostics else []
        missing = self._design_size() - (self.observation_count - self._turbo_active_start)
        if missing > 0:
            requested = min(count, max(0, missing - len(pending_physical)))
            configs, diagnostics = self._warmup_suggestions(requested, pending_physical) if requested else ([], [])
            self._restart_remaining = max(0, self._restart_remaining - len(configs))
        else:
            configs, diagnostics = self._model_suggestions(count, pending_physical)
        reserved = self.space.configs_to_tensor(configs)
        self._pending = self._merge_physical(self._pending, reserved)
        self._proposal_count += len(configs)
        return (configs, diagnostics) if return_diagnostics else configs

    # Merge external reservations with internal pending points
    def _merge_pending(self, external: torch.Tensor) -> torch.Tensor:
        return self._merge_physical(self._pending, external)

    # Merge physical rows while removing exact typed duplicates
    def _merge_physical(self, first: torch.Tensor, second: torch.Tensor) -> torch.Tensor:
        if len(first) == 0:
            return second.clone()
        if len(second) == 0:
            return first.clone()
        combined = self.space.canonical_unit(self.space.to_unit(torch.cat((first, second), dim=0)))
        selected = []
        references = None
        for row in combined:
            candidate = row.unsqueeze(0)
            distance = self.space.minimum_topology_distance(candidate, references)
            if float(distance[0].cpu()) <= self.duplicate_tolerance:
                continue
            selected.append(row)
            references = candidate if references is None else torch.cat((references, candidate), dim=0)
        unique = torch.stack(selected)
        return self.space.from_unit(unique, round_integers=True)

    # Draw distinct Sobol warm-up suggestions
    def _warmup_suggestions(self, count: int, pending_physical: torch.Tensor) -> tuple[list[dict], list[dict]]:
        references = self._reference_units(pending_physical)
        selected = []
        diagnostics = []
        attempt = 0
        while len(selected) < count and attempt < 8:
            candidates = self.space.sobol(max(2 * count, 32), seed=self.seed + self._proposal_count + 104729 * attempt)
            distance = self.space.minimum_topology_distance(candidates, references)
            for index in torch.argsort(distance, descending=True).tolist():
                candidate = candidates[index : index + 1]
                current_distance = self.space.minimum_topology_distance(candidate, references)
                if float(current_distance[0].cpu()) <= self.duplicate_tolerance:
                    continue
                selected.append(candidate[0])
                references = candidate if references is None else torch.cat((references, candidate), dim=0)
                diagnostics.append(
                    {"acquisition": "sobol_warmup", "proposal_index": self._proposal_count + len(selected) - 1}
                )
                if len(selected) == count:
                    break
            attempt += 1
        if len(selected) < count:
            raise RuntimeError("ICEBO could not find enough distinct warm-up configurations")
        physical = self.space.from_unit(torch.stack(selected), round_integers=True)
        return self.space.tensor_to_configs(physical), diagnostics

    # Compute observed and pending unit coordinates
    def _reference_units(self, pending_physical: torch.Tensor) -> torch.Tensor | None:
        observed = self.space.to_unit(self._physical)
        pending = self.space.to_unit(pending_physical)
        excluded = self.space.to_unit(self._excluded)
        if len(observed) + len(pending) + len(excluded) == 0:
            return None
        return torch.unique(torch.cat((observed, pending, excluded), dim=0), dim=0)

    # Select batches greedily with pending-conditioned qLogNEI [REFERENCE: Wilson et al., NeurIPS 2018, arXiv:1805.10196]
    def _model_suggestions(self, count: int, pending_physical: torch.Tensor) -> tuple[list[dict], list[dict]]:
        self._fit_if_needed()
        pending_unit = self.space.to_unit(pending_physical)
        references = self._reference_units(pending_physical)
        if self.trust_region:
            center = self._trust_center(self._surrogate)
            lower, upper = turbo_bounds(center, self._surrogate.lengthscale, self.space.geometry, self._turbo.length)
        else:
            lower = torch.zeros(self.space.dimension, device=self.device, dtype=self.dtype)
            upper = torch.ones(self.space.dimension, device=self.device, dtype=self.dtype)
        maximizer = QLogNEIOptimizer(
            self.space, settings=self._acquisition_settings(), seed=self.seed + self._proposal_count
        )
        selected, diagnostics = [], []
        # Condition each greedy addition on the pending points already selected
        for batch_index in range(count):
            candidate, common = maximizer.maximize(
                self._surrogate, q=1, pending=pending_unit, references=references,
                lower=lower, upper=upper, iteration=self._proposal_count + batch_index,
            )
            selected.append(candidate)
            pending_unit = torch.cat((pending_unit, candidate), dim=0)
            references = torch.cat((references, candidate), dim=0)
            diagnostics.append(
                {
                    **common,
                    "batch_index": batch_index,
                    "proposal_index": self._proposal_count + batch_index,
                    "device": str(self.device),
                    "dtype": str(self.dtype),
                    "surrogate": self._surrogate.diagnostics()["surrogate"],
                    "trust_region_length": self._turbo.length if self.trust_region else None,
                }
            )
        physical = self.space.from_unit(torch.cat(selected), round_integers=True)
        return self.space.tensor_to_configs(physical), diagnostics

    # Compute the posterior-mean TuRBO center for noisy observations
    def _trust_center(self, surrogate: ExactGaussianProcess) -> torch.Tensor:
        with torch.no_grad():
            mean = surrogate.predict_mean(surrogate.X).flatten()
            index = int(torch.argmin(mean).cpu())
            self._turbo.best_value = float(self._transform.inverse(mean[index]).cpu())
            self._turbo_center = surrogate.X[index : index + 1].detach().clone()
        return self.space.to_unit(self._turbo_center)[0]

    # Build the complete qLogNEI numerical configuration
    def _acquisition_settings(self) -> QLogNEISettings:
        return QLogNEISettings(**self.acquisition_settings, duplicate_tolerance=self.duplicate_tolerance)

    # Fit the exact GP to observations from the active trust region run
    def _fit_if_needed(self) -> None:
        if not self._dirty and self._surrogate is not None:
            self._apply_turbo_update()
            return
        active = self._turbo_active_start if self.trust_region else 0
        if self.observation_count <= active:
            raise RuntimeError("ICEBO requires observations from the active trust region for prediction")
        data = torch.cat((self._physical[active:], self._values[active:, None], self._errors[active:, None]), dim=1)
        order = torch.arange(len(data), device=self.device)
        for column in reversed(range(data.shape[1])):
            order = order[torch.argsort(data[order, column], stable=True)]
        data = data[order]
        checksum = hashlib.sha256(data.detach().cpu().numpy().tobytes())
        checksum.update(json.dumps(self.space.metadata(), sort_keys=True).encode("ascii"))
        digest = checksum.hexdigest()
        unit = self.space.to_unit(data[:, :-2])
        inputs = self.space.surrogate_inputs_from_unit(unit)
        transform = StandardizeOutput()
        transform.fit(data[:, -2])
        target = transform.forward(data[:, -2])
        transformed_errors = data[:, -1] / transform.scale
        surrogate = ExactGaussianProcess(
            inputs,
            target,
            Y_err=transformed_errors,
            geometry=self.space.geometry,
            settings=self.gp_settings,
            device=self.device,
            dtype=self.dtype,
            seed=self.seed + self.observation_count,
        )
        saved = self._fit_state
        compatible = (saved is not None and saved["statistics"]["training_start"] == active
                      and saved["statistics"]["settings"] == self.gp_settings)
        if compatible:
            surrogate._restore_parameters({name: torch.as_tensor(value, device=self.device, dtype=self.dtype)
                                           for name, value in saved["parameters"].items()})
        if compatible and saved["data_digest"] == digest:
            fit_stats = {**saved["statistics"], "elapsed_seconds": 0.0}
            with torch.no_grad():
                surrogate._refresh()
        else:
            fit_stats = surrogate.fit(warm_start=compatible)
        surrogate.freeze_for_inference()
        self._surrogate = surrogate
        self._transform = transform
        self._fit_stats = {**fit_stats, **surrogate.diagnostics(), "training_start": active}
        self._fit_state = {
            "data_digest": digest,
            "parameters": {name: value.detach().cpu().tolist() for name, value in surrogate.named_parameters()},
            "statistics": {key: copy.deepcopy(value) for key, value in self._fit_stats.items() if key != "elapsed_seconds"},
        }
        self._dirty = False
        self._apply_turbo_update()

    # Predict objective moments for arbitrary configurations
    def predict(self, configs: dict | Sequence[dict]) -> dict:
        if self.observation_count < 2:
            raise RuntimeError("ICEBO requires at least two observations for prediction")
        self._fit_if_needed()
        physical = self.space.configs_to_tensor(configs)
        unit = self.space.to_unit(physical)
        with torch.no_grad():
            posterior = self._surrogate.posterior(self.space.surrogate_inputs_from_unit(unit))
            mean = self._transform.inverse(posterior.mean)
            epistemic = self._transform.variance_inverse(posterior.epistemic_variance)
            aleatoric = self._transform.variance_inverse(posterior.aleatoric_variance)
            if self.direction == "maximize":
                mean = -mean
        return {
            "mean": mean.detach().cpu().numpy(),
            "epistemic_std": torch.sqrt(epistemic).detach().cpu().numpy(),
            "aleatoric_std": torch.sqrt(aleatoric).detach().cpu().numpy(),
        }

    # Compute JSON-safe optimizer and fitted-surrogate diagnostics
    def diagnostics(self) -> dict:
        return {
            "name": "icebo",
            "acquisition": "qlognei",
            "direction": self.direction,
            "device": str(self.device),
            "dtype": str(self.dtype),
            "observations": self.observation_count,
            "pending": len(self._pending),
            "proposal_count": self._proposal_count,
            "warmup": self.warmup,
            "trust_region": self.trust_region,
            "turbo": self._turbo.diagnostics(),
            "restart_design_remaining": self._restart_remaining,
            "turbo_active_start": self._turbo_active_start,
            "settings": copy.deepcopy(self.settings),
            "space": self.space.metadata(),
            "output_transform": None if self._transform is None else self._transform.metadata(),
            "fit": copy.deepcopy(self._fit_stats),
        }

    # Compute replayable non-model search state
    def state_dict(self) -> dict:
        return {
            "schema_version": 1,
            "dimension": self.space.dimension,
            "proposal_count": self._proposal_count,
            "restart_design_remaining": self._restart_remaining,
            "turbo_active_start": self._turbo_active_start,
            "turbo": self._turbo.diagnostics(),
            "turbo_batches": sorted(self._turbo_batches),
            "model_state": copy.deepcopy(self._fit_state),
            "turbo_center": None if self._turbo_center is None else self._turbo_center.detach().cpu().tolist(),
        }

    # Restore replayable search state after rebuilding observations
    def load_state_dict(self, state: dict) -> None:
        if int(state.get("schema_version", -1)) != 1:
            raise ValueError("Unsupported ICEBO state schema")
        if int(state.get("dimension", -1)) != self.space.dimension:
            raise ValueError("ICEBO state dimension does not match the parameter space")
        self._proposal_count = max(0, int(state.get("proposal_count", 0)))
        self._restart_remaining = max(0, int(state.get("restart_design_remaining", 0)))
        self._turbo_active_start = max(0, int(state.get("turbo_active_start", 0)))
        self._turbo_batches = set(state.get("turbo_batches", []))
        self._fit_state = copy.deepcopy(state.get("model_state"))
        center = state.get("turbo_center")
        self._turbo_center = None if center is None else torch.as_tensor(center, device=self.device, dtype=self.dtype)
        self._dirty = True
        turbo = state.get("turbo", {})
        if int(turbo.get("dimension", self.space.dimension)) != self.space.dimension:
            raise ValueError("ICEBO TuRBO state dimension does not match the parameter space")
        for name in (
            "batch_size",
            "length",
            "length_min",
            "length_max",
            "success_counter",
            "failure_counter",
            "success_tolerance",
            "best_value",
            "restart_triggered",
            "restart_count",
        ):
            if name in turbo and turbo[name] is not None:
                setattr(self._turbo, name, turbo[name])
        if turbo.get("best_value") is None:
            self._turbo.best_value = math.inf

    # Restore the deterministic proposal sequence counter during history replay
    def set_replay_proposal_count(self, count: int) -> None:
        self._proposal_count = max(self._proposal_count, int(count))
