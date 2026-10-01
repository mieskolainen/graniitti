# TuRBO state transitions and anisotropic trust region construction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from dataclasses import dataclass

import torch

from core.tune.parameters.kernel import ProductManifoldGeometry


@dataclass
class TurboState:
    """Canonical TuRBO-1 state for scalar minimization"""

    dimension: int
    length_initial: float
    length_min: float
    length_max: float
    success_tolerance: int
    failure_base: float
    failure_dimension_scale: float
    success_length_multiplier: float
    failure_length_multiplier: float
    improvement_relative_tolerance: float
    batch_size: int = 1
    length: float = math.nan
    success_counter: int = 0
    failure_counter: int = 0
    best_value: float = math.inf
    restart_triggered: bool = False
    restart_count: int = 0

    # Build canonical TuRBO state from the authoritative nested settings
    @classmethod
    def from_settings(cls, dimension: int, settings: dict) -> TurboState:
        return cls(
            dimension=int(dimension),
            length_initial=float(settings["length_initial"]),
            length=float(settings["length_initial"]),
            length_min=float(settings["length_min"]),
            length_max=float(settings["length_max"]),
            success_tolerance=int(settings["success_tolerance"]),
            failure_base=float(settings["failure_base"]),
            failure_dimension_scale=float(settings["failure_dimension_scale"]),
            success_length_multiplier=float(settings["success_length_multiplier"]),
            failure_length_multiplier=float(settings["failure_length_multiplier"]),
            improvement_relative_tolerance=float(settings["improvement_relative_tolerance"]),
        )

    # Compute the canonical batch-aware TuRBO failure tolerance
    @property
    def failure_tolerance(self) -> int:
        return int(
            math.ceil(
                max(
                    self.failure_base / max(1, self.batch_size),
                    self.failure_dimension_scale * self.dimension / max(1, self.batch_size),
                )
            )
        )

    # Update trust region counters from one completed evaluation batch
    def update(self, values: torch.Tensor) -> None:
        finite = values[torch.isfinite(values)]
        if len(finite) == 0:
            return
        batch_best = float(torch.min(finite).detach().cpu())
        tolerance = self.improvement_relative_tolerance * abs(self.best_value)
        improved = not math.isfinite(self.best_value) or batch_best < self.best_value - tolerance
        if improved:
            self.success_counter += 1
            self.failure_counter = 0
        else:
            self.success_counter = 0
            self.failure_counter += 1
        self.best_value = min(self.best_value, batch_best)
        if self.success_counter >= self.success_tolerance:
            self.length = min(self.success_length_multiplier * self.length, self.length_max)
            self.success_counter = 0
        elif self.failure_counter >= self.failure_tolerance:
            self.length *= self.failure_length_multiplier
            self.failure_counter = 0
        self.restart_triggered = self.length < self.length_min

    # Reset TuRBO state before a new global Sobol restart design
    def restart(self) -> None:
        self.length = self.length_initial
        self.success_counter = 0
        self.failure_counter = 0
        self.best_value = math.inf
        self.restart_triggered = False
        self.restart_count += 1

    # Compute a JSON safe TuRBO state summary
    def diagnostics(self) -> dict:
        return {
            "dimension": self.dimension,
            "batch_size": self.batch_size,
            "length": self.length,
            "length_initial": self.length_initial,
            "length_min": self.length_min,
            "length_max": self.length_max,
            "success_counter": self.success_counter,
            "failure_counter": self.failure_counter,
            "success_tolerance": self.success_tolerance,
            "failure_tolerance": self.failure_tolerance,
            "failure_base": self.failure_base,
            "failure_dimension_scale": self.failure_dimension_scale,
            "success_length_multiplier": self.success_length_multiplier,
            "failure_length_multiplier": self.failure_length_multiplier,
            "improvement_relative_tolerance": self.improvement_relative_tolerance,
            "best_value": self.best_value if math.isfinite(self.best_value) else None,
            "restart_triggered": self.restart_triggered,
            "restart_count": self.restart_count,
        }


# Construct the canonical length scale weighted TuRBO hyperrectangle
def turbo_bounds(
    center: torch.Tensor, lengthscale: torch.Tensor, geometry: ProductManifoldGeometry, length: float
) -> tuple[torch.Tensor, torch.Tensor]:
    weights = geometry.coordinate_weights(lengthscale)
    lower = center - 0.5 * float(length) * weights
    upper = center + 0.5 * float(length) * weights
    periodic = torch.zeros(geometry.dimension, dtype=torch.bool, device=center.device)
    periodic[geometry.periodic_indices] = True
    lower = torch.where(periodic, lower, torch.clamp(lower, min=0.0, max=1.0))
    upper = torch.where(periodic, upper, torch.clamp(upper, min=0.0, max=1.0))
    return lower, upper
