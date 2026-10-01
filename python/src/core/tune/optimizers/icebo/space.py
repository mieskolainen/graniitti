# Typed parameter space for the standalone ICEBO optimizer
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from collections.abc import Sequence

import torch

from core.tune.parameters.kernel import ProductManifoldGeometry
from core.tune.parameters.space import is_integer_bound, normalize_param_space, typed_config
from core.tune.parameters.topology import indexed_parameter_topology


class IceboSpace:
    """Finite mixed numeric parameter space with differentiable topology features"""

    # Initialize ordered physical bounds and Torch representations
    def __init__(
        self, bounds: dict, *, parameter_topology: dict | None, device: torch.device, dtype: torch.dtype
    ) -> None:
        self.bounds = normalize_param_space(bounds)
        if not self.bounds:
            raise ValueError("ICEBO requires at least one parameter")
        self.names = sorted(self.bounds)
        self.parameter_topology = parameter_topology or {}
        indexed_parameter_topology(self.names, self.parameter_topology)
        limits = [[float(self.bounds[name]["lower"]), float(self.bounds[name]["upper"])] for name in self.names]
        self.limits = torch.tensor(limits, device=device, dtype=dtype)
        self.integer_mask = torch.tensor(
            [is_integer_bound(self.bounds[name]) for name in self.names], device=device, dtype=torch.bool
        )
        groups = list((self.parameter_topology or {}).get("groups", []))
        periodic_names = set()
        for group in groups:
            kind = str(group.get("kind"))
            parameters = [str(name) for name in group["parameters"]]
            if kind in {"phase", "phase_projective", "polar_components"}:
                periodic_names.add(parameters[0])
            elif kind == "polar_projective":
                periodic_names.add(parameters[1])
            elif kind == "polar_sphere":
                periodic_names.update((parameters[1], parameters[-1]))
            elif kind == "sphere":
                periodic_names.add(parameters[-1])
        self.periodic_mask = torch.tensor(
            [name in periodic_names for name in self.names], device=device, dtype=torch.bool
        )
        self.direction_groups = []
        for group in groups:
            kind = str(group.get("kind"))
            parameters = [str(name) for name in group["parameters"]]
            if kind in {"projective", "sphere"}:
                direction_kind = kind
                direction = parameters
            elif kind in {"polar_projective", "cartesian_projective", "phase_projective"}:
                direction_kind = "projective"
                direction = parameters[1:] if kind == "phase_projective" else parameters[2:]
            elif kind == "polar_sphere":
                direction_kind = "sphere"
                direction = parameters[2:]
            elif kind == "radial_projective":
                direction_kind = "projective"
                direction = parameters[1:]
            else:
                continue
            if direction:
                self.direction_groups.append(
                    {"kind": direction_kind, "indices": [self.names.index(name) for name in direction]}
                )
        self.device = device
        self.dtype = dtype
        self.geometry = ProductManifoldGeometry.from_metadata(
            self.names, tuple((float(row[0]), float(row[1])) for row in limits), self.parameter_topology
        ).to(device=device, dtype=dtype)

    # Compute the number of physical optimizer coordinates
    @property
    def dimension(self) -> int:
        return len(self.names)

    # Convert one or more configuration dictionaries to physical rows
    def configs_to_tensor(self, configs: dict | Sequence[dict]) -> torch.Tensor:
        rows = [configs] if isinstance(configs, dict) else list(configs)
        if not rows:
            return torch.empty((0, self.dimension), device=self.device, dtype=self.dtype)
        typed = [typed_config(row, self.bounds) for row in rows]
        return torch.tensor(
            [[float(row[name]) for name in self.names] for row in typed], device=self.device, dtype=self.dtype
        )

    # Convert physical rows to typed configuration dictionaries
    def tensor_to_configs(self, values: torch.Tensor) -> list[dict]:
        rows = torch.atleast_2d(values).detach().cpu().numpy()
        output = []
        for row in rows:
            output.append(typed_config({name: float(row[index]) for index, name in enumerate(self.names)}, self.bounds))
        return output

    # Map physical rows into the unit hypercube
    def to_unit(self, physical: torch.Tensor) -> torch.Tensor:
        lower = self.limits[:, 0]
        width = self.limits[:, 1] - lower
        return (physical.to(device=self.device, dtype=self.dtype) - lower) / width

    # Map unit-hypercube rows into physical coordinates
    def from_unit(self, unit: torch.Tensor, *, round_integers: bool) -> torch.Tensor:
        rows = unit.to(device=self.device, dtype=self.dtype)
        rows = torch.where(self.periodic_mask, torch.remainder(rows, 1.0), torch.clamp(rows, 0.0, 1.0))
        physical = self.limits[:, 0] + rows * (self.limits[:, 1] - self.limits[:, 0])
        if round_integers and bool(torch.any(self.integer_mask)):
            physical = torch.where(self.integer_mask, torch.round(physical), physical)
        return torch.maximum(torch.minimum(physical, self.limits[:, 1]), self.limits[:, 0])

    # Compute differentiable physical inputs for the direct topology kernel
    def surrogate_inputs_from_unit(self, unit: torch.Tensor) -> torch.Tensor:
        return self.from_unit(unit, round_integers=False)

    # Canonicalize unit coordinates after integer rounding
    def canonical_unit(self, unit: torch.Tensor) -> torch.Tensor:
        return self.to_unit(self.from_unit(unit, round_integers=True)).clamp(0.0, 1.0)

    # Encode unit directions into their spherical or projective angle coordinates
    def _direction_angles(self, vectors: torch.Tensor, kind: str) -> torch.Tensor:
        values = vectors
        if kind == "projective":
            nonzero = torch.abs(values) > torch.finfo(values.dtype).eps
            pivot = torch.argmax(nonzero.to(torch.int64), dim=1)
            sign = torch.sign(values.gather(1, pivot[:, None]))
            values = values * torch.where(sign == 0.0, torch.ones_like(sign), sign)
        angles = []
        count = values.shape[1] - (2 if kind == "sphere" else 1)
        for index in range(count):
            radius = torch.sqrt(values[:, 0].square() + values[:, index + 1 :].square().sum(dim=1)).clamp_min(
                torch.finfo(values.dtype).tiny
            )
            ratio = torch.clamp(values[:, index + 1] / radius, -1.0, 1.0)
            angles.append(torch.asin(ratio))
        if kind == "sphere":
            angles.append(torch.atan2(values[:, -1], values[:, 0]))
        return torch.stack(angles, dim=1)

    # Draw deterministic scrambled Sobol points with uniform direction groups
    def sobol(self, count: int, *, seed: int, topology_uniform: bool = True) -> torch.Tensor:
        if count < 1:
            return torch.empty((0, self.dimension), device=self.device, dtype=self.dtype)
        extra = len(self.direction_groups) if topology_uniform else 0
        engine = torch.quasirandom.SobolEngine(self.dimension + extra, scramble=True, seed=int(seed))
        latent = engine.draw(int(count), dtype=self.dtype).to(device=self.device)
        points = latent[:, : self.dimension]
        if topology_uniform:
            epsilon = torch.finfo(self.dtype).eps
            for offset, group in enumerate(self.direction_groups):
                indices = group["indices"]
                uniforms = torch.cat(
                    (points[:, indices], latent[:, self.dimension + offset : self.dimension + offset + 1]), dim=1
                )
                arguments = torch.clamp(2.0 * uniforms - 1.0, -1.0 + epsilon, 1.0 - epsilon)
                normals = math.sqrt(2.0) * torch.erfinv(arguments)
                norm = torch.linalg.vector_norm(normals, dim=1, keepdim=True).clamp_min(
                    torch.finfo(self.dtype).tiny
                )
                vectors = normals / norm
                angles = self._direction_angles(vectors, group["kind"])
                limits = self.limits[indices]
                points[:, indices] = (angles - limits[:, 0]) / (limits[:, 1] - limits[:, 0])
        return self.canonical_unit(points)

    # Compute normalized feature distances to a reference set
    def minimum_topology_distance(self, candidates: torch.Tensor, references: torch.Tensor | None) -> torch.Tensor:
        if references is None or len(references) == 0:
            return torch.full((len(candidates),), torch.inf, device=self.device, dtype=self.dtype)
        first = self.surrogate_inputs_from_unit(candidates)
        second = self.surrogate_inputs_from_unit(references)
        lengthscale = torch.ones(self.geometry.effective_dim, device=self.device, dtype=self.dtype)
        squared = self.geometry.squared_distance(first, second, lengthscale)
        return torch.sqrt(torch.clamp(torch.min(squared, dim=1).values, min=0.0))

    # Compute a JSON-safe description of the optimizer space
    def metadata(self) -> dict:
        return {
            "bounds": {name: dict(self.bounds[name]) for name in self.names},
            "dimension": self.dimension,
            "kernel_dimension": self.geometry.effective_dim,
            "integer_parameters": [name for name in self.names if is_integer_bound(self.bounds[name])],
            "parameter_topology": self.parameter_topology,
        }


# Convert a Torch dtype request into the selected arithmetic precision
def resolve_icebo_dtype(device: torch.device, dtype: str | torch.dtype) -> torch.dtype:
    if isinstance(dtype, torch.dtype):
        return dtype
    name = str(dtype).lower()
    if name == "auto":
        return torch.float32 if device.type == "cuda" else torch.float64
    if name == "float32":
        return torch.float32
    if name == "float64":
        return torch.float64
    raise ValueError("ICEBO dtype must be auto, float32 or float64")


# Select an available Torch device from an explicit or automatic request
def resolve_icebo_device(device: str | torch.device) -> torch.device:
    if isinstance(device, torch.device):
        selected = device
    elif str(device).lower() == "auto":
        selected = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        selected = torch.device(str(device))
    if selected.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("ICEBO CUDA device requested but CUDA is unavailable")
    return selected
