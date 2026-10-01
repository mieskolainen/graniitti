# Shared product manifold geometry for tuning surrogate kernels
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math

import torch
import torch.nn as nn


# Decode hyperspherical angles into normalized real vectors
def projective_vectors(angles: torch.Tensor) -> torch.Tensor:
    remaining = torch.ones((*angles.shape[:-1], 1), dtype=angles.dtype, device=angles.device)
    components = []
    for index in range(angles.shape[-1]):
        angle = angles[..., [index]]
        components.append(remaining * torch.sin(angle))
        remaining = remaining * torch.cos(angle)
    return torch.cat([remaining, *components], dim=-1)


# Decode first orthant hyperspherical angles into square roots of simplex probabilities
def simplex_amplitudes(angles: torch.Tensor) -> torch.Tensor:
    remaining = torch.ones((*angles.shape[:-1], 1), dtype=angles.dtype, device=angles.device)
    components = []
    for index in range(angles.shape[-1]):
        angle = angles[..., [index]]
        components.append(remaining * torch.cos(angle))
        remaining = remaining * torch.sin(angle)
    return torch.cat([*components, remaining], dim=-1)


class ProductManifoldGeometry(nn.Module):
    """Direct metric for ordinary, circular, spherical, projective and simplex parameters"""

    # Initialize an indexed product manifold without expanding input coordinates
    def __init__(
        self, dimension: int, topology_groups: tuple[dict, ...], numeric_bounds: tuple[tuple[float, float], ...]
    ) -> None:
        super().__init__()
        self.dimension = int(dimension)
        self.topology_groups = tuple(
            {
                "kind": str(group["kind"]),
                "indices": tuple(int(index) for index in group["indices"]),
                **(
                    {"period": float(group.get("period", math.nan))}
                    if str(group["kind"])
                    in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"}
                    else {}
                ),
                **({"weights": tuple(float(value) for value in group["weights"])} if "weights" in group else {}),
            }
            for group in topology_groups
        )
        self.numeric_bounds = tuple((float(bound[0]), float(bound[1])) for bound in numeric_bounds)
        self._validate()
        grouped = {index for group in self.topology_groups for index in group["indices"]}
        ordinary = [index for index in range(self.dimension) if index not in grouped]
        periodic = []
        for group in self.topology_groups:
            if group["kind"] in {"phase", "phase_projective", "polar_components"}:
                periodic.append(group["indices"][0])
            elif group["kind"] == "polar_projective":
                periodic.append(group["indices"][1])
            elif group["kind"] == "polar_sphere":
                periodic.extend((group["indices"][1], group["indices"][-1]))
            elif group["kind"] == "sphere":
                periodic.append(group["indices"][-1])
        self.register_buffer("ordinary_indices", torch.tensor(ordinary, dtype=torch.long))
        self.register_buffer("periodic_indices", torch.tensor(periodic, dtype=torch.long))
        self.register_buffer(
            "ordinary_lower", torch.tensor([self.numeric_bounds[index][0] for index in ordinary], dtype=torch.float64)
        )
        self.register_buffer(
            "ordinary_width",
            torch.tensor(
                [self.numeric_bounds[index][1] - self.numeric_bounds[index][0] for index in ordinary], dtype=torch.float64
            ),
        )
        self.effective_dim = len(ordinary) + len(self.topology_groups)

    # Build a geometry from parameter names and shared topology metadata
    @classmethod
    def from_metadata(
        cls, param_names: list[str], numeric_bounds: tuple[tuple[float, float], ...], topology: dict | None
    ) -> ProductManifoldGeometry:
        names = [str(name) for name in param_names]
        positions = {name: index for index, name in enumerate(names)}
        groups = []
        for group in (topology or {}).get("groups", []):
            parameters = [str(name) for name in group.get("parameters", [])]
            groups.append(
                {
                    "kind": str(group.get("kind")),
                    "indices": tuple(positions[name] for name in parameters),
                    **(
                        {"period": group.get("period")}
                        if str(group.get("kind"))
                        in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"}
                        else {}
                    ),
                    **({"weights": group["weights"]} if "weights" in group else {}),
                }
            )
        return cls(len(names), tuple(groups), numeric_bounds)

    # Validate the complete indexed manifold layout
    def _validate(self) -> None:
        if len(self.numeric_bounds) != self.dimension:
            raise ValueError("Topology kernel bounds do not match the parameter dimension")
        if any(not math.isfinite(lower) or not math.isfinite(upper) or upper <= lower
               for lower, upper in self.numeric_bounds):
            raise ValueError("Topology kernel bounds must be finite with positive width")
        indices = []
        for group in self.topology_groups:
            kind = group["kind"]
            selected = group["indices"]
            if kind not in {
                "phase",
                "projective",
                "sphere",
                "simplex",
                "polar_projective",
                "cartesian_projective",
                "phase_projective",
                "polar_sphere",
                "polar_components",
                "radial_projective",
            }:
                raise ValueError(f'Unsupported topology kernel group "{kind}"')
            if kind == "phase" and len(selected) != 1:
                raise ValueError("Phase topology requires exactly one parameter")
            if kind in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"} and (
                not math.isfinite(group["period"]) or group["period"] <= 0.0
            ):
                raise ValueError("Phase topology requires a finite positive period")
            if kind in {"projective", "sphere", "simplex"} and not selected:
                raise ValueError(f"{kind} topology requires at least one parameter")
            if kind in {"polar_projective", "cartesian_projective"} and len(selected) < 2:
                raise ValueError(f"{kind} topology requires coefficient coordinates")
            if kind in {"polar_sphere", "polar_components"} and len(selected) < 2:
                raise ValueError(f"{kind} topology requires a phase and coupling coordinates")
            if kind == "radial_projective" and len(selected) < 2:
                raise ValueError("radial_projective topology requires a norm and direction")
            if kind == "phase_projective" and not selected:
                raise ValueError("phase_projective topology requires a phase")
            if "weights" in group:
                direction_size = {
                    "projective": len(selected) + 1,
                    "polar_projective": len(selected) - 1,
                    "cartesian_projective": len(selected) - 1,
                    "phase_projective": len(selected),
                    "radial_projective": len(selected),
                }.get(kind)
                if (
                    direction_size is None
                    or len(group["weights"]) != direction_size
                    or any(not math.isfinite(value) or value <= 0.0 for value in group["weights"])
                ):
                    raise ValueError("Topology kernel completion weights do not match the direction")
            if any(index < 0 or index >= self.dimension for index in selected):
                raise ValueError("Topology kernel index is outside the parameter space")
            indices.extend(selected)
        if len(indices) != len(set(indices)):
            raise ValueError("Topology kernel groups overlap")

    # Normalize ungrouped coordinates to the unit interval
    def _ordinary_values(self, values: torch.Tensor) -> torch.Tensor:
        selected = values.index_select(-1, self.ordinary_indices)
        return (selected - self.ordinary_lower.to(values)) / self.ordinary_width.to(values)

    # Compute one topology group distance for aligned or pairwise rows
    def _group_distance(self, first: torch.Tensor, second: torch.Tensor, group: dict, diag: bool) -> torch.Tensor:
        indices = list(group["indices"])
        if group["kind"] == "phase":
            scale = math.pi / group["period"]
            if diag:
                delta = first[..., indices[0]] - second[..., indices[0]]
            else:
                delta = first[..., :, indices[0]].unsqueeze(-1) - second[..., :, indices[0]].unsqueeze(-2)
            return 4.0 * torch.sin(scale * delta).square()

        if group["kind"] == "projective":
            decoded_first = projective_vectors(first[..., indices])
            decoded_second = projective_vectors(second[..., indices])
            decoded_first = self._weighted_direction(decoded_first, group)
            decoded_second = self._weighted_direction(decoded_second, group)
            inner = (
                (decoded_first * decoded_second).sum(dim=-1)
                if diag
                else torch.matmul(decoded_first, decoded_second.transpose(-1, -2))
            )
            return (1.0 - inner.square()).clamp_min(0.0)

        if group["kind"] in {"polar_projective", "cartesian_projective", "phase_projective"}:
            return self._complex_projective_distance(first, second, group, diag)

        if group["kind"] in {
            "polar_sphere",
            "polar_components",
            "radial_projective",
        }:
            first_features = self._coupling_features(first, group)
            second_features = self._coupling_features(second, group)
            if diag:
                return (first_features - second_features).square().sum(dim=-1)
            difference = first_features.unsqueeze(-2) - second_features.unsqueeze(-3)
            return difference.square().sum(dim=-1)

        if group["kind"] == "sphere":
            decoded_first = projective_vectors(first[..., indices])
            decoded_second = projective_vectors(second[..., indices])
            inner = (
                (decoded_first * decoded_second).sum(dim=-1)
                if diag
                else torch.matmul(decoded_first, decoded_second.transpose(-1, -2))
            )
            return (2.0 * (1.0 - inner)).clamp_min(0.0)

        decoded_first = simplex_amplitudes(first[..., indices])
        decoded_second = simplex_amplitudes(second[..., indices])
        if diag:
            return 0.5 * (decoded_first - decoded_second).square().sum(dim=-1)
        difference = decoded_first.unsqueeze(-2) - decoded_second.unsqueeze(-3)
        return 0.5 * difference.square().sum(dim=-1)

    # Decode one coupled or quadratically sewn coupling group into real features
    def _coupling_features(self, values: torch.Tensor, group: dict) -> torch.Tensor:
        indices = list(group["indices"])
        selected = values[..., indices]
        kind = group["kind"]
        if kind == "polar_sphere":
            scale = max(abs(value) for value in self.numeric_bounds[indices[0]])
            phase = 2.0 * math.pi * selected[..., 1] / group["period"]
            coefficient = torch.complex(selected[..., 0] / scale, torch.zeros_like(selected[..., 0])) * torch.exp(
                torch.complex(torch.zeros_like(phase), phase)
            )
            vector = coefficient.unsqueeze(-1) * projective_vectors(selected[..., 2:])
            return torch.cat((vector.real, vector.imag), dim=-1)
        if kind == "polar_components":
            scale = max(abs(value) for index in indices[1:] for value in self.numeric_bounds[index])
            phase = 2.0 * math.pi * selected[..., 0] / group["period"]
            coefficient = torch.exp(torch.complex(torch.zeros_like(phase), phase))
            vector = coefficient.unsqueeze(-1) * selected[..., 1:] / scale
            return torch.cat((vector.real, vector.imag), dim=-1)
        if kind == "radial_projective":
            scale = max(abs(value) for value in self.numeric_bounds[indices[0]])
            coefficient = selected[..., 0] / scale
            direction = projective_vectors(selected[..., 1:])
        if "weights" in group:
            weights = torch.as_tensor(group["weights"], dtype=direction.dtype, device=direction.device)
            direction = torch.sqrt(weights) * direction
            direction = direction / torch.linalg.vector_norm(direction, dim=-1, keepdim=True)
        projector = direction.unsqueeze(-1) * direction.unsqueeze(-2)
        triangle = torch.triu_indices(direction.shape[-1], direction.shape[-1], device=direction.device)
        packed = projector[..., triangle[0], triangle[1]]
        factors = torch.where(
            triangle[0] == triangle[1],
            torch.ones_like(triangle[0], dtype=direction.dtype),
            torch.full_like(triangle[0], math.sqrt(2.0), dtype=direction.dtype),
        )
        return coefficient.unsqueeze(-1) * packed * factors

    # Compute distance between the physical complex signed-real coupling vectors
    def _complex_projective_distance(
        self, first: torch.Tensor, second: torch.Tensor, group: dict, diag: bool
    ) -> torch.Tensor:
        indices = list(group["indices"])
        kind = group["kind"]
        if kind == "polar_projective":
            scale = max(abs(value) for value in self.numeric_bounds[indices[0]])
            phase_first = 2.0 * math.pi * first[..., indices[1]] / group["period"]
            phase_second = 2.0 * math.pi * second[..., indices[1]] / group["period"]
            coefficient_first = torch.complex(
                first[..., indices[0]] / scale, torch.zeros_like(first[..., indices[0]])
            ) * torch.exp(torch.complex(torch.zeros_like(phase_first), phase_first))
            coefficient_second = torch.complex(
                second[..., indices[0]] / scale, torch.zeros_like(second[..., indices[0]])
            ) * torch.exp(torch.complex(torch.zeros_like(phase_second), phase_second))
            angle_indices = indices[2:]
        elif kind == "cartesian_projective":
            scale = max(abs(value) for index in indices[:2] for value in self.numeric_bounds[index])
            coefficient_first = torch.complex(first[..., indices[0]], first[..., indices[1]]) / scale
            coefficient_second = torch.complex(second[..., indices[0]], second[..., indices[1]]) / scale
            angle_indices = indices[2:]
        else:
            phase_first = 2.0 * math.pi * first[..., indices[0]] / group["period"]
            phase_second = 2.0 * math.pi * second[..., indices[0]] / group["period"]
            coefficient_first = torch.exp(
                torch.complex(torch.zeros_like(phase_first), phase_first)
            )
            coefficient_second = torch.exp(
                torch.complex(torch.zeros_like(phase_second), phase_second)
            )
            angle_indices = indices[1:]

        direction_first = projective_vectors(first[..., angle_indices])
        direction_second = projective_vectors(second[..., angle_indices])
        direction_first = self._weighted_direction(direction_first, group)
        direction_second = self._weighted_direction(direction_second, group)
        inner = (
            (direction_first * direction_second).sum(dim=-1)
            if diag
            else torch.matmul(direction_first, direction_second.transpose(-1, -2))
        )
        norm_first = coefficient_first.abs().square()
        norm_second = coefficient_second.abs().square()
        if diag:
            cross = coefficient_first.real * coefficient_second.real + coefficient_first.imag * coefficient_second.imag
        else:
            norm_first = norm_first.unsqueeze(-1)
            norm_second = norm_second.unsqueeze(-2)
            cross = (
                coefficient_first.real.unsqueeze(-1) * coefficient_second.real.unsqueeze(-2)
                + coefficient_first.imag.unsqueeze(-1) * coefficient_second.imag.unsqueeze(-2)
            )
        return (norm_first + norm_second - 2.0 * inner * cross).clamp_min(0.0)

    # Map compact directions into the normalized completed Frobenius basis
    def _weighted_direction(self, direction: torch.Tensor, group: dict) -> torch.Tensor:
        if "weights" not in group:
            return direction
        weights = torch.as_tensor(group["weights"], dtype=direction.dtype, device=direction.device)
        completed = torch.sqrt(weights) * direction
        return completed / torch.linalg.vector_norm(completed, dim=-1, keepdim=True)

    # Compute squared product manifold distances
    def squared_distance(
        self, first: torch.Tensor, second: torch.Tensor, lengthscale: torch.Tensor, *, diag: bool = False
    ) -> torch.Tensor:
        if lengthscale.numel() != self.effective_dim:
            raise ValueError("Topology kernel length scales do not match its effective dimension")
        scale = lengthscale.reshape(-1)
        ordinary_count = self.ordinary_indices.numel()
        if diag:
            ordinary = self._ordinary_values(first) - self._ordinary_values(second)
            distance = ((ordinary / scale[:ordinary_count]) ** 2).sum(dim=-1)
        else:
            ordinary = self._ordinary_values(first).unsqueeze(-2) - self._ordinary_values(second).unsqueeze(-3)
            distance = ((ordinary / scale[None, None, :ordinary_count]) ** 2).sum(dim=-1)

        for group_index, group in enumerate(self.topology_groups):
            component = self._group_distance(first, second, group, diag)
            distance = distance + component / scale[ordinary_count + group_index].square()
        return distance

    # Expand group length scales to physical coordinate trust region weights
    def coordinate_weights(self, lengthscale: torch.Tensor) -> torch.Tensor:
        scale = lengthscale.reshape(-1)
        weights = torch.ones(self.dimension, dtype=scale.dtype, device=scale.device)
        ordinary_count = self.ordinary_indices.numel()
        if ordinary_count:
            weights[self.ordinary_indices] = scale[:ordinary_count]
        for group_index, group in enumerate(self.topology_groups):
            weights[list(group["indices"])] = scale[ordinary_count + group_index]
        weights = weights / torch.mean(weights)
        geometric_mean = torch.exp(torch.mean(torch.log(torch.clamp(weights, min=1.0e-12))))
        return weights / geometric_mean

    # Evaluate a Matern five halves covariance from the shared direct distance
    def matern52(self, first: torch.Tensor, second: torch.Tensor, lengthscale: torch.Tensor) -> torch.Tensor:
        squared = self.squared_distance(first, second, lengthscale)
        distance = torch.sqrt(torch.clamp(squared, min=torch.finfo(first.dtype).tiny))
        scaled = math.sqrt(5.0) * distance
        return (1.0 + scaled + scaled.square() / 3.0) * torch.exp(-scaled)
