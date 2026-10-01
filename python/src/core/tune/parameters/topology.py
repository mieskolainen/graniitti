# Shared topology-aware surrogate feature construction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>

from __future__ import annotations

import math

import numpy as np
import torch

from core.tune.parameters import tools
from core.tune.parameters.space import to_unit


# Validate topology groups and index them by optimizer parameter name
def indexed_parameter_topology(param_names: list[str], topology: dict | None) -> dict:
    if not topology:
        return {}
    if int(topology.get("schema_version", -1)) != 1:
        raise ValueError("Unsupported parameter-topology schema")

    indexed = {}
    known = set(str(name) for name in param_names)
    for group in topology.get("groups", []):
        kind = str(group.get("kind"))
        parameters = [str(name) for name in group.get("parameters", [])]
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
        } or not parameters:
            raise ValueError(f"Invalid parameter-topology group: {group}")
        if any(name not in known for name in parameters):
            raise ValueError(f"Parameter-topology group is incompatible with the parameter space: {group}")
        if any(name in indexed for name in parameters):
            raise ValueError(f"Parameter-topology groups overlap: {group}")
        period = group.get("period")
        if kind in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"} and (
            isinstance(period, bool)
            or not isinstance(period, (int, float))
            or not math.isfinite(float(period))
            or float(period) <= 0.0
        ):
            raise ValueError(f"Phase topology requires a finite positive period: {group}")
        if kind in {"polar_sphere", "polar_components"} and len(parameters) < 2:
            raise ValueError(f"{kind} topology requires a phase and coupling coordinates: {group}")
        if kind == "radial_projective" and len(parameters) < 2:
            raise ValueError(f"{kind} topology requires a norm and direction: {group}")
        entry = {"base": str(group.get("base", parameters[0])), "kind": kind, "parameters": parameters}
        if kind in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"}:
            entry["period"] = float(period)
        if "weights" in group:
            weights = [float(value) for value in group["weights"]]
            direction_size = {
                "projective": len(parameters) + 1,
                "polar_projective": len(parameters) - 1,
                "cartesian_projective": len(parameters) - 1,
                "phase_projective": len(parameters),
                "radial_projective": len(parameters),
            }.get(kind)
            if direction_size is None or len(weights) != direction_size or any(
                not math.isfinite(value) or value <= 0.0 for value in weights
            ):
                raise ValueError(f"Invalid parameter-topology completion weights: {group}")
            entry["weights"] = weights
        for name in parameters:
            indexed[name] = entry
    return indexed


# Compute ordered plain-parameter and topology-group index blocks
def _feature_plan(param_names: list[str], topology: dict | None) -> list[tuple[list[int], dict | None]]:
    names = [str(name) for name in param_names]
    indexed = indexed_parameter_topology(names, topology)
    positions = {name: index for index, name in enumerate(names)}
    emitted = set()
    plan = []
    for index, name in enumerate(names):
        group = indexed.get(name)
        if group is None:
            plan.append(([index], None))
            continue
        base = group["base"]
        if base in emitted:
            continue
        emitted.add(base)
        plan.append(([positions[item] for item in group["parameters"]], group))
    return plan


# Compute topology-aware surrogate feature names in emitted order
def model_feature_names(param_names: list[str], topology: dict | None) -> list[str]:
    names = [str(name) for name in param_names]
    output = []
    for indices, group in _feature_plan(names, topology):
        if group is None:
            output.append(names[indices[0]])
            continue
        base = group["base"]
        if group["kind"] == "phase":
            output.extend([f"{base}:cos", f"{base}:sin"])
        elif group["kind"] in {"polar_projective", "cartesian_projective", "phase_projective", "polar_sphere", "polar_components"}:
            vector_size = len(indices) - 1
            if group["kind"] == "phase_projective":
                vector_size += 1
            output.extend(
                [
                    *[f"{base}:re[{index}]" for index in range(vector_size)],
                    *[f"{base}:im[{index}]" for index in range(vector_size)],
                ]
            )
        elif group["kind"] == "radial_projective":
            output.extend(tools.projective_embedding_feature_names(base, len(indices)))
        elif group["kind"] == "projective":
            output.extend(tools.projective_embedding_feature_names(base, len(indices) + 1))
        elif group["kind"] == "sphere":
            output.extend(tools.spherical_embedding_feature_names(base, len(indices) + 1))
        else:
            output.extend(f"{base}[{index}]" for index in range(len(indices) + 1))
    return output


# Decode one topology group into its surrogate feature vector
def _decode_group(group: dict, values: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    kind = group["kind"]
    if kind == "phase":
        phase = 2.0 * math.pi * float(values[0]) / group["period"]
        return np.asarray([np.cos(phase), np.sin(phase)], dtype=np.float64)
    if kind == "projective":
        if "weights" in group:
            direction = np.asarray(tools.hemisphere_vector_from_angles(values.tolist()), dtype=np.float64)
            weights = np.asarray(group["weights"], dtype=np.float64)
            direction = np.sqrt(weights) * direction / np.sqrt(np.sum(weights * direction**2))
            return np.asarray(tools.projective_embedding_from_vector(direction.tolist()), dtype=np.float64)
        return np.asarray(tools.projective_embedding_from_angles(values.tolist()), dtype=np.float64)
    if kind in {"polar_projective", "cartesian_projective", "phase_projective"}:
        if kind == "polar_projective":
            scale = max(abs(float(bounds[0, 0])), abs(float(bounds[0, 1])))
            if scale <= 0.0:
                raise ValueError("polar_projective requires a nonzero norm scale")
            phase = 2.0 * math.pi * float(values[1]) / group["period"]
            coefficient = (float(values[0]) / scale) * np.exp(1j * phase)
            angles = values[2:]
        elif kind == "cartesian_projective":
            scale = max(abs(float(bounds[:2].min())), abs(float(bounds[:2].max())))
            if scale <= 0.0:
                raise ValueError("cartesian_projective requires a nonzero component scale")
            coefficient = complex(float(values[0]), float(values[1])) / scale
            angles = values[2:]
        else:
            phase = 2.0 * math.pi * float(values[0]) / group["period"]
            coefficient = np.exp(1j * phase)
            angles = values[1:]
        direction = np.asarray(tools.hemisphere_vector_from_angles(angles.tolist()), dtype=np.float64)
        if "weights" in group:
            weights = np.asarray(group["weights"], dtype=np.float64)
            direction = np.sqrt(weights) * direction / np.sqrt(np.sum(weights * direction**2))
        return np.concatenate((coefficient.real * direction, coefficient.imag * direction))
    if kind == "polar_sphere":
        scale = max(abs(float(bounds[0, 0])), abs(float(bounds[0, 1])))
        phase = 2.0 * math.pi * float(values[1]) / group["period"]
        coefficient = (float(values[0]) / scale) * np.exp(1j * phase)
        direction = np.asarray(tools.spherical_vector_from_angles(values[2:].tolist()), dtype=np.float64)
        return np.concatenate((coefficient.real * direction, coefficient.imag * direction))
    if kind == "polar_components":
        scale = float(np.max(np.abs(bounds[1:])))
        phase = 2.0 * math.pi * float(values[0]) / group["period"]
        vector = (values[1:] / scale) * np.exp(1j * phase)
        return np.concatenate((vector.real, vector.imag))
    if kind == "radial_projective":
        scale = float(np.max(np.abs(bounds[0])))
        direction = np.asarray(tools.hemisphere_vector_from_angles(values[1:].tolist()), dtype=np.float64)
        if "weights" in group:
            weights = np.asarray(group["weights"], dtype=np.float64)
            direction = np.sqrt(weights) * direction / np.sqrt(np.sum(weights * direction**2))
        projector = np.asarray(tools.projective_embedding_from_vector(direction.tolist()), dtype=np.float64)
        return (float(values[0]) / scale) * projector
    if kind == "sphere":
        return np.asarray(tools.spherical_embedding_from_angles(values.tolist()), dtype=np.float64)
    if kind == "simplex":
        return np.asarray(tools.simplex_probabilities_from_angles(values.tolist()), dtype=np.float64)
    raise ValueError(f'Unknown parameter-topology kind "{kind}"')


# Map physical parameters to topology-aware surrogate coordinates
def model_features(X: np.ndarray, bounds: np.ndarray, param_names: list[str], topology: dict | None) -> np.ndarray:
    values = np.asarray(X, dtype=np.float64)
    single_row = values.ndim == 1
    rows = np.atleast_2d(values)
    if rows.ndim != 2 or rows.shape[1] != len(param_names):
        raise ValueError("model_features: input dimension does not match parameter names")
    if not topology:
        return to_unit(values, bounds)

    unit = to_unit(rows, bounds)
    output_columns = []
    for indices, group in _feature_plan(param_names, topology):
        if group is None:
            output_columns.append(unit[:, indices[0]])
            continue
        decoded = np.asarray(
            [_decode_group(group, row[indices], np.asarray(bounds)[indices]) for row in rows],
            dtype=np.float64,
        )
        output_columns.extend(decoded[:, column] for column in range(decoded.shape[1]))
    output = np.column_stack(output_columns)
    return output[0] if single_row else output


# Compute inclusive and exclusive cumulative products for batched angles
def _cumulative_products(values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    ones = torch.ones((values.shape[0], 1), dtype=values.dtype, device=values.device)
    inclusive = torch.cumprod(values, dim=1)
    return inclusive, torch.cat((ones, inclusive[:, :-1]), dim=1)


# Decode signed projective angles into sign-invariant Torch features
def _projective_embedding_torch(angles: torch.Tensor) -> torch.Tensor:
    cosine_prefix, exclusive_cosine = _cumulative_products(torch.cos(angles))
    components = exclusive_cosine * torch.sin(angles)
    vector = torch.cat((cosine_prefix[:, -1:], components), dim=1)
    projector = vector[:, :, None] * vector[:, None, :]
    triangle = torch.triu_indices(vector.shape[1], vector.shape[1], device=vector.device)
    features = projector[:, triangle[0], triangle[1]]
    scales = torch.where(
        triangle[0] == triangle[1],
        torch.ones_like(triangle[0], dtype=vector.dtype),
        torch.full_like(triangle[0], np.sqrt(2.0), dtype=vector.dtype),
    )
    return features * scales


# Decode hemisphere angles without applying an independent antipodal quotient
def _hemisphere_vector_torch(angles: torch.Tensor) -> torch.Tensor:
    if angles.shape[1] == 0:
        return torch.ones((angles.shape[0], 1), dtype=angles.dtype, device=angles.device)
    cosine_prefix, exclusive_cosine = _cumulative_products(torch.cos(angles))
    components = exclusive_cosine * torch.sin(angles)
    return torch.cat((cosine_prefix[:, -1:], components), dim=1)


# Decode oriented spherical angles into unit vector Torch features
def _spherical_embedding_torch(angles: torch.Tensor) -> torch.Tensor:
    cosine_prefix, exclusive_cosine = _cumulative_products(torch.cos(angles))
    components = exclusive_cosine * torch.sin(angles)
    return torch.cat((cosine_prefix[:, -1:], components), dim=1)


# Decode first-orthant simplex angles into Torch probabilities
def _simplex_probabilities_torch(angles: torch.Tensor) -> torch.Tensor:
    sine_prefix, exclusive_sine = _cumulative_products(torch.sin(angles))
    amplitudes = torch.cat((exclusive_sine * torch.cos(angles), sine_prefix[:, -1:]), dim=1)
    return amplitudes.square()


# Map physical Torch parameters to differentiable surrogate coordinates
def model_features_torch(
    X: torch.Tensor, bounds: np.ndarray | torch.Tensor, param_names: list[str], topology: dict | None
) -> torch.Tensor:
    values = X
    single_row = values.ndim == 1
    rows = values.unsqueeze(0) if single_row else values
    if rows.ndim != 2 or rows.shape[1] != len(param_names):
        raise ValueError("model_features_torch: input dimension does not match parameter names")
    limits = torch.as_tensor(bounds, dtype=rows.dtype, device=rows.device)
    if limits.shape != (len(param_names), 2):
        raise ValueError("model_features_torch: bounds have inconsistent dimensions")
    unit = (rows - limits[:, 0]) / (limits[:, 1] - limits[:, 0])
    if not topology:
        return unit[0] if single_row else unit

    output_columns = []
    for indices, group in _feature_plan(param_names, topology):
        if group is None:
            output_columns.append(unit[:, indices])
            continue
        angles = rows[:, indices]
        if group["kind"] == "phase":
            phase = 2.0 * math.pi * angles[:, 0] / group["period"]
            decoded = torch.stack((torch.cos(phase), torch.sin(phase)), dim=1)
        elif group["kind"] == "projective":
            if "weights" in group:
                direction = _hemisphere_vector_torch(angles)
                weights = torch.as_tensor(group["weights"], dtype=rows.dtype, device=rows.device)
                direction = torch.sqrt(weights) * direction / torch.sqrt((weights * direction.square()).sum(dim=1))[:, None]
                projector = direction[:, :, None] * direction[:, None, :]
                triangle = torch.triu_indices(direction.shape[1], direction.shape[1], device=direction.device)
                decoded = projector[:, triangle[0], triangle[1]]
                decoded = decoded * torch.where(
                    triangle[0] == triangle[1],
                    torch.ones_like(triangle[0], dtype=direction.dtype),
                    torch.full_like(triangle[0], np.sqrt(2.0), dtype=direction.dtype),
                )
            else:
                decoded = _projective_embedding_torch(angles)
        elif group["kind"] in {"polar_projective", "cartesian_projective", "phase_projective"}:
            if group["kind"] == "polar_projective":
                scale = torch.max(torch.abs(limits[indices[0]]))
                phase = 2.0 * math.pi * angles[:, 1] / group["period"]
                coefficient = torch.complex(angles[:, 0] / scale, torch.zeros_like(angles[:, 0])) * torch.exp(
                    1j * phase
                )
                direction = _hemisphere_vector_torch(angles[:, 2:])
            elif group["kind"] == "cartesian_projective":
                scale = torch.max(torch.abs(limits[indices[:2]]))
                coefficient = torch.complex(angles[:, 0], angles[:, 1]) / scale
                direction = _hemisphere_vector_torch(angles[:, 2:])
            else:
                phase = 2.0 * math.pi * angles[:, 0] / group["period"]
                coefficient = torch.exp(torch.complex(torch.zeros_like(phase), phase))
                direction = _hemisphere_vector_torch(angles[:, 1:])
            if "weights" in group:
                weights = torch.as_tensor(group["weights"], dtype=rows.dtype, device=rows.device)
                direction = torch.sqrt(weights) * direction / torch.sqrt((weights * direction.square()).sum(dim=1))[:, None]
            decoded = torch.cat(
                (coefficient.real[:, None] * direction, coefficient.imag[:, None] * direction), dim=1
            )
        elif group["kind"] == "polar_sphere":
            scale = torch.max(torch.abs(limits[indices[0]]))
            phase = 2.0 * math.pi * angles[:, 1] / group["period"]
            coefficient = torch.complex(angles[:, 0] / scale, torch.zeros_like(angles[:, 0])) * torch.exp(
                1j * phase
            )
            direction = _spherical_embedding_torch(angles[:, 2:])
            decoded = torch.cat(
                (coefficient.real[:, None] * direction, coefficient.imag[:, None] * direction), dim=1
            )
        elif group["kind"] == "polar_components":
            scale = torch.max(torch.abs(limits[indices[1:]]))
            phase = 2.0 * math.pi * angles[:, 0] / group["period"]
            coefficient = torch.exp(torch.complex(torch.zeros_like(phase), phase))
            vector = (angles[:, 1:] / scale) * coefficient[:, None]
            decoded = torch.cat((vector.real, vector.imag), dim=1)
        elif group["kind"] == "radial_projective":
            scale = torch.max(torch.abs(limits[indices[0]]))
            direction = _hemisphere_vector_torch(angles[:, 1:])
            if "weights" in group:
                weights = torch.as_tensor(group["weights"], dtype=rows.dtype, device=rows.device)
                direction = torch.sqrt(weights) * direction
                direction = direction / torch.linalg.vector_norm(direction, dim=1, keepdim=True)
            projector = direction[:, :, None] * direction[:, None, :]
            triangle = torch.triu_indices(direction.shape[1], direction.shape[1], device=direction.device)
            decoded = projector[:, triangle[0], triangle[1]]
            decoded = decoded * torch.where(
                triangle[0] == triangle[1],
                torch.ones_like(triangle[0], dtype=direction.dtype),
                torch.full_like(triangle[0], np.sqrt(2.0), dtype=direction.dtype),
            )
            decoded = (angles[:, :1] / scale) * decoded
        elif group["kind"] == "sphere":
            decoded = _spherical_embedding_torch(angles)
        elif group["kind"] == "simplex":
            decoded = _simplex_probabilities_torch(angles)
        else:
            raise ValueError(f'Unknown parameter-topology kind "{group["kind"]}"')
        output_columns.append(decoded)
    output = torch.cat(output_columns, dim=1)
    return output[0] if single_row else output
