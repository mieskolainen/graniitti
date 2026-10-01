# GPU accelerated bounded L-BFGS optimization and local covariance
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from numbers import Real

import numpy as np
import torch


@dataclass(frozen=True)
class Minimum:
    """Store the bounded optimizer minimum summary"""

    fval: float
    is_valid: bool
    message: str


@dataclass
class Result:
    """Expose fitted values and covariance to shared post-fit diagnostics"""

    values: dict[str, float]
    covariance: np.ndarray
    errors: dict[str, float]
    fmin: Minimum


# Validate the configurable step scale and Armijo backtracking controls
def validate_line_search(learning_rate: float, line_search_decay: float, line_search_armijo: float) -> None:
    for name, value in (("learning_rate", learning_rate), ("line_search_decay", line_search_decay),
                        ("line_search_armijo", line_search_armijo)):
        if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"L-BFGS {name} must be finite and positive")
        if name != "learning_rate" and value >= 1.0:
            raise ValueError(f"L-BFGS {name} must be less than one")


# Evaluate one scalar objective and optionally defer its gradient until acceptance
def value_gradient(
    objective: Callable[[torch.Tensor], torch.Tensor], values: torch.Tensor, *, defer: bool = False
) -> tuple[torch.Tensor, torch.Tensor | Callable[[], torch.Tensor]]:
    with torch.enable_grad():
        point = values.detach().clone().requires_grad_(True)
        objective_value = objective(point)

    # Differentiate the retained forward evaluation only when requested
    def derivative():
        return torch.autograd.grad(objective_value, point)[0].detach()

    return objective_value.detach(), derivative if defer else derivative()


# Validate and broadcast finite, one-sided or unbounded optimization coordinates
def _bounds(start: torch.Tensor, bounds: tuple | None) -> tuple[torch.Tensor, torch.Tensor]:
    lower, upper = (-torch.inf, torch.inf) if bounds is None else bounds
    lower, upper = (torch.broadcast_to(torch.as_tensor(value, dtype=start.dtype, device=start.device), start.shape)
                    for value in (lower, upper))
    if bool(torch.any(torch.isnan(lower) | torch.isnan(upper) | (lower > upper) | lower.isposinf() | upper.isneginf())):
        raise ValueError("L-BFGS bounds must be ordered and admit finite coordinates")
    return lower, upper


# Select coordinates whose descent directions are feasible within the bounds
def _free_mask(values: torch.Tensor, gradient: torch.Tensor, bounds: tuple = (0.0, 1.0)) -> torch.Tensor:
    tolerance = 16.0 * torch.finfo(values.dtype).eps * values.abs().clamp_min(1.0)
    at_lower = values <= bounds[0] + tolerance
    at_upper = values >= bounds[1] - tolerance
    blocked = (at_lower & (gradient > 0.0)) | (at_upper & (gradient < 0.0))
    return ~blocked


# Remove gradient components whose descent directions leave the unit box
def projected_gradient(values: torch.Tensor, gradient: torch.Tensor, bounds: tuple = (0.0, 1.0)) -> torch.Tensor:
    return torch.where(_free_mask(values, gradient, bounds), gradient, torch.zeros_like(gradient))


# Restrict curvature directions to the free coordinates and the feasible tangent cone
def _project_direction(values: torch.Tensor, gradient: torch.Tensor, direction: torch.Tensor,
                       bounds: tuple = (0.0, 1.0)) -> torch.Tensor:
    feasible = -projected_gradient(values, -direction, bounds)
    return torch.where(_free_mask(values, gradient, bounds), feasible, torch.zeros_like(direction))


# Apply the L-BFGS inverse Hessian recursion on the selected Torch device
def _lbfgs_direction(
    gradient: torch.Tensor, step_history: list[torch.Tensor], gradient_history: list[torch.Tensor]
) -> torch.Tensor:
    valid = torch.ones(1, dtype=torch.bool, device=gradient.device)
    history = [(step[None], change[None], valid)
               for step, change in zip(step_history, gradient_history, strict=True)]
    return _batch_lbfgs_direction(gradient[None], history)[0]


# Limit scalar or batched search directions to the specified bounds
def _box_step_limit(values: torch.Tensor, direction: torch.Tensor, bounds: tuple = (0.0, 1.0)) -> torch.Tensor:
    distance = torch.where(direction > 0.0, bounds[1] - values, bounds[0] - values)
    moving = direction.abs() > 0.0
    limits = distance / torch.where(moving, direction, torch.ones_like(direction))
    return torch.where(moving, limits, torch.inf).amin(dim=-1)


# Backtrack one Armijo line search entirely with Torch evaluations
def _bounded_line_search(
    objective: Callable[[torch.Tensor], torch.Tensor] | None,
    values: torch.Tensor,
    objective_value: torch.Tensor,
    gradient: torch.Tensor,
    direction: torch.Tensor,
    max_steps: int,
    learning_rate: float,
    line_search_decay: float,
    line_search_armijo: float,
    evaluator: Callable[[torch.Tensor], tuple[torch.Tensor, torch.Tensor | Callable[[], torch.Tensor]]] | None = None,
    bounds: tuple = (0.0, 1.0),
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, float, int] | None:
    slope = torch.dot(gradient, direction)
    if not torch.isfinite(slope) or float(slope.item()) >= 0.0:
        return None
    step_size = min(learning_rate, float(_box_step_limit(values, direction, bounds)))
    if step_size <= 0.0:
        return None
    for evaluation in range(1, int(max_steps) + 1):
        candidate = torch.clamp(values + step_size * direction, *bounds)
        candidate_value, candidate_gradient = (
            value_gradient(objective, candidate, defer=True) if evaluator is None else evaluator(candidate)
        )
        armijo_limit = objective_value + line_search_armijo * step_size * slope
        if torch.isfinite(candidate_value) and bool(candidate_value <= armijo_limit):
            candidate_gradient = candidate_gradient() if callable(candidate_gradient) else candidate_gradient
            if bool(torch.all(torch.isfinite(candidate_gradient))):
                return candidate, candidate_value, candidate_gradient, step_size, evaluation
        del candidate_gradient
        step_size *= line_search_decay
    return None


# Retain one stable curvature pair in the finite L-BFGS history
def _update_history(
    step_history: list[torch.Tensor],
    gradient_history: list[torch.Tensor],
    step: torch.Tensor,
    gradient_change: torch.Tensor,
    history_size: int,
) -> None:
    curvature = torch.dot(step, gradient_change)
    curvature_floor = (
        torch.finfo(step.dtype).eps * torch.linalg.vector_norm(step) * torch.linalg.vector_norm(gradient_change)
    )
    if not bool(torch.isfinite(curvature)) or float(curvature.item()) <= float(curvature_floor.item()):
        return
    step_history.append(step)
    gradient_history.append(gradient_change)
    if len(step_history) > int(history_size):
        step_history.pop(0)
        gradient_history.pop(0)


# Minimize one differentiable objective with bounded or unbounded Torch L-BFGS
def minimize(
    objective: Callable[[torch.Tensor], torch.Tensor] | None,
    start: torch.Tensor,
    maxiter: int = 200,
    history_size: int = 20,
    gradient_tolerance: float = 1.0e-7,
    relative_gradient_tolerance: float = 0.0,
    line_search_steps: int = 30,
    retry_failed_search: bool = True,
    progress_callback: Callable[[int, float], None] | None = None,
    evaluator: Callable[[torch.Tensor], tuple[torch.Tensor, torch.Tensor | Callable[[], torch.Tensor]]] | None = None,
    *,
    bounds: tuple | None = (0.0, 1.0),
    function_tolerance: float = 0.0,
    learning_rate: float = 1.0,
    line_search_decay: float = 0.5,
    line_search_armijo: float = 1.0e-4,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, dict]:
    validate_line_search(learning_rate, line_search_decay, line_search_armijo)
    if not np.isfinite(function_tolerance) or function_tolerance < 0.0:
        raise ValueError("L-BFGS function tolerance must be finite and non-negative")
    if gradient_tolerance < 0.0 or relative_gradient_tolerance < 0.0:
        raise ValueError("L-BFGS gradient tolerances must be non-negative")
    bounds = _bounds(start, bounds)
    values = torch.clamp(start.detach().clone(), *bounds)
    if objective is None and evaluator is None:
        raise ValueError("L-BFGS requires an objective or a value and gradient evaluator")
    function_evaluations = gradient_evaluations = 0

    # Count actual forward and backward evaluations including rejected searches
    def evaluate(point):
        nonlocal function_evaluations, gradient_evaluations
        function_evaluations += 1
        value, derivative = value_gradient(objective, point, defer=True) if evaluator is None else evaluator(point)
        if not callable(derivative):
            gradient_evaluations += 1
            return value, derivative

        # Materialize a deferred derivative after the objective passes the line search
        def gradient():
            nonlocal gradient_evaluations
            gradient_evaluations += 1
            return derivative()

        return value, gradient

    objective_value, gradient = evaluate(values)
    gradient = gradient() if callable(gradient) else gradient
    step_history: list[torch.Tensor] = []
    gradient_history: list[torch.Tensor] = []
    initial_projected_norm = float(torch.max(torch.abs(projected_gradient(values, gradient, bounds))).item())
    effective_gradient_tolerance = max(
        float(gradient_tolerance), float(relative_gradient_tolerance) * max(1.0, initial_projected_norm)
    )
    history_resets = 0
    message = "maximum iterations reached"
    success = False
    iteration = 0
    for iteration in range(1, int(maxiter) + 1):
        completed_iterations = iteration
        if not bool(torch.isfinite(objective_value)) or not bool(torch.all(torch.isfinite(gradient))):
            break
        projected = projected_gradient(values, gradient, bounds)
        projected_norm = float(torch.max(torch.abs(projected)).item())
        if projected_norm <= effective_gradient_tolerance:
            message = "projected gradient tolerance reached"
            success = True
            break
        direction = _lbfgs_direction(projected, step_history, gradient_history)
        direction = _project_direction(values, gradient, direction, bounds)
        if float(torch.dot(gradient, direction).item()) >= 0.0:
            direction = -projected
        line_result = _bounded_line_search(
            objective=objective,
            values=values,
            objective_value=objective_value,
            gradient=gradient,
            direction=direction,
            max_steps=line_search_steps,
            learning_rate=learning_rate,
            line_search_decay=line_search_decay,
            line_search_armijo=line_search_armijo,
            bounds=bounds,
            evaluator=evaluate,
        )
        if line_result is None and step_history and retry_failed_search:
            step_history.clear()
            gradient_history.clear()
            history_resets += 1
            line_result = _bounded_line_search(
                objective=objective,
                values=values,
                objective_value=objective_value,
                gradient=gradient,
                direction=-projected,
                max_steps=line_search_steps,
                learning_rate=learning_rate,
                line_search_decay=line_search_decay,
                line_search_armijo=line_search_armijo,
                bounds=bounds,
                evaluator=evaluate,
            )
        if line_result is None:
            message = (
                "bounded line search failed after L-BFGS history reset"
                if history_resets > 0
                else "bounded line search failed"
            )
            break
        candidate, candidate_value, candidate_gradient, _, _ = line_result
        small_change = function_tolerance > 0.0 and bool(
            abs(candidate_value - objective_value)
            <= function_tolerance * max(1.0, abs(float(objective_value)), abs(float(candidate_value)))
        )
        # Rebuild curvature when the free coordinate subspace changes
        if bool(torch.any(_free_mask(values, gradient, bounds) != _free_mask(candidate, candidate_gradient, bounds))):
            step_history.clear()
            gradient_history.clear()
            history_resets += 1
        else:
            _update_history(
                step_history=step_history,
                gradient_history=gradient_history,
                step=candidate - values,
                gradient_change=candidate_gradient - gradient,
                history_size=history_size,
            )
        values = candidate
        objective_value = candidate_value
        gradient = candidate_gradient
        if progress_callback is not None:
            progress_callback(iteration, float(objective_value.detach().cpu()))
        if small_change:
            success, message = True, "relative objective tolerance reached"
            break
    projected = projected_gradient(values, gradient, bounds)
    if not bool(torch.isfinite(objective_value)) or not bool(torch.all(torch.isfinite(gradient))):
        message = "non-finite objective or gradient"
    elif float(torch.max(torch.abs(projected)).item()) <= effective_gradient_tolerance:
        success = True
        message = "projected gradient tolerance reached"
    diagnostics = {
        "success": success,
        "message": message,
        "iterations": completed_iterations if iteration else 0,
        "function_evaluations": function_evaluations,
        "gradient_evaluations": gradient_evaluations,
        "max_abs_unit_gradient": float(torch.max(torch.abs(gradient)).item()),
        "max_abs_projected_unit_gradient": float(torch.max(torch.abs(projected)).item()),
        "initial_max_abs_projected_unit_gradient": initial_projected_norm,
        "effective_gradient_tolerance": effective_gradient_tolerance,
        "relative_gradient_tolerance": float(relative_gradient_tolerance),
        "device": str(values.device),
        "dtype": str(values.dtype),
        "history_size": len(step_history),
        "history_resets": history_resets,
    }
    return values, objective_value, gradient, diagnostics


# Evaluate batched objectives and optionally defer their row-wise gradients
def batch_value_gradient(
    objective: Callable[[torch.Tensor], torch.Tensor], values: torch.Tensor, *, defer: bool = False
) -> tuple[torch.Tensor, torch.Tensor | Callable[[], torch.Tensor]]:
    with torch.enable_grad():
        points = values.detach().clone().requires_grad_(True)
        objective_values = objective(points).reshape(-1)
        if len(objective_values) != len(points):
            raise ValueError("Batched objective must return one value per row")
        finite_sum = torch.sum(
            torch.where(torch.isfinite(objective_values), objective_values, torch.zeros_like(objective_values))
        )

    # Differentiate finite rows only when the line search accepts a candidate
    def derivative():
        return torch.autograd.grad(finite_sum, points)[0].detach()

    return objective_values.detach(), derivative if defer else derivative()


# Compute one dot product per row
def _row_dot(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    return torch.sum(left * right, dim=1)


# Apply independent L-BFGS recursions to all active batch rows
def _batch_lbfgs_direction(
    gradient: torch.Tensor, history: list[tuple[torch.Tensor, torch.Tensor, torch.Tensor]]
) -> torch.Tensor:
    if not history:
        return -gradient
    # Apply the compact inverse BFGS form [REFERENCE: Byrd, Nocedal and Schnabel, Math. Programming 63 (1994) 129-156]
    steps, changes, valid = (torch.stack(items, dim=1) for items in zip(*history, strict=True))
    steps = steps.masked_fill(~valid[..., None], 0.0)
    changes = changes.masked_fill(~valid[..., None], 0.0)
    curvature = (steps * changes).sum(dim=-1)
    scales = curvature / changes.square().sum(dim=-1)
    usable = valid & torch.isfinite(scales) & (scales > 0.0)
    last = torch.where(usable, torch.arange(len(history), device=gradient.device), -1).amax(dim=1)
    scale = torch.where(last[:, None] >= 0, scales.gather(1, last.clamp_min(0)[:, None]), 1.0)
    upper = (steps @ changes.transpose(-1, -2)).triu()
    upper.diagonal(dim1=-2, dim2=-1).copy_(torch.where(valid, curvature, 1.0))
    alpha = torch.linalg.solve_triangular(upper, steps @ gradient[..., None], upper=True)
    vector = scale[..., None] * (gradient[..., None] - changes.transpose(-1, -2) @ alpha)
    correction = torch.linalg.solve_triangular(
        upper.transpose(-1, -2), curvature[..., None] * alpha - changes @ vector, upper=False
    )
    return -(vector + steps.transpose(-1, -2) @ correction).squeeze(-1)


# Run independent Armijo searches with one batched objective evaluation per step
def _batch_line_search(
    objective: Callable[[torch.Tensor], torch.Tensor],
    values: torch.Tensor,
    objective_values: torch.Tensor,
    gradient: torch.Tensor,
    direction: torch.Tensor,
    active: torch.Tensor,
    max_steps: int,
    learning_rate: float,
    line_search_decay: float,
    line_search_armijo: float,
    bounds: tuple = (0.0, 1.0),
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, int]:
    slope = _row_dot(gradient, direction)
    step_size = torch.minimum(torch.full_like(slope, learning_rate), _box_step_limit(values, direction, bounds))
    pending = active & torch.isfinite(slope) & (slope < 0.0) & (step_size > 0.0)
    accepted = torch.zeros_like(active)
    output_values = values.clone()
    output_objective = objective_values.clone()
    output_gradient = gradient.clone()
    evaluations = 0
    for _ in range(int(max_steps)):
        if not bool(torch.any(pending)):
            break
        candidate = torch.clamp(values + step_size[:, None] * direction, *bounds)
        candidate = torch.where(pending[:, None], candidate, values)
        candidate_objective, candidate_gradient = batch_value_gradient(objective, candidate, defer=True)
        evaluations += 1
        armijo_limit = objective_values + line_search_armijo * step_size * slope
        accept = (
            pending
            & torch.isfinite(candidate_objective)
            & (candidate_objective <= armijo_limit)
        )
        if bool(torch.any(accept)):
            candidate_gradient = candidate_gradient()
            accept &= torch.all(torch.isfinite(candidate_gradient), dim=1)
        else:
            del candidate_gradient
            step_size = torch.where(pending, line_search_decay * step_size, step_size)
            continue
        output_values = torch.where(accept[:, None], candidate, output_values)
        output_objective = torch.where(accept, candidate_objective, output_objective)
        output_gradient = torch.where(accept[:, None], candidate_gradient, output_gradient)
        accepted |= accept
        pending &= ~accept
        step_size = torch.where(pending, line_search_decay * step_size, step_size)
    return output_values, output_objective, output_gradient, accepted, evaluations


# Retain one batched L-BFGS curvature update with a row validity mask
def _batch_history_update(
    history: list[tuple[torch.Tensor, torch.Tensor, torch.Tensor]],
    step: torch.Tensor,
    gradient_change: torch.Tensor,
    accepted: torch.Tensor,
    history_size: int,
) -> None:
    curvature = _row_dot(step, gradient_change)
    curvature_floor = (
        torch.finfo(step.dtype).eps
        * torch.linalg.vector_norm(step, dim=1)
        * torch.linalg.vector_norm(gradient_change, dim=1)
    )
    valid = accepted & torch.isfinite(curvature) & (curvature > curvature_floor)
    history.append((step, gradient_change, valid))
    if len(history) > int(history_size):
        history.pop(0)


# Invalidate prior curvature pairs for selected batch rows
def _batch_history_reset(history: list[tuple[torch.Tensor, torch.Tensor, torch.Tensor]], reset: torch.Tensor) -> None:
    for index, (step, change, valid) in enumerate(history):
        history[index] = (step, change, valid & ~reset)


# Minimize independent objectives in parallel with batched bounded Torch L-BFGS
def minimize_batch(
    objective: Callable[[torch.Tensor], torch.Tensor],
    start: torch.Tensor,
    maxiter: int = 200,
    history_size: int = 20,
    gradient_tolerance: float = 1.0e-7,
    relative_gradient_tolerance: float = 0.0,
    line_search_steps: int = 30,
    retry_failed_search: bool = True,
    progress_callback: Callable[[int, int, float], None] | None = None,
    *,
    bounds: tuple | None = (0.0, 1.0),
    function_tolerance: float = 0.0,
    learning_rate: float = 1.0,
    line_search_decay: float = 0.5,
    line_search_armijo: float = 1.0e-4,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, dict]:
    validate_line_search(learning_rate, line_search_decay, line_search_armijo)
    if not np.isfinite(function_tolerance) or function_tolerance < 0.0:
        raise ValueError("L-BFGS function tolerance must be finite and non-negative")
    if start.ndim != 2 or start.shape[0] < 1 or start.shape[1] < 1:
        raise ValueError("Batched L-BFGS start must be a non-empty matrix")
    if gradient_tolerance < 0.0 or relative_gradient_tolerance < 0.0:
        raise ValueError("Batched L-BFGS gradient tolerances must be non-negative")
    bounds = _bounds(start, bounds)
    values = torch.clamp(start.detach().clone(), *bounds)
    objective_values, gradient = batch_value_gradient(objective, values)
    initial_projected_norm = torch.max(torch.abs(projected_gradient(values, gradient, bounds)), dim=1).values
    effective_gradient_tolerance = torch.maximum(
        torch.full_like(initial_projected_norm, float(gradient_tolerance)),
        float(relative_gradient_tolerance) * torch.clamp(initial_projected_norm, min=1.0),
    )
    finite = torch.isfinite(objective_values) & torch.all(torch.isfinite(gradient), dim=1)
    active = finite.clone()
    success = torch.zeros_like(active)
    function_converged = torch.zeros_like(active)
    line_search_failed = torch.zeros_like(active)
    iterations = torch.zeros(len(values), dtype=torch.int64, device=values.device)
    history: list[tuple[torch.Tensor, torch.Tensor, torch.Tensor]] = []
    batch_evaluations = 1
    line_search_retries = 0
    for iteration in range(1, int(maxiter) + 1):
        projected = projected_gradient(values, gradient, bounds)
        projected_norm = torch.max(torch.abs(projected), dim=1).values
        converged = active & (projected_norm <= effective_gradient_tolerance)
        success |= converged
        active &= ~converged
        if not bool(torch.any(active)):
            break
        direction = _batch_lbfgs_direction(projected, history)
        direction = _project_direction(values, gradient, direction, bounds)
        slope = _row_dot(gradient, direction)
        steepest = active & (~torch.isfinite(slope) | (slope >= 0.0))
        direction = torch.where(steepest[:, None], -projected, direction)
        new_values, new_objective, new_gradient, accepted, evaluations = _batch_line_search(
            objective=objective,
            values=values,
            objective_values=objective_values,
            gradient=gradient,
            direction=direction,
            active=active,
            max_steps=line_search_steps,
            learning_rate=learning_rate,
            line_search_decay=line_search_decay,
            line_search_armijo=line_search_armijo,
            bounds=bounds,
        )
        batch_evaluations += evaluations
        retry = active & ~accepted
        if retry_failed_search and bool(torch.any(retry)):
            _batch_history_reset(history, retry)
            (retry_values, retry_objective, retry_gradient, retry_accepted, retry_evaluations) = _batch_line_search(
                objective=objective,
                values=values,
                objective_values=objective_values,
                gradient=gradient,
                direction=-projected,
                active=retry,
                max_steps=line_search_steps,
                learning_rate=learning_rate,
                line_search_decay=line_search_decay,
                line_search_armijo=line_search_armijo,
                bounds=bounds,
            )
            batch_evaluations += retry_evaluations
            line_search_retries += int(torch.count_nonzero(retry).item())
            new_values = torch.where(retry_accepted[:, None], retry_values, new_values)
            new_objective = torch.where(retry_accepted, retry_objective, new_objective)
            new_gradient = torch.where(retry_accepted[:, None], retry_gradient, new_gradient)
            accepted |= retry_accepted
        failed = active & ~accepted
        line_search_failed |= failed
        active &= ~failed
        if function_tolerance > 0.0:
            scale = torch.maximum(torch.ones_like(objective_values), torch.maximum(objective_values.abs(), new_objective.abs()))
            small_change = accepted & ((new_objective - objective_values).abs() <= function_tolerance * scale)
            function_converged |= small_change
            success |= small_change
            active &= ~small_change
        iterations += accepted.to(torch.int64)
        # Rebuild curvature independently for rows whose free coordinates change
        changed = torch.any(_free_mask(values, gradient, bounds) != _free_mask(new_values, new_gradient, bounds), dim=1)
        _batch_history_reset(history, changed)
        _batch_history_update(
            history=history,
            step=new_values - values,
            gradient_change=new_gradient - gradient,
            accepted=accepted & ~changed,
            history_size=history_size,
        )
        values = new_values
        objective_values = new_objective
        gradient = new_gradient
        if progress_callback is not None:
            finite_objectives = objective_values[torch.isfinite(objective_values)]
            best_objective = (
                float(torch.min(finite_objectives).detach().cpu()) if len(finite_objectives) > 0 else float("nan")
            )
            progress_callback(iteration, int(torch.count_nonzero(active).item()), best_objective)
    projected = projected_gradient(values, gradient, bounds)
    projected_norm = torch.max(torch.abs(projected), dim=1).values
    converged = active & (projected_norm <= effective_gradient_tolerance)
    success |= converged
    messages = np.full(len(values), "maximum iterations reached", dtype=object)
    messages[~finite.detach().cpu().numpy()] = "non-finite initial objective or gradient"
    messages[line_search_failed.detach().cpu().numpy()] = "bounded line search failed"
    messages[success.detach().cpu().numpy()] = "projected gradient tolerance reached"
    messages[function_converged.detach().cpu().numpy()] = "relative objective tolerance reached"
    diagnostics = {
        "success": success,
        "messages": messages.tolist(),
        "iterations": iterations,
        "batch_evaluations": batch_evaluations,
        "row_evaluations": batch_evaluations * len(values),
        "max_abs_projected_unit_gradient": projected_norm,
        "initial_max_abs_projected_unit_gradient": initial_projected_norm,
        "effective_gradient_tolerance": effective_gradient_tolerance,
        "relative_gradient_tolerance": float(relative_gradient_tolerance),
        "line_search_retries": line_search_retries,
        "device": str(values.device),
        "dtype": str(values.dtype),
        "history_size": len(history),
    }
    return values, objective_values, gradient, diagnostics


# Compute the exact physical Hessian and its local likelihood covariance
def hessian_covariance(
    objective: Callable[[torch.Tensor], torch.Tensor],
    values: torch.Tensor,
    errordef: float = 1.0,
    free_mask: torch.Tensor | None = None,
    *,
    vectorize: bool = True,
) -> tuple[torch.Tensor, torch.Tensor, dict]:
    with torch.enable_grad():
        point = values.detach().clone().requires_grad_(True)
        try:
            hessian = torch.autograd.functional.hessian(objective, point, vectorize=vectorize)
        except (NotImplementedError, RuntimeError):
            if not vectorize:
                raise
            hessian = torch.autograd.functional.hessian(objective, point, vectorize=False)
    hessian = 0.5 * (hessian + hessian.T)
    full_eigenvalues = torch.linalg.eigvalsh(hessian)
    if free_mask is None:
        free_mask = torch.ones(len(values), dtype=torch.bool, device=values.device)
    else:
        free_mask = torch.as_tensor(free_mask, dtype=torch.bool, device=values.device)
    if free_mask.shape != values.shape:
        raise ValueError("Hessian free-parameter mask has inconsistent dimensions")
    free_indices = torch.where(free_mask)[0]
    covariance = torch.full_like(hessian, torch.nan)
    machine_floor = torch.as_tensor(torch.finfo(values.dtype).eps, dtype=values.dtype, device=values.device)
    if len(free_indices) > 0:
        free_hessian = hessian.index_select(0, free_indices).index_select(1, free_indices)
        free_eigenvalues, free_eigenvectors = torch.linalg.eigh(free_hessian)
        largest = torch.max(torch.abs(free_eigenvalues))
        eigenvalue_floor = torch.maximum(largest * 1.0e-10, machine_floor)
        regularized = torch.clamp(free_eigenvalues, min=eigenvalue_floor)
        free_covariance = (2.0 * float(errordef)) * (
            (free_eigenvectors * torch.reciprocal(regularized).unsqueeze(0)) @ free_eigenvectors.T
        )
        covariance[free_indices[:, None], free_indices[None, :]] = free_covariance
        regularized_count = int(torch.count_nonzero(free_eigenvalues < eigenvalue_floor).detach().cpu())
        free_positive_definite = bool(torch.all(free_eigenvalues > 0.0).detach().cpu())
        free_minimum = float(torch.min(free_eigenvalues).detach().cpu())
        free_maximum = float(torch.max(free_eigenvalues).detach().cpu())
    else:
        eigenvalue_floor = machine_floor
        regularized_count = 0
        free_positive_definite = True
        free_minimum = None
        free_maximum = None
    diagnostics = {
        "minimum_eigenvalue": float(torch.min(full_eigenvalues).detach().cpu()),
        "maximum_eigenvalue": float(torch.max(full_eigenvalues).detach().cpu()),
        "eigenvalue_floor": float(eigenvalue_floor.detach().cpu()),
        "regularized_eigenvalues": regularized_count,
        "positive_definite": bool(torch.all(full_eigenvalues > 0.0).detach().cpu()),
        "free_minimum_eigenvalue": free_minimum,
        "free_maximum_eigenvalue": free_maximum,
        "free_positive_definite": free_positive_definite,
        "free_dimension": int(len(free_indices)),
        "active_dimension": int(len(values) - len(free_indices)),
        "device": str(values.device),
        "dtype": str(values.dtype),
    }
    return hessian.detach(), covariance.detach(), diagnostics


# Build a shared fitted result after the Torch optimization and covariance steps
def build_result(
    param_names: Sequence[str],
    values: torch.Tensor,
    objective_value: torch.Tensor,
    covariance: torch.Tensor,
    success: bool,
    message: str,
) -> Result:
    fitted = values.detach().cpu().numpy().astype(np.float64, copy=False)
    covariance_np = covariance.detach().cpu().numpy().astype(np.float64, copy=False)
    errors = np.sqrt(np.clip(np.diag(covariance_np), 0.0, np.inf))
    return Result(
        values={str(name): float(fitted[index]) for index, name in enumerate(param_names)},
        covariance=covariance_np,
        errors={str(name): float(errors[index]) for index, name in enumerate(param_names)},
        fmin=Minimum(fval=float(objective_value.detach().cpu()), is_valid=bool(success), message=str(message)),
    )
