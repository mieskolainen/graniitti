# Differentiable parameter maps and full covariance propagation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from collections.abc import Callable, Mapping, Sequence

import numpy as np
import torch
from scipy.optimize import minimize


# Compute standard deviations and correlations, leaving undefined entries unavailable
def covariance_summary(covariance):
    covariance = np.asarray(covariance, dtype=float)
    variance = np.diag(covariance)
    errors = np.sqrt(np.where(np.isfinite(variance) & (variance >= 0.0), variance, np.nan))
    denominator = np.outer(errors, errors)
    correlation = np.divide(covariance, denominator, out=np.full_like(covariance, np.nan),
                            where=np.isfinite(covariance) & np.isfinite(denominator) & (denominator > 0.0))
    return errors, correlation


# Keep named decoder outputs and their derivatives in the same order
class ParameterTransform:
    # Bind a differentiable named decoder at one reference input
    def __init__(self, names: Sequence[str], decode: Callable, reference, periods: Mapping | None = None):
        self.input_names = tuple(names)
        self.decode = decode
        point = torch.as_tensor(reference, dtype=torch.float64)
        self.names = tuple(decode(dict(zip(self.input_names, point, strict=True))))
        if not self.names or len(set(self.input_names)) != len(self.input_names):
            raise ValueError("Parameter transform requires nonempty outputs and distinct input names")
        self.periods = np.array([float((periods or {}).get(name, 0.0)) for name in self.names])
        if not np.all(np.isfinite(self.periods)) or np.any(self.periods < 0.0):
            raise ValueError("Parameter periods must be finite and nonnegative")

    # Evaluate the original decoder without detaching tensor values
    def __call__(self, point):
        point = torch.as_tensor(point, dtype=torch.float64)
        if point.shape != (len(self.input_names),):
            raise ValueError("Parameter transform input has inconsistent dimensions")
        values = self.decode(dict(zip(self.input_names, point, strict=True)))
        if tuple(values) != self.names:
            raise ValueError("Parameter transform changed its output names")
        return torch.stack([torch.as_tensor(values[name], dtype=point.dtype, device=point.device) for name in self.names])

    # Compute the local displacement in each output coordinate, including periodic coordinates
    def difference(self, values, reference):
        values = torch.as_tensor(values, dtype=torch.float64)
        reference = torch.as_tensor(reference, dtype=values.dtype, device=values.device)
        delta = values - reference
        periods = torch.as_tensor(self.periods, dtype=delta.dtype, device=delta.device)
        scale = torch.where(periods > 0.0, periods, 1.0)
        angle = delta * (2.0 * torch.pi) / scale
        wrapped = torch.atan2(torch.sin(angle), torch.cos(angle)) * scale / (2.0 * torch.pi)
        return torch.where(periods > 0.0, wrapped, delta)

    # Compute output values and the full autograd Jacobian
    def linearize(self, point):
        point = torch.as_tensor(point, dtype=torch.float64)
        values = self(point)
        jacobian = torch.autograd.functional.jacobian(self, point)
        return values.detach().cpu().numpy(), jacobian.detach().cpu().numpy()

    # Propagate the full input covariance with the local autograd Jacobian
    def propagate(self, point, covariance):
        values, jacobian = self.linearize(point)
        covariance = torch.as_tensor(covariance, dtype=torch.float64).detach().cpu().numpy()
        if covariance.shape != (len(self.input_names), len(self.input_names)):
            raise ValueError("Parameter transform covariance has inconsistent dimensions")
        active = np.abs(jacobian) > 0.0
        missing = active.astype(int) @ (~np.isfinite(covariance)).astype(int) @ active.T.astype(int) > 0
        transformed = jacobian @ np.where(np.isfinite(covariance), covariance, 0.0) @ jacobian.T
        transformed[missing] = np.nan
        errors, correlation = covariance_summary(transformed)
        return {"names": list(self.names), "values": values, "jacobian": jacobian,
                "covariance": transformed, "errors": errors, "correlation": correlation,
                "method": "autograd_delta", "periods": self.periods.copy()}

    # Transform samples exactly through the decoder instead of linearizing them
    def samples(self, points):
        with torch.no_grad():
            return torch.stack([self(point) for point in points]).cpu().numpy()

    # Constrain one decoded coordinate within the covariance support and original input bounds
    def constrain(self, point, covariance, bounds, index, displacement, *, maxiter, objective=None):
        point = torch.as_tensor(point, dtype=torch.float64).detach().cpu()
        failure = {"point": point.numpy(), "success": False, "iterations": 0, "evaluations": 0}
        displacement = torch.as_tensor(displacement, dtype=point.dtype).detach().cpu()
        covariance = torch.as_tensor(covariance, dtype=point.dtype).detach().cpu()
        known = torch.isfinite(torch.diag(covariance))
        if not torch.all(torch.isfinite(covariance[known][:, known])):
            return {**failure, "message": "Incomplete covariance"}
        covariance = torch.where(known[:, None] & known[None, :], covariance, 0.0)
        bounds = torch.as_tensor(bounds, dtype=point.dtype).detach().cpu()
        eigenvalues, eigenvectors = torch.linalg.eigh(covariance)
        tolerance = torch.finfo(point.dtype).eps * len(point) * eigenvalues.abs().max()
        active = eigenvalues > tolerance
        if torch.any(eigenvalues < -tolerance):
            return {**failure, "message": "Covariance is not positive semidefinite"}
        if not torch.any(active):
            return {**failure, "message": "No covariance support"}
        inverse = eigenvectors[:, active] / torch.sqrt(eigenvalues[active])
        null = eigenvectors[:, ~active].T
        width = bounds[:, 1] - bounds[:, 0]
        width = torch.where(width > 0.0, width, 1.0)
        unit_bounds = np.column_stack((np.zeros(len(point)), ((bounds[:, 1] - bounds[:, 0]) / width).numpy()))
        reference = self(point).detach()

        # Evaluate either the covariance distance or the supplied differentiable likelihood
        def value_gradient(unit):
            unit = torch.as_tensor(unit, dtype=point.dtype).requires_grad_()
            x = bounds[:, 0] + width * unit
            if objective is None:
                value = 0.5 * ((x - point) @ inverse).square().sum()
                gradient, = torch.autograd.grad(value, unit)
                return value.item(), gradient.detach().numpy()
            value, gradient = objective(x.detach().numpy())
            return float(value), (width * torch.as_tensor(gradient, dtype=point.dtype)).numpy()

        # Evaluate the physical equality while preserving the covariance support
        def constraints(unit):
            x = bounds[:, 0] + width * unit
            residual = self.difference(self(x), reference)[index] - displacement
            return torch.cat((residual.reshape(-1), null @ (x - point)))

        # Compute constraint values and derivatives through the same decoder graph
        def constraint_values(unit):
            return constraints(torch.as_tensor(unit, dtype=point.dtype)).detach().numpy()

        # Compute the complete autograd constraint Jacobian
        def constraint_jacobian(unit):
            return torch.autograd.functional.jacobian(constraints, torch.as_tensor(unit, dtype=point.dtype)).numpy()

        initial = ((point - bounds[:, 0]) / width).numpy()
        result = minimize(value_gradient, initial, jac=True, method="SLSQP", bounds=unit_bounds,
                          options={"maxiter": maxiter},
                          constraints={"type": "eq", "fun": constraint_values, "jac": constraint_jacobian})
        return {"point": (bounds[:, 0] + width * torch.as_tensor(result.x)).numpy(),
                "success": bool(result.success), "message": str(result.message), "value": float(result.fun),
                "iterations": int(result.nit), "evaluations": int(result.nfev)}
