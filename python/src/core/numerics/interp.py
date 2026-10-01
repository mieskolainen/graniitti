# Differentiable polynomial interpolation on Chebyshev grids
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import torch


# Store cardinal polynomials with extra resolution near the parameter bounds
class ChebyshevGrid:
    # Construct a finite Chebyshev axis including both parameter bounds
    def __init__(self, lower, upper, count):
        if not np.isfinite([lower, upper]).all() or upper <= lower or type(count) is not int or count < 2:
            raise ValueError("Chebyshev interpolation requires increasing bounds and at least two nodes")
        self.lower, self.upper = float(lower), float(upper)
        unit = -np.cos(np.linspace(0.0, np.pi, count))
        self.nodes = torch.as_tensor(lower + (upper - lower) * (unit + 1.0) / 2.0)
        self.coefficients = torch.as_tensor(np.polynomial.chebyshev.chebfit(unit, np.eye(count), count - 1))

    # Evaluate cardinal polynomials with Clenshaw recursion, including derivatives at nodes
    def __call__(self, value):
        value = torch.as_tensor(value, dtype=self.nodes.dtype)
        unit = (2.0 * (value - self.lower) / (self.upper - self.lower) - 1.0)[..., None]
        coefficients = self.coefficients.to(value.device)
        first = torch.zeros((*value.shape, len(self.nodes)), dtype=value.dtype, device=value.device)
        second = torch.zeros_like(first)
        for coefficient in reversed(coefficients[1:]):
            first, second = coefficient + 2.0 * unit * first - second, first
        weights = coefficients[0] + unit * first - second
        valid = (value >= self.lower) & (value <= self.upper)
        return torch.where(valid[..., None], weights, torch.nan)

    # Interpolate the residual after extracting an exponential dependence
    def exponential(self, value, scale):
        scale = torch.as_tensor(scale, dtype=self.nodes.dtype)
        value = torch.as_tensor(value, dtype=self.nodes.dtype, device=scale.device)
        weights = self(value)
        return weights * torch.exp((self.nodes.to(scale.device) - value[..., None]) * scale[..., None])
