# Chebyshev interpolation values and autograd derivatives at and between grid nodes
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import torch
from core.numerics.interp import ChebyshevGrid


# Compare polynomial values and derivatives with NumPy, including interior grid nodes
def test_chebyshev_grid_autograd():
    grid = ChebyshevGrid(0.0, 2.5, 7)
    polynomial = np.polynomial.Chebyshev([0.2, 0.7, -0.3, 0.0, 0.1, 0.03], domain=[0.0, 2.5])
    nodes = grid.nodes.numpy()
    values = polynomial(nodes)
    points = torch.tensor([0.17, nodes[1], 0.5, nodes[3], 1.4, 1.8], dtype=torch.float64, requires_grad=True)
    prediction = grid(points) @ torch.as_tensor(values)
    derivative, = torch.autograd.grad(prediction.sum(), points)
    np.testing.assert_allclose(prediction.detach(), polynomial(points.detach().numpy()), atol=1e-12)
    np.testing.assert_allclose(derivative, polynomial.deriv()(points.detach().numpy()), atol=1e-12)
    assert torch.autograd.gradcheck(grid, (points,))
    assert torch.autograd.gradgradcheck(grid, (points,))
    torch.testing.assert_close(grid(torch.as_tensor(nodes)), torch.eye(len(nodes), dtype=torch.float64))
    assert torch.isnan(grid(torch.tensor(-0.1))).all()


# Preserve exponential tails and their derivatives at and between native grid nodes
def test_exponential_grid_autograd():
    grid = ChebyshevGrid(0.0, 1.5, 5)
    scale = torch.tensor([0.0, 4.0, 16.0], dtype=torch.float64)
    samples = torch.exp(-scale[:, None] * grid.nodes) * (1 + grid.nodes.square())

    # Interpolate a polynomial residual multiplying a known exponential
    def evaluate(value):
        return (grid.exponential(value, scale) * samples).sum(-1)

    for index, value in enumerate([*grid.nodes, torch.tensor(0.37, dtype=torch.float64)]):
        expected = torch.exp(-scale * value) * (1 + value.square())
        torch.testing.assert_close(evaluate(value), expected, atol=1e-12, rtol=1e-12)
        point = value.detach().clone().requires_grad_()
        derivative, = torch.autograd.grad(evaluate(point).sum(), point)
        expected_derivative = (torch.exp(-scale * value) * (2 * value - scale * (1 + value.square()))).sum()
        torch.testing.assert_close(derivative, expected_derivative, atol=1e-12, rtol=1e-12)
        if index not in (0, len(grid.nodes) - 1):
            assert torch.autograd.gradcheck(evaluate, (point,))
            assert torch.autograd.gradgradcheck(evaluate, (point,))
