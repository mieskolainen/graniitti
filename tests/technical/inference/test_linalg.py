# Tests for shared triangular solves and their derivatives
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest
import torch
from core.numerics.linalg import solve_triangular


# Check noncontiguous batched right hand sides against native Torch values and derivatives
@pytest.mark.parametrize("batch", [(), (3,), (2, 3)])
@pytest.mark.parametrize("upper", [False, True])
def test_shared_triangular_solve(batch, upper):
    generator = torch.Generator().manual_seed(19)
    matrix = torch.randn(4, 4, generator=generator, dtype=torch.float64)
    factor = torch.linalg.cholesky(matrix @ matrix.T + torch.eye(4))
    factor = (factor.T if upper else factor).detach().requires_grad_()
    values = torch.randn(*batch, 2, 4, generator=generator, dtype=torch.float64).transpose(-1, -2).requires_grad_()
    actual = solve_triangular(factor, values, upper=upper)
    expected = torch.linalg.solve_triangular(factor, values, upper=upper)
    torch.testing.assert_close(actual, expected)
    first = torch.autograd.grad(actual.square().sum(), (factor, values), create_graph=True)
    reference = torch.autograd.grad(expected.square().sum(), (factor, values), create_graph=True)
    for derivative, target in zip(first, reference, strict=True):
        torch.testing.assert_close(derivative, target)
    second = torch.autograd.grad(sum(value.square().sum() for value in first), (factor, values))
    reference_second = torch.autograd.grad(sum(value.square().sum() for value in reference), (factor, values))
    for derivative, target in zip(second, reference_second, strict=True):
        torch.testing.assert_close(derivative, target)
