# Tests for algebraic smooth positive parts and maxima
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
import torch
from core.numerics.smooth import fatmax, log_fatplus


# Check value bounds, positive gradients and finite extreme tails
@pytest.mark.parametrize("dtype", [torch.float32, torch.float64])
def test_log_fatplus_bounds_and_tails(dtype):
    values = torch.tensor([-1e30, -1e6, -1.0, 0.0, 1.0, 1e6, 1e30], dtype=dtype, requires_grad=True)
    tau = 1e-6
    result = log_fatplus(values, tau=tau)
    gradient, = torch.autograd.grad(result.sum(), values)
    assert torch.isfinite(result).all() and torch.isfinite(gradient).all()
    assert (gradient > 0.0).all()
    torch.testing.assert_close(result[-2:], values[-2:].log())
    assert result[3].exp().item() == pytest.approx(tau * (math.log(2.0) + 0.1), rel=1e-6)
    assert result[0].item() == pytest.approx(math.log(0.1 * tau**3) - 2 * math.log(1e30), rel=1e-6)


# Check the analytic smoothing bound and shift invariance of the maximum
@pytest.mark.parametrize("q", [1, 2, 16])
def test_fatmax_bounds_and_gradients(q):
    values = torch.linspace(-2.0, 1.0, 3 * q, dtype=torch.float64).reshape(3, q).requires_grad_()
    tau = 0.01
    result = fatmax(values, tau=tau)
    maximum = values.amax(dim=-1)
    assert torch.all(result >= maximum - 1e-14)
    assert torch.all(result <= maximum + tau * math.log(q) + 1e-14)
    torch.testing.assert_close(fatmax(values + 10.0, tau=tau), result + 10.0)
    gradient, = torch.autograd.grad(result.sum(), values)
    assert (gradient > 0.0).all()
    torch.testing.assert_close(gradient.sum(dim=-1), torch.ones_like(result))
    assert torch.autograd.gradcheck(lambda x: fatmax(x, tau=tau), (values,))


# Check the soft positive part through the transition and both tails
def test_log_fatplus_gradient():
    values = torch.tensor([-40.0, -2.0, 0.0, 2.0, 40.0], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(log_fatplus, (values,))
