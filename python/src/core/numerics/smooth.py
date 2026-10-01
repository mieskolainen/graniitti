# Smooth positive part and maximum with algebraic tails
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import torch


# Compute the log of a softplus plus Cauchy tail [REFERENCE: arXiv:2310.20708]
def log_fatplus(values: torch.Tensor, tau: float = 1.0) -> torch.Tensor:
    scaled = values / tau
    cutoff = math.log(torch.finfo(values.dtype).eps)
    log_softplus = torch.where(
        scaled < cutoff, scaled, torch.nn.functional.softplus(scaled.clamp_min(cutoff)).log()
    )
    log_cauchy = math.log(0.1) - 2.0 * torch.hypot(scaled, torch.ones_like(scaled)).log()
    return math.log(tau) + torch.logaddexp(log_softplus, log_cauchy)


# Compute the smooth maximum with quadratic Pareto tails [REFERENCE: arXiv:2310.20708]
def fatmax(values: torch.Tensor, dim: int = -1, tau: float = 1.0) -> torch.Tensor:
    maximum = values.amax(dim=dim, keepdim=True)
    distance = (maximum - values) / tau
    log_tail = math.log(2.0) - 2.0 * torch.hypot(distance + 1.0, torch.ones_like(distance)).log()
    return maximum.squeeze(dim) + tau * torch.logsumexp(log_tail, dim=dim)
