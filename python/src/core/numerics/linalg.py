# Torch linear solves with a shared factor across independent batches
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import torch


# Solve all right hand sides together without expanding the shared triangular factor
def solve_triangular(factor: torch.Tensor, values: torch.Tensor, *, upper: bool = False) -> torch.Tensor:
    right = values.movedim(-2, 0)
    solved = torch.linalg.solve_triangular(factor, right.reshape(factor.shape[0], -1), upper=upper)
    return solved.reshape(right.shape).movedim(0, -2)
