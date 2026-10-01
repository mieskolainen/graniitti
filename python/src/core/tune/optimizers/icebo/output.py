# Affine objective standardization for icebo
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math

import torch


class StandardizeOutput:
    """Exact affine standardization with variance propagation"""

    # Initialize an unfitted affine transformation
    def __init__(self) -> None:
        self.location = 0.0
        self.scale = 1.0

    # Fit sample mean and population standard deviation
    def fit(self, values: torch.Tensor) -> None:
        if values.numel() == 0 or not bool(torch.all(torch.isfinite(values))):
            raise ValueError("ICEBO output standardization requires finite observations")
        location = torch.mean(values)
        scale = torch.std(values, unbiased=False)
        minimum = math.sqrt(torch.finfo(values.dtype).eps)
        self.location = float(location.detach().cpu())
        self.scale = max(float(scale.detach().cpu()), minimum)

    # Map physical objectives into GP units
    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return (values - self.location) / self.scale

    # Map physical variances into GP units
    def variance_forward(self, variances: torch.Tensor) -> torch.Tensor:
        return variances / self.scale**2

    # Map GP means into physical objective units
    def inverse(self, values: torch.Tensor) -> torch.Tensor:
        return self.location + self.scale * values

    # Map GP variances into physical objective units
    def variance_inverse(self, variances: torch.Tensor) -> torch.Tensor:
        return self.scale**2 * variances

    # Compute a JSON safe transformation description
    def metadata(self) -> dict:
        return {"kind": "standardize", "location": self.location, "scale": self.scale}
