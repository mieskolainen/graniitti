# Shared NumPy and Torch arrays for differentiable histogram and amplitude equations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import sys

import numpy as np


# Select the numerical library without copying or detaching a tensor
def namespace(value):
    torch = sys.modules.get("torch")
    return torch if torch is not None and isinstance(value, getattr(torch, "Tensor", ())) else np


# Convert constants to the prediction device while preserving existing autograd graphs
def asarray(value, *, like=None, dtype=None):
    if namespace(value) is not np or namespace(like) is not np:
        reference = like if namespace(like) is not np else value
        torch = namespace(reference)
        dtype = {float: torch.float64, complex: torch.complex128, bool: torch.bool, int: torch.int64}.get(dtype, dtype)
        return torch.as_tensor(value, dtype=reference.dtype if dtype is None else dtype, device=reference.device)
    return np.asarray(value, dtype=dtype)


# Detach completed predictions when constructing numerical reports and plots
def to_numpy(value):
    return value.detach().cpu().numpy() if namespace(value) is not np else np.asarray(value)


# Stack scalar predictions without detaching any fitted coordinate
def stack(values):
    reference = next((value for value in values if namespace(value) is not np), None)
    return (
        np.asarray(values) if reference is None else namespace(reference).stack([asarray(value, like=reference) for value in values])
    )


# Concatenate predictions with the same numerical library and device
def concatenate(values):
    reference = next((value for value in values if namespace(value) is not np), None)
    arrays = [asarray(value, like=reference) for value in values]
    return namespace(reference).concatenate(arrays)


# Compute a nonnegative square root with a finite zero derivative at empty MC bins
def sqrt(value):
    xp = namespace(value)
    positive = value > 0.0
    return xp.where(positive, xp.sqrt(xp.where(positive, value, 1.0)), 0.0)[()]


# Combine uncertainties without overflow or an undefined autograd derivative at the origin
def hypot(first, second):
    first = asarray(first, like=second)
    xp = namespace(first)
    second = asarray(second, like=first)
    zero = (first == 0.0) & (second == 0.0)
    return xp.where(zero, 0.0, xp.hypot(xp.where(zero, 1.0, first), second))[()]
