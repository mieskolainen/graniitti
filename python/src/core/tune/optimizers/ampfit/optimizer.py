# Bounded Torch L-BFGS proposals from worker amplitude values and autograd gradients
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import torch

from core.numerics import lbfgsb
from core.tune.parameters.space import is_integer_bound, sample_config


# Suspend deterministic L-BFGS replay until a worker evaluates the next coordinate
class EvaluationRequired(Exception):
    # Preserve the physical optimizer coordinates requested by the line search
    def __init__(self, config):
        self.config = config
        super().__init__("Amplitude value and gradient require a worker trial")


# Rebuild bounded L-BFGS restarts from immutable completed worker evaluations
class AmplitudeSearch:
    # Resolve continuous bounds and deterministic restart coordinates once
    def __init__(self, *, bounds, initial_points, settings, seed):
        if any(is_integer_bound(bound) for bound in bounds.values()):
            raise ValueError("ampfit requires continuous active parameters")
        self.names = sorted(bounds)
        self.lower = torch.tensor([bounds[name]["lower"] for name in self.names], dtype=torch.float64)
        self.width = torch.tensor(
            [bounds[name]["upper"] - bounds[name]["lower"] for name in self.names], dtype=torch.float64)
        excluded = {"bank", "covariance", "covariance_interval", "starts", "start_relative_range"}
        self.settings = {key: value for key, value in settings.items() if key not in excluded}
        starts = settings["starts"]
        if not isinstance(starts, int) or starts < 1:
            raise ValueError("ampfit starts must be a positive integer")
        relative = settings["start_relative_range"]
        if initial_points is None:
            raise ValueError("ampfit relative starts require initial parameter values")
        # Perturb active continuous coordinates, retaining discrete choices in the model decoder
        nearby = {name: {**bound,
                         "lower": max(bound["lower"], initial_points[name] - relative * abs(initial_points[name])),
                         "upper": min(bound["upper"], initial_points[name] + relative * abs(initial_points[name]))}
                  for name, bound in bounds.items()}
        rng = np.random.default_rng(seed)
        self.starts = [dict(initial_points)] + [
            sample_config(nearby, uniform=rng.uniform, integer=rng.integers) for _ in range(starts - 1)]
        self.diagnostics = []
        self.finished = False

    # Compute a stable identity in the same physical optimizer coordinates as trial history
    def key(self, config):
        return tuple(float(config[name]) for name in self.names)

    # Propose at most one unevaluated line search point from each independent restart
    def ask(self, *, records, failures, pending, count, cost):
        completed = {self.key(record["config"]): record for record in records}
        failed = {self.key(record["config"]) for record in failures}
        reserved = {self.key(config) for config in pending}
        proposals = []
        self.diagnostics = []
        waiting = False

        # Read actual worker values and gradients without evaluating physics on the head
        def evaluate(unit):
            values = self.lower + self.width * unit
            config = dict(zip(self.names, values.tolist(), strict=True))
            key = self.key(config)
            if key in failed:
                return unit.new_tensor(torch.inf), torch.zeros_like(unit)
            record = completed.get(key)
            if record is None:
                raise EvaluationRequired(config)
            derivative = record.get("gradient")
            if not isinstance(derivative, dict) or set(derivative) != set(self.names):
                raise ValueError("ampfit history requires an autograd gradient for every active parameter")
            gradient = unit.new_tensor([derivative[name] for name in self.names]) * self.width
            value = unit.new_tensor(record["metrics"][cost])
            if not torch.isfinite(value) or not torch.all(torch.isfinite(gradient)):
                raise ValueError("Completed amplitude value and gradient must be finite")
            return value, gradient

        for index, start in enumerate(self.starts):
            unit = (self.lower.new_tensor(self.key(start)) - self.lower) / self.width
            try:
                _, value, _, diagnostics = lbfgsb.minimize(None, unit, evaluator=evaluate, **self.settings)
                if not torch.isfinite(value):
                    diagnostics.update(success=False, message="starting amplitude evaluation failed")
                self.diagnostics.append({"start": index, **diagnostics})
            except EvaluationRequired as request:
                key = self.key(request.config)
                waiting = True
                if key not in reserved and len(proposals) < count:
                    initial = index == 0 and not records and not failures
                    proposals.append((request.config, {"kind": "initial" if initial else "ampfit", "start": index}))
                    reserved.add(key)
        self.finished = not waiting
        return proposals
