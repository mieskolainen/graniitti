# Joint logarithmic noisy expected improvement for icebo
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from dataclasses import dataclass

import torch

from core.inference.gp import ExactGaussianProcess, robust_cholesky
from core.numerics import lbfgsb
from core.numerics.linalg import solve_triangular
from core.numerics.smooth import fatmax, log_fatplus
from core.tune.optimizers.icebo.space import IceboSpace


@dataclass(frozen=True)
class QLogNEISettings:
    """Complete numerical settings for joint qLogNEI maximization"""

    mc_samples: int
    raw_samples: int
    restarts: int
    gradient_steps: int
    optimization_batch_size: int
    mc_cholesky_jitter: float
    lbfgs_tolerance_change: float
    lbfgs_tolerance_grad: float
    tau_max: float
    tau_relu: float
    duplicate_tolerance: float
    prune_baseline: bool


# Transform uniform Sobol points into standard Normal samples
def sobol_normal(count: int, dimension: int, *, seed: int, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
    engine = torch.quasirandom.SobolEngine(dimension=max(1, int(dimension)), scramble=True, seed=int(seed))
    uniform = engine.draw(int(count), dtype=dtype).to(device=device)
    epsilon = torch.finfo(dtype).eps
    uniform = torch.clamp(uniform, min=epsilon, max=1.0 - epsilon)
    return math.sqrt(2.0) * torch.erfinv(2.0 * uniform - 1.0)


class QLogNoisyExpectedImprovement:
    """Reparameterized joint qLogNEI with exact pending conditioning"""

    # Initialize fixed QMC samples over observed and pending reference points
    def __init__(
        self, surrogate: ExactGaussianProcess, *, pending: torch.Tensor, q: int, settings: QLogNEISettings, seed: int
    ) -> None:
        self.surrogate = surrogate
        self.settings = settings
        self.q = int(q)
        reference = surrogate.X
        if len(pending):
            reference = torch.cat((reference, pending), dim=0)
        with torch.no_grad():
            mean, covariance = surrogate.posterior_covariance(reference)
            cholesky, _ = robust_cholesky(
                covariance.detach(),
                attempts=surrogate.settings["cholesky_attempts"],
                initial_jitter=surrogate.settings["initial_jitter"],
                jitter_multiplier=surrogate.settings["jitter_multiplier"],
            )
            if settings.prune_baseline:
                normals = sobol_normal(settings.mc_samples, len(reference), seed=seed + 1,
                                       device=surrogate.device, dtype=surrogate.dtype)
                samples = mean + normals @ cholesky.T
                selected = samples[:, : len(surrogate.X)].argmin(dim=-1).unique()
                selected = torch.cat((selected, torch.arange(len(surrogate.X), len(reference), device=surrogate.device)))
                reference, mean = reference[selected], mean[selected]
                covariance = covariance[selected][:, selected]
                cholesky, _ = robust_cholesky(
                    covariance, attempts=surrogate.settings["cholesky_attempts"],
                    initial_jitter=surrogate.settings["initial_jitter"],
                    jitter_multiplier=surrogate.settings["jitter_multiplier"],
                )
            self.reference, self.reference_mean, self.reference_cholesky = reference, mean, cholesky
            self.reference_solve = torch.cholesky_solve(surrogate.kernel(surrogate.X, reference), surrogate.cholesky)
        normals = sobol_normal(
            settings.mc_samples, len(reference) + self.q, seed=seed, device=surrogate.device, dtype=surrogate.dtype
        )
        self.reference_normals = normals[:, : len(reference)]
        self.candidate_normals = normals[:, len(reference) :]
        self.reference_samples = self.reference_mean + (self.reference_normals @ self.reference_cholesky.T)
        self.reference_best = torch.min(self.reference_samples, dim=-1).values

    # Evaluate exact conditional candidate samples for many proposed batches
    def _conditional_samples(self, candidates: torch.Tensor) -> torch.Tensor:
        mean, covariance = self.surrogate.posterior_covariance(candidates)
        cross = self.surrogate.kernel(candidates, self.reference)
        cross = cross - self.surrogate.kernel(candidates, self.surrogate.X) @ self.reference_solve
        right = cross.transpose(-1, -2)
        weights = solve_triangular(self.reference_cholesky, right)
        conditional = covariance - weights.transpose(-1, -2) @ weights
        identity = torch.eye(self.q, device=conditional.device, dtype=conditional.dtype)
        conditional = conditional + self.settings.mc_cholesky_jitter * identity
        conditional_factor, _ = robust_cholesky(
            conditional,
            attempts=self.surrogate.settings["cholesky_attempts"],
            initial_jitter=self.surrogate.settings["initial_jitter"],
            jitter_multiplier=self.surrogate.settings["jitter_multiplier"],
        )
        reference_shift = torch.einsum("sn,bnq->bsq", self.reference_normals, weights)
        innovation = torch.einsum("sq,bkq->bsk", self.candidate_normals, conditional_factor)
        return mean[:, None, :] + reference_shift + innovation

    # Compute one qLogNEI value per candidate batch
    def __call__(self, candidates: torch.Tensor) -> torch.Tensor:
        points = candidates.unsqueeze(0) if candidates.ndim == 2 else candidates
        samples = self._conditional_samples(points)
        improvement = self.reference_best[None, :, None] - samples
        log_improvement = log_fatplus(improvement, tau=self.settings.tau_relu)
        log_improvement = fatmax(log_improvement, dim=-1, tau=self.settings.tau_max)
        return torch.logsumexp(log_improvement, dim=-1) - math.log(self.settings.mc_samples)


class QLogNEIOptimizer:
    """Deterministic GPU batched multistart maximizer for joint qLogNEI"""

    # Initialize one acquisition optimizer
    def __init__(self, space: IceboSpace, *, settings: QLogNEISettings, seed: int) -> None:
        self.space = space
        self.settings = settings
        self.seed = int(seed)

    # Draw complete q batches from one scrambled Sobol design
    def _raw_batches(self, q: int, lower: torch.Tensor, upper: torch.Tensor, iteration: int) -> torch.Tensor:
        engine = torch.quasirandom.SobolEngine(
            q * self.space.dimension, scramble=True, seed=self.seed + 7919 * int(iteration) + 1
        )
        base = engine.draw(self.settings.raw_samples, dtype=self.space.dtype).to(self.space.device)
        base = base.reshape(self.settings.raw_samples, q, self.space.dimension)
        raw = lower + (upper - lower) * base
        canonical = self.space.canonical_unit(raw)
        return torch.where(self.space.periodic_mask, raw, canonical)

    # Evaluate acquisition batches with bounded GPU memory
    def _scores(self, acquisition: QLogNoisyExpectedImprovement, batches: torch.Tensor) -> torch.Tensor:
        scores = []
        step = self.settings.optimization_batch_size
        for start in range(0, len(batches), step):
            physical = self.space.surrogate_inputs_from_unit(batches[start : start + step])
            scores.append(acquisition(physical))
        return torch.cat(scores)

    # Refine independent restarts with the shared bounded Torch L-BFGS solver
    def _refine(
        self, seeds: torch.Tensor, acquisition: QLogNoisyExpectedImprovement, lower: torch.Tensor, upper: torch.Tensor
    ) -> torch.Tensor:
        if self.settings.gradient_steps == 0 or bool(self.space.integer_mask.all()):
            return seeds.detach().clone()
        refined = []
        span = upper - lower
        free = ~self.space.integer_mask & (span > 0.0)
        scale = torch.where(free, span, torch.ones_like(span))
        for batch in seeds.split(self.settings.optimization_batch_size):
            # Map local unit coordinates to the trust region and preserve integer values
            def decode(values, batch=batch):
                return torch.where(free, lower + span * values.reshape_as(batch), batch)

            # Minimize the negative acquisition independently for each restart
            def objective(values):
                return -acquisition(self.space.surrogate_inputs_from_unit(decode(values)))

            candidates, _, _, _ = lbfgsb.minimize_batch(
                objective, ((batch - lower) / scale).flatten(start_dim=1),
                maxiter=self.settings.gradient_steps, gradient_tolerance=self.settings.lbfgs_tolerance_grad,
                function_tolerance=self.settings.lbfgs_tolerance_change,
            )
            refined.append(decode(candidates).detach())
        return torch.cat(refined)

    # Maximize one joint qLogNEI acquisition
    def maximize(
        self,
        surrogate: ExactGaussianProcess,
        *,
        q: int,
        pending: torch.Tensor,
        references: torch.Tensor | None,
        lower: torch.Tensor,
        upper: torch.Tensor,
        iteration: int,
    ) -> tuple[torch.Tensor, dict]:
        acquisition = QLogNoisyExpectedImprovement(
            surrogate,
            pending=self.space.surrogate_inputs_from_unit(pending),
            q=q,
            settings=self.settings,
            seed=self.seed + 3571 * int(iteration) + 17,
        )
        raw = self._raw_batches(q, lower, upper, iteration)
        with torch.no_grad():
            raw_scores = self._scores(acquisition, raw)
            # Preserve the best start and sample Boltzmann restarts [REFERENCE: botorch.optim.initializers.initialize_q_batch]
            valid = torch.nonzero(torch.isfinite(raw_scores), as_tuple=True)[0]
            if not len(valid):
                raise RuntimeError("ICEBO qLogNEI search found no finite acquisition values")
            values = raw_scores[valid]
            restart_count = min(self.settings.restarts, len(valid))
            scale = values.std(unbiased=False).clamp_min(torch.finfo(values.dtype).eps)
            weights = torch.softmax((values - values.max()) / scale, dim=0).clamp_min(torch.finfo(values.dtype).tiny)
            generator = torch.Generator(device=raw.device).manual_seed(self.seed + int(iteration))
            chosen = torch.multinomial(weights, restart_count, replacement=False, generator=generator)
            best = torch.argmax(values)
            if not bool(torch.any(chosen == best)):
                chosen[-1] = best
            seeds = raw[valid[chosen]]
        refined = self._refine(seeds, acquisition, lower, upper)
        combined = torch.cat((refined, raw), dim=0)
        with torch.no_grad():
            scores = torch.cat((self._scores(acquisition, refined), raw_scores))
            finite = torch.nonzero(torch.isfinite(scores), as_tuple=True)[0]
            order = finite[torch.argsort(scores[finite], descending=True)]
            for index in order.tolist():
                selected = self.space.canonical_unit(combined[index])
                if self._is_valid_batch(selected, references):
                    return selected, {
                        "acquisition": "qlognei",
                        "log_acquisition": float(scores[index].detach().cpu()),
                        "mc_samples": self.settings.mc_samples,
                        "raw_batches": len(raw),
                        "gradient_restarts": restart_count,
                        "joint_batch_size": q,
                    }
        raise RuntimeError("ICEBO qLogNEI search exhausted the finite parameter space")

    # Check completed, pending and within-batch topology duplicates
    def _is_valid_batch(self, candidates: torch.Tensor, references: torch.Tensor | None) -> bool:
        active = references
        for row in candidates:
            point = row.unsqueeze(0)
            distance = self.space.minimum_topology_distance(point, active)
            if float(distance[0].detach().cpu()) <= self.settings.duplicate_tolerance:
                return False
            active = point if active is None else torch.cat((active, point), dim=0)
        return True
