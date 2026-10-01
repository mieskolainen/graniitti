# HEBO wrappers with phase topology treatment
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
import random

import numpy as np
import pandas as pd
import torch
from gpytorch.constraints import Positive
from gpytorch.kernels import Kernel, ScaleKernel
from gpytorch.priors import GammaPrior
from hebo.models.base_model import BaseModel
from hebo.models.gp.gp import GP
from hebo.models.gp.svgp import SVGP
from hebo.models.scalers import TorchIdentityScaler
from hebo.optimizers.hebo import HEBO

from core.io import logger as log
from core.tune.optimizers.hebo.config import load_hebo_config, validate_hebo_config
from core.tune.parameters.kernel import ProductManifoldGeometry

TOPOLOGY_MODEL_NAME = "icetune_topology_gp"
HYBRID_MODEL_NAME = "icetune_hybrid_gp"
logger = log.get_logger(__name__)


# Compute common exact and sparse HEBO model settings
def _hybrid_model_config(settings: dict) -> dict:
    return {**settings["gp"], "sparse_config": dict(settings["sparse_gp"])}


# Extract surrogate-aware manifold groups from topology metadata
def surrogate_parameter_groups(parameter_topology: dict | None) -> list[dict]:
    if not parameter_topology:
        return []
    if int(parameter_topology.get("schema_version", -1)) != 1:
        raise ValueError("Unsupported HEBO parameter-topology schema")

    groups = []
    seen = set()
    for group in parameter_topology.get("groups", []):
        kind = str(group.get("kind"))
        if kind not in {
            "phase",
            "projective",
            "sphere",
            "simplex",
            "polar_projective",
            "cartesian_projective",
            "phase_projective",
            "polar_sphere",
            "polar_components",
            "radial_projective",
        }:
            continue
        parameters = [str(name) for name in group.get("parameters", [])]
        if kind == "phase" and len(parameters) != 1:
            raise ValueError(f"HEBO phase topology requires one parameter per group: {group}")
        period = group.get("period")
        if kind in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"} and (
            isinstance(period, bool)
            or not isinstance(period, (int, float))
            or not math.isfinite(float(period))
            or float(period) <= 0.0
        ):
            raise ValueError(f"HEBO phase topology requires a finite positive period: {group}")
        if kind in {"projective", "sphere", "simplex"} and not parameters:
            raise ValueError(f"HEBO {kind} topology requires at least one parameter: {group}")
        if kind in {"polar_projective", "cartesian_projective"} and len(parameters) < 2:
            raise ValueError(f"HEBO {kind} topology requires coefficient coordinates: {group}")
        if kind in {"polar_sphere", "polar_components"} and len(parameters) < 2:
            raise ValueError(f"HEBO {kind} topology requires phase and coupling coordinates: {group}")
        if kind == "radial_projective" and len(parameters) < 2:
            raise ValueError(f"HEBO {kind} topology requires a norm and direction: {group}")
        if kind == "phase_projective" and not parameters:
            raise ValueError(f"HEBO {kind} topology requires a phase coordinate: {group}")
        if any(name in seen for name in parameters):
            raise ValueError("HEBO parameter-topology groups overlap")
        seen.update(parameters)
        entry = {"kind": kind, "base": str(group.get("base", parameters[0])), "parameters": parameters}
        if kind in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"}:
            entry["period"] = float(period)
        if "weights" in group:
            entry["weights"] = [float(value) for value in group["weights"]]
        groups.append(entry)
    return groups


class TopologyMaternKernel(Kernel):
    """Matern kernel using direct circular and constrained vector chordal distances."""

    # Initialize groupwise length scales without expanding the HEBO inputs
    def __init__(
        self,
        num_cont: int,
        topology_groups: tuple[dict, ...],
        numeric_bounds: tuple[tuple[float, float], ...],
        *,
        nu: float = 1.5,
        lengthscale_scale: float = 1.0,
        ard: bool = True,
    ):
        super().__init__()
        self.geometry = ProductManifoldGeometry(num_cont, topology_groups, numeric_bounds)
        self.effective_dim = self.geometry.effective_dim
        self.matern52 = nu > 2.0
        if self.effective_dim == 0:
            raise ValueError("HEBO topology kernel requires at least one parameter group")
        constraint = Positive()
        initial_scale = torch.full((1, self.effective_dim if ard else 1), lengthscale_scale * math.sqrt(self.effective_dim))
        parameter = torch.nn.Parameter(constraint.inverse_transform(initial_scale))
        self.register_parameter("raw_group_lengthscale", parameter)
        self.register_constraint("raw_group_lengthscale", constraint)

    # Compute positive length scales for ordinary coordinates and topology groups
    @property
    def group_lengthscale(self) -> torch.Tensor:
        return self.raw_group_lengthscale_constraint.transform(self.raw_group_lengthscale).expand(1, self.effective_dim)

    # Evaluate the selected radial Matern kernel on the product manifold distance
    def forward(self, x1: torch.Tensor, x2: torch.Tensor, diag: bool = False, **params):
        squared = self.geometry.squared_distance(x1, x2, self.group_lengthscale, diag=diag)
        distance = torch.sqrt(squared.clamp_min(1.0e-30))
        if self.matern52:
            scaled = math.sqrt(5.0) * distance
            return (1.0 + scaled + scaled.square() / 3.0) * torch.exp(-scaled)
        scaled = math.sqrt(3.0) * distance
        return (1.0 + scaled) * torch.exp(-scaled)


class HybridGP(BaseModel):
    """HEBO GP which becomes sparse above a fixed observation count."""

    support_grad = True

    # Initialize the exact model and retain the sparse model settings
    def __init__(self, num_cont: int, num_enum: int, num_out: int, **conf):
        defaults = (
            _hybrid_model_config(load_hebo_config())
            if any(key not in conf for key in ("max_exact_trials", "sparse_config"))
            else {}
        )
        conf = {**defaults, **conf}
        super().__init__(num_cont, num_enum, num_out, **conf)
        self.max_exact_trials = int(conf.pop("max_exact_trials"))
        if self.max_exact_trials < 0:
            raise ValueError("HEBO exact model limit must be non-negative")
        self.synthetic_trials = int(conf.pop("synthetic_trials", 0))
        if self.synthetic_trials < 0:
            raise ValueError("HEBO synthetic trial count must be non-negative")
        self.sparse_config = dict(conf.pop("sparse_config"))
        self.identity_scaler = bool(conf.pop("identity_scaler", False))
        self.model_config = conf
        self.model = self._model(GP)

    # Construct one exact or sparse HEBO model with common settings
    def _model(self, model_type):
        conf = dict(self.model_config)
        if model_type is SVGP:
            conf.update(self.sparse_config)
        model = model_type(self.num_cont, self.num_enum, self.num_out, **conf)
        if self.identity_scaler:
            model.xscaler = TorchIdentityScaler()
        return model

    # Fit the exact or sparse model selected by the observation count
    def fit(self, Xc: torch.Tensor, Xe: torch.Tensor, y: torch.Tensor):
        fitted = max(0, int(torch.isfinite(y).all(dim=1).sum()) - self.synthetic_trials)
        model_type = GP if fitted <= self.max_exact_trials else SVGP
        if not isinstance(self.model, model_type):
            self.model = self._model(model_type)
        return self.model.fit(Xc, Xe, y)

    # Predict through the selected model
    def predict(self, Xc: torch.Tensor, Xe: torch.Tensor):
        return self.model.predict(Xc, Xe)

    # Draw predictive samples through the selected model
    def sample_y(self, Xc: torch.Tensor, Xe: torch.Tensor, n_samples: int = 1):
        return self.model.sample_y(Xc, Xe, n_samples=n_samples)

    # Compute the fitted observation noise variance
    @property
    def noise(self):
        return self.model.noise


class TopologyGP(HybridGP):
    """Hybrid GP adapter with a circular and real-projective kernel."""

    # Initialize the fixed-bound manifold kernel and wrapped HEBO GP
    def __init__(self, num_cont: int, num_enum: int, num_out: int, **conf):
        if num_enum != 0:
            raise ValueError("TopologyGP supports numeric HEBO spaces only")
        self.topology_groups = tuple(
            {**group, "indices": tuple(int(index) for index in group["indices"])}
            for group in conf.pop("topology_groups")
        )
        self.numeric_bounds = tuple((float(bound[0]), float(bound[1])) for bound in conf.pop("numeric_bounds"))
        self._validate_layout(num_cont)
        topology = conf.pop("topology_config", None)
        if topology is None:
            topology = load_hebo_config()["topology"]
        conf["kern"] = ScaleKernel(
            TopologyMaternKernel(
                num_cont,
                self.topology_groups,
                self.numeric_bounds,
                nu=topology["nu"],
                lengthscale_scale=topology["lengthscale_scale"],
                ard=conf.get("ard_kernel", True),
            ),
            outputscale_prior=GammaPrior(topology["outputscale_prior_shape"], topology["outputscale_prior_rate"]),
        )
        conf["identity_scaler"] = True
        super().__init__(num_cont, num_enum, num_out, **conf)

    # Validate topology indices and fixed numerical bounds
    def _validate_layout(self, num_cont: int) -> None:
        if len(self.numeric_bounds) != num_cont:
            raise ValueError("HEBO numerical bounds do not match the model dimension")
        if any(upper <= lower for lower, upper in self.numeric_bounds):
            raise ValueError("HEBO numerical bounds must have positive width")
        all_indices = []
        for group in self.topology_groups:
            if group["kind"] not in {
                "phase",
                "projective",
                "sphere",
                "simplex",
                "polar_projective",
                "cartesian_projective",
                "phase_projective",
                "polar_sphere",
                "polar_components",
                "radial_projective",
            }:
                raise ValueError(f'Unsupported HEBO topology kind "{group["kind"]}"')
            indices = group["indices"]
            if group["kind"] == "phase" and len(indices) != 1:
                raise ValueError("HEBO phase topology requires exactly one index")
            if group["kind"] in {"projective", "sphere", "simplex"} and not indices:
                raise ValueError(f"HEBO {group['kind']} topology requires at least one index")
            if group["kind"] in {"polar_projective", "cartesian_projective"} and len(indices) < 2:
                raise ValueError(f"HEBO {group['kind']} topology requires coefficient coordinates")
            if group["kind"] in {"polar_sphere", "polar_components"} and len(indices) < 2:
                raise ValueError(f"HEBO {group['kind']} topology requires phase and coupling coordinates")
            if group["kind"] == "radial_projective" and len(indices) < 2:
                raise ValueError("HEBO radial_projective topology requires a norm and direction")
            if group["kind"] == "phase_projective" and not indices:
                raise ValueError("HEBO phase_projective topology requires a phase coordinate")
            if any(index < 0 or index >= num_cont for index in indices):
                raise ValueError("HEBO topology index is outside the numerical parameter space")
            all_indices.extend(indices)
        if len(all_indices) != len(set(all_indices)):
            raise ValueError("HEBO topology indices overlap")


class SeededEvolutionOpt:
    """Small deterministic form of HEBO's Pymoo acquisition optimizer."""

    # Initialize the wrapped HEBO acquisition and fixed Pymoo seed
    def __init__(self, space, acquisition, *, seed: int, pop: int = 100, iterations: int = 100):
        self.space = space
        self.acquisition = acquisition
        self.seed = int(seed)
        self.pop = int(pop)
        self.iterations = int(iterations)

    # Optimize one acquisition without Pymoo drawing entropy from the operating system
    def optimize(self, initial, fix_input=None) -> pd.DataFrame:
        from hebo.acq_optimizers.evolution_optimizer import BOProblem, get_init_pop
        from pymoo.algorithms.moo.nsga2 import NSGA2
        from pymoo.core.mixed import MixedVariableDuplicateElimination, MixedVariableGA, MixedVariableMating
        from pymoo.optimize import minimize

        problem = BOProblem(self.acquisition, self.space, fix_input)
        population = get_init_pop(self.space, self.pop, initial, True)
        if self.acquisition.num_obj == 1:
            algorithm = MixedVariableGA(pop_size=self.pop, sampling=population)
        else:
            duplicate = MixedVariableDuplicateElimination()
            algorithm = NSGA2(
                pop_size=self.pop,
                sampling=population,
                mating=MixedVariableMating(eliminate_duplicates=duplicate),
                eliminate_duplicates=duplicate,
            )
        result = minimize(problem, algorithm, ("n_gen", self.iterations), seed=self.seed, verbose=False)
        values = result.X if result.X is not None else [point.X for point in result.pop]
        if isinstance(values, dict):
            values = [values]
        if isinstance(values, np.ndarray):
            values = values.tolist()
        numeric = pd.DataFrame(values)[self.space.para_names].to_numpy(dtype=float, copy=True)
        continuous = torch.from_numpy(numeric[:, : self.space.num_numeric])
        categorical = torch.from_numpy(numeric[:, self.space.num_numeric :])
        output = self.space.inverse_transform(continuous, categorical)
        if fix_input is not None:
            for key, value in fix_input.items():
                output[key] = value
        return output


# Transform completed costs without hiding surrogate fitting errors
def _transform_costs(optimizer) -> torch.Tensor:
    from sklearn.preprocessing import power_transform

    raw = torch.as_tensor(optimizer.y, dtype=torch.float32).clone()
    deviation = float(optimizer.y.std())
    if optimizer.settings["output"]["transform"] == "standardize" or deviation <= np.finfo(float).tiny:
        return raw
    methods = ("yeo-johnson",) if optimizer.y.min() <= 0.0 else ("box-cox", "yeo-johnson")
    for method in methods:
        try:
            values = torch.as_tensor(power_transform(optimizer.y / deviation, method=method), dtype=torch.float32)
        except (ValueError, FloatingPointError, OverflowError):
            continue
        if torch.isfinite(values).all() and values.std() >= 0.5:
            return values
    logger.warning("HEBO power transformation failed, using standardized costs")
    return raw


# Select exploitation first so a partial scheduler batch retains the best mean
def _select_recommendations(model, space, recommendations: pd.DataFrame, count: int) -> pd.DataFrame:
    with torch.no_grad():
        mean, variance = model.predict(*space.transform(recommendations))
        best = int(mean.reshape(-1).argmin())
        uncertain = int(variance.reshape(-1).argmax())
    selected = [best]
    if count > 1 and uncertain != best:
        selected.append(uncertain)
    remaining = [index for index in range(len(recommendations)) if index not in selected]
    selected.extend(np.random.choice(remaining, count - len(selected), replace=False).tolist())
    return recommendations.iloc[selected].copy()


# Compute deterministic HEBO suggestions from one fully rebuilt observation state
def suggest_seeded(optimizer, *, n_suggestions: int, seed: int) -> pd.DataFrame:
    from hebo.acquisitions.acq import MACE
    from hebo.models.model_factory import get_model

    count = int(n_suggestions)
    if count < 1:
        raise ValueError("HEBO requires at least one suggestion")
    if optimizer.acq_cls is not MACE and count != 1:
        raise RuntimeError("Parallel optimization is supported only for MACE acquisition")
    if optimizer.X.shape[0] < max(2, optimizer.rand_sample):
        return optimizer.quasi_sample(count)

    continuous, categorical = optimizer.space.transform(optimizer.X)
    values = _transform_costs(optimizer)
    model = get_model(
        optimizer.model_name, optimizer.space.num_numeric, optimizer.space.num_categorical, 1, **optimizer.model_config
    )
    model.fit(continuous, categorical, values)

    best = optimizer.X.iloc[[optimizer.get_best_id(None)]]
    predicted, _ = model.predict(*optimizer.space.transform(best))
    best_mean = predicted.detach().numpy().squeeze()
    completed = optimizer.X.shape[0] - int(optimizer.model_config.get("synthetic_trials", 0))
    iteration = max(1, completed // count)
    config = optimizer.settings["acquisition"]
    kappa = np.sqrt(
        config["kappa_scale"]
        * 2.0
        * ((2.0 + optimizer.X.shape[1] / 2.0) * np.log(iteration) + np.log(np.pi**2 / config["delta"]))
    )
    acquisition = optimizer.acq_cls(model, best_y=best_mean, kappa=kappa, eps=config["eps"])
    recommendations = SeededEvolutionOpt(
        optimizer.space, acquisition, seed=seed, pop=config["population"], iterations=config["generations"]
    ).optimize(best)
    recommendations = recommendations.drop_duplicates()
    recommendations = recommendations[optimizer.check_unique(recommendations)]
    for _ in range(4):
        if recommendations.shape[0] >= count:
            break
        random_rows = optimizer.quasi_sample(count - recommendations.shape[0])
        random_rows = random_rows[optimizer.check_unique(random_rows)]
        recommendations = pd.concat([recommendations, random_rows], ignore_index=True)
    if recommendations.shape[0] < count:
        recommendations = pd.concat(
            [recommendations, optimizer.quasi_sample(count - recommendations.shape[0])], ignore_index=True
        )

    return _select_recommendations(model, optimizer.space, recommendations, count)


class SeededHEBO(HEBO):
    """HEBO with an explicit deterministic acquisition seed."""

    # Compute one deterministic adaptive or Sobol suggestion block
    def suggest(self, n_suggestions=1, fix_input=None):
        if fix_input is not None:
            raise ValueError("Seeded HEBO fixed inputs are not used by icetune")
        seed = int(getattr(self, "acquisition_seed", self.scramble_seed or 0))
        python_state = random.getstate()
        numpy_state = np.random.get_state()
        torch_state = torch.random.get_rng_state()
        try:
            random.seed(seed)
            np.random.seed(seed % (2**32))
            torch.manual_seed(seed)
            return suggest_seeded(self, n_suggestions=int(n_suggestions), seed=seed)
        finally:
            random.setstate(python_state)
            np.random.set_state(numpy_state)
            torch.random.set_rng_state(torch_state)


# Build the HEBO topology model configuration from its numeric design space
def _topology_model_config(design, parameter_topology: dict | None, settings: dict) -> dict | None:
    groups = surrogate_parameter_groups(parameter_topology)
    if not groups or not settings["topology"]["enabled"]:
        return None
    numeric_names = list(design.numeric_names)
    group_names = {name for group in groups for name in group["parameters"]}
    missing = sorted(group_names - set(numeric_names))
    if missing:
        raise ValueError(f"HEBO topology parameters are absent from the numerical design space: {missing}")
    topology_groups = [
        {
            "kind": group["kind"],
            "base": group["base"],
            "indices": [numeric_names.index(name) for name in group["parameters"]],
            **(
                {"period": group["period"]}
                if group["kind"]
                in {"phase", "polar_projective", "phase_projective", "polar_sphere", "polar_components"}
                else {}
            ),
            **({"weights": group["weights"]} if "weights" in group else {}),
        }
        for group in groups
    ]
    numeric_bounds = [[float(design.paras[name].lb), float(design.paras[name].ub)] for name in numeric_names]
    return {
        **_hybrid_model_config(settings),
        "topology_groups": topology_groups,
        "numeric_bounds": numeric_bounds,
        "topology_config": dict(settings["topology"]),
    }


# Register the local exact and sparse HEBO models
def _register_models() -> None:
    from hebo.models import model_factory

    for name, model in ((HYBRID_MODEL_NAME, HybridGP), (TOPOLOGY_MODEL_NAME, TopologyGP)):
        registered = model_factory.model_dict.get(name)
        if registered not in {None, model}:
            raise RuntimeError(f'HEBO model name "{name}" is already registered')
        model_factory.model_dict[name] = model


# Construct standard HEBO or the phase-topology variant
def create_hebo(
    design, *, parameter_topology: dict | None, rand_sample: int, scramble_seed: int, settings: dict | None = None
):
    _register_models()
    settings = load_hebo_config() if settings is None else validate_hebo_config(settings)
    model_config = _topology_model_config(design, parameter_topology, settings)
    optimizer = SeededHEBO(
        design,
        model_name=TOPOLOGY_MODEL_NAME if model_config is not None else HYBRID_MODEL_NAME,
        model_config=model_config if model_config is not None else _hybrid_model_config(settings),
        rand_sample=rand_sample,
        scramble_seed=scramble_seed,
    )
    optimizer.settings = settings
    return optimizer
