# Simulator driver interface for fitting and validation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import pathlib
from abc import ABC, abstractmethod
from collections.abc import Callable


class SimulatorDriver(ABC):
    """Abstract base class for simulation or black-box optimization drivers."""

    DEFAULT_TUNE = None
    AMPLITUDE_FIT = False

    # Compute the stable command-line and summary identifier
    @classmethod
    @abstractmethod
    def driver_name(cls) -> str: ...


    # Compute whether a metadata-free summary belongs to this driver
    @classmethod
    @abstractmethod
    def matches_summary(cls, summary: dict) -> bool: ...

    # Add driver-owned command-line settings
    @classmethod
    def add_cli_arguments(cls, parser) -> None:
        return None

    # Validate driver-owned command-line settings
    @classmethod
    def validate_cli_arguments(cls, parser, args) -> None:
        return None

    # Compute the default external library directory
    @classmethod
    def default_library_path(cls) -> str | None:
        return None

    # Build the driver-owned run steering payload
    def build_run_steering(self, args) -> dict:
        return {"tune_default": args.tune_default, "tunesetup_name": args.tunesetup}

    # Build or apply the driver runtime environment
    @abstractmethod
    def runtime_environment(
        self, *, cdir: str, libdir: str, python_version: str, apply: bool = False
    ) -> dict[str, str]: ...

    # Build the selected tuning card in the driver context
    @classmethod
    def build_tunesetup(cls, *, config, cdir, tune_default):
        raise NotImplementedError("The simulator must build its JSON tuning definition")

    # Freeze driver inputs before saving a tuning definition or launching jobs
    def prepare_tunesetup(self, *, tunesetup, args) -> None:
        return None

    # Prepare driver-owned immutable state for repeated trial evaluations
    def prepare_trial_runtime(self, param: dict) -> None:
        return None

    # Describe the physical parameter map for a simulator with direct input coordinates
    def parameter_transform(self, names, reference, *, cdir: str, metadata: dict):
        from core.stats.transform import ParameterTransform

        return ParameterTransform(names, dict, reference)

    # Identify shared inputs which must retain their absolute source paths
    def shared_runtime_files(self, param: dict) -> list[str]:
        return []

    # List driver inputs that must accompany an uploaded worker runtime
    def runtime_files(self, param: dict) -> list[str]:
        return list(param.get("aux_param_space", {}).get("runtime_files", {}))

    # Compute a driver-owned identifier for trial-local outputs
    def trial_tunename(
        self, *, trial_id: str, node_id: str | None = None, pid: int | None = None, purpose: str = "trial"
    ) -> str:
        logical_id = trial_id if purpose == "trial" else f"{trial_id}-{purpose}"
        from core.tune.core import unique_trial_identifier

        return unique_trial_identifier(prefix=self.TRIAL_PREFIX, trial_id=logical_id, node_id=node_id, pid=pid)

    # Compute bootstrap datacards for a selected backend
    def prepare_init_datacards(self, datacards: list[dict]) -> list[dict]:
        return copy.deepcopy(datacards)

    # Initialize driver data
    def init_data(self, *args, **kwargs):
        raise NotImplementedError(f"{self.__class__.__name__} does not implement init_data()")

    # Build driver-owned callbacks for one selected backend process
    @abstractmethod
    def initialize_backend(self, *args, **kwargs) -> dict: ...

    @abstractmethod
    def initialize(self, *args, **kwargs):
        """Driver initialization."""

    # Compute a content fingerprint for Ray restart validation
    @abstractmethod
    def physics_fingerprint(self, *, param: dict, tunesetup) -> str: ...

    @abstractmethod
    def get_initial_param(
        self, param_space: dict, aux_param_space: dict, cdir: str, tune_default: str | None = None
    ) -> dict:
        """Return default parameters for the search space."""

    @abstractmethod
    def create_steering_card(
        self, param_space: dict, tunename: str, cdir: str | None = None, tune_default: str | None = None
    ) -> dict:
        """Stage driver-specific steering outputs for one evaluation."""

    # Resolve the common optimized config and publication target
    def _resolve_push_inputs(self, summary: dict, target_path: str, cdir: str) -> tuple[dict, pathlib.Path]:
        config = summary.get("config")
        if not isinstance(config, dict) or not config:
            raise ValueError(f"{self.driver_name()} push summary has no parameters")
        target = pathlib.Path(target_path).expanduser()
        if not target.is_absolute():
            target = pathlib.Path(cdir) / target
        return config, target.resolve()

    # Publish fitted parameters only after the caller approves the prepared rows
    def push_parameters(
        self,
        *,
        summary: dict,
        target_path: str,
        cdir: str,
        options: dict,
        confirm: Callable[[list[tuple[str, object, object]]], bool],
    ) -> list[tuple[str, object, object]]:
        raise NotImplementedError(f"{self.__class__.__name__} does not implement push_parameters()")

    @abstractmethod
    def compute(self, *args, **kwargs) -> dict:
        """Run the driver for one evaluation."""

    @abstractmethod
    def evaluate_trial_outputs(self, *, config: dict, param: dict, trial_id: str, tunename: str) -> dict:
        """Evaluate one trial and return a normalized output bundle."""

    @abstractmethod
    def render_trial_figures_to_dir(
        self, *, outputs: dict, param: dict, summary_payload: dict, output_dir: str, summary_file: str | None = None
    ) -> dict:
        """Render driver-specific visualizations into a caller-provided directory."""

    def cleanup_trial_outputs(self, *args, **kwargs) -> None:  # noqa: B027
        """Clean up driver-specific temporary outputs after one trial."""
