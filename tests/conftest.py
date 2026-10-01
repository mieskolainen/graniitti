# Pytest configuration file
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import pathlib
import sys

import pyjson5 as json5
import pytest

# Keep generated Condor workspaces out of test collection
collect_ignore = ["condor/runs"]

# Make tests.* helper imports stable under pytest importlib collection
PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tests.physics.validation._registry import LITERATURE_VALIDATIONS, VALIDATION_PREFIX  # noqa: E402


# Parse boolean-like pytest options to a stable integer 0/1 convention
def parse_bool_int(value):
    if isinstance(value, bool):
        return int(value)
    text = str(value).strip().lower()
    if text in ("1", "true", "yes", "on"):
        return 1
    if text in ("0", "false", "no", "off"):
        return 0
    raise pytest.UsageError(f"Expected boolean integer 0/1 or true/false, got {value!r}")


# Parse a positive event-count override
def parse_positive_int(value):
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise pytest.UsageError(f"Expected a positive integer, got {value!r}") from exc
    if result <= 0:
        raise pytest.UsageError(f"Expected a positive integer, got {value!r}")
    return result


# Register project-specific pytest command-line options
def pytest_addoption(parser):
    parser.addoption(
        "--test-output-dir",
        type=str,
        default=None,
        help=(
            "root directory for pytest-driven test results, overriding GRANIITTI_TEST_OUTPUT_DIR"
        ),
    )
    parser.addoption(
        "--LOOPSCREEN",
        type=parse_bool_int,
        default=None,
        help="override Pomeron-loop screening in generator-driven physics tests",
    )
    parser.addoption(
        "--WEIGHTED",
        type=parse_bool_int,
        default=1,
        help="select weighted or unweighted generation in generator-driven physics tests",
    )
    parser.addoption(
        "--NEVENTS",
        type=parse_positive_int,
        default=None,
        help="override generated events per process in generator-driven physics tests",
    )
    parser.addoption("--GENERATE", type=parse_bool_int, default=1)
    parser.addoption("--ANALYZE", type=parse_bool_int, default=1)
    parser.addoption("--MODELPARAM", type=str, default="TUNE0")
    parser.addoption(
        "--run-integration",
        action="store_true",
        default=False,
        help="run technical generator tests below tests/technical/integration and tests/technical/basic",
    )
    parser.addoption(
        "--run-physics",
        action="store_true",
        default=False,
        help="run generator-driven tests below tests/physics",
    )


# Apply the command-line test-output root before test modules are collected
def pytest_configure(config):
    output_dir = config.getoption("--test-output-dir")
    if output_dir:
        os.environ["GRANIITTI_TEST_OUTPUT_DIR"] = output_dir


# Compute the configured root directory for pytest-driven test results
@pytest.fixture(scope="session")
def test_output_dir():
    from tests.technical.support.output import get_test_output_root

    path = get_test_output_root()
    path.mkdir(parents=True, exist_ok=True)
    return path


# Mark generator-driven integration and physics tests and keep them opt-in
def pytest_collection_modifyitems(config, items):
    run_integration = config.getoption("--run-integration")
    run_physics = config.getoption("--run-physics")
    skip_integration = pytest.mark.skip(reason="integration tests require --run-integration")
    skip_physics = pytest.mark.skip(reason="physics tests require --run-physics")
    integration_root = PROJECT_ROOT / "tests" / "technical" / "integration"
    basic_root = PROJECT_ROOT / "tests" / "technical" / "basic"
    physics_root = PROJECT_ROOT / "tests" / "physics"
    physics_unit_root = physics_root / "unit"
    for item in items:
        relative_path = item.path.resolve().relative_to(PROJECT_ROOT).as_posix()
        try:
            is_integration = item.path.is_relative_to(integration_root) or item.path.is_relative_to(
                basic_root
            )
            is_physics = item.path.is_relative_to(physics_root) and not item.path.is_relative_to(
                physics_unit_root
            )
        except AttributeError:
            is_integration = str(item.fspath).startswith(str(integration_root)) or str(
                item.fspath
            ).startswith(str(basic_root))
            is_physics = str(item.fspath).startswith(str(physics_root)) and not str(
                item.fspath
            ).startswith(str(physics_unit_root))
        if is_integration:
            item.add_marker(pytest.mark.integration)
            if not run_integration:
                item.add_marker(skip_integration)
        if is_physics:
            item.add_marker(pytest.mark.physics)
            validation_path = relative_path
            literature = False
            closure = relative_path == f"{VALIDATION_PREFIX}/test_symbolic.py"
            if relative_path == f"{VALIDATION_PREFIX}/test_icepacks.py":
                callspec = getattr(item, "callspec", None)
                dataset_path = None if callspec is None else callspec.params.get("dataset_path")
                if dataset_path is not None:
                    validation_path = (
                        pathlib.Path(dataset_path)
                        .resolve()
                        .relative_to(PROJECT_ROOT)
                        .as_posix()
                    )
                    dataset = json5.loads(pathlib.Path(dataset_path).read_text(encoding="utf-8"))
                    literature = any(entry.get("data", True) for entry in dataset["sets"])
                    closure = dataset.get("validation", {}).get("mc_reference") is not None
                elif item.originalname == "test_integrated_xs_table":
                    literature = True
            if literature or validation_path in LITERATURE_VALIDATIONS:
                item.add_marker(pytest.mark.literature)
            if closure:
                item.add_marker(pytest.mark.closure)
            if not run_physics:
                item.add_marker(skip_physics)


# Compute the unified generator steering for generator-driven physics tests
@pytest.fixture(scope="session")
def physics_screening(request):
    from tests.technical.support.screening import PhysicsScreening

    return PhysicsScreening(
        enabled=bool(request.config.getoption("--LOOPSCREEN")),
        weighted=bool(request.config.getoption("--WEIGHTED")),
        event_override=request.config.getoption("--NEVENTS"),
    )


# Compute the unified event-count selector for generator-driven physics tests
@pytest.fixture(scope="session")
def physics_events(request):
    from tests.technical.support.event_count import PhysicsEventCount

    return PhysicsEventCount(override=request.config.getoption("--NEVENTS"))


# Compute the event count used by older process tests
@pytest.fixture(scope="session")
def NEVENTS(physics_events):
    return physics_events.select(50000)


# Compute the weighted-generation mode used by older process tests
@pytest.fixture(scope="session")
def WEIGHTED(request):
    return request.config.getoption("--WEIGHTED")


# Compute the generation switch used by older process tests
@pytest.fixture(scope="session")
def GENERATE(request):
    value = request.config.option.GENERATE
    if value is None:
        pytest.skip()
    return value


# Compute the analysis switch used by older process tests
@pytest.fixture(scope="session")
def ANALYZE(request):
    value = request.config.option.ANALYZE
    if value is None:
        pytest.skip()
    return value


# Compute the model tune used by older process tests
@pytest.fixture(scope="session")
def MODELPARAM(request):
    value = request.config.option.MODELPARAM
    if value is None:
        pytest.skip()
    return value


# Construct tuning domains from the repository study settings and source model
@pytest.fixture
def tuning_context():
    from core.io.serialize import load_json_file
    from core.tune.drivers.graniitti.tunesetup.domains import context

    from submit import CAMPAIGN_DIR
    settings = load_json_file(CAMPAIGN_DIR / "tunecards/graniitti/_defaults.json")
    with context(cdir=PROJECT_ROOT, model_path=PROJECT_ROOT / "modeldata/TUNE0", settings=settings):
        yield
