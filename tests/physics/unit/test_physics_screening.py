# Unit tests for unified physics-test CLI steering
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest

from tests.technical.support.screening import PhysicsScreening


# Add screening and weighting to arguments accepted by the iceplot helper
@pytest.mark.parametrize(("enabled", "expected"), [(False, "0"), (True, "1")])
def test_generator_args_add_selection(enabled, expected):
    selected = PhysicsScreening(enabled, weighted=False).generator_args(["-i", "gencard/test.json"])
    assert selected == ["-i", "gencard/test.json", "-l", expected, "-w", "0"]


# Replace existing generator selection and event-count values
def test_generator_args_replace_selection():
    selected = PhysicsScreening(True, weighted=False, event_override=123).generator_args(
        ["-i", "gencard/test.json", "-l", "false", "-w", "true", "-n", "5"]
    )
    assert selected == ["-i", "gencard/test.json", "-l", "1", "-w", "0", "-n", "123"]


# Apply shared CLI selection only to a full argv that invokes gr
def test_argv_selects_only_gr_commands():
    screening = PhysicsScreening(True, weighted=True, event_override=17)
    assert screening.argv(["./bin/gr", "-i", "gencard/test.json"]) == [
        "./bin/gr",
        "-i",
        "gencard/test.json",
        "-l",
        "1",
        "-w",
        "1",
        "-n",
        "17",
    ]
    assert screening.argv(["python", "python/src/core/iceplot.py"]) == ["python", "python/src/core/iceplot.py"]


# Apply shared CLI selection only to shell commands that invoke gr
def test_shell_selects_only_gr_commands():
    screening = PhysicsScreening(False, weighted=False, event_override=19)
    assert screening.shell("./bin/gr -i gencard/test.json") == (
        "./bin/gr -i gencard/test.json -l 0 -w 0 -n 19"
    )
    assert screening.shell("./bin/gr -l true -w true -n 4 -i gencard/test.json") == (
        "./bin/gr -l 0 -w 0 -n 19 -i gencard/test.json"
    )
    assert screening.shell("python -m core.iceplot") == "python -m core.iceplot"


# Reject malformed or ambiguous screening options
@pytest.mark.parametrize(
    "args",
    [
        ["-l"],
        ["-l", "invalid"],
        ["-l", "true", "-l", "false"],
    ],
)
def test_gen_args_reject_invalid_screening(args):
    with pytest.raises(ValueError):
        PhysicsScreening(True).generator_args(args)


# Reject malformed or ambiguous shell screening options
@pytest.mark.parametrize(
    "command",
    [
        "./bin/gr -l",
        "./bin/gr -l invalid",
        "./bin/gr -l true -l false",
    ],
)
def test_shell_reject_invalid_screening(command):
    with pytest.raises(ValueError):
        PhysicsScreening(True).shell(command)


# Reject malformed weighting and event-count options in shell commands
@pytest.mark.parametrize(
    "command",
    [
        "./bin/gr -w",
        "./bin/gr -w invalid",
        "./bin/gr -w 1 -w 0",
        "./bin/gr -n",
        "./bin/gr -n invalid",
        "./bin/gr -n 1 -n 2",
    ],
)
def test_shell_reject_invalid_shared_options(command):
    with pytest.raises(ValueError):
        PhysicsScreening(True, event_override=10).shell(command)


# Reject malformed weighting and event-count options when overriding them
@pytest.mark.parametrize(
    "args",
    [
        ["-w"],
        ["-w", "invalid"],
        ["-w", "1", "-w", "0"],
        ["-n"],
        ["-n", "invalid"],
        ["-n", "1", "-n", "2"],
    ],
)
def test_gen_args_reject_invalid_shared_options(args):
    with pytest.raises(ValueError):
        PhysicsScreening(True, event_override=10).generator_args(args)


# Preserve local event counts when no CLI override was requested
def test_gen_args_local_event_count():
    selected = PhysicsScreening(False, weighted=True).generator_args(
        ["-n", "0", "-i", "gencard/test.json"]
    )
    assert selected == ["-n", "0", "-i", "gencard/test.json", "-l", "0", "-w", "1"]
