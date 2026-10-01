# Shared CLI steering for generator driven physics tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""Shared CLI steering for generator-driven physics tests."""

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

GRANIITTI_SHELL_PATTERN = re.compile(r"(?:^|\s)(?:\S*/)?bin/gr(?:\s|$)")
BOOLEAN_VALUES = frozenset(("true", "false", "0", "1"))
SHELL_PATTERNS = {
    "-l": re.compile(r"(?<!\S)-l(?:\s+(\S+))?"),
    "-w": re.compile(r"(?<!\S)-w(?:\s+(\S+))?"),
    "-n": re.compile(r"(?<!\S)-n(?:\s+(\S+))?"),
}


# Add or replace one validated option in generator arguments
def replace_generator_option(args, option, value, valid_value):
    command = list(args)
    indices = [index for index, item in enumerate(command) if item == option]
    if len(indices) > 1:
        raise ValueError(f"Generator command contains multiple {option} options")
    if not indices:
        command.extend((option, value))
        return command

    index = indices[0]
    if index + 1 >= len(command):
        raise ValueError(f"Generator command has {option} without a value")
    if not valid_value(command[index + 1]):
        raise ValueError(f"Generator command has invalid {option} value {command[index + 1]!r}")
    command[index + 1] = value
    return command


# Add or replace one validated option in a generator shell command
def replace_shell_option(command, option, value, valid_value):
    pattern = SHELL_PATTERNS[option]
    matches = list(pattern.finditer(command))
    if len(matches) > 1:
        raise ValueError(f"Generator command contains multiple {option} options")
    if not matches:
        return f"{command} {option} {value}"

    current = matches[0].group(1)
    if current is None or not valid_value(current):
        raise ValueError(f"Generator command has invalid {option} value {current!r}")
    return pattern.sub(f"{option} {value}", command, count=1)


# Accept the generator's numeric and textual boolean spellings
def is_boolean_value(value):
    return value in BOOLEAN_VALUES


# Accept non-negative integer generator event counts
def is_event_count(value):
    try:
        return int(value) >= 0 and str(value).strip() == str(int(value))
    except (TypeError, ValueError):
        return False


@dataclass(frozen=True)
class PhysicsScreening:
    """Apply one immutable CLI selection to physics-test commands"""

    enabled: bool
    weighted: bool = True
    event_override: int | None = None

    # Validate the optional shared event-count override
    def __post_init__(self):
        if self.event_override is None:
            return
        if (
            not isinstance(self.event_override, int)
            or isinstance(self.event_override, bool)
            or self.event_override <= 0
        ):
            raise ValueError(
                "Physics event-count override must be a positive integer, "
                f"got {self.event_override!r}"
            )

    # Compute the generator-compatible numeric screening value
    @property
    def value(self) -> str:
        return str(int(self.enabled))

    # Compute the generator-compatible numeric weighting value
    @property
    def weighted_value(self) -> str:
        return str(int(self.weighted))

    # Compute a collision-resistant label for the shared physics selection
    @property
    def suffix(self) -> str:
        return f"loopscreen_{self.value}_weighted_{self.weighted_value}"

    # Add or replace shared CLI options in generator arguments
    def generator_args(self, args):
        command = replace_generator_option(args, "-l", self.value, is_boolean_value)
        command = replace_generator_option(command, "-w", self.weighted_value, is_boolean_value)
        if self.event_override is not None:
            command = replace_generator_option(
                command, "-n", str(self.event_override), is_event_count
            )
        return command

    # Add or replace shared CLI options only when argv invokes gr
    def argv(self, args):
        command = list(args)
        if not command or Path(command[0]).name != "gr":
            return command
        return [command[0], *self.generator_args(command[1:])]

    # Add or replace shared CLI options in one generator shell command
    def shell(self, command):
        if not GRANIITTI_SHELL_PATTERN.search(command):
            return command
        selected = replace_shell_option(command, "-l", self.value, is_boolean_value)
        selected = replace_shell_option(selected, "-w", self.weighted_value, is_boolean_value)
        if self.event_override is not None:
            selected = replace_shell_option(
                selected, "-n", str(self.event_override), is_event_count
            )
        return selected

    # Execute generator argv with the shared CLI selection
    def run(self, args, **kwargs):
        return subprocess.run(self.argv(args), **kwargs)

    # Run one generator with the shared CLI selection through iceplot support
    def run_graniitti(self, args):
        from tests.technical.support.iceplot import run_graniitti

        return run_graniitti(self.generator_args(args))

    # Execute shell commands with the shared CLI selection through test support
    def execute(self, commands, expect="[gr: done]"):
        from tests.technical.support.commands import execute

        selected = (
            [self.shell(command) for command in commands]
            if isinstance(commands, list)
            else self.shell(commands)
        )
        return execute(selected, expect=expect)
