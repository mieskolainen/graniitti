# Shared event count selection for generator driven physics tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
"""Shared event-count selection for generator-driven physics tests."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PhysicsEventCount:
    """Apply one optional command-line event-count override"""

    override: int | None

    # Validate the optional command-line override once
    def __post_init__(self):
        if self.override is None:
            return
        if (
            not isinstance(self.override, int)
            or isinstance(self.override, bool)
            or self.override <= 0
        ):
            raise ValueError(
                f"Physics event-count override must be a positive integer, got {self.override!r}"
            )

    # Select the override or retain the physics test's positive default
    def select(self, default):
        if not isinstance(default, int) or isinstance(default, bool) or default <= 0:
            raise ValueError(
                f"Physics event-count default must be a positive integer, got {default!r}"
            )
        return default if self.override is None else self.override
