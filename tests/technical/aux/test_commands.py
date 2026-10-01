# Tests for command helpers and executable test layout
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest

from tests.technical.support.commands import execute


# Accept a successful command only when its completion marker is present
def test_execute_requires_completion_marker():
    assert execute("printf '[gr: done]\\n'") == "[gr: done]\n"
    with pytest.raises(AssertionError, match="missing"):
        execute("printf 'silent success\\n'")


# Propagate a nonzero shell exit even if it prints the expected marker
def test_execute_rejects_nonzero_exit():
    with pytest.raises(AssertionError, match="exit code 7"):
        execute("printf '[gr: done]\\n'; exit 7")
