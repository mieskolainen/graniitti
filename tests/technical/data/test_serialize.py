# Numerical JSON serialization for persisted fit outputs
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json

import numpy as np
import pytest
from core.io.serialize import json_default, json_safe


# Convert extended precision scalars and arrays without recursive scalar conversion
def test_extended_precision_json():
    values = {
        "scalar": np.longdouble("1.25"),
        "array": np.array([1.25, np.inf, np.nan], dtype=np.longdouble),
        "complex": np.clongdouble(1 + 2j),
    }
    converted = json_safe(values)
    assert converted == {"scalar": 1.25, "array": [1.25, None, None], "complex": "(1+2j)"}
    assert json.loads(json.dumps(converted, allow_nan=False)) == converted
    assert json.loads(json.dumps(values["scalar"], default=json_default, allow_nan=False)) == 1.25


# Preserve the selected nonfinite representation for extended precision values
@pytest.mark.parametrize("value", [np.inf, -np.inf, np.nan])
def test_extended_precision_nonfinite_json(value):
    assert json_safe(np.longdouble(value), nonfinite="string") == str(value)
