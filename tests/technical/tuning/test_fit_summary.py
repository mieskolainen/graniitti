# Tests for shared best fit summary contracts
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.tune import summary as fit_summary


# Check the shared best-fit writer preserves values and optional uncertainties
def test_build_best_fit_contract():
    payload = fit_summary.build_best_fit(
        values={"x": 0.25, "n": 2},
        uncertainties={"x": 0.05},
        objective_name="Z",
        objective_value=12.5,
        source="surrogate",
        uncertainty_method="torch_autograd_hessian",
    )

    assert payload == {
        "schema_version": 1,
        "source": "surrogate",
        "objective": {"name": "Z", "value": 12.5},
        "uncertainty_method": "torch_autograd_hessian",
        "parameters": {
            "x": {"value": 0.25, "uncertainty": 0.05},
            "n": {"value": 2.0, "uncertainty": None},
        },
    }
    assert fit_summary.best_fit_values(payload) == {"x": 0.25, "n": 2.0}


# Check unavailable numerical diagnostics serialize as null
def test_build_best_fit_null_unavailable_diagnostics():
    payload = fit_summary.build_best_fit(
        values={"x": 0.25},
        uncertainties={"x": math.nan},
        objective_name="chi2",
        objective_value=math.inf,
        source="realized",
    )

    assert payload["parameters"]["x"]["uncertainty"] is None
    assert payload["objective"]["value"] is None
    assert payload["uncertainty_method"] is None


# Check corrupt common blocks cannot reach a steering-card writer
@pytest.mark.parametrize(
    "payload",
    [
        {},
        {"schema_version": 2, "parameters": {"x": {"value": 0.5}}},
        {"schema_version": 1, "parameters": {}},
        {"schema_version": 1, "parameters": {"x": {"uncertainty": 0.1}}},
        {"schema_version": 1, "parameters": {"x": {"value": math.nan}}},
    ],
)
def test_best_fit_values_invalid_contract(payload):
    with pytest.raises(ValueError):
        fit_summary.best_fit_values(payload)
