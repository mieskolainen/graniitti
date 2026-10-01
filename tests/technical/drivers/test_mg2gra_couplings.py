# Unit tests for MadGraph coupling conversion invariants
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest

from develop.MG2GRA.modules.mg5_couplings import alpha_zero_definition


# Compute minimal generated independent couplings for converter tests
def parameter_source() -> str:
    return """
void Parameters_model::setIndependentCouplings() {
  GC_3 = -(mdl_q * mdl_complexi);
  GC_5 = mdl_q__exp__2 * mdl_complexi;
  GC_QCD = mdl_G * mdl_complexi;
}
"""


# Verify alpha(0) generation adapts only charge-dependent model couplings
def test_alpha_qed_zero_replacement_charge() -> None:
    definition = alpha_zero_definition(
        parameter_source(),
        "Parameters_model",
        "mdl_q",
        "mdl_q__exp__2",
        137.03599908,
    )

    assert definition is not None
    assert "mdl_q__exp__2_NEW = mdl_q_NEW * mdl_q_NEW;" in definition
    assert "GC_3 = -(mdl_q_NEW * mdl_complexi);" in definition
    assert "GC_5 = mdl_q__exp__2_NEW * mdl_complexi;" in definition
    assert "GC_QCD" not in definition


# Verify a model without configured charge couplings fails loudly
def test_alpha_qed_zero_requires_charge_dependency() -> None:
    source = parameter_source().replace("mdl_q", "mdl_unused")
    with pytest.raises(RuntimeError, match="No generated coupling depends"):
        alpha_zero_definition(
            source,
            "Parameters_model",
            "mdl_q",
            "mdl_q__exp__2",
            137.03599908,
        )


# Verify models without an alpha(0) request require no generated adapter
def test_alpha_qed_zero_optional() -> None:
    assert (
        alpha_zero_definition(
            parameter_source(), "Parameters_model", None, None, 137.03599908
        )
        is None
    )


# Verify comments, literals and nested initializers cannot forge coupling use
def test_alpha_qed_zero_cpp_lexing() -> None:
    source = r'''
void Parameters_model::setIndependentCouplings() {
  // GC_COMMENT = mdl_q;
  GC_TEXT = label("mdl_q; not a coupling");
  GC_COMPLEX = std::complex<double>{
      mdl_q,
      0.0};
  GC_LAMBDA = [=]() { return mdl_q__exp__2; }();
  /* GC_BLOCK = mdl_q__exp__2; */
}
'''
    definition = alpha_zero_definition(
        source,
        "Parameters_model",
        "mdl_q",
        "mdl_q__exp__2",
        137.03599908,
    )

    assert definition is not None
    assert "GC_COMMENT" not in definition
    assert "GC_BLOCK" not in definition
    assert "GC_TEXT" not in definition
    assert 'label("mdl_q; not a coupling")' not in definition
    assert "GC_COMPLEX = std::complex<double>{ mdl_q_NEW, 0.0};" in definition
    assert "GC_LAMBDA = [=]() { return mdl_q__exp__2_NEW; }();" in definition


# Follow charge aliases through higher powers and mixed QCD couplings
def test_alpha_qed_zero_follows_param_dependencies() -> None:
    source = parameter_source() + """
void Parameters_model::setIndependentParameters(SLHAReader &slha) {
  mdl_q = sqrt(4.0 * M_PI / inverse);
  q4 = pow(mdl_q__exp__2, 2);
  vertex = q4 / vev;
  unused = mdl_q * mass;
}
void Parameters_model::setDependentParameters() {
  mixed = vertex * mdl_G;
}
void Parameters_model::setDependentCouplings() {
  GC_MIX = mixed * mdl_complexi;
}
"""
    definition = alpha_zero_definition(source, "Parameters_model", "mdl_q", "mdl_q__exp__2", 137.03599908)
    assert "q4_NEW = pow(mdl_q__exp__2_NEW, 2)" in definition
    assert "vertex_NEW = q4_NEW / vev" in definition
    assert "mixed_NEW = vertex_NEW * mdl_G" in definition
    assert "GC_MIX = mixed_NEW * mdl_complexi" in definition
    assert "unused" not in definition
