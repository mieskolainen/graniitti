# ZEUS tagged elastic phi photoproduction selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common import exclusive_jpsi

# [REFERENCE: arXiv:hep-ex/9910038, Sections 5.2, 6.2 and 10.2]
cut_param = {"W": [85.0, 105.0], "Q2_MAX": 0.01, "THETA_E_MAX": 0.003, "M": [0.99, 1.06], "ABS_T_MAX": 3.0}


# Apply the positron tagger scattering angle limit relative to the incident beam
# [REFERENCE: arXiv:hep-ex/9910038, Section 6.2]
def tagged(event):
    incoming = exclusive_jpsi.unique_momentum(event, 11, True)
    outgoing = exclusive_jpsi.unique_momentum(event, 11, False)
    return incoming.angle(outgoing) < event.cut_param["THETA_E_MAX"]


# Apply the published elastic phase space and positron tagger selection
def cut_func(event):
    return exclusive_jpsi.accepted(event) and tagged(event)
