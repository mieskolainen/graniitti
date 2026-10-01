# ZEUS tagged proton dissociative phi photoproduction selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from ..._common import dissociative
from ..cuts import cut_param, tagged

# [REFERENCE: arXiv:hep-ex/9910038, Table 6]
cut_param = {**cut_param, "MY2_W2_MAX": 0.1}


# Apply the published target mass fraction and positron tagger selection
def cut_func(event):
    return dissociative.accepted(event) and tagged(event)
