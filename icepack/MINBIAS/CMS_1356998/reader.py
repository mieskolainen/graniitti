# Disabled CMS diffraction measurement awaiting HEPData JSON
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


# Require a supported HEPData JSON measurement before enabling this input
def read(filename, rebin_factor=None, region=None, dataset=None):
    raise ValueError("Measurement unavailable: a supported HEPData JSON table is required")
