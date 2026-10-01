# Disabled measurement inputs awaiting supported HEPData JSON tables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


# Require a supported HEPData JSON measurement before enabling this input
def energy():
    raise ValueError('Measurement unavailable: a supported HEPData JSON input is required')


# Require a supported HEPData JSON measurement before enabling this input
def slope():
    raise ValueError('Measurement unavailable: a supported HEPData JSON input is required')


# Require a supported HEPData JSON measurement before enabling this input
def slope_tmax():
    raise ValueError('Measurement unavailable: a supported HEPData JSON input is required')
