# Shared MC parameter fitting and optimizer package
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


# Compute the compact hash used in icetune paths and status messages
def short_id(digest: str) -> str:
    return digest[:16]
