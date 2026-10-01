# GRANIITTI tuning-card dataset construction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


# Build one process-independent icetune datacard
def card(
    path: str,
    *,
    nevents: int,
    loopscreen: bool,
    xsmode: str,
    weighted: bool = True,
    force_density: bool = False,
    swap_process: str | None = None,
    integrator: str | None = None,
    **steering: object,
) -> dict[str, object]:
    if not isinstance(force_density, bool):
        raise TypeError("ICETUNE force_density must be boolean")
    output = {
        "datacard": path,
        "nevents": nevents,
        "weighted": weighted,
        "loopscreen": loopscreen,
        "xsmode": xsmode,
        **steering,
    }
    if force_density:
        output["force_density"] = True
    if swap_process is not None:
        output["swap_process"] = swap_process
    if integrator is not None:
        value = str(integrator).upper()
        if value not in {"VEGAS", "NEUROJAC"}:
            raise ValueError("ICETUNE integrator must be VEGAS or NEUROJAC")
        output["integrator"] = value
    return output
