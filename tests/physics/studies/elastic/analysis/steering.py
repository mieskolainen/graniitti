# Read exact eikonal outputs from explicitly steered elastic scans
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re
from pathlib import Path

from core.tune.drivers.graniitti.eikonal import inspect_matrix_output

MODELS = {"single": 1, "double": 2, "triple": 3}


# Read only cache paths reported by completed scans of the requested model
def scan_outputs(scan_dir: Path, model: str, beam: str) -> list:
    files = {}
    beams = ("pp", "ppbar") if beam == "auto" else (beam,)
    for name in beams:
        path = scan_dir / name / "scan.log"
        text = path.read_text()
        if "[xscan: done]" not in text:
            raise ValueError(f"Elastic scan is incomplete: {path}")
        paths = re.findall(r"^(?:Loaded|Saved) matrix eikonal cache: (.+)$", text, re.MULTILINE)
        if not paths:
            raise ValueError(f"Elastic scan reported no matrix outputs: {path}")
        for value in paths:
            source = inspect_matrix_output(Path(value.strip()))
            expected_beam = (2212, 2212 if name == "pp" else -2212)
            if source.nchannels != MODELS[model] or (source.beam1, source.beam2) != expected_beam:
                raise ValueError(f"Elastic scan output disagrees with {model} {name}: {source.path}")
            files[source.path] = source
    return list(files.values())
