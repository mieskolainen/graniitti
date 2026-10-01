# Amplitude bank storage and initialization failures through the real filesystem
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import errno
import json
import resource
import shutil
import signal
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import torch
from core.tune.drivers.graniitti.ampfit.amplitude import Basis, combine_basis, stage_bank
from core.tune.init import run_initialization


# Restrict real file writes while restoring the process limits and signal handler
@contextmanager
def file_limit(size):
    limits = resource.getrlimit(resource.RLIMIT_FSIZE)
    handler = signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
    try:
        resource.setrlimit(resource.RLIMIT_FSIZE, (size, limits[1]))
        yield
    finally:
        resource.setrlimit(resource.RLIMIT_FSIZE, limits)
        signal.signal(signal.SIGXFSZ, handler)


# Prepare interleaved complex components with identical event kinematics
@pytest.fixture
def basis(tmp_path):
    values = np.arange(4 * 64 * 2).reshape(4, 64, 2)
    amplitudes = values + 1j * (values + 1)
    intensity = (np.abs(amplitudes)**2).sum(axis=2)
    kinematics = np.arange(64 * 5, dtype=np.float64).reshape(64, 5)
    parts = []
    for number, indices in enumerate(([2, 0], [], [3, 1])):
        part = tmp_path / f"part_{number}.bin"
        amplitudes[indices].tofile(part)
        intensity[indices].tofile(str(part) + ".intensity")
        kinematics.tofile(str(part) + ".kinematics")
        Path(str(part) + ".json").write_text(json.dumps({
            "shape": [len(indices), 64, 2], "layout": "component,event,helicity,re_im",
            "indices": indices, "components": [f"component_{index}" for index in indices],
            "failures": [0] * len(indices),
        }))
        parts.append(part)
    return parts, amplitudes, intensity, kinematics


# Preserve every complex amplitude and its component order through a merge
def test_combine_basis(tmp_path, basis):
    parts, amplitudes, intensity, kinematics = basis
    output = tmp_path / "amplitudes.bin"
    combine_basis(parts, output)
    np.testing.assert_array_equal(np.fromfile(output, dtype=np.complex128).reshape(amplitudes.shape), amplitudes)
    np.testing.assert_array_equal(np.fromfile(str(output) + ".intensity").reshape(intensity.shape), intensity)
    np.testing.assert_array_equal(np.fromfile(str(output) + ".kinematics").reshape(kinematics.shape), kinematics)
    metadata = json.loads(Path(str(output) + ".json").read_text())
    assert metadata["indices"] == list(range(4))
    assert metadata["components"] == [f"component_{index}" for index in range(4)]
    assert not list(tmp_path.glob("*.partial"))


# Keep failed writes out of complete bank files and recover from the same shards
@pytest.mark.parametrize("existing", [False, True])
def test_combine_write_failure(tmp_path, basis, existing):
    parts, amplitudes, _, _ = basis
    output = tmp_path / "amplitudes.bin"
    if existing:
        combine_basis(parts, output)
    before = {path.name: path.read_bytes() for path in tmp_path.glob("amplitudes.bin*")}
    with file_limit(1024), pytest.raises(OSError) as caught:
        combine_basis(parts, output)
    assert caught.value.errno == errno.EFBIG
    assert any(str(output) in note for note in caught.value.__notes__)
    assert {path.name: path.read_bytes() for path in tmp_path.glob("amplitudes.bin*")
            if not path.name.endswith(".partial")} == before
    assert Path(str(output) + ".partial").stat().st_size == 1024
    combine_basis(parts, output)
    np.testing.assert_array_equal(np.fromfile(output, dtype=np.complex128).reshape(amplitudes.shape), amplitudes)


# Validate all shard sizes before changing an existing complete bank
def test_combine_incomplete_intensity(tmp_path, basis):
    parts, _, _, _ = basis
    output = tmp_path / "amplitudes.bin"
    combine_basis(parts, output)
    before = {path.name: path.read_bytes() for path in tmp_path.glob("amplitudes.bin*")}
    Path(str(parts[-1]) + ".intensity").write_bytes(b"")
    with pytest.raises(ValueError, match="file size"):
        combine_basis(parts, output)
    assert {path.name: path.read_bytes() for path in tmp_path.glob("amplitudes.bin*")} == before


# Preserve the merge exception when initialization cannot publish its failure record
def test_init_write_failure(tmp_path, basis):
    parts, _, _, _ = basis
    output = tmp_path / "amplitudes.bin"
    args = SimpleNamespace(backend="ray", phase="init", coord_dir=str(tmp_path / "coord"))

    # Exhaust permitted storage only after the running state has been published
    def bootstrap():
        resource.setrlimit(resource.RLIMIT_FSIZE, (0, resource.getrlimit(resource.RLIMIT_FSIZE)[1]))
        combine_basis(parts, output)

    with file_limit(4096), pytest.raises(OSError) as caught:
        run_initialization(args=args, bootstrapper=bootstrap)
    assert caught.value.errno == errno.EFBIG
    assert any("Amplitude bank merge" in note for note in caught.value.__notes__)
    assert any("Could not write initialization failure" in note for note in caught.value.__notes__)
    assert json.loads((Path(args.coord_dir) / "init.json").read_text())["status"] == "running"
    assert not output.exists()


# Read and stage interleaved shards without duplicating them in the completed shared bank
def test_stage_shards(tmp_path, basis):
    parts, amplitudes, _, _ = basis
    shared = tmp_path / "shared"
    (shared / "parts").mkdir(parents=True)
    for part in parts:
        for suffix in ("", ".json", ".intensity", ".kinematics"):
            shutil.copyfile(str(part) + suffix, shared / "parts" / (part.name + suffix))
    combined = tmp_path / "combined/amplitudes.bin"
    combine_basis(parts, combined)
    metadata = json.loads(Path(str(combined) + ".json").read_text())
    metadata["parts"] = [f"parts/{part.name}" for part in parts]
    (shared / "amplitudes.bin.json").write_text(json.dumps(metadata))
    (shared / "finalized.pkl").write_bytes(b"completed predictions")
    before = {path: path.read_bytes() for path in shared.rglob("*") if path.is_file()}
    mapped = Basis(shared, metadata)
    np.testing.assert_array_equal(mapped[:].numpy(), amplitudes)
    np.testing.assert_array_equal(mapped[0].numpy(), amplitudes[0])
    np.testing.assert_array_equal(mapped[1:3, 7:13].numpy(), amplitudes[1:3, 7:13])
    assert mapped[0:0, 7:13].shape == (0, 6, 2)
    weights = torch.ones(4, dtype=torch.complex128, requires_grad=True)
    actual = torch.einsum("c,ceh->eh", weights, mapped[:])
    expected = torch.einsum("c,ceh->eh", weights, torch.from_numpy(amplitudes))
    torch.testing.assert_close(actual, expected, rtol=0, atol=0)
    torch.testing.assert_close(torch.autograd.grad(actual.abs().square().sum(), weights)[0],
                              torch.autograd.grad(expected.abs().square().sum(), weights)[0], rtol=0, atol=0)
    staged = tmp_path / "staged"
    stage_bank(shared, staged, predictions=False)
    assert not (staged / "parts").exists()
    assert not (staged / "finalized.pkl").exists()
    assert not (shared / "amplitudes.bin").exists()
    assert before == {path: path.read_bytes() for path in before}
    staged_metadata = json.loads((staged / "amplitudes.bin.json").read_text())
    np.testing.assert_array_equal(Basis(staged, staged_metadata)[:].numpy(), amplitudes)
    assert not list(staged.rglob("*.bin"))
    assert staged_metadata["directory"] == str(shared)
    # A second local staging reads the same shared bytes after an interrupted job
    stage_bank(shared, tmp_path / "retry")
    assert (tmp_path / "retry/finalized.pkl").read_bytes() == b"completed predictions"


# Reject incomplete or overlapping mapped components before evaluating fit amplitudes
@pytest.mark.parametrize("damage", ["missing", "overlap", "truncated"])
def test_mapped_basis_invalid(tmp_path, basis, damage):
    parts, _, _, _ = basis
    output = tmp_path / "merged/amplitudes.bin"
    combine_basis(parts, output)
    metadata = json.loads(Path(str(output) + ".json").read_text())
    metadata["parts"] = [part.name for part in parts]
    if damage == "missing":
        metadata["parts"].pop()
    elif damage == "overlap":
        metadata["parts"].append(parts[0].name)
    else:
        parts[-1].write_bytes(b"")
    with pytest.raises(ValueError, match="components|dimensions"):
        Basis(tmp_path, metadata)
