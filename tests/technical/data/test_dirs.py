# Shared directory creation and nested JSON output
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import pytest
from core.io.files import ensure_dir
from core.io.serialize import load_json_file, write_json_file


# Allow concurrent writers to prepare the same missing directory tree
def test_shared_directory(tmp_path):
    """Create nested folders concurrently and reject files in their place"""
    target = tmp_path / "figs" / "nested plots"
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(ensure_dir, [target] * 16))
    assert target.is_dir()
    filename = target / "record.json"
    filename.write_text("{}")
    with pytest.raises(NotADirectoryError):
        ensure_dir(filename)
    with pytest.raises(NotADirectoryError):
        ensure_dir(filename / "nested")


# Keep directory allocation exclusive when multiple workers request one name
def test_exclusive_directory(tmp_path):
    """Allow exactly one worker to allocate an exclusive directory"""
    target = tmp_path / "missing" / "trial"
    with ThreadPoolExecutor(max_workers=4) as pool:
        tasks = [pool.submit(ensure_dir, target, exist_ok=False, mode=0o700) for _ in range(8)]
    assert sum(task.exception() is None for task in tasks) == 1
    assert all(task.exception() is None or isinstance(task.exception(), FileExistsError) for task in tasks)
    assert target.is_dir()


# Let the shared JSON writer prepare output parents without caller setup
def test_nested_json(tmp_path):
    """Write and read JSON beneath initially missing folders"""
    filename = tmp_path / "output" / "nested" / "record.json"
    payload = {"values": [1, 2, 3]}
    write_json_file(filename, payload)
    assert load_json_file(filename) == payload


# Always attempt directory creation even when retries are disabled
def test_directory_without_retries(tmp_path):
    """Create a real output directory with the retry count set to zero"""
    target = tmp_path / "runs" / "nested"
    subprocess.run(
        [sys.executable, "-c", "import sys; from core.io.files import ensure_dir; ensure_dir(sys.argv[1])", str(target)],
        env={**os.environ, "ICETUNE_IO_DIR_RETRIES": "0"}, check=True,
    )
    assert target.is_dir()
