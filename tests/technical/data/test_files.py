# Check independent dated output directories under concurrent creation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from concurrent.futures import ThreadPoolExecutor

from core.io.files import dated_directory


# Concurrent runs preserve separate outputs even when started in the same second
def test_dated_directories(tmp_path):
    with ThreadPoolExecutor(max_workers=4) as pool:
        paths = list(pool.map(lambda _: dated_directory(tmp_path, prefix="icepacks."), range(16)))
    assert len(set(paths)) == len(paths)
    for index, path in enumerate(paths):
        (path / "result.txt").write_text(str(index))
    assert [int((path / "result.txt").read_text()) for path in paths] == list(range(len(paths)))
