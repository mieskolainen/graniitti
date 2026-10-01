# Tests for the icetune bootstrap content cache
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from core.tune import cache as icetune_cache
from core.tune import short_id

FINGERPRINT = "a" * 64


# Preserve literal process names when extracting VEGAS grids from the bootstrap
@pytest.mark.parametrize("tar_options", ["", "--wildcards"])
def test_bootstrap_stages_literal_process_names(tmp_path, monkeypatch, tar_options):
    from core.tune.drivers.graniitti.driver import GraniittiDriver
    from core.tune.drivers.graniitti.runtime import allowed_bootstrap_path

    monkeypatch.setenv("TAR_OPTIONS", tar_options)
    source = tmp_path / "source"
    relative = Path("vgrid/CMS_2752118_pipi__GP[RES+CON]<F>__INTEGRATOR_VEGAS__TUNE0.vgrid")
    output = source / relative
    output.parent.mkdir(parents=True)
    output.write_bytes(b"VEGAS grid bytes")
    temporary = tmp_path / "publish"
    temporary.mkdir()
    bootstrap = icetune_cache.publish_content_archive(
        cache_base_url=str(tmp_path / "cache"), kind="graniitti", fingerprint=FINGERPRINT,
        root=source, records=icetune_cache.file_records([output], root=source, allowed=allowed_bootstrap_path),
        temporary_dir=temporary, schema_version=GraniittiDriver.BOOTSTRAP_SCHEMA_VERSION, label="GRANIITTI",
    )
    root = tmp_path / "worker"
    driver = GraniittiDriver()
    driver.stage_bootstrap(cdir=str(root), bootstrap=bootstrap)
    assert (root / relative).read_bytes() == output.read_bytes()


# Reuse shared amplitude references under another run name and worker checkout
def test_amp_bootstrap_rebases_run_worker(tmp_path):
    from core.tune.drivers.graniitti.ampfit.amplitude import bank_files
    from core.tune.drivers.graniitti.driver import GraniittiDriver
    from core.tune.drivers.graniitti.runtime import allowed_bootstrap_path

    source = tmp_path / "source"
    relative = Path("runs/icetune/initial/results/amplitude/0/nominal")
    payload = bytes(range(256))
    output = source / relative / "amplitudes.bin"
    output.parent.mkdir(parents=True)
    output.write_bytes(payload)
    for name in ("request.json", "steering.json", "amplitudes.bin.json"):
        (output.parent / name).write_text("{}")
    (output.parent / "events.hepmc3").touch()
    files = bank_files(source / relative.parents[1])
    temporary = tmp_path / "publish"
    temporary.mkdir()
    bootstrap = icetune_cache.publish_content_archive(
        cache_base_url=str(tmp_path / "cache"), kind="graniitti", fingerprint=FINGERPRINT,
        root=source, records=icetune_cache.file_records(files, root=source, allowed=allowed_bootstrap_path),
        temporary_dir=temporary, schema_version=GraniittiDriver.BOOTSTRAP_SCHEMA_VERSION, label="GRANIITTI",
    )
    for name in ("head", "worker"):
        root = tmp_path / name
        previous = root / "sudakov"
        previous.mkdir(parents=True)
        (previous / ".gitignore").write_text("*\n!.gitignore\n")
        (previous / "other.json").write_text("{\"cache\": 1}\n")
        driver = GraniittiDriver()
        driver.stage_bootstrap(cdir=str(root), bootstrap=bootstrap)
        saved = list(root.glob("sudakov.*._old"))
        assert len(saved) == 1
        assert (saved[0] / "other.json").read_text() == "{\"cache\": 1}\n"
        assert (previous / ".gitignore").read_text() == "*\n!.gitignore\n"
        assert not (previous / "other.json").exists()
        directory = driver.amplitude_directory(cdir=str(root), run_name="resumed")
        metadata = json.loads((directory / "0/nominal/amplitudes.bin.json").read_text())
        assert (Path(metadata["directory"]) / "amplitudes.bin").read_bytes() == payload
        assert not list(directory.rglob("*.bin"))
        driver.stage_bootstrap(cdir=str(root), bootstrap=bootstrap)
        files = driver.runtime_files({"cdir": str(root), "run_name": "resumed",
                                      "optimization": {"optimizer": "ampfit"}})
        assert set(files) == {(relative / name).as_posix() for name in
                              ("request.json", "steering.json", "events.hepmc3", "amplitudes.bin.json")}


# Check immutable publication copies bytes and permits only one concurrent writer
def test_immutable_publication_source_snapshot(tmp_path):
    sources = [tmp_path / "first.dat", tmp_path / "second.dat"]
    for source in sources:
        source.write_bytes(source.name.encode("ascii"))
    target = tmp_path / "cache" / "output.dat"
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(icetune_cache.publish_file_immutable, source, str(target)) for source in sources]
        installed = [future.result() for future in futures]
    assert installed.count(True) == 1
    expected = sources[installed.index(True)].read_bytes()
    for source in sources:
        source.write_bytes(b"scratch file reused")
    assert target.read_bytes() == expected


# Compute cache records for one generated GRANIITTI file
def _records(root: Path) -> list[dict]:
    output = root / "vgrid" / "grid.dat"
    return icetune_cache.file_records(
        [output],
        root=root,
        allowed=lambda relative: relative.parts[0] == "vgrid",
    )


# Check the input manifest commits an archive stored by its own checksum
def test_content_archive_output_checksum(tmp_path):
    root = tmp_path / "root"
    output = root / "vgrid" / "grid.dat"
    output.parent.mkdir(parents=True)
    output.write_text("first\n", encoding="utf-8")
    temporary_dir = tmp_path / "temporary"
    temporary_dir.mkdir()
    cache = tmp_path / "cache"

    manifest = icetune_cache.publish_content_archive(
        cache_base_url=str(cache),
        kind="graniitti",
        fingerprint=FINGERPRINT,
        root=root,
        records=_records(root),
        temporary_dir=temporary_dir,
        schema_version=1,
        label="GRANIITTI",
    )

    archive_sha256 = manifest["archive_sha256"]
    assert manifest["archive_url"] == str(
        cache / "graniitti" / "objects" / f"{short_id(archive_sha256)}.tar.zst"
    )
    assert (cache / "graniitti" / f"{short_id(FINGERPRINT)}.json").is_file()
    assert Path(manifest["archive_url"]).is_file()
    assert icetune_cache.sha256_file(manifest["archive_url"]) == archive_sha256
    with pytest.raises(icetune_cache.PermanentConfigurationError, match="identity mismatch"):
        icetune_cache.reuse_content_archive(
            cache_base_url=str(cache), kind="graniitti", fingerprint=short_id(FINGERPRINT) + "b" * 48,
            temporary_dir=temporary_dir, schema_version=1, label="GRANIITTI")


# Check orphan archives and concurrent differing publishers cannot cause collisions
def test_content_manifest_has_one_winner(tmp_path):
    root = tmp_path / "root"
    output = root / "vgrid" / "grid.dat"
    output.parent.mkdir(parents=True)
    output.write_text("first\n", encoding="utf-8")
    cache = tmp_path / "cache"
    cache_kind = cache / "graniitti"
    cache_kind.mkdir(parents=True)
    orphan = cache_kind / f"{short_id(FINGERPRINT)}.tar.zst"
    orphan.write_bytes(b"incomplete earlier publication")
    orphan_bytes = b"uncommitted content object"
    orphan_sha256 = hashlib.sha256(orphan_bytes).hexdigest()
    content_orphan = cache_kind / "objects" / f"{short_id(orphan_sha256)}.tar.zst"
    content_orphan.parent.mkdir()
    content_orphan.write_bytes(orphan_bytes)
    first_temporary = tmp_path / "first"
    first_temporary.mkdir()

    first = icetune_cache.publish_content_archive(
        cache_base_url=str(cache),
        kind="graniitti",
        fingerprint=FINGERPRINT,
        root=root,
        records=_records(root),
        temporary_dir=first_temporary,
        schema_version=1,
        label="GRANIITTI",
    )
    first_archive = Path(first["archive_url"])
    first_bytes = first_archive.read_bytes()

    output.write_text("second\n", encoding="utf-8")
    winner = icetune_cache.publish_content_archive(
        cache_base_url=str(cache),
        kind="graniitti",
        fingerprint=FINGERPRINT,
        root=root,
        records=_records(root),
        temporary_dir=first_temporary,
        schema_version=1,
        label="GRANIITTI",
    )

    assert winner == first
    assert first_archive.read_bytes() == first_bytes
    assert len(list((cache_kind / "objects").glob("*.tar.zst"))) == 2
    assert orphan.read_bytes() == b"incomplete earlier publication"
    assert content_orphan.read_bytes() == orphan_bytes


# Check previous damaged or missing cache content can be rebuilt
@pytest.mark.parametrize("reuse", [False, True])
@pytest.mark.parametrize("orphan", [False, True])
@pytest.mark.parametrize("damage", ["checksum", "size", "missing"])
def test_previous_content_can_be_rebuilt(tmp_path, reuse, orphan, damage):
    root = tmp_path / "root"
    output = root / "vgrid" / "grid.dat"
    output.parent.mkdir(parents=True)
    output.write_bytes(b"grid data")
    temporary_dir = tmp_path / "temporary"
    temporary_dir.mkdir()
    options = dict(
        cache_base_url=str(tmp_path / "cache"),
        kind="graniitti",
        fingerprint=FINGERPRINT,
        root=root,
        records=_records(root),
        temporary_dir=temporary_dir,
        schema_version=1,
        label="GRANIITTI",
    )
    manifest = icetune_cache.publish_content_archive(**options)
    archive = Path(manifest["archive_url"])
    if orphan:
        (tmp_path / "cache" / "graniitti" / f"{short_id(FINGERPRINT)}.json").rename(temporary_dir / "manifest._old")
    damaged = b"x" * (manifest["archive_size"] if damage == "checksum" else 1)
    if damage == "missing":
        archive.rename(temporary_dir / "missing._old")
    else:
        archive.write_bytes(damaged)

    if reuse:
        assert icetune_cache.reuse_content_archive(
            **{key: value for key, value in options.items() if key not in {"root", "records"}}
        ) is None
        if not orphan:
            output.write_bytes(b"regenerated grid data")
            options["records"] = _records(root)
    rebuilt = icetune_cache.publish_content_archive(**options)
    assert rebuilt["fingerprint"] == manifest["fingerprint"]
    assert (rebuilt != manifest) == (reuse and not orphan)
    assert icetune_cache.sha256_file(rebuilt["archive_url"]) == rebuilt["archive_sha256"]
    if damage != "missing":
        previous = list(archive.parent.glob(f"{archive.name}.*._old"))
        assert len(previous) == 1
        assert previous[0].read_bytes() == damaged
