# Download complete versioned HEPData records as unmodified JSON responses
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

"""
Usage: python HEPData/download.py [PHOTOPROD ...] [--output tmp/hepdata]
"""

import argparse
import json
import re
import sys
from pathlib import Path
from shutil import copy2
from tempfile import NamedTemporaryFile

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


# Discover the requested versioned record directories without a separate record list
def records(root: Path, selections: list[str]) -> list[Path]:
    found = set()
    for selection in selections or [str(root)]:
        path = Path(selection)
        path = (path if path.exists() else root / path).resolve()
        if not path.is_relative_to(root) or not path.is_dir():
            raise ValueError(f"Not a directory under {root}: {selection}")
        matches = [
            p
            for p in [path, *path.rglob("HEPData-*-json")]
            if p.is_dir() and re.fullmatch(r"HEPData-ins\d+-v[1-9]\d*-json", p.name)
        ]
        if not matches:
            raise ValueError(f"No versioned HEPData records in {selection}")
        found.update(matches)
    return sorted(found)


# Preserve existing table filenames by DOI and use HEPData names for new tables
def table_names(source: Path, tables: list[dict]) -> list[str]:
    existing = {}
    for path in source.glob("*.json"):
        doi = json.loads(path.read_bytes())["doi"]
        if doi in existing:
            raise ValueError(f"Duplicate local table DOI {doi} in {source}")
        existing[doi] = path.name
    names = [existing.get(t["doi"], t["processed_name"] + ".json") for t in tables]
    if len(set(names)) != len(names) or len({t["doi"] for t in tables}) != len(tables):
        raise ValueError(f"Duplicate table names or DOIs in {source}")
    if any(Path(name).name != name for name in names):
        raise ValueError(f"Invalid table filename in {source}")
    return names


# Reject incomplete tables and incorrect versions before writing any response
def check_table(raw: bytes, doi: str) -> None:
    data = json.loads(raw)
    if (
        not isinstance(data, dict)
        or data.get("doi") != doi
        or not isinstance(data.get("headers"), list)
        or not isinstance(data.get("values"), list)
    ):
        raise ValueError(f"Response is not the complete JSON table {doi}")


# Save response bytes atomically and retain every replaced file as ._old
def save(path: Path, raw: bytes, temporary: Path) -> str:
    if path.exists() and path.read_bytes() == raw:
        return "unchanged"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=temporary, prefix="hepdata-", suffix=".json", delete=False) as stream:
        stream.write(raw)
        staged = Path(stream.name)
    backup = None
    if path.exists():
        backup = path.with_name(path.name + "._old")
        while backup.exists():
            backup = backup.with_name(backup.name + "._old")
        copy2(path, backup)
    staged.replace(path)
    return "replaced" if backup is not None else "downloaded"


# Fetch every table from the selected record version using the published data identifiers
def download(session: requests.Session, source: Path, destination: Path, temporary: Path, timeout: float) -> int:
    record, version = re.fullmatch(r"HEPData-(ins\d+)-v([1-9]\d*)-json", source.name).groups()
    url = f"https://www.hepdata.net/record/{record}?version={version}&format=json"
    response = session.get(url, timeout=timeout)
    response.raise_for_status()
    data = response.json()
    if data["version"] != int(version) or str(data["record"]["inspire_id"]) != record[3:]:
        raise ValueError(f"Incorrect record or version from {url}")
    tables = data["data_tables"]
    if not tables:
        raise ValueError(f"No tables listed at {url}")
    names = table_names(source, tables)
    print(f"{record} v{version}: {len(tables)} tables -> {destination}", flush=True)
    failures = 0
    for table, name in zip(tables, names, strict=True):
        url = f"https://www.hepdata.net/record/data/{data['recid']}/{table['id']}/{version}/"
        try:
            response = session.get(url, timeout=timeout)
            response.raise_for_status()
            check_table(response.content, table["doi"])
            status = save(destination / name, response.content, temporary)
            print(f"  {status}: {name}", flush=True)
        except (requests.RequestException, OSError, ValueError) as error:
            failures += 1
            print(f"  FAILED {url}: {error}", file=sys.stderr, flush=True)
    return failures


# Download all repository records by default or selected categories and record directories
def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="HEPData categories or record directories (default: all)")
    parser.add_argument("--output", type=Path, default=root, help="output root, preserving the category directories")
    parser.add_argument("--list", action="store_true", help="list selected record versions without downloading")
    parser.add_argument("--timeout", type=float, default=30.0, help="HTTP timeout in seconds (default: 30)")
    parser.add_argument("--search", nargs="+", help="search HEPData publication metadata without downloading tables")
    args = parser.parse_args()
    if not 0 < args.timeout < float("inf"):
        parser.error("--timeout must be finite and positive")
    try:
        selected = [] if args.search else records(root, args.paths)
    except ValueError as error:
        parser.error(str(error))
    if args.list:
        for source in selected:
            print(source.relative_to(root))
        return 0
    failures = 0
    with requests.Session() as session:
        retry = Retry(total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"])
        session.mount("https://", HTTPAdapter(max_retries=retry))
        for query in args.search or []:
            try:
                response = session.get("https://www.hepdata.net/search/",
                                       params={"q": query, "format": "json", "size": 100}, timeout=args.timeout)
                response.raise_for_status()
                data = response.json()
                print(json.dumps({"query": query, "total": data["total"], "results": [
                    {key: row.get(key) for key in ("id", "inspire_id", "title", "doi", "hepdata_doi")}
                    for row in data["results"]]}, indent=2), flush=True)
            except (requests.RequestException, ValueError, KeyError) as error:
                failures += 1
                print(f"FAILED search {query}: {error}", file=sys.stderr, flush=True)
        for source in selected:
            destination = args.output.resolve() / source.relative_to(root)
            try:
                failures += download(session, source, destination, root.parent / "tmp", args.timeout)
            except (requests.RequestException, OSError, ValueError, KeyError) as error:
                failures += 1
                print(f"FAILED {source.name}: {error}", file=sys.stderr, flush=True)
    print(f"{len(selected)} records processed, {failures} failures", flush=True)
    return int(failures > 0)


if __name__ == "__main__":
    sys.exit(main())
