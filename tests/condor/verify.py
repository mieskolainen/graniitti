# Verify all Condor test results and write final reports
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import json
import math
import pathlib
import xml.etree.ElementTree as element_tree


# Parse the final DAG node invocation
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=pathlib.Path)
    parser.add_argument("dag_status", type=int)
    parser.add_argument("failed_count", type=int)
    return parser.parse_args()


# Read one indexed array manifest
def read_manifest(path: pathlib.Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    seen_indices: set[str] = set()
    seen_items: set[str] = set()
    with path.open(encoding="ascii") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split(maxsplit=1)
            if len(fields) != 2 or not fields[0].isdigit():
                raise RuntimeError(f"invalid manifest line {path}:{line_number}")
            index, item = fields
            if index in seen_indices or item in seen_items:
                raise RuntimeError(f"duplicate manifest entry {path}:{line_number}")
            seen_indices.add(index)
            seen_items.add(item)
            entries.append((index, item))
    if not entries:
        raise RuntimeError(f"empty manifest: {path}")
    return entries


# Compute every expected result grouped by category
def expected_results(run_dir: pathlib.Path) -> tuple[dict[str, list[tuple[str, str]]], list[str]]:
    config = json.loads((run_dir / "config.json").read_text(encoding="ascii"))
    expected = {"build": [("000000", "build")]} if config["build"] else {}
    errors: list[str] = []
    for category in config["categories"]:
        path = run_dir / "manifests" / f"{category}.tsv"
        try:
            expected[category] = read_manifest(path)
        except (OSError, RuntimeError) as exc:
            expected[category] = []
            errors.append(str(exc))
    return expected, errors


# Load and validate one result record
def load_result(
    path: pathlib.Path,
    category: str,
    index: str,
    item: str,
    members: list[str],
) -> dict:
    with path.open(encoding="ascii") as handle:
        result = json.load(handle)
    identity = (result.get("category"), result.get("index"), result.get("item"))
    if identity != (category, index, item):
        raise RuntimeError(f"result identity mismatch: {path}")
    if result.get("members") != members:
        raise RuntimeError(f"result member mismatch: {path}")
    if result.get("outcome") not in ("PASS", "FAIL"):
        raise RuntimeError(f"invalid result outcome: {path}")
    return result


# Compute the concrete tests represented by one array entry
def expected_members(run_dir: pathlib.Path, category: str, item: str) -> list[str]:
    if category != "python":
        return [item]
    shard_path = (run_dir / item).resolve()
    task_root = (run_dir / "tasks").resolve()
    if not shard_path.is_relative_to(task_root):
        raise RuntimeError(f"Python shard is outside the task directory: {item}")
    with shard_path.open(encoding="ascii") as handle:
        shard = json.load(handle)
    members = shard.get("members") if isinstance(shard, dict) else None
    if (
        not isinstance(members, list)
        or not members
        or not all(isinstance(member, str) and member for member in members)
    ):
        raise RuntimeError(f"invalid Python shard: {shard_path}")
    return members


# Collect suite result counts and failure details
def collect_results(
    run_dir: pathlib.Path,
    expected: dict[str, list[tuple[str, str]]],
) -> tuple[dict[str, dict], list[dict], list[str]]:
    suites: dict[str, dict] = {}
    results: list[dict] = []
    errors: list[str] = []
    for category in expected:
        suite = {"expected": 0, "failed": 0, "missing": 0, "passed": 0}
        for index, item in expected[category]:
            try:
                members = expected_members(run_dir, category, item)
            except (OSError, ValueError, RuntimeError) as exc:
                suite["failed"] += 1
                errors.append(str(exc))
                continue
            member_count = len(members)
            suite["expected"] += member_count
            path = run_dir / "results" / category / f"{index}.json"
            if not path.is_file():
                suite["missing"] += member_count
                errors.append(f"missing result [{category}:{index}] {item}")
                continue
            try:
                result = load_result(path, category, index, item, members)
            except (OSError, ValueError, RuntimeError) as exc:
                suite["failed"] += member_count
                errors.append(str(exc))
                continue
            results.append(result)
            if result["outcome"] == "PASS":
                suite["passed"] += member_count
            else:
                suite["failed"] += member_count
                errors.append(
                    f"failed result [{category}:{index}] {item} exit={result['exit_code']}"
                )
        suites[category] = suite
    return suites, results, errors


# Aggregate Python JUnit counts from available reports
def collect_junit(run_dir: pathlib.Path) -> tuple[dict[str, float | int], list[str]]:
    totals: dict[str, float | int] = {
        "errors": 0,
        "failures": 0,
        "skipped": 0,
        "tests": 0,
        "time": 0.0,
    }
    errors: list[str] = []
    for path in sorted((run_dir / "junit").glob("python.*.xml")):
        try:
            root = element_tree.parse(path).getroot()
        except (OSError, element_tree.ParseError) as exc:
            errors.append(f"invalid JUnit report {path}: {exc}")
            continue
        suites = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
        for suite in suites:
            for key in ("errors", "failures", "skipped", "tests"):
                totals[key] += int(suite.attrib.get(key, 0))
            totals["time"] += float(suite.attrib.get("time", 0.0))
    return totals, errors


# Build reusable duration history from test results and timing plugins
def collect_timings(
    run_dir: pathlib.Path, results: list[dict]
) -> tuple[dict[str, float], list[str]]:
    timings: dict[str, float] = {}
    errors: list[str] = []
    for result in results:
        category = result["category"]
        if category in ("cpp", "study"):
            timings[f"{category}::{result['item']}"] = float(result["duration_seconds"])
    for path in sorted((run_dir / "junit").glob("python.*.timings.json")):
        try:
            with path.open(encoding="ascii") as handle:
                shard_timings = json.load(handle)
            for node, value in shard_timings.items():
                duration = float(value)
                if duration < 0.0 or not math.isfinite(duration):
                    raise ValueError(f"invalid duration {value!r}")
                timings[f"python::{node}"] = duration
        except (OSError, ValueError, AttributeError) as exc:
            errors.append(f"invalid timing report {path}: {exc}")
    return timings, errors


# Render the human readable report
def render_report(payload: dict) -> str:
    lines = [
        "GRANIITTI Condor test report",
        f"Run: {payload['config'].get('run_name', '<unknown>')}",
        f"Git revision: {payload['config'].get('git_revision') or '<unknown>'}",
        f"DAG status: {payload['dag_status']}",
        f"Failed DAG nodes: {payload['failed_count']}",
        "",
        "Suite       Expected  Passed  Failed  Missing",
    ]
    for category, suite in payload["suites"].items():
        lines.append(
            f"{category:<11} {suite['expected']:>8} {suite['passed']:>7} "
            f"{suite['failed']:>7} {suite['missing']:>8}"
        )
    junit = payload["junit"]
    lines.extend(
        (
            "",
            "Python test cases: "
            f"{junit['tests']} total, {junit['failures']} failed, {junit['errors']} errors, "
            f"{junit['skipped']} skipped, {junit['time']:.3f} seconds",
        )
    )
    if payload["errors"]:
        lines.extend(("", "Failures and missing results:"))
        lines.extend(f"  {error}" for error in payload["errors"])
    lines.extend(("", f"Overall: {payload['outcome']}"))
    return "\n".join(lines) + "\n"


# Verify the complete run and make the verifier exit status authoritative
def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    with (run_dir / "config.json").open(encoding="ascii") as handle:
        config = json.load(handle)
    expected, manifest_errors = expected_results(run_dir)
    suites, results, result_errors = collect_results(run_dir, expected)
    junit, junit_errors = collect_junit(run_dir)
    timings, timing_errors = collect_timings(run_dir, results)
    errors = manifest_errors + result_errors + junit_errors + timing_errors
    if args.dag_status != 0:
        errors.append(f"DAG status is {args.dag_status}")
    if args.failed_count != 0:
        errors.append(f"DAG reports {args.failed_count} failed nodes")
    payload = {
        "config": config,
        "dag_status": args.dag_status,
        "errors": errors,
        "failed_count": args.failed_count,
        "junit": junit,
        "outcome": "PASS" if not errors else "FAIL",
        "results": results,
        "suites": suites,
    }
    report = render_report(payload)
    (run_dir / "report.txt").write_text(report, encoding="ascii")
    (run_dir / "report.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    (run_dir / "timings.json").write_text(
        json.dumps(timings, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    print(report, end="")
    return 0 if payload["outcome"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
