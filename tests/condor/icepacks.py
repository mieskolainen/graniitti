# Distributed measurement icepacks and their complete result table
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import csv
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

from core.io.files import dated_directory
from core.io.steering import load_dataset
from core.stats.validation import check_physics, comparison_samples

from tests.condor.job_runner import ROOT, write_json


# Discover measurement cards by default or all cards below a selected icepack folder
def datasets(name=None):
    base = (ROOT / "icepack").resolve()
    selected = (base / name).resolve() if name else base
    if name and (Path(name).is_absolute() or selected == base
                 or not selected.is_relative_to(base) or not selected.is_dir()):
        raise ValueError(f"Unknown icepack folder: {name}")
    entries = []
    for path in sorted(selected.rglob("dataset.json")):
        if any(part.startswith("_") or "._old" in part for part in path.relative_to(base).parts):
            continue
        card, _ = load_dataset(str(path), cdir=ROOT)
        measurement = "measurement" in card.get("validation", {})
        if name or measurement:
            entries.append(dict(entry=str(path.parent.relative_to(base)), active=card["active"], measurement=measurement))
    if not entries:
        raise ValueError("No icepacks selected")
    return entries


# Keep generator outputs and caches separate even when two submissions use the same dataset
def workspace(run, index):
    work = dated_directory(run / "work", prefix=f"{index}.")
    for name in ("HEPData", "modeldata", "MG5cards", "python", "install", "tests"):
        (work / name).symlink_to(ROOT / name, target_is_directory=(ROOT / name).is_dir())
    shutil.copytree(ROOT / "icepack", work / "icepack", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy2(ROOT / "VERSION.json", work / "VERSION.json")
    (work / "bin").mkdir()
    for name in ("gr", "pythia_lhe_hadronize"):
        if (ROOT / "bin" / name).is_file():
            shutil.copy2(ROOT / "bin" / name, work / "bin" / name)
    for name in ("eikonal", "nuclear", "sudakov", "vgrid"):
        shutil.copytree(ROOT / name, work / name) if (ROOT / name).is_dir() else (work / name).mkdir()
    return work


# Run the normal generator and iceplot path and always retain its exit status
def worker(run, index):
    config = json.loads((run / "config.json").read_text())
    item = config["entries"][index]
    start = time.monotonic()
    result = dict(entry=item["entry"], status="RUNNING", failures=[], comparisons=0)
    write_json(run / "results" / f"{index}.json", result)
    result["status"] = "ERROR"
    try:
        work = workspace(run, index)
        output = run / "outputs" / str(index)
        command = ["bash", "icepack/run.sh", "--output-dir", str(output), *config["steering"], item["entry"]]
        environment = {key: value for key, value in os.environ.items()
                       if key not in {"NEVENTS", "LOOPSCREEN", "WEIGHTED", "INTEGRATOR", "DENSITY", "ICEPACK_REPORT_DIR"}}
        with (run / "logs" / f"{index}.log").open("w") as log:
            code = subprocess.run(command, cwd=work, env=environment, stdout=log, stderr=subprocess.STDOUT).returncode
        path = output / item["entry"].replace("/", "__") / "reports" / (item["entry"].replace("/", "__") + ".json")
        result.update(exit_code=code, log=f"logs/{index}.log", report=str(path.relative_to(run)))
        report = json.loads(path.read_text())
        result["comparisons"] = len(list(comparison_samples(report)))
        if item["measurement"] and "measurement" not in report.get("validation", {}):
            result["failures"].append("Measurement criteria missing from report")
        else:
            try:
                check_physics(report)
            except ValueError as exc:
                result.update(status="FAIL", failures=str(exc).splitlines())
            else:
                result["failures"].extend(report["measurement"]["failures"])
                if result["failures"]:
                    result["status"] = "FAIL"
                elif code == 0:
                    result["status"] = "PASS" if item["measurement"] else "DONE"
        if code and not result["failures"]:
            result["failures"].append(f"Run exited with status {code}")
    except Exception as exc:
        result["failures"].append(f"{type(exc).__name__}: {exc}")
    result["seconds"] = round(time.monotonic() - start, 4)
    write_json(run / "results" / f"{index}.json", result)
    return int(result["status"] not in {"PASS", "DONE"})


# Include every expected dataset so failed or absent jobs cannot disappear from the table
def summary(run, final=False):
    config = json.loads((run / "config.json").read_text())
    if (run / "report.json").is_file():
        final = final or json.loads((run / "report.json").read_text()).get("final", False)
    rows = []
    for index, item in enumerate(config["entries"]):
        path = run / "results" / f"{index}.json"
        row = dict(entry=item["entry"], status="MISSING" if final else "PENDING", comparisons=0, seconds="", failures=[])
        if not item["active"]:
            row.update(status="SKIP", failures=["Inactive dataset card"])
        elif path.is_file():
            try:
                result = json.loads(path.read_text())
                if result["entry"] != item["entry"]:
                    raise ValueError("Result does not match dataset")
                row.update(result)
                if final and row["status"] == "RUNNING":
                    row.update(status="MISSING", failures=["Worker did not finish"])
            except (ValueError, KeyError) as exc:
                row.update(status="ERROR", failures=[str(exc)])
        rows.append(row)
    write_json(run / "report.json", dict(final=final, runs=rows))
    columns = ("Dataset", "Result", "Comparisons", "Seconds", "Reason", "Report / log")
    values = [(r["entry"], r["status"], r["comparisons"], r["seconds"], " / ".join(r["failures"]),
               " / ".join(r.get(key, "") for key in ("report", "log"))) for r in rows]
    with (run / "report.csv").open("w", newline="") as stream:
        csv.writer(stream).writerows([columns, *values])
    table = "\n".join("| " + " | ".join(str(value).replace("|", "\\|").replace("\n", " ") for value in row) + " |"
                      for row in [columns, ["---"] * len(columns), *values]) + "\n"
    (run / "report.md").write_text(table)
    print(table, end="")
    return int(any(row["status"] not in {"PASS", "DONE", "SKIP"} for row in rows))


# Parse measurement steering once before any worker starts
def steering():
    arguments = []
    for name in ("NEVENTS", "LOOPSCREEN", "WEIGHTED", "INTEGRATOR"):
        value = os.environ.get(name, "VEGAS" if name == "INTEGRATOR" else "auto")
        if value != "auto":
            valid = value.isdigit() and int(value) > 0 if name == "NEVENTS" else value in (
                {"VEGAS", "NEUROJAC"} if name == "INTEGRATOR" else {"0", "1", "true", "false"})
            if not valid:
                raise ValueError(f"Invalid {name}={value}")
            arguments.extend([f"--{name}", value])
    return arguments


# Execute one measurement job or produce its final result table
def main():
    parser = argparse.ArgumentParser(description="Execute an icepack or refresh its measurement table")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--run", nargs=2, metavar=("DIRECTORY", "INDEX"))
    action.add_argument("--report", type=Path)
    parser.add_argument("--final", action="store_true")
    args = parser.parse_args()
    if args.run:
        return worker(Path(args.run[0]).resolve(), int(args.run[1]))
    return summary(args.report.resolve(), args.final)


if __name__ == "__main__":
    raise SystemExit(main())
