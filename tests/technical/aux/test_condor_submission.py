# Test Condor DAG scheduling and manifest generation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import os
import shlex
import stat
import subprocess
import sys
from pathlib import Path

import pyjson5
import pytest

from tests.condor import job_runner as runner
from tests.condor import submit
from tests.condor import verify as verifier

ICEPACKS = sorted(
    path.name for path in (Path(__file__).resolve().parents[3] / "icepack").iterdir()
    if path.is_dir() and not path.name.startswith("_") and "._old" not in path.name
)


# Prepare every suite through the real launcher from outside the checkout
@pytest.mark.parametrize(("suite", "name"), [
    *((suite, None) for suite in ("all", "cpp", "pytest", "studies", "icepacks")),
    ("icepacks", "UPC/PHOTOPROD/ALICE_2658375"),
    ("icepacks", "UPC/PHOTOPROD/ALICE_2658375/jpsi_incoherent/"),
    ("icepacks", "GAMMA/continuum"),
])
def test_launchers(suite, name, tmp_path):
    root = Path(__file__).resolve().parents[3]
    output = tmp_path / "runs"
    subprocess.run(["bash", str(root / "tests/condor/submit.sh"), suite, "--dry-run",
                    "--output-root", str(output), *(["--name", name] if name else [])],
                   cwd=tmp_path, capture_output=True, text=True, check=True)
    run, = output.iterdir()
    config = json.loads((run / "config.json").read_text())
    assert config["suite"] == suite
    dag = (run / "condor/run.dag").read_text()
    assert "FINAL " in dag
    assert ("JOB BUILD " in dag) == (suite in {"all", "cpp"})
    for path in (run / "condor").glob("*.sub"):
        fields = dict(line.split(" = ", 1) for line in path.read_text().splitlines() if " = " in line)
        settings = submit.resource_settings()
        resources = {**settings["defaults"], **settings["jobs"][path.stem]}
        assert all(fields[key] == str(value) for key, value in resources.items())
        executable = Path(fields["executable"])
        assert executable.is_file() and executable.stat().st_mode & stat.S_IXUSR
        assert shlex.split(fields["arguments"].strip('"'))[0] == str(root)
        assert Path(fields["output"]).is_relative_to(run / "logs")
    if suite == "icepacks":
        if name:
            from tests.condor.icepacks import datasets
            selected = {item["entry"] for item in datasets(Path(name).parts[0])
                        if Path(item["entry"]).is_relative_to(name)}
            assert {item["entry"] for item in config["entries"]} == selected
        expected = {str(i) for i, item in enumerate(config["entries"]) if item["active"]}
        assert {line.split()[1][1:] for line in dag.splitlines() if line.startswith("JOB ")} == expected
        rows = json.loads((run / "report.json").read_text())["runs"]
        assert len(rows) == len(config["entries"])
        assert all(row["status"] == ("PENDING" if item["active"] else "SKIP")
                   for row, item in zip(rows, config["entries"], strict=True))
    else:
        categories = {"all": {"cpp", "python", "study"}, "cpp": {"cpp"},
                      "pytest": {"python"}, "studies": {"study"}}[suite]
        assert set(config["categories"]) == categories
        tasks = [json.loads(path.read_text()) for path in (run / "tasks").glob("*.json")]
        assert {task["category"] for task in tasks} == categories - {"cpp"}
        assert sum(config["array_limits"]) == config["max_jobs"]
        members = [member for task in tasks if task["category"] == "python" for member in task["members"]]
        assert len(members) == len(set(members)) == config["python_tests"]
        if "python" in categories:
            assert sorted(members) == sorted(json.loads((run / "manifests/python_nodes.json").read_text()))


# Discover each selected family without dropping nested or inactive dataset cards
@pytest.mark.parametrize("pack", [None, *ICEPACKS])
def test_icepack_discovery(pack):
    from tests.condor.icepacks import datasets
    entries = datasets(pack)
    assert entries and len({item["entry"] for item in entries}) == len(entries)
    root = Path(__file__).resolve().parents[3]
    for item in entries:
        card = pyjson5.loads((root / "icepack" / item["entry"] / "dataset.json").read_text())
        assert item["active"] == card["active"]
        assert item["entry"].split("/")[0] == pack if pack else item["measurement"]


# Reject missing, helper, and outside icepack folders
@pytest.mark.parametrize("name", ["UNKNOWN", "_common", "../tests", ".", "/",
                                  "UPC/PHOTOPROD/UNKNOWN", "GAMMA/continuum/_common"])
def test_icepack_invalid_folder(name):
    from tests.condor.icepacks import datasets
    with pytest.raises(ValueError, match="icepack"):
        datasets(name)


# Keep all expected runs in the table and propagate failed or missing jobs to the final status
def test_icepack_complete_results(tmp_path):
    from tests.condor.icepacks import summary
    from tests.condor.job_runner import write_json
    entries = [dict(entry=name, active=name != 'inactive') for name in ('passed', 'failed', 'error', 'running', 'missing', 'inactive')]
    write_json(tmp_path / 'config.json', dict(entries=entries))
    for index, status in enumerate(('PASS', 'FAIL', 'ERROR', 'RUNNING')):
        write_json(tmp_path / 'results' / f'{index}.json', dict(entry=entries[index]['entry'], status=status,
                   comparisons=1, seconds=1, failures=[] if status == 'PASS' else ['Comparison or run failed']))
    assert summary(tmp_path) == 1
    assert json.loads((tmp_path / 'report.json').read_text())['runs'][4]['status'] == 'PENDING'
    assert summary(tmp_path, final=True) == 1
    rows = json.loads((tmp_path / 'report.json').read_text())['runs']
    assert [row['status'] for row in rows] == ['PASS', 'FAIL', 'ERROR', 'MISSING', 'MISSING', 'SKIP']
    for index, item in enumerate(entries[:-1]):
        write_json(tmp_path / 'results' / f'{index}.json', dict(entry=item['entry'], status='PASS', failures=[]))
    assert summary(tmp_path, final=True) == 0
    subprocess.run(["bash", str(submit.ROOT / "tests/condor/submit.sh"), "icepacks",
                    "--report", str(tmp_path)], cwd=tmp_path, capture_output=True, check=True)
    assert json.loads((tmp_path / "report.json").read_text())["final"]


# Write one synthetic verifier result record
def write_result(
    run_dir: Path,
    category: str,
    index: str,
    item: str,
    members: list[str],
) -> None:
    path = run_dir / "results" / category / f"{index}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "category": category,
        "duration_seconds": 1.0,
        "exit_code": 0,
        "index": index,
        "item": item,
        "members": members,
        "outcome": "PASS",
    }
    path.write_text(json.dumps(payload) + "\n", encoding="ascii")


# Require longest first Python sharding to cover nodes once and balance load
def test_python_shards_cover_nodes_once_balance():
    nodes = ["test.py::a", "test.py::b", "test.py::c", "test.py::d"]
    timings = {
        "python::test.py::a": 60.0,
        "python::test.py::b": 30.0,
        "python::test.py::c": 20.0,
        "python::test.py::d": 10.0,
    }

    shards = submit.shard_python_tests(nodes, timings, 2)

    assert sorted(node for shard in shards for node in shard["nodes"]) == nodes
    assert [shard["estimated_seconds"] for shard in shards] == [60.0, 60.0]


# Keep publication generation and its aggregate table in one pytest process
@pytest.mark.parametrize("with_table", [False, True])
def test_integrated_icepacks_share_table_shard(with_table):
    from tests.physics.validation.test_icepacks import ICEPACK_ROOT, INTEGRATED_DATASETS

    prefix = "tests/physics/validation/test_icepacks.py::"
    datasets = [
        f"{prefix}test_icepack_event_flow[{path.parent.relative_to(ICEPACK_ROOT).as_posix()}]"
        for path in INTEGRATED_DATASETS
    ]
    aggregate = f"{prefix}test_integrated_xs_table"
    nodes = [*datasets, "test.py::independent"]
    if with_table:
        nodes.append(aggregate)
    shards = submit.shard_python_tests(nodes, {}, len(nodes))
    assert sorted(node for shard in shards for node in shard["nodes"]) == sorted(nodes)
    assert all(shard["nodes"] for shard in shards)
    if with_table:
        combined = next(shard["nodes"] for shard in shards if aggregate in shard["nodes"])
        assert combined == sorted(datasets) + [aggregate]
        assert len(shards) == 2
    else:
        assert all(len(shard["nodes"]) == 1 for shard in shards)


# Accept supported compiler parallelism and reject values outside the limits
def test_build_jobs_are_limited():
    assert submit.parse_args([]).build_jobs == submit.resource_settings()["jobs"]["build"]["request_cpus"]
    for suite in ("all", "cpp"):
        for jobs in (1, 8):
            assert submit.parse_args([suite, "--build-jobs", str(jobs)]).build_jobs == jobs
        for jobs in (0, 9):
            with pytest.raises(SystemExit) as error:
                submit.parse_args([suite, "--build-jobs", str(jobs)])
            assert error.value.code == 2


# Report build initialization failure to DAGMan as well as the final verifier
def test_build_failure_stops_dependents(tmp_path):
    root = Path(__file__).resolve().parents[3]
    (tmp_path / "config.json").write_text("{}\n", encoding="ascii")
    result = subprocess.run(
        [sys.executable, str(root / "tests/condor/job_runner.py"),
         "build", "000000", "build", str(tmp_path)],
        cwd=root, capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0
    status = json.loads((tmp_path / "results/build/000000.json").read_text())
    assert status["outcome"] == "FAIL"
    assert status["exit_code"] == result.returncode


# Require every canonical study launcher to be discovered recursively
def test_study_discovery_recursive_name_driven(tmp_path):
    first = tmp_path / "tests" / "physics" / "studies" / "family" / "first" / "run.sh"
    second = tmp_path / "tests" / "physics" / "studies" / "family" / "second" / "run.sh"
    ignored = (
        tmp_path / "tests" / "physics" / "studies" / "family" / "second" / "analyze.sh"
    )
    for path in (first, second, ignored):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    assert submit.discover_studies(tmp_path) == [
        "tests/physics/studies/family/first/run.sh",
        "tests/physics/studies/family/second/run.sh",
    ]


# Execute a local study through the real Condor workspace and subprocess helpers
def test_study_action_and_event_count(tmp_path, monkeypatch):
    source = tmp_path / "source"
    study = source / "tests" / "physics" / "studies" / "sample" / "run.sh"
    study.parent.mkdir(parents=True)
    for name in ("vgrid", "eikonal", "sudakov"):
        (source / name).mkdir()
    run_dir = source / "runs/tests/run"
    workspace = run_dir / "work/study/000001"
    study.write_text(
        "#!/usr/bin/env bash\nset -eu\n"
        'test "$NEVENTS" = 37\n'
        'test "$STUDY_ACTION" = both\n'
        'test -z "${EVENTS+x}"\n'
        'test -d "$GRANIITTI_TEST_OUTPUT_DIR"\n'
        'test -d "$MPLCONFIGDIR"\n'
        f'test "$PWD" = {shlex.quote(str(workspace))}\n',
        encoding="ascii",
    )
    (source / "runs").mkdir()
    monkeypatch.setattr(runner, "ROOT", source)
    for nevents, expected in ((37, 0), (36, 1)):
        status, _, _ = runner.run_study(
            "tests/physics/studies/sample/run.sh", "000001", {"study_events": nevents},
            run_dir, {**os.environ, "EVENTS": "stale"},
        )
        assert status == expected
        assert not (workspace / "runs").exists()


# Require Condor to copy studies and link the remaining test tree
def test_study_workspace_flat_physics_studies(tmp_path):
    source = tmp_path / "source"
    studies = source / "physics" / "studies"
    studies.mkdir(parents=True)
    (studies / "run.sh").write_text("#!/usr/bin/env bash\n", encoding="ascii")
    unit = source / "physics" / "unit"
    unit.mkdir()
    technical = source / "technical"
    technical.mkdir()
    destination = tmp_path / "destination"

    runner.link_study_tests(source, destination)

    copied = destination / "physics" / "studies"
    assert copied.is_dir() and not copied.is_symlink()
    assert (copied / "run.sh").is_file()
    assert (destination / "physics" / "unit").is_symlink()
    assert (destination / "technical").is_symlink()


# Require work arrays to honor one exact global materialization budget
def test_dag_array_limits_sum_worker_budget(tmp_path):
    dag_path = tmp_path / "test.dag"

    submit.write_dag(dag_path, tmp_path, [67, 67, 66])
    text = dag_path.read_text(encoding="ascii")

    assert text.count("JOB WORK_") == 3
    assert text.count("MAX_MATERIALIZE=") == 3
    assert 'MAX_MATERIALIZE="67"' in text
    assert 'MAX_MATERIALIZE="66"' in text
    assert "FINAL VERIFIER" in text
    assert "PARENT BUILD CHILD WORK_0 WORK_1 WORK_2" in text


# Require mixed task pools to balance estimated duration across task types
def test_work_manifests_balance_mixed_tasks(tmp_path):
    task_dir = tmp_path / "tasks"
    manifest_dir = tmp_path / "manifests"
    task_dir.mkdir()
    manifest_dir.mkdir()
    categories = ("cpp", "python", "study", "cpp", "python", "study")
    durations = (9.0, 8.0, 7.0, 6.0, 5.0, 4.0)
    for index, (category, duration) in enumerate(zip(categories, durations, strict=True)):
        item = f"{category}_{index}"
        payload = {
            "category": category,
            "estimated_seconds": duration,
            "item": item,
            "members": [item],
        }
        (task_dir / f"task.{index}.json").write_text(json.dumps(payload) + "\n", encoding="ascii")

    runner.write_work_manifests(tmp_path)

    pool_sizes = [
        len((manifest_dir / f"work_{index}.tsv").read_text(encoding="ascii").splitlines())
        for index in range(3)
    ]
    assert pool_sizes == [2, 2, 2]
    assert len((manifest_dir / "cpp.tsv").read_text(encoding="ascii").splitlines()) == 2
    assert len((manifest_dir / "python.tsv").read_text(encoding="ascii").splitlines()) == 2
    assert len((manifest_dir / "study.tsv").read_text(encoding="ascii").splitlines()) == 2


# Require the verifier to count concrete tests inside a successful shard
@pytest.mark.parametrize("categories", [["cpp"], ["python"], ["study"], ["cpp", "python", "study"]])
def test_verifier_counts_selected_tests(tmp_path, categories):
    build = "cpp" in categories
    runner.write_json(tmp_path / "config.json", {"build": build, "categories": categories})
    for directory in ("junit", "manifests", "results", "tasks"):
        (tmp_path / directory).mkdir()
    python_item = "tasks/python.000000.json"
    python_members = ["test_sample.py::test_a", "test_sample.py::test_b"]
    (tmp_path / python_item).write_text(
        json.dumps({"members": python_members}) + "\n", encoding="ascii"
    )
    manifests = {
        "cpp": "000001 cpp.test_sample\n",
        "python": f"000002 {python_item}\n",
        "study": "000003 tests/physics/studies/sample/run.sh\n",
    }
    for category, text in manifests.items():
        (tmp_path / "manifests" / f"{category}.tsv").write_text(text, encoding="ascii")
    write_result(tmp_path, "build", "000000", "build", ["build"])
    write_result(tmp_path, "cpp", "000001", "cpp.test_sample", ["cpp.test_sample"])
    write_result(tmp_path, "python", "000002", python_item, python_members)
    study = "tests/physics/studies/sample/run.sh"
    write_result(tmp_path, "study", "000003", study, [study])

    expected, manifest_errors = verifier.expected_results(tmp_path)
    suites, results, result_errors = verifier.collect_results(tmp_path, expected)

    assert manifest_errors == []
    assert result_errors == []
    assert len(results) == len(categories) + int(build)
    assert set(suites) == set(categories) | ({"build"} if build else set())
    if "python" in categories:
        assert suites["python"] == {"expected": 2, "failed": 0, "missing": 0, "passed": 2}

    category = categories[-1]
    index, _ = expected[category][0]
    path = tmp_path / "results" / category / f"{index}.json"
    path.rename(path.with_name(path.name + "._old"))
    suites, _, errors = verifier.collect_results(tmp_path, expected)
    assert suites[category]["missing"] == (2 if category == "python" else 1)
    assert errors


# Execute a real Python shard and propagate success or collection failure through the final verifier
@pytest.mark.parametrize("test", ["test_build_jobs_are_limited", "test_does_not_exist"])
def test_python_worker_and_verifier(test, tmp_path):
    node = f"tests/technical/aux/test_condor_submission.py::{test}"
    item = "tasks/python.000000.json"
    runner.write_json(tmp_path / "config.json", {"build": False, "categories": ["python"], "loopscreen": 1})
    runner.write_json(tmp_path / item, {"category": "python", "item": item, "members": [node]})
    (tmp_path / "manifests").mkdir()
    runner.write_manifest_records(tmp_path / "manifests/python.tsv", [("000000", item)])
    command = ["bash", str(submit.ROOT / "tests/condor/run_job.sh"), str(submit.ROOT)]
    environment = {**os.environ, "GRANIITTI_PYTEST_COLLECTION": str(tmp_path / "collection.json")}
    subprocess.run([*command, "job_runner", "work", "000000", item, str(tmp_path)],
                   cwd=tmp_path, env=environment, capture_output=True, text=True, check=True)
    record = json.loads((tmp_path / "results/python/000000.json").read_text())
    passed = test == "test_build_jobs_are_limited"
    assert record["outcome"] == ("PASS" if passed else "FAIL")
    assert (record["exit_code"] == 0) == passed
    assert record["members"] == [node]
    assert not (tmp_path / "collection.json").exists()
    timings = json.loads((tmp_path / "junit/python.000000.timings.json").read_text())
    assert set(timings) == ({node} if passed else set())
    assert all(seconds >= 0 for seconds in timings.values())
    verified = subprocess.run([*command, "verify", str(tmp_path), "0", "0"],
                              cwd=tmp_path, capture_output=True, text=True, check=False)
    assert verified.returncode == int(not passed)
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["outcome"] == record["outcome"]
    assert report["suites"]["python"]["expected"] == 1


# Reject suite-specific options and keep explicit Condor arguments separate from test options
def test_submission_options():
    args = submit.parse_args(["cpp", "--max-jobs", "8", "--", "-maxidle", "2"])
    assert args.suite == "cpp" and args.max_jobs == 8
    assert args.condor_args == ["-maxidle", "2"]
    assert submit.parse_args([]).output_root == submit.ROOT / "runs/tests"
    assert submit.parse_args(["icepacks"]).output_root == submit.ROOT / "figs/icepack"
    assert submit.parse_args(["icepacks", "--output-root", "custom"]).output_root == Path("custom")
    for options in (["cpp", "--local"], ["pytest", "--name", "UPC"], ["studies", "--report", "run"],
                    ["icepacks", "--local", "--validate"], ["cpp", "--unknown"]):
        with pytest.raises(SystemExit) as error:
            submit.parse_args(options)
        assert error.value.code == 2
