# Fast technical process-construction tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

# Run with: pytest tests/technical/basic/test_processes_fast.py --run-integration -s -rP


import json
import math
import re
import subprocess
from collections import Counter
from pathlib import Path

from tests.technical.support.hepmc import read_events
from tests.technical.support.iceplot import run_iceplot
from tests.technical.support.output import PhysicsOutput

OUTPUT = PhysicsOutput("processes_fast")


# Read event weights from fully reconstructed HepMC3 events
def read_hepmc3_weights(path):
    return [event.weights()[0] for event in read_events(path)]


# Extract particle momenta, parent identities and color flow from HepMC3
def read_hepmc3_particles(path):
    return [{
        'id': p.id(), 'parent': p.parents()[0].id() if len(p.parents()) == 1 else 0,
        'pdg': p.pid(), 'status': p.status(), 'mass': p.generated_mass(),
        'px': p.momentum().px(), 'py': p.momentum().py(), 'pz': p.momentum().pz(),
        'energy': p.momentum().e(),
        'flow1': int(p.attribute_as_string('flow1') or 0),
        'flow2': int(p.attribute_as_string('flow2') or 0),
    } for event in read_events(path) for p in event.particles()]


# Read the incoming hard-process PDG ids from one HepMC3 event
def read_hepmc3_pdf_ids(path):
    return tuple(map(int, read_events(path)[0].attribute_as_string("GenPdfInfo").split()[:2]))


# Compute direct stable partons from the central-system vertex
def read_hepmc3_central_partons(particles):
    system_ids = {particle["id"] for particle in particles if particle["pdg"] == 90}
    return [
        particle
        for particle in particles
        if particle["status"] == 1
        and particle["parent"] in system_ids
        and (1 <= abs(particle["pdg"]) <= 5 or particle["pdg"] == 21)
    ]


# Check that every stored color line closes on two stable particles
def assert_hepmc3_color_closure(particles):
    tags = [
        particle[key]
        for particle in particles
        if particle["status"] == 1
        for key in ("flow1", "flow2")
        if particle[key] != 0
    ]
    assert tags
    assert all(count == 2 for count in Counter(tags).values())


# Count fully decoded HepMC3 events
def count_hepmc3_events(path):
    return len(read_events(path))


# ------------------------------------------------------------------------
# INPUT-OUTPUT


def test_input_output(physics_screening):
    """
    Test input-output procedures
    """

    cmd = []
    cmd.append("./bin/gr -i gencard/test.json -h 0 -n 0 -p 'GP[CON]<F> -> pi+ pi-' -o GP_CON")
    cmd.append(
        "./bin/gr -i gencard/test.json -h 0 -n 0 -p 'GP[RES+CON]<F> -> pi+ pi-' -o GP_RES_CON"
    )
    cmd.append(
        "./bin/gr -i gencard/test.json -d ./vgrid/GP_CON.vgrid -h 0 -n 1000 -p 'GP[CON]<F> -> pi+ pi-' -o GP_CON -f hepmc3"
    )
    cmd.append(
        "./bin/gr -i gencard/test.json -d ./vgrid/GP_RES_CON.vgrid -h 0 -n 1000 -p 'GP[RES+CON]<F> -> pi+ pi-' -o GP_RES_CON -f hepmc3"
    )
    physics_screening.execute(cmd)

    run_iceplot(
        hepmc3_tags=["GP_CON", "GP_RES_CON"],
        plot_tag="processes_fast_input_output",
        labels=["continuum", "resonance plus continuum"],
        pid=[[211, -211], [211, -211]],
        output=OUTPUT,
    )


# ------------------------------------------------------------------------
# GENERAL


def test_iceplot(physics_screening):
    """
    Test iceplot
    """

    physics_screening.execute("./bin/gr -i icepack/SOFTCEP/STAR_1792394/pipi/gencard.json -n 10000")
    physics_screening.execute("./bin/gr -i icepack/SOFTCEP/CMS_2752118/pipi/gencard.json -n 10000")

    run_iceplot(
        hepmc3_tags=["STAR_1792394_pipi"],
        plot_tag="processes_fast_star",
        labels=["GRANIITTI"],
        analysis="icepack/SOFTCEP/STAR_1792394/pipi",
        output=OUTPUT,
    )
    run_iceplot(
        hepmc3_tags=["CMS_2752118_pipi"],
        plot_tag="processes_fast_cms",
        labels=["GRANIITTI"],
        analysis="icepack/SOFTCEP/CMS_2752118/pipi",
        output=OUTPUT,
    )


# Validate generated and reloaded integration grids
def test_vgrid(physics_screening):
    """
    Test reading pre-computed MC grid
    """
    cmd = []

    command = "./bin/gr -i gencard/test.json -h 0 -n 0 -p 'GP[CON]<F> -> pi+ pi-' -o vgrid_test"
    cmd.append(f"{command}")
    cmd.append(f"{command} -d ./vgrid/vgrid_test.vgrid")
    physics_screening.execute(cmd)
    grid_path = Path("vgrid/vgrid_test.vgrid")
    assert grid_path.is_file()
    with grid_path.open(encoding="utf-8") as stream:
        sigma = float(json.load(stream)["STAT"]["sigma"])
    assert math.isfinite(sigma)
    assert sigma > 0.0


# Validate requested output names
def test_output(physics_screening):
    """
    Test output name
    """
    cmd = []
    outputs = ["output_A", "output_B"]
    for output in outputs:
        cmd.append(
            f"./bin/gr -i gencard/test.json -h 0 -n 0 -p 'GP[CON]<F> -> pi+ pi-' -o {output}"
        )
    physics_screening.execute(cmd)
    for output in outputs:
        assert (Path("vgrid") / f"{output}.vgrid").is_file()


# Validate each supported event output format
def test_format(physics_screening):
    """
    Test output formats
    """
    cmd = []
    formats = ["hepmc2", "hepmc3", "hepevt"]
    for output_format in formats:
        cmd.append(
            f"./bin/gr -i gencard/test.json -h 0 -n 2 -p 'GP[CON]<F> -> pi+ pi-' -f {output_format} -o format_test_{output_format}"
        )
    physics_screening.execute(cmd)
    for output_format in formats:
        output_path = Path("output") / f"format_test_{output_format}.{output_format}"
        assert output_path.is_file()
        assert output_path.stat().st_size > 0


# Compare finite physical rates from independent single and multithreaded VEGAS runs
def test_cores(physics_screening):
    rates = []
    for cores in [1, 2, 4]:
        physics_screening.execute(
            f"./bin/gr -i gencard/test.json -g VEGAS -h 0 -n 0 -p 'GP[CON]<F> -> pi+ pi-' "
            f"-c {cores} -r {27183 + 7919 * cores} -o cores_test_{cores}"
        )
        stat = json.loads(Path(f"vgrid/cores_test_{cores}.vgrid").read_text())["STAT"]
        rate, error = stat["sigma"], stat["sigma_err"]
        assert math.isfinite(rate) and rate > 0.0
        assert math.isfinite(error) and 0.0 < error <= 0.01 * rate
        rates.append((rate, error))
    for rate, error in rates[1:]:
        assert abs(rate - rates[0][0]) <= 5.0 * math.hypot(error, rates[0][1])


# Validate the requested generated event count
def test_nevents(physics_screening, physics_events):
    """
    Test number of events
    """
    events = physics_events.select(100)
    output = f"nevents_test_{events}_{physics_screening.suffix}"
    command = (
        "./bin/gr -i gencard/test.json -h 0 -f hepmc3 "
        f"-p 'GP[CON]<F> -> pi+ pi-' -n {events} -o {output}"
    )
    physics_screening.execute(command)
    assert count_hepmc3_events(Path("output") / f"{output}.hepmc3") == events


# Validate the globally selected HepMC3 weighting mode
def test_weighted(physics_screening, physics_events):
    """
    Test the selected event generation mode
    """
    events = physics_events.select(100)
    output = f"weighted_test_{events}_{physics_screening.suffix}"
    command = (
        "./bin/gr -i gencard/test.json -h 0 -f hepmc3 "
        f"-p 'GP[CON]<F> -> pi+ pi-' -n {events} -o {output}"
    )
    physics_screening.execute(command)
    weights = read_hepmc3_weights(Path("output") / f"{output}.hepmc3")
    assert len(weights) == events
    if physics_screening.weighted:
        assert all(math.isfinite(value) and value > 0.0 for value in weights)
        if events > 1:
            assert len(set(weights)) > 1
    else:
        assert weights == [1.0] * events


# Validate the explicit model-tune override
def test_modelparam(physics_screening):
    """
    Test changing modeltune
    """
    physics_screening.execute(
        "./bin/gr -i gencard/test.json -m TUNE0 -h 0 -n 0 -p 'GP[CON]<F> -> pi+ pi-' -o modelparam_test_TUNE0"
    )
    assert Path("vgrid/modelparam_test_TUNE0.vgrid").is_file()


# ---------------------------------------------------------
# PROCESS


def test_process(physics_screening):
    cmd = []
    for proc in [
        "MP[RES+CON]<F> -> pi+ pi-",
        "MP[CON]<C> -> K+ K-",
        "GP[CON]<F> -> pi+ pi-",
        "GP[RES+CON]<F> -> pi+ pi- @RES{f2_1270:1}",
        "MP[CON]<C> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}",
        "MP[RES+CON]<F> -> pi+ pi- @RES{f0_500:0,rho_770:1,f0_980:1,f2_1270:1} @R[f0_980]{M:0.98,W:0.065}",
    ]:
        cmd.append(f"./bin/gr -i gencard/test.json -h 0 -n 0 -p '{proc}' -o process_test")
    physics_screening.execute(cmd)


# Generate physical on-shell partons through each generic selector process family
def test_parton_mass_shell(physics_screening):
    # Read the bottom mass from the fixed columns of the generator's PDG table
    pdg_rows = Path("modeldata/mass_width_2026.mcd").read_text(encoding="utf-8").splitlines()
    bottom_mass = next(float(row[33:51]) for row in pdg_rows if row[:8].strip() == "5")
    common = (
        "./bin/gr -w true -h 0 -n 1 -c 1 -f hepmc3 -g VEGAS "
        "--set 'INTEGRATOR.min_samples=100' "
        "--set 'INTEGRATOR.max_samples=5000' "
        "--set 'INTEGRATOR.precision=1.0' "
        "--set 'INTEGRATOR.VEGAS.ncall=1000' "
        "--set 'INTEGRATOR.VEGAS.rounds=2' "
        "--set 'NUMERICS.json:NUMERICS_VEGAS.automatic_convergence=false' "
        "--set 'FIDCUTS.active=false'"
    )
    gluon_pomeron = "--set 'GENERAL.json:PARAM_HARDPOMERON.parton_flavours=[21]'"
    cases = (
        {
            "card": "icepack/DURHAM/_common/gencard.json",
            "process": "gg[QCD]<F> -> j ~j @j={b}",
            "tag": "generic_parton_durham_b",
            "extra": "-q MMHT2014lo68cl",
            "count": 2,
            "pdgs": {5},
            "mass": 4.7,
            "pdf_ids": (21, 21),
        },
        {
            "card": "icepack/DURHAM/_common/gencard.json",
            "process": "gg[QCD]<F> -> j ~j @j={g}",
            "tag": "generic_parton_durham_g",
            "extra": "-q MMHT2014lo68cl",
            "count": 2,
            "pdgs": {21},
            "mass": 0.0,
            "pdf_ids": (21, 21),
        },
        {
            "card": "icepack/GAMMA/mumu_with_pythia/gencard.json",
            "process": "yy[EPA]<F> -> j ~j @j={b,g}",
            "tag": "generic_parton_epa",
            "extra": "",
            "count": 2,
            "pdgs": {5},
            "mass": bottom_mass,
            "pdf_ids": (22, 22),
        },
        {
            "card": "icepack/HARDPOM/jj/gencard_IPIP_jj.json",
            "process": "yy[jj]<F> -> j ~j @j={u}",
            "tag": "generic_parton_yy_jj_u",
            "extra": "",
            "count": 2,
            "pdgs": {2},
            "mass": 0.0,
            "pdf_ids": (22, 22),
        },
        {
            "card": "icepack/HARDPOM/z/gencard_IPp_ZJ_mumu_u.json",
            "process": "IPp[Zj]<F> -> Z > {mu+ mu-} j @j={u}",
            "tag": "generic_parton_ipp",
            "extra": "",
            "count": 1,
            "pdgs": {2},
            "mass": 0.0,
            "pdf_ids": None,
        },
        {
            "card": "icepack/HARDPOM/z/gencard_yy_ZJJ_mumu_uu.json",
            "process": "yy[Zjj]<F> -> Z > {mu+ mu-} j ~j @j={u}",
            "tag": "generic_parton_yy_zjj_u",
            "extra": "",
            "count": 2,
            "pdgs": {2},
            "mass": 0.0,
            "pdf_ids": (22, 22),
        },
        {
            "card": "icepack/HARDPOM/z/gencard_yy_ZJJ_mumu_uu.json",
            "process": "yy[Zjj]<F> -> Z > {mu+ mu-} j ~j @j={u,c,g}",
            "tag": "generic_parton_yy_zjj_ucg",
            "extra": "",
            "count": 2,
            "pdgs": {2, 4},
            "mass": 0.0,
            "pdf_ids": (22, 22),
        },
        {
            "card": "icepack/HARDPOM/jj/gencard_IPIP_jj.json",
            "process": "IPIP[jj]<F> -> j ~j @j={u}",
            "tag": "generic_parton_ipip_u",
            "extra": f"{gluon_pomeron}",
            "count": 2,
            "pdgs": {2},
            "mass": 0.0,
            "pdf_ids": (21, 21),
        },
        {
            "card": "icepack/HARDPOM/jj/gencard_IPIP_jj.json",
            "process": "IPIP[jj]<F> -> j ~j @j={u,c}",
            "tag": "generic_parton_ipip_uc",
            "extra": f"{gluon_pomeron}",
            "count": 2,
            "pdgs": {2, 4},
            "mass": 0.0,
            "pdf_ids": (21, 21),
        },
    )
    commands = [
        f"{common} -i {case['card']} -p '{case['process']}' -o {case['tag']} {case['extra']}"
        for case in cases
    ]
    outputs = physics_screening.execute(commands)

    for output, case in zip(outputs, cases, strict=True):
        match = re.search(r"Fiducial cross section:\s+\[([0-9.Ee+-]+)", output)
        assert match is not None
        assert math.isfinite(float(match.group(1)))
        assert float(match.group(1)) > 0.0

        path = Path("output") / f"{case['tag']}.hepmc3"
        weights = read_hepmc3_weights(path)
        assert len(weights) == 1
        assert math.isfinite(weights[0]) and weights[0] > 0.0

        particles = read_hepmc3_particles(path)
        assert all(abs(particle["pdg"]) != 89 for particle in particles)
        partons = read_hepmc3_central_partons(particles)
        assert len(partons) == case["count"]
        assert all(abs(particle["pdg"]) in case["pdgs"] for particle in partons)
        assert all(particle["flow1"] != 0 or particle["flow2"] != 0 for particle in partons)
        assert_hepmc3_color_closure(particles)
        if case["pdf_ids"] is not None:
            assert read_hepmc3_pdf_ids(path) == case["pdf_ids"]

        for particle in partons:
            momentum2 = particle["px"] ** 2 + particle["py"] ** 2 + particle["pz"] ** 2
            mass2 = particle["energy"] ** 2 - momentum2
            scale = max(particle["energy"] ** 2, momentum2, case["mass"] ** 2, 1.0)
            assert math.isclose(particle["mass"], case["mass"], rel_tol=1e-9, abs_tol=1e-4)
            assert abs(mass2 - case["mass"] ** 2) <= 1e-10 * scale


def test_cascade_c_phase_space_processes(physics_screening):
    """
    Test <C> cascaded vector-pair phase space for continuum models
    """

    cascade = "rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}"
    common = (
        "./bin/gr -i gencard/test.json -h 0 -n 0 -g VEGAS -c 1 "
        "--set 'INTEGRATOR.min_samples=1000' --set 'INTEGRATOR.max_samples=10000' "
        "--set 'INTEGRATOR.VEGAS.ncall=1000' --set 'INTEGRATOR.VEGAS.rounds=2' "
        "--set 'NUMERICS.json:NUMERICS_VEGAS.automatic_convergence=false' "
        "--set 'INTEGRATOR.precision=1.0'"
    )
    cases = [
        ("MP", "CON", "false"),
        ("MP", "CON", "true"),
        ("XP", "CON", "false"),
        ("XP", "CON", "true"),
        ("GP", "CON", "false"),
        ("GP", "CON", "true"),
        ("TP", "CON", "false"),
        ("TP", "CON", "true"),
    ]

    cmd = []
    for family, channel, decay_sym in cases:
        output = f"cascade_c_{family.lower()}_{channel.lower()}_sym_{decay_sym}"
        cmd.append(
            f"{common} -p '{family}[{channel}]<C> -> {cascade}' -o {output} "
            f"--set 'GENERAL.json:PARAM_SPIN.DECAY_SYM={decay_sym}'"
        )
    physics_screening.execute(cmd)


def test_mcon_cascade_accepts_decay_sym_switch(physics_screening):
    """
    Test that MP vector-cascade generation accepts both DECAY_SYM modes
    for both <F> and <C> phase-space mappings
    """

    cascade = "rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}"

    for phase_space in ("F", "C"):
        common = (
            "./bin/gr -i gencard/test.json -h 0 -n 0 -g VEGAS -c 1 "
            f"-p 'MP[CON]<{phase_space}> -> {cascade}' "
            "--set 'INTEGRATOR.min_samples=1000' --set 'INTEGRATOR.max_samples=10000' "
            "--set 'INTEGRATOR.VEGAS.ncall=1000' --set 'INTEGRATOR.VEGAS.rounds=2' "
            "--set 'NUMERICS.json:NUMERICS_VEGAS.automatic_convergence=false' "
            "--set 'INTEGRATOR.precision=1.0'"
        )

        for decay_sym in ("false", "true"):
            output = subprocess.check_output(
                physics_screening.shell(
                    f"{common} -o cascade_{phase_space.lower()}_con_compare_{decay_sym} "
                    f"--set 'GENERAL.json:PARAM_SPIN.DECAY_SYM={decay_sym}'"
                ),
                shell=True,
                text=True,
                stderr=subprocess.STDOUT,
            )
            match = re.search(r"Fiducial cross section:\s+\[([0-9.Ee+-]+)", output)
            assert match is not None
            assert float(match.group(1)) > 0.0


def test_energy(physics_screening):
    cmd = []
    for energy in ["7000", "13000"]:
        cmd.append(
            f"./bin/gr -i gencard/test.json -h 0 -n 0 -p 'GP[CON]<C> -> pi+ pi-' -e {energy} -o energy_test_{energy}"
        )
    physics_screening.execute(cmd)


def test_nstars(physics_screening):
    """
    Test forward excitation (N*)
    """
    cmd = []
    for excite in [0, 1, 2]:
        cmd.append(
            f"./bin/gr -i gencard/test.json -h 0 -n 0 -p 'GP[CON]<F> -> pi+ pi-' -s {excite} -o excite_test_{excite}"
        )
    physics_screening.execute(cmd)


# Validate embedded histogram production
def test_hist(physics_screening):
    """
    Test fast embedded histograms
    """
    output = f"hist_test_{physics_screening.suffix}"
    command = f"./bin/gr -i gencard/test.json -h 1 -n 0 -p 'GP[CON]<C> -> pi+ pi-' -o {output}"
    physics_screening.execute(command)
    histogram_path = Path("output") / f"{output}.hfast"
    assert histogram_path.is_file()
    with histogram_path.open(encoding="utf-8") as stream:
        histograms = json.load(stream)["h1"]
    assert histograms
    filled_histograms = 0
    for histogram in histograms.values():
        weights = histogram["weights"]
        filled_histograms += int(histogram["fills"]) > 0
        assert all(math.isfinite(float(value)) for value in weights)
    assert filled_histograms > 0


# Validate random-seed reproducibility and separation
def test_rndseed(physics_screening):
    """
    Test random seeding
    """
    cases = [("123_a", 123), ("123_b", 123), ("456", 456)]
    cmd = []
    for tag, seed in cases:
        cmd.append(
            f"./bin/gr -i gencard/test.json -h 0 -n 100 -c 1 -f hepmc3 -p 'GP[CON]<C> -> pi+ pi-' -r {seed} -o rndseed_test_{tag}"
        )
    physics_screening.execute(cmd)
    weights = {
        tag: read_hepmc3_weights(Path("output") / f"rndseed_test_{tag}.hepmc3") for tag, _ in cases
    }
    assert weights["123_a"] == weights["123_b"]
    assert weights["123_a"] != weights["456"]


# ------------------------------------------------------------------------
# On-the-flight parameter syntax


# Exercise valid combinations of process steering and parameter overrides
def test_otf(physics_screening):
    cases = [
        ("MP[CON]", "@FLATAMP:1"),
        ("MP[CON]", "@FLATMASS2:true"),
        ("MP[CON]", "@OFFSHELL:7"),
        ("MP[CON]", "@PDG[995]{M:350.0, W:5.0}"),
        ("MP[RES]", "@RES{f0_980:1} @R[f0_980]{M:0.990, W:0.065}"),
        ("MP[RES]", "@RES{f0_980:1, f2_1270:1}"),
        ("MP[RES]", "@RES{f2_1270:1} @SPINGEN:true"),
        ("MP[RES]", "@RES{f2_1270:1} @SPINDEC:true"),
        ("MP[RES]", "@RES{f2_1270:1} @MP_FRAME:CS"),
        ("MP[RES]", "@RES{f2_1270:1} @MP_FRAME:HX"),
        ("MP[RES]", "@RES{f2_1270:1} @MP_FRAME:CM @R[f2_1270]{JZ0:1.0, JZ1:0.0, JZ2:0.0}"),
        ("MP[RES]", "@RES{f2_1270:1} @R[f2_1270]{JZ0:0.5, JZ1:0.0, JZ2:0.5}"),
        ("GP[RES]", "@RES{f2_1270:1} @MMAX:4"),
        ("TP[RES]", "@RES{f0_980:1} @R[f0_980]{g0:1.0, g1:0.2}"),
        ("TP[RES]", "@RES{f2_1270:1} @R[f2_1270]{g0:1.0, g1:0.1, g2:0.4, g3:0.3, g4:0.2, g5:0.1, g6:0.0}"),
    ]
    for index, (process, override) in enumerate(cases):
        output = f"otf_{index}_{physics_screening.suffix}"
        physics_screening.execute(
            f"./bin/gr -i gencard/test.json -h 0 -n 100 -f hepmc3 "
            f"-p '{process}<F> -> pi+ pi- {override}' -o {output}"
        )
        weights = read_hepmc3_weights(Path("output") / f"{output}.hepmc3")
        assert len(weights) == (physics_screening.event_override or 100)
        assert all(math.isfinite(weight) and weight > 0.0 for weight in weights)


# ------------------------------------------------------------------------
# PROCESSES


def test_PP(physics_screening):
    """
    Test soft processes
    """
    proclist = ["MP[CON]", "MP[RES+CON]", "GP[RES]", "MP[RES]"]

    cmd = []
    for proc in proclist:
        for fs in ["pi+ pi-", "K+ K-"]:
            cmd.append(
                f"./bin/gr -i gencard/test.json -h 0 -n 100 -p '{proc}<F> -> {fs}' -o test_{proc}"
            )
    cmd.append(
        "./bin/gr -i gencard/test.json -h 0 -n 100 "
        "-p 'MP[RES]<F> -> pi+ pi- @RES{rho_770_odd:1}' -o test_MP_RES_odd_rho"
    )
    cmd.append(
        "./bin/gr -i gencard/test.json -h 0 -n 100 "
        "-p 'MP[RES]<F> -> K+ K- @RES{phi_1020_odd:1}' -o test_MP_RES_odd_phi"
    )
    physics_screening.execute(cmd)


def test_photoproduction(physics_screening):
    """
    # Test gamma-Pomeron channels selected inside PP resonance processes
    """
    proclist = ["TP[RES]", "MP[RES]"]
    RES = "@RES{rho_770:1}"

    cmd = []
    for proc in proclist:
        for fs in ["pi+ pi-"]:
            cmd.append(
                f"./bin/gr -i gencard/test.json -h 0 -n 100 -p '{proc}<F> -> {fs} {RES}' -o test_photoproduction"
            )
    physics_screening.execute(cmd)


# Keep all allowed LS rows, with only the selected L=0, S=2 operator active
def test_regge_res_gamma_gamma_channels(physics_screening):
    """
    Test explicit gamma-gamma channels inside PP resonance processes
    """
    gp_block = {
        "basis": "g_ls",
        "Lambda": 1.0,
        "g_ls": [[2, 0, 0.0, 0.0], [0, 2, None, 0.0], [2, 2, 0.0, 0.0], [4, 2, 0.0, 0.0]],
        "CP": [True, True],
    }
    xp_block = {
        "basis": "g_ls",
        "Lambda": 1.0,
        "g_ls": [[2, 0, 0.0, 0.0], [0, 2, None, 0.0], [2, 2, 0.0, 0.0], [4, 2, 0.0, 0.0]],
        "CP": [True, True],
    }
    res_block = {
        "basis": "auto_min_L",
        "Lambda": 1.0,
        "g": [None, 0.0],
        "polarization": {"mode": "none"},
        "CP": [True, True],
    }
    cases = [
        ("MP", 'PARAM_RES.MODELS.MP["[22,22]"]', res_block),
        (
            "XP",
            'PARAM_RES.MODELS.XP["[22,22]"]',
            xp_block,
        ),
        (
            "GP",
            'PARAM_RES.MODELS.GP["[22,22]"]',
            gp_block,
        ),
    ]

    cmd = []
    for model, block_path, block in cases:
        block_json = json.dumps(block, separators=(",", ":"))
        cmd.append(
            "./bin/gr -i gencard/test.json -h 0 -n 20 "
            f"-p '{model}[RES]<F> -> pi+ pi- @RES{{f2_1270_yy:1}}' "
            f"-o test_pp_{model.lower()}_gamma_gamma "
            f"--set 'f2_1270_yy.json:{block_path}={block_json}'"
        )
    physics_screening.execute(cmd)
    for model, _, _ in cases:
        path = Path("output") / f"test_pp_{model.lower()}_gamma_gamma.hepmc3"
        weights = read_hepmc3_weights(path)
        assert len(weights) == (physics_screening.event_override or 20)
        assert all(math.isfinite(weight) and weight > 0.0 for weight in weights)


def test_yy_res_is_removed(physics_screening):
    """
    Test that the removed yy[RES] process fails loudly
    """
    command = (
        "./bin/gr -i gencard/test.json -h 0 -n 0 "
        "-p 'yy[RES]<F> -> pi+ pi- @RES{f0_980:1}' -o test_yy_res_removed"
    )
    output = subprocess.run(
        physics_screening.shell(command),
        shell=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    assert output.returncode != 0


def test_yy_DZ(physics_screening):
    """
    Test collinear gamma-gamma process routes
    """
    cmd = []
    for initial, card in [
        ("yy_DZ", "gencard_dz.json"),
        ("yy_LUX", "gencard_lux.json"),
    ]:
        for fs in ["e+ e-", "mu+ mu-", "tau+ tau-"]:
            cmd.append(
                "./bin/gr -i "
                f"icepack/GAMMA/continuum/_common/{card} -h 0 -n 100 "
                f"-p '{initial}[EPA]<P> -> {fs}' -o test_{initial}_EPA"
            )
        for process, final_state in [
            ("jj", "u u~"),
            ("WW", "W+ > {mu+ vm} W- > {e- ve~}"),
            ("Zjj", "Z > {mu+ mu-} u u~"),
        ]:
            cmd.append(
                "./bin/gr -i "
                f"icepack/GAMMA/continuum/_common/{card} -h 0 -n 0 "
                f"-p '{initial}[{process}]<P> -> {final_state}' "
                f"-o test_{initial}_{process}"
            )
    cmd.append(
        "./bin/gr -i icepack/GAMMA/continuum/_common/gencard_dz.json "
        "-h 0 -n 100 -p 'yy_DZ[FLUX]<P> -> 22 22' -o test_yy_DZ_FLUX"
    )
    physics_screening.execute(cmd)


def test_yy(physics_screening):
    """
    Test gamma-gamma processes
    """
    proclist = ["yy[EPA]", "yy[FLUX]", "yy[QED]"]

    # Higgs
    physics_screening.execute(
        [
            "./bin/gr -i icepack/GAMMA/resonances/higgs/gencard.json -h 0 -n 100 -p 'yy[Higgs]<F> &> b b~' -o test_Higgs_bbar"
        ]
    )

    # Others
    cmd = []
    for proc in proclist:
        for fs in ["e+ e-", "mu+ mu-"]:
            cmd.append(
                f"./bin/gr -i gencard/test.json -h 0 -n 100 -p '{proc}<F> -> {fs}' -o test_yy"
            )
    physics_screening.execute(cmd)
