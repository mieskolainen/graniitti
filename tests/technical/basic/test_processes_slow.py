# Slow technical process-construction tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#
# Run with: pytest tests/technical/basic/test_processes_slow.py --run-integration -s -rP


import shlex

from core.io.steering import build_generation_plan, load_dataset

# ------------------------------------------------------------------------


def test_PP(physics_screening):
    """
    Test soft processes
    """
    proclist = ["TP[RES]", "TP[RES+CON]"]

    cmd = []
    for proc in proclist:
        for fs in ["pi+ pi-", "K+ K-"]:
            cmd.append(
                f"./bin/gr -i gencard/test.json -h 0 -n 100 -p '{proc}<F> -> {fs} @RES{{rho_770:1, f0_980:1, f2_1270:1}}' -o test_PP_{proc}"
            )
    physics_screening.execute(cmd)


# Exercise the unit-amplitude gluon flux in its configured hard-process phase space
def test_gg_flux(physics_screening):
    dataset, path = load_dataset("icepack/DURHAM/flux/dataset.json", cdir=".")
    sample = build_generation_plan(dataset, dataset_path=path, cdir=".", output_prefix="test_gg_flux")["samples"][0]
    overrides = " ".join(f"--set {shlex.quote(value)}" for value in sample["overrides"])
    physics_screening.execute(
        f"./bin/gr -i {shlex.quote(sample['gencard'])} -h 0 -n 1 -o {sample['output']} {overrides}"
    )


# Exercise Durham QCD partons, meson pairs and charmonium channels
def test_gg(physics_screening):
    """
    Test Durham QCD processes
    """
    cmd = []

    # Exercise every generated Durham MG5 flavor plus both generic jet aliases
    proclist = ["gg[QCD]<C> -> g g"]
    proclist.extend(f"gg[QCD]<C> -> {flavor} {flavor}~" for flavor in ("u", "d", "s", "c", "b"))
    proclist.append("gg[QCD]<C> -> j j~")
    for index, proc in enumerate(proclist):
        cmd.append(
            "./bin/gr -i icepack/DURHAM/_common/gencard.json "
            f"-h 0 -n 1 -p '{proc}' -o test_gg_{index}"
        )
        if proc == "gg[QCD]<C> -> j j~":
            cmd[-1] += " --set 'INTEGRATOR.precision=0.05'"

    proclist = [f"gg[QCD]<F> -> {flavor} {flavor}~ g" for flavor in ("u", "d", "s", "c", "b")]
    proclist.extend(("gg[QCD]<F> -> j j~ g", "gg[QCD]<F> -> g g g"))
    for index, proc in enumerate(proclist):
        cmd.append(
            "./bin/gr -i icepack/DURHAM/_common/gencard.json "
            f"-h 0 -n 1 -p '{proc}' -o test_gg_threebody_{index}"
        )

    proclist = ["gg[MM]<F> -> pi+ pi-", "gg[MM]<F> -> K+ K-"]
    for proc in proclist:
        cmd.append(
            f"./bin/gr -i icepack/DURHAM/_common/gencard.json -h 0 -n 100 -p '{proc}' -o test_MMbar"
        )

    proclist = ["gg[chic(0)]<F> &> pi+ pi-", "gg[chic(0)]<F> &> K+ K-"]
    for proc in proclist:
        cmd.append(
            f"./bin/gr -i icepack/DURHAM/_common/gencard.json -h 0 -n 100 -p '{proc}' -o test_chic0"
        )
    proclist = [
        (
            "icepack/DURHAM/_common/gencard.json",
            "gg[chic(1)]<F> &> phi(1020)0 phi(1020)0",
            "test_chic1",
        ),
        (
            "icepack/DURHAM/_common/gencard.json",
            "gg[chic(2)]<F> &> phi(1020)0 phi(1020)0",
            "test_chic2",
        ),
    ]
    for card, proc, output in proclist:
        cmd.append(f"./bin/gr -i {card} -h 0 -n 100 -p '{proc}' -o {output}")

    physics_screening.execute(cmd)


def test_mg5_process_families(physics_screening):
    """
    Test each automated multi-subprocess MG5 family through the generator
    """
    processes = (
        (
            "icepack/HARDPOM/z/gencard_yy_ZJJ_mumu_uu.json",
            "yy[Zjj]<F> -> Z > {mu+ mu-} u u~",
            "yy_zjj",
        ),
        (
            "icepack/HARDPOM/z/gencard_IPIP_Z_mumu.json",
            "IPIP[Z]<F> -> mu+ mu-",
            "ipip_z",
        ),
        (
            "icepack/HARDPOM/z/gencard_IPp_ZJ_mumu_u.json",
            "IPp[Zj]<F> -> Z > {mu+ mu-} u",
            "ipp_zj",
        ),
        (
            "icepack/HARDPOM/jj/gencard_IPIP_jj.json",
            "IPIP[jj]<F> -> u u~",
            "ipip_jj",
        ),
        (
            "icepack/HARDPOM/jj/gencard_IPIP_jj.json",
            "IPp[W]<F> -> mu+ vm",
            "ipp_w",
        ),
        (
            "icepack/HARDPOM/jj/gencard_IPIP_jj.json",
            "yy[jj]<F> -> u u~",
            "yy_jj",
        ),
        (
            "icepack/HARDPOM/jj/gencard_IPIP_jj.json",
            "yy[WW]<F> -> W+ > {mu+ vm} W- > {e- ve~}",
            "yy_ww_mue",
        ),
    )
    cmd = []
    for card, process, tag in processes:
        cmd.append(
            f"./bin/gr -i {card} -h 0 -n 1 -p '{process}' "
            f"-o test_mg5_family_{tag} "
        )
        # Generate the W region where the proton PDF and DPDF share five-flavour evolution
        if tag == "ipp_w":
            cmd[-1] += "--set 'GENCUTS.<F>.M=[60,120]'"
    physics_screening.execute(cmd)


# Validate the globally selected Pomeron-loop mode
def test_loopscreen(physics_screening):
    command = "./bin/gr -i gencard/test.json -h 0 -n 0 -p 'GP[CON]<C> -> pi+ pi-'"
    physics_screening.execute(command)


def test_lhapdf(physics_screening):
    """
    Test setting lhapdf
    """
    cmd = []
    for pdf in ["CT10nlo", "MMHT2014lo68cl"]:
        cmd.append(
            f"./bin/gr -i icepack/DURHAM/_common/gencard.json -h 0 -n 0 -p 'gg[QCD]<C> -> g g' -o gg2gg -f hepmc3 -q {pdf}"
        )
    physics_screening.execute(cmd)
