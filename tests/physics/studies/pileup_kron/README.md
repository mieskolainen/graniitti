# Pileup jet comparison

Compare `Jet`, `JetPUPPI` and `JetGraphKron` reconstruction with Delphes. Run from the repository root with Delphes installed:

```bash
conda activate graniitti
source install/setenv.sh
RUN_DELPHES=1 DELPHES_DIR=/path/to/Delphes NEVENTS=300 PU=50 bash tests/physics/studies/pileup_kron/run.sh
```

Events, detector simulation and jet comparisons are saved under `output/graph_pu50/` by default. `OUTPUT_DIR` changes the directory. `INCLUDE_GRAPH_KRON=0` omits the graph comparison. See [run_delphes.sh](run_delphes.sh) for steering variables and `python tests/physics/studies/pileup_kron/analyze_delphes.py --help` for analysis options.

`GenHardJet` is the hard scattering reference, clustered from all visible stable particles with `Particle.IsPU == 0` before applying jet axis acceptance. `JetGraphKron` uses EFlow candidates and vertices independently of the other reconstructed jets.

Changes to generation settings, cards or driver sources cause the detector simulation to be repeated. Previous outputs are retained as backups.
