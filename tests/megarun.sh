#!/bin/sh
#
# GLOBAL benchmark test of the program
#
# Run with: source ./tests/megarun.sh > megarun.out


make -j4 TEST=TRUE

read -p "megarun: Generate events (or only analyze)? [y/n] " -n 1 -r
echo # New line


if [[ $REPLY =~ ^[Yy]$ ]]
then
CMD=y
else
CMD=n
fi


# Not set
if [ -z "$EVENTS" ]
then
	EVENTS=5000
fi


# Processes and plotting
yes $CMD | source ./tests/run_screening/run.sh
yes $CMD | source ./tests/run_cdf_single/run.sh
yes $CMD | source ./tests/run_alice_multi/run.sh
yes $CMD | source ./tests/run_alice_single/run.sh
yes $CMD | source ./tests/run_cms_multi/run.sh
yes $CMD | source ./tests/run_excitation/run.sh
yes $CMD | source ./tests/run_rhorho_phiphi/run.sh
yes $CMD | source ./tests/run_odderon_phi/run.sh

# Spin density generation and control
yes $CMD | source ./tests/run_JW_polarization/run.sh
yes $CMD | source ./tests/run_JW_frames/run.sh

# kt-EPA vs full QED yy in low-mass and high-mass domain
yes $CMD | source ./tests/run_cms_mumu_low_high/run.sh

# Gamma-gamma etc
yes $CMD | source ./tests/run_luxtest/run.sh
yes $CMD | source ./tests/run_flux/run.sh

# Eikonal/elastic scattering
yes $CMD | source ./tests/run_elastic/run.sh

# Minimum bias
yes $CMD | source ./tests/run_minbias/run.sh

# Durham QCD
yes $CMD | source ./tests/run_durham_chic0/run.sh
yes $CMD | source ./tests/run_durham_mmbar/run.sh

# Spherical harmonic expansion
yes $CMD | source ./tests/run_cms_harmonic/run.sh
yes $CMD | source ./tests/run_alice_harmonic/run.sh

# Tensor Pomeron
yes $CMD | source ./tests/run_tensor0_multi/run.sh
yes $CMD | source ./tests/run_tensor2_multi/run.sh
yes $CMD | source ./tests/run_tensor_spectrum/run.sh


# Catch2 (C++) driven tests
./bin/testbench0
./bin/testbench1
./bin/testbench2
./bin/testbench3


# Pytest driven tests (fast)
pytest ./tests/testbench_global_fast.py -s
pytest ./tests/testbench_integrators.py -s
pytest ./tests/testbench_vgrid.py -s


# Pytest driven tests (slow)
pytest ./tests/testbench_cepdata.py -s
pytest ./tests/testbench_global_slow.py -s
pytest ./tests/testbench_exloop.py -s


echo "[megarun.sh: done]"
