# !/bin/sh
#
# Visualize events with mcviz and test hepmc record structure validity


# Generate events
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS17_2pi_E0.json -o test0 -f hepmc2 -p "PP[CON]<F> -> pi+ pi-" -s 0 -n 10
./tests/run_mcviz/mcviz/mcv ./output/test0.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

./bin/gr -i ./tests/LHC_TEVATRON_RHIC/STAR18_2pi.json -o test0B -f hepmc2 -p "PP[CON]<F> -> pi+ pi-" -s 0 -n 10
./tests/run_mcviz/mcviz/mcv ./output/test0B.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

# Generate events
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS17_2pi_E0.json -o test1 -f hepmc2 -s 1 -p "PP[CON]<F> -> pi+ pi-" -s 1 -n 10
./tests/run_mcviz/mcviz/mcv ./output/test1.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

./bin/gr -i ./tests/LHC_TEVATRON_RHIC/STAR18_2pi.json -o test1B -f hepmc2 -s 1 -p "PP[CON]<F> -> pi+ pi-" -s 1 -n 10
./tests/run_mcviz/mcviz/mcv ./output/test1B.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

# Generate events
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS17_2pi_E0.json -o test2 -f hepmc2 -s 2 -p "PP[CON]<F> -> pi+ pi-" -s 2 -n 10
./tests/run_mcviz/mcviz/mcv ./output/test2.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

./bin/gr -i ./tests/LHC_TEVATRON_RHIC/STAR18_2pi.json -o test2B -f hepmc2 -s 2 -p "PP[CON]<F> -> pi+ pi-" -s 2 -n 10
./tests/run_mcviz/mcviz/mcv ./output/test2B.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg


# ----------------------------------------


# Generate events
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS17_2pi_E0.json -o test4 -f hepmc2 -p "X[EL]<Q>" -l true -n 10
./tests/run_mcviz/mcviz/mcv ./output/test4.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

./bin/gr -i ./tests/LHC_TEVATRON_RHIC/STAR18_2pi.json -o test4B -f hepmc2 -p "X[EL]<Q>" -l true -n 10
./tests/run_mcviz/mcviz/mcv ./output/test4B.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

# Generate events
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS17_2pi_E0.json -o test5 -f hepmc2 -p "X[SD]<Q>" -l true -n 10
./tests/run_mcviz/mcviz/mcv ./output/test5.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

./bin/gr -i ./tests/LHC_TEVATRON_RHIC/STAR18_2pi.json -o test5B -f hepmc2 -p "X[SD]<Q>" -l true -n 10
./tests/run_mcviz/mcviz/mcv ./output/test5B.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

# Generate events
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS17_2pi_E0.json -o test6 -f hepmc2 -p "X[DD]<Q>" -l true -n 10
./tests/run_mcviz/mcviz/mcv ./output/test6.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

./bin/gr -i ./tests/LHC_TEVATRON_RHIC/STAR18_2pi.json -o test6B -f hepmc2 -p "X[DD]<Q>" -l true -n 10
./tests/run_mcviz/mcviz/mcv ./output/test6B.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

# Generate events
./bin/gr -i ./tests/LHC_TEVATRON_RHIC/CMS17_2pi_E0.json -o test7 -f hepmc2 -p "X[ND]<Q>" -l true -n 10
./tests/run_mcviz/mcviz/mcv ./output/test7.hepmc2:0 -sFancyLines --verbose && firefox mcviz.svg

