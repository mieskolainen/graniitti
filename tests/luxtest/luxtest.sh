./bin/gr -i ./tests/luxtest/LUX.json -w true -n 50000

# Scale factor 3 x for three lepton flavors (we generate a sample only for one flavor)
./bin/analyze -i LUX -g 11 -n 2 -l 'l^{+}l^{-} / 13 TeV ' -M 5000 -Y 1.5 -P 2.0 -u fb -S 3 #-X 1000
