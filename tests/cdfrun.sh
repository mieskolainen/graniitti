make && ./bin/gr -i ./tests/experiment/CDF_2pi0.json -p 1

make -j4 ROOT=TRUE && ./bin/analyzer -i CDF_2pi0 -g 22 -n 4 -l \
	'2#pi^{0}#rightarrow2(#gamma#gamma)' -M 2.5 -Y 1.5 -P 1.5 -u ub
