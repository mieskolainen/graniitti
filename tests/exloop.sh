#!/bin/sh
#
# Loop over fiducial measurements and construct latex array
#
# Run with: bash ./tests/exloop.sh


read -p "run: Screening on? [y/n] " -n 1 -r
echo # New line

if [[ $REPLY =~ ^[Yy]$ ]]
then
SCREENING=1
else
SCREENING=0
fi

echo ''
echo 'Hint: run cat /tmp/GRANIITTI.out or cat /tmp/GRANIITTI_S2.out for the progress ...'
echo ''
echo ''

# Print table
echo '\begin{table}'
echo '\begin{center}'
echo '\begin{tabular}{|lllr|rrc|}'
echo '\hline'
echo ' & & & \tiny{MEASUREMENT} & & \tiny{GRANIITTI} & \\'
echo '\hline'
echo ' & $\sqrt{s}$ & \tiny{PHASE SPACE CUTS} & value $\pm$ stat $\pm$ syst & $\sigma_{S^2}$ & $\sigma_0$ & $\langle S^2 \rangle$ \\'
echo '\hline'

# Loop over steering cards
FILES_A=./tests/processes/*.json
FILES_B=./tests/LHC_TEVATRON_RHIC/*.json

for F in $FILES_A $FILES_B
do

	# Extract CUTS + strip of extension
	CARD=$(basename $F)
	CARD="${CARD%.*}"

	#echo $CARD

	# ====================================================================
	CUTS="-"
	CHANNEL="-"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="-"
	# ====================================================================
	# Simulation without experimental reference


	if   [ "$CARD" == "ALICE19_2pi_E0" ]
	then
	CUTS="$|\eta| < 0.9, p_t > 0.15$ GeV"
	CHANNEL="$\pi^+\pi^-_{EL}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "ALICE19_2pi_E1" ]
	then
	CUTS="$|\eta| < 0.9, p_t > 0.15, M_{1} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{SD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "ALICE19_2pi_E2" ]
	then
	CUTS="$|\eta| < 0.9, p_t > 0.15, M_{1,2} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{DD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "ALICE19_2pi_PWA" ]
	then
	CUTS="$|Y_X| < 0.9$"
	CHANNEL="$\pi^+\pi^-$ "
	BIBTEX="-"
	MEASUREMENT=(- - -) # 31 0.5 2
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "STAR18_2pi" ]
	then
	CUTS="$|\eta| < 0.7, p_t > 0.2$ GeV"
	CHANNEL="\$\pi^+\pi^-$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="0.2"

	elif [ "$CARD" == "gg2gg" ]
	then
	CUTS="$|y| < 2.5, p_t > 20$ GeV"
	CHANNEL="\$gg$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=9
	SQRTS="13"
	
	elif [ "$CARD" == "WW7TeV" ]
	then
	CUTS="Full \$4\pi$"
	CHANNEL="\$W^+W^-$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=15
	SQRTS="7"
	
	# ====================================================================
	# Simulations with LHC/Tevatron measurement reference


	# https://arxiv.org/abs/1502.01391
	elif [ "$CARD" == "CDF14_2pi" ]
	then
	CUTS="CDF"
	CHANNEL="$\pi^+\pi^-$"
	BIBTEX="aaltonen2015measurement"
	MEASUREMENT=(x x x) # Integrated not given in the paper (one should integrate histograms)
	SQRTS="1.96"

	# https://arxiv.org/abs/1706.08310
	elif [ "$CARD" == "CMS17_2pi_E0" ]
	then
	CUTS="CMS"
	CHANNEL="$\pi^+\pi^-_{EL}$"
	BIBTEX="khachatryan2017exclusive"
	MEASUREMENT=(26.5 0.3 5.12) # 26.5 +/- 0.3 (stat) +/- 5.0 (syst) +/- 1.1 (lumi) microbarns.
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "CMS17_2pi_E1" ]
	then
	CUTS="$|y| < 2, p_t > 0.2, M_{1} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{SD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "CMS17_2pi_E2" ]
	then
	CUTS="$|y| < 2, p_t > 0.2, M_{1,2} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{DD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	# http://cds.cern.ch/record/2679648/files/FSQ-16-006-pas.pdf
	elif [ "$CARD" == "CMS19_2pi_E0" ]
	then
	CUTS="CMS"
	CHANNEL="$\pi^+\pi^-_{EL}$"
	BIBTEX="CMS-PAS-FSQ-16-006"
	MEASUREMENT=(19.0 0.6 3.2)
	UNIT=6
	SQRTS="13"
	
	elif [ "$CARD" == "CMS19_2pi_E1" ]
	then
	CUTS="$|\eta| < 2.4, p_t > 0.2, M_{1} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{SD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="13"
	
	elif [ "$CARD" == "CMS19_2pi_E2" ]
	then
	CUTS="$|\eta| < 2.4, p_t > 0.2, M_{1,2} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{DD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="13"
	
	# https://discoverycenter.nbi.ku.dk/teaching/thesis_page/MasterEmilBolsFinal.pdf
	elif [ "$CARD" == "ATLAS17_2pi" ]
	then
	CUTS="ATLAS (Thesis)"
	CHANNEL="$\pi^+\pi^-$"
	BIBTEX="Bols:2288372"
	MEASUREMENT=(18.75 0.048 0.770)
	UNIT=6
	SQRTS="13"
	
#	elif [ "$CARD" == "ATLAS17_4pi_0" ]
#	then
#	CUTS="ATLAS (Thesis)"
#	CHANNEL="\tiny{$\pi^+\pi^-\pi^+\pi^-$}"
#	BIBTEX="Bols:2288372"
#	MEASUREMENT=(3.575 0.065 0.338)
#	UNIT=6
#	SQRTS="13"

#	elif [ "$CARD" == "ATLAS17_4pi_1_2f0500" ]
#	then
#	CUTS="$|\eta| < 2.5, p_t > 0.1$ GeV"
#	CHANNEL="\tiny{\$2(f_0 \rightarrow \pi^+\pi^-)$}"
#	BIBTEX="-"
#	MEASUREMENT=(- - -)
#	UNIT=6
#	SQRTS="13"

	# https://arxiv.org/pdf/hep-ex/0611040.pdf
	elif [ "$CARD" == "CDF07_ee" ]
	then
	CUTS="CDF"
	CHANNEL="\$e^+e^-$"
	BIBTEX="abulencia2007observation"
	MEASUREMENT=(1.6 0.5 0.3)
	UNIT=12
	SQRTS="1.96"
	
	# https://arxiv.org/pdf/1112.0858.pdf
	elif [ "$CARD" == "CDF11_ee" ]
	then
	CUTS="CDF"
	CHANNEL="\$e^+e^-$"
	BIBTEX="aaltonen2012observation"
	MEASUREMENT=(2.88 0.57 0.63)
	UNIT=12
	SQRTS="1.96"
	
	# https://arxiv.org/abs/1111.5536
	elif [ "$CARD" == "CMS11_mumu" ]
	then
	CUTS="CMS"
	CHANNEL="$\mu^+\mu^-$"
	BIBTEX="chatrchyan2012exclusive"
	MEASUREMENT=(3.38 0.58 0.21)
	UNIT=12
	SQRTS="7"
	
	# https://arxiv.org/abs/1506.07098
	elif [ "$CARD" == "ATLAS15_ee" ]
	then
	CUTS="ATLAS"
	CHANNEL="\$e^+e^-$"
	BIBTEX="atlas2015measurement"
	MEASUREMENT=(0.428 0.035 0.018)
	UNIT=12
	SQRTS="7"
	
	# https://arxiv.org/abs/1506.07098
	elif [ "$CARD" == "ATLAS15_mumu" ]
	then
	CUTS="ATLAS"
	CHANNEL="$\mu^+\mu^-$"	
	BIBTEX="atlas2015measurement"
	MEASUREMENT=(0.628 0.032 0.021)
	UNIT=12
	SQRTS="7"
	
	# https://arxiv.org/abs/1708.04053
	elif [ "$CARD" == "ATLAS17_mumu" ]
	then
	CUTS="ATLAS"
	CHANNEL="$\mu^+\mu^-$"
	BIBTEX="aaboud2018measurement"
	MEASUREMENT=(3.12 0.07 0.14)
	UNIT=12
	SQRTS="13"
	
	# Not in our list
	else
	continue
	fi

	# ====================================================================
	if [ $UNIT == 3 ]
	then
	BARN="mb"
	fi
	if [ $UNIT == 6 ]
	then
	BARN="$\mu\$b"
	fi
	if [ $UNIT == 9 ]
	then
	BARN="nb"
	fi
	if [ $UNIT == 12 ]
	then
	BARN="pb"
	fi
	if [ $UNIT == 15 ]
	then
	BARN="fb"
	fi

	# ====================================================================
	# Generate
	./bin/gr -i $F -l false -h 0 -n 0 > /tmp/GRANIITTI.out

	if [ $SCREENING == 1 ]
	then
	./bin/gr -i $F -l true  -h 0 -n 0 > /tmp/GRANIITTI_S2.out		
	fi
	
	# ====================================================================
	# Extract integrated cross section
	XS=$(grep -E 'barn' /tmp/GRANIITTI.out)
	XS=${XS#*[}
	XS=${XS%]*}
	XS=${XS%+-*} # Remove stat errors

	# Extract screened integrated cross section
	XS_S2=0.0
	if [ "$SCREENING" == 1 ]
	then
	XS_S2=$(grep -E 'barn' /tmp/GRANIITTI_S2.out)
	XS_S2=${XS_S2#*[}
	XS_S2=${XS_S2%]*}
	XS_S2=${XS_S2%+-*} # Remove stat errors
	fi

	# ====================================================================
	# Add +- characters
	printf -v X_STAT_SYST ' $\pm$ %s' "${MEASUREMENT[@]}" # yields "value +- stat +- syst"
	X_STAT_SYST=${X_STAT_SYST:6}                          # remove the leading ' +- '

	# Format
	XS=$(python -c    "print('%0.3g' % ($XS*1E$UNIT))")
	XS_S2=$(python -c "print('%0.3g' % ($XS_S2*1E$UNIT))")	

	# Calculate division
	S2=$(python -c "print('%0.2f' % ($XS_S2/$XS))")

	# Do we have a measurement to cite
	if [ $BIBTEX == "-" ]
	then
	CUTS_AND_CITE="\tiny{$CUTS}"
	X_STAT_SYST= # no measurement
	else
	CUTS_AND_CITE="$CUTS \cite{$BIBTEX}"
	fi

	echo $CHANNEL ' & ' $SQRTS ' & ' $CUTS_AND_CITE ' & ' $X_STAT_SYST ' & ' $XS_S2 ' & ' $XS ' ' $BARN ' & ' $S2 '\\'
	#printf '%s & %s & %s & %s \n' $CARD $X_STAT_SYST $XS_S2 $XS

done

echo '\hline'
echo '\end{tabular}'
echo '\caption{Measurements versus GRANIITTI.}'
echo '\label{table:xstable}'
echo '\end{center}'
echo '\end{table}'
echo ''

# Analyze
#  // <https://indico.cern.ch/event/713101/contributions/3102315/attachments/1705771/2748440/Diffraction2018_RafalSikora.pdf>
#  input.push_back(MEASUREMENT(PATH + "STAR18_2pi.json" 0 0 0));

