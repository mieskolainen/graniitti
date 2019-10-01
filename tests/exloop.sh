#!/bin/sh
#
# Loop over fiducial measurements and construct latex array
#
# Run with: source ./tests/exloop.sh


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

echo '\hline \\'
echo 'Fiducial measurement' '&' '\sqrt{s} TeV' '&' 'xs +- stat +- syst' '&' 'screened xs' '&' 'bare xs' '&' '<S^2>' '\\'
echo '\hline \\'
echo ''

# Loop over steering cards
FILES=./tests/LHC_TEVATRON_RHIC/*.json
for F in $FILES
do

	# Extract filename + strip of extension
	CARD=$(basename $F)
	CARD="${CARD%.*}"
	
	# ====================================================================
	NAME=$CARD
	BIBTEX="-"
	MEASUREMENT=(0.0 0.0 0.0)
	UNIT=6
	SQRTS="-"
	# ====================================================================

	# https://arxiv.org/abs/1502.01391
	if [ "$CARD" == "CDF14_2pi" ]
	then
	NAME="CDF $\pi^+\pi^-$"
	BIBTEX="aaltonen2015measurement"
	MEASUREMENT=(0 0 0)
	SQRTS="1.96"
	fi
	# http://cds.cern.ch/record/2679648/files/FSQ-16-006-pas.pdf
	if [ "$CARD" == "CMS19_2pi" ]
	then
	NAME="CMS $\pi^+\pi^-$"
	BIBTEX="CMS-PAS-FSQ-16-006"
	MEASUREMENT=(19.0 0.6 3.2)
	UNIT=6
	SQRTS="13"
	fi
	#  <ALICE 7 TeV PWA |y(pi+pi-)| < 0.9
	if [ "$CARD" == "ALICE19_2pi_PWA" ]
	then
	NAME="ALICE $\pi^+\pi^-$"
	BIBTEX="-"
	MEASUREMENT=(31 0.5 2)
	UNIT=6
	SQRTS="7"
	fi
	# https://discoverycenter.nbi.ku.dk/teaching/thesis_page/MasterEmilBolsFinal.pdf
	if [ "$CARD" == "ATLAS17_2pi" ]
	then
	NAME="ATLAS $\pi^+\pi^-$"
	BIBTEX="Bols:2288372"
	MEASUREMENT=(18.75 0.048 0.770)
	UNIT=6
	SQRTS="13"
	fi
	# https://arxiv.org/pdf/hep-ex/0611040.pdf
	if [ "$CARD" == "CDF07_ee" ]
	then
	NAME="CDF \$e^+e^-$"
	BIBTEX="abulencia2007observation"
	MEASUREMENT=(1.6 0.5 0.3)
	UNIT=12
	SQRTS="1.96"
	fi
	# https://arxiv.org/pdf/1112.0858.pdf
	if [ "$CARD" == "CDF11_ee" ]
	then
	NAME="CDF \$e^+e^-$"
	BIBTEX="aaltonen2012observation"
	MEASUREMENT=(2.88 0.57 0.63)
	UNIT=12
	SQRTS="1.96"
	fi
	# https://arxiv.org/abs/1111.5536
	if [ "$CARD" == "CMS11_mumu" ]
	then
	NAME="CMS $\mu^+\mu^-$"
	BIBTEX="cms2011exclusive"
	MEASUREMENT=(3.38 0.58 0.21)
	UNIT=12
	SQRTS="7"
	fi
	# https://arxiv.org/abs/1506.07098
	if [ "$CARD" == "ATLAS15_ee" ]
	then
	NAME="ATLAS \$e^+e^-$"
	BIBTEX="atlas2015measurement"
	MEASUREMENT=(0.428 0.035 0.018)
	UNIT=12
	SQRTS="7"
	fi
	# https://arxiv.org/abs/1506.07098
	if [ "$CARD" == "ATLAS15_mumu" ]
	then
	NAME="ATLAS $\mu^+\mu^-$"	
	BIBTEX="atlas2015measurement"
	MEASUREMENT=(0.628 0.032 0.021)
	UNIT=12
	SQRTS="7"
	fi
	# https://arxiv.org/abs/1708.04053
	if [ "$CARD" == "ATLAS17_mumu" ]
	then
	NAME="ATLAS $\mu^+\mu^-$"
	BIBTEX="aaboud2018measurement"
	MEASUREMENT=(3.12 0.07 0.14)
	UNIT=12
	SQRTS="13"
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

	# Generate
	./bin/gr -i ./tests/LHC_TEVATRON_RHIC/$CARD.json -l false -n 0 > /tmp/GRANIITTI.out

	if [ $SCREENING == 1 ]
	then
	./bin/gr -i ./tests/LHC_TEVATRON_RHIC/$CARD.json -l true -n 0 > /tmp/GRANIITTI_S2.out		
	fi

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

	# Add +- characters
	printf -v ids_d ' $\pm$ %s' "${MEASUREMENT[@]}" # yields "value +- stat +- syst"
	ids_d=${ids_d:6}                             # remove the leading ' +- '


	XS_S2=$(python -c "print('%0.2g' % ($XS_S2*1E$UNIT))")	
	XS=$(python -c "print('%0.2g'    % ($XS*1E$UNIT))")

	# Calculate division
	S2=$(python -c "print('%0.2f' % ($XS_S2/$XS))")

	# Print out
	echo $NAME ' ' '\cite{'$BIBTEX'}' ' & ' $SQRTS ' & ' $ids_d ' ' $BARN ' & ' $XS_S2 ' & ' $XS ' & ' $S2 '\\'
	#printf '%s & %s & %s & %s \n' $CARD $ids_d $XS_S2 $XS

done

echo ''
echo 'exloop.sh [DONE]'


# Analyze
#  // <https://indico.cern.ch/event/713101/contributions/3102315/attachments/1705771/2748440/Diffraction2018_RafalSikora.pdf>
#  input.push_back(MEASUREMENT(PATH + "STAR18_2pi.json" 0 0 0));

