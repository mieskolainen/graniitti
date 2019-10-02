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
echo 'Channel' '&' '\sqrt{s} TeV' '&' 'Measurement' '&' 'xs +- stat +- syst' '&' 'screened xs' '&' 'bare xs' '&' '<S^2>' '\\'
echo '\hline \\'
echo ''

# Loop over steering cards
FILES_A=./tests/processes/*.json
FILES_B=./tests/LHC_TEVATRON_RHIC/*.json

echo ${FILES[@]}
echo ''

for F in $FILES_A $FILES_B
do

	# Extract name + strip of extension
	CARD=$(basename $F)
	CARD="${CARD%.*}"

	#echo $CARD

	# ====================================================================
	NAME="-"
	CHANNEL="-"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="-"
	# ====================================================================
	
	if [ "$CARD" == "ALICE19_2pi" ]
	then
	NAME="$|\eta| < 0.9, p_t > 0.15$ GeV"
	CHANNEL="$\pi^+\pi^-$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "ALICE19_2pi_E1" ]
	then
	NAME="$|\eta| < 0.9, p_t > 0.15$ GeV, $M_{1} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{SD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "ALICE19_2pi_E2" ]
	then
	NAME="$|\eta| < 0.9, p_t > 0.15$ GeV, $M_{1,2} < 5$ GeV"
	CHANNEL="$\pi^+\pi^-_{DD}$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "ALICE19_2pi_PWA" ]
	then
	NAME="$|Y_X| < 0.9$"
	CHANNEL="$\pi^+\pi^-$ "
	BIBTEX="-"
	MEASUREMENT=(- - -) # 31 0.5 2
	UNIT=6
	SQRTS="7"

	elif [ "$CARD" == "STAR18_2pi" ]
	then
	NAME="$|\eta| < 0.7, p_t > 0.2$ GeV"
	CHANNEL="\$\pi^+\pi^-$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=6
	SQRTS="0.2"

	elif [ "$CARD" == "gg2gg" ]
	then
	NAME="$|y| < 2.5, p_t > 20$ GeV"
	CHANNEL="\$gg$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=9
	SQRTS="13"
	
	elif [ "$CARD" == "WW7TeV" ]
	then
	NAME="Full $\4\pi$"
	CHANNEL="\$W^+W^-$"
	BIBTEX="-"
	MEASUREMENT=(- - -)
	UNIT=15
	SQRTS="7"
	
	# ====================================================================

	# LHC/Tevatron measurements

	# https://arxiv.org/abs/1502.01391
	elif [ "$CARD" == "CDF14_2pi" ]
	then
	NAME="CDF"
	CHANNEL="$\pi^+\pi^-$"
	BIBTEX="aaltonen2015measurement"
	MEASUREMENT=(- - -) # Integrated not given in the paper (one should integrate histograms)
	SQRTS="1.96"
	
	# http://cds.cern.ch/record/2679648/files/FSQ-16-006-pas.pdf
	elif [ "$CARD" == "CMS19_2pi" ]
	then
	NAME="CMS"
	CHANNEL="$\pi^+\pi^-$"
	BIBTEX="CMS-PAS-FSQ-16-006"
	MEASUREMENT=(19.0 0.6 3.2)
	UNIT=6
	SQRTS="13"
	
	# https://discoverycenter.nbi.ku.dk/teaching/thesis_page/MasterEmilBolsFinal.pdf
	elif [ "$CARD" == "ATLAS17_2pi" ]
	then
	NAME="ATLAS"
	CHANNEL="$\pi^+\pi^-$"
	BIBTEX="Bols:2288372"
	MEASUREMENT=(18.75 0.048 0.770)
	UNIT=6
	SQRTS="13"
	
	# https://arxiv.org/pdf/hep-ex/0611040.pdf
	elif [ "$CARD" == "CDF07_ee" ]
	then
	NAME="CDF"
	CHANNEL="\$e^+e^-$"
	BIBTEX="abulencia2007observation"
	MEASUREMENT=(1.6 0.5 0.3)
	UNIT=12
	SQRTS="1.96"
	
	# https://arxiv.org/pdf/1112.0858.pdf
	elif [ "$CARD" == "CDF11_ee" ]
	then
	NAME="CDF"
	CHANNEL="\$e^+e^-$"
	BIBTEX="aaltonen2012observation"
	MEASUREMENT=(2.88 0.57 0.63)
	UNIT=12
	SQRTS="1.96"
	
	# https://arxiv.org/abs/1111.5536
	elif [ "$CARD" == "CMS11_mumu" ]
	then
	NAME="CMS"
	CHANNEL="$\mu^+\mu^-$"
	BIBTEX="chatrchyan2012exclusive"
	MEASUREMENT=(3.38 0.58 0.21)
	UNIT=12
	SQRTS="7"
	
	# https://arxiv.org/abs/1506.07098
	elif [ "$CARD" == "ATLAS15_ee" ]
	then
	NAME="ATLAS"
	CHANNEL="\$e^+e^-$"
	BIBTEX="atlas2015measurement"
	MEASUREMENT=(0.428 0.035 0.018)
	UNIT=12
	SQRTS="7"
	
	# https://arxiv.org/abs/1506.07098
	elif [ "$CARD" == "ATLAS15_mumu" ]
	then
	NAME="ATLAS"
	CHANNEL="$\mu^+\mu^-$"	
	BIBTEX="atlas2015measurement"
	MEASUREMENT=(0.628 0.032 0.021)
	UNIT=12
	SQRTS="7"
	
	# https://arxiv.org/abs/1708.04053
	elif [ "$CARD" == "ATLAS17_mumu" ]
	then
	NAME="ATLAS"
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
	./bin/gr -i $F -l false -n 0 > /tmp/GRANIITTI.out

	if [ $SCREENING == 1 ]
	then
	./bin/gr -i $F -l true -n 0 > /tmp/GRANIITTI_S2.out		
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
	XS_S2=$(python -c "print('%0.3g' % ($XS_S2*1E$UNIT))")	
	XS=$(python -c "print('%0.3g'    % ($XS*1E$UNIT))")

	# Calculate division
	S2=$(python -c "print('%0.2f' % ($XS_S2/$XS))")

	# Print out
	echo $CHANNEL ' & ' $SQRTS ' & ' $NAME '\cite{'$BIBTEX'}' ' & ' $X_STAT_SYST ' & ' $XS_S2 ' & ' $XS ' ' $BARN ' & ' $S2 '\\'
	#printf '%s & %s & %s & %s \n' $CARD $X_STAT_SYST $XS_S2 $XS

done

echo ''
echo 'exloop.sh [DONE]'


# Analyze
#  // <https://indico.cern.ch/event/713101/contributions/3102315/attachments/1705771/2748440/Diffraction2018_RafalSikora.pdf>
#  input.push_back(MEASUREMENT(PATH + "STAR18_2pi.json" 0 0 0));

