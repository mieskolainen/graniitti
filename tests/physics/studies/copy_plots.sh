#!/usr/bin/env bash
#
# Copy pdfs from the simulations and analysis

set -euo pipefail

# Select the two explicit copy targets
: "${FIGURE_TARGET:?Set FIGURE_TARGET to the analysis figure directory}"
: "${MC_FIGURE_TARGET:?Set MC_FIGURE_TARGET to the Monte Carlo figure directory}"
F0="${FIGURE_TARGET}"
F="${MC_FIGURE_TARGET}"
mkdir -p "$F0" "$F"

FOLDER=ATLAS+ATLAS_NEG+ATLAS_POS
for FILE in h1_S_M hP_S_M_Pt; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_ATLASPIPI.pdf"
done

FOLDER=continuum+continuum_screened
for FILE in h1_PP_dpt_logy h1_PP_t1_logy; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_KK.pdf"
done

FOLDER=ALICE_2pi+ALICE_2K+ALICE_ppbar
for FILE in h1_S_M_logy; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_PIPIKKPPBAR.pdf"
done

FOLDER=2pi_excite_0+2pi_excite_1+2pi_excite_2
for FILE in hP_S_M_Pt; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_EXCITE.pdf"
done

FOLDER=f0_980+rho_JZ0+rho_JZ1+f2_JZ0+f2_JZ1+f2_JZ2
for FILE in h1_PP_dpt_logy h1_2B_acop_logy h1_2B_diffrap_logy h1_PP_dphi_logy; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_JW.pdf"
done
for FILE in h1_costheta_CM_logy h1_costheta_CS_logy h1_costheta_HX_logy h1_costheta_GJ_logy; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_JW.pdf"
done
for FILE in h1_phi_CM_logy h1_phi_CS_logy h1_phi_HX_logy h1_phi_GJ_logy; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_JW.pdf"
done


## Tensor Pomeron

FOLDER=f2_0+f2_1+f2_2+f2_3+f2_4+f2_5+f2_6
for FILE in h1_costheta_GJ_logy h1_PP_t1_logy; do
	cp "./figs/$FOLDER/$FILE.pdf" "$F/${FILE}_TENSOR.pdf"
done


# Python analysis
for BEAM in pp ppbar; do
    mkdir -p "$F0/elastic/$BEAM"
    cp -R "./figs/elastic/$BEAM/." "$F0/elastic/$BEAM/"
done
#
cp ./tests/physics/studies/cross_section_scan/analysis/figs/xsec.pdf "$F0/"
cp ./tests/physics/studies/cross_section_scan/analysis/figs/xsecratios.pdf "$F0/"
#
cp ./tests/physics/studies/total_cross_sections/analysis/figs/sd.pdf "$F0/"
cp ./tests/physics/studies/total_cross_sections/analysis/figs/dd.pdf "$F0/"
cp ./tests/physics/studies/total_cross_sections/analysis/figs/total.pdf "$F0/"
cp ./tests/physics/studies/total_cross_sections/analysis/figs/ratio.pdf "$F0/"
#
mkdir -p "$F0/sudakov"
cp -R ./tests/physics/studies/sudakov/figs/*.json "$F0/sudakov/"

echo "[copy_plots.sh: done]"
