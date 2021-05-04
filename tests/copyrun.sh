#/bin/sh
#
# Copy pdfs from the simulations and analysis

# Copy target folder
F0=/home/user/Dropbox/GRANIITTI_tex/figs
F=/home/user/Dropbox/GRANIITTI_tex/mcfigs

#yes $CMD | source ./tests/run_atlas_multi/run.sh
FOLDER=ATLAS+ATLAS_NEG+ATLAS_POS
for FILE in h1_S_M hP_S_M_Pt; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_ATLASPIPI.pdf
done

#yes $CMD | source ./tests/run_screening/run.sh
FOLDER=continuum+continuum_screened
for FILE in h1_PP_dpt_logy h1_PP_t1_logy; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_KK.pdf
done

#yes $CMD | source ./tests/run_alice_multi/run.sh
FOLDER=ALICE_2pi+ALICE_2K+ALICE_ppbar
for FILE in h1_S_M_logy; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_PIPIKKPPBAR.pdf
done

#yes $CMD | source ./tests/run_excitation/run.sh
FOLDER=2pi_excite_0+2pi_excite_1+2pi_excite_2
for FILE in hP_S_M_Pt; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_EXCITE.pdf
done

#yes $CMD | source ./tests/run_JW_polarization/run.sh
FOLDER=f0_980+rho_JZ0+rho_JZ1+f2_JZ0+f2_JZ1+f2_JZ2
for FILE in h1_PP_dpt_logy h1_2B_acop_logy h1_2B_diffrap_logy h1_PP_dphi_logy; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_JW.pdf
done
for FILE in h1_costheta_CM_logy h1_costheta_CS_logy h1_costheta_HX_logy h1_costheta_GJ_logy; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_JW.pdf
done
for FILE in h1_phi_CM_logy h1_phi_CS_logy h1_phi_HX_logy h1_phi_GJ_logy; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_JW.pdf
done

#yes $CMD | source ./tests/run_JW_frames/run.sh

## Tensor Pomeron
#yes $CMD | source ./tests/run_tensor0_multi/run.sh

#yes $CMD | source ./tests/run_tensor2_multi/run.sh
FOLDER=f2_0+f2_1+f2_2+f2_3+f2_4+f2_5+f2_6
for FILE in h1_costheta_GJ_logy h1_PP_t1_logy; do
	cp ./figs/$FOLDER/$FILE.pdf $F/$FILE\_TENSOR.pdf
done

#yes $CMD | source ./tests/run_tensor_spectrum/run.sh

## Spherical harmonic expansion
#yes $CMD | source ./tests/run_cms_harmonic/run.sh
FOLDER=harmonicfit/SH_2pi_J0_CMS___SH_2pi_J0_CMS+SH_2pi_CMS
for SUBFOLDER in OBS_0_CM OBS_0_CS OBS_0_HX OBS_0_PG OBS_0_GJ OBS_0_AH; do
for FILE in "h_Response" "h_{Moments}[MPP]<fid>" "h_{Moments}[MPP]<fla>"; do
	cp ./figs/$FOLDER/$SUBFOLDER/$FILE.pdf $F\harmonic/$SUBFOLDER\_$FILE.pdf
done
done

#yes $CMD | source ./tests/run_alice_harmonic/run.sh

# MATLAB analysis
cp ./tests/analysis_eikonal/figs/t_space_dsigmadt.pdf $F0/
#
cp ./tests/analysis_scan_xs/figs/xsec.pdf $F0/
cp ./tests/analysis_scan_xs/figs/xsecratios.pdf $F0/
#
cp ./tests/analysis_xsec/figs/sd.pdf $F0/
cp ./tests/analysis_xsec/figs/dd.pdf $F0/
cp ./tests/analysis_xsec/figs/total.pdf $F0/
cp ./tests/analysis_xsec/figs/ratio.pdf $F0/
#
cp ./tests/analysis_sudakov/figs/sudakov_1D_mu.pdf $F0/
cp ./tests/analysis_sudakov/figs/sudakov_1D_Q2.pdf $F0/

echo "[copyrun.sh: done]"
