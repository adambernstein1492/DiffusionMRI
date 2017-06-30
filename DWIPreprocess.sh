#! /bin/bash
dwiDICOMpath=$1
RPEDICOMpath=$2
SavePath=$3

# Create Output Directory
mkdir -p $SavePath
mkdir $SavePath/MIFs


### CONVERT TO MRTRIX FORMAT ##################################################
mrconvert $dwiDICOMpath $SavePath/MIFs/DWI.mif
mrconvert $RPEDICOMpath $SavePath/MIFs/RPE.mif
cd ${SavePath}/MIFs
###############################################################################



### ISOLATE B=0 IMAGES ########################################################
dwiextract -bzero DWI.mif APB0.mif
mrcat APB0.mif RPE.mif B0s.mif
###############################################################################


### RUN TOPUP AND EDDY ########################################################
dwipreproc DWI.mif DWI_top_eddy.mif -rpe_header -se_epi B0s.mif
###############################################################################



### RUN BRAIN EXTRACTION ######################################################
dwi2mask DWI_top_eddy.mif mask.mif
###############################################################################



### RUN BIAS FIELD CORRECTION #################################################
#dwibiascorrect -ants -mask mask.mif DWI_top_eddy.mif DWI_Preprocessed.mif
mv DWI_top_eddy.mif DWI_Preprocessed.mif
###############################################################################



### SAVE NIFTI VERSIONS #######################################################
mkdir ${SavePath}/NIFTIs
mrconvert mask.mif ${SavePath}/NIFTIs/mask.nii
mrconvert -export_grad_fsl ${SavePath}/NIFTIs/bvec ${SavePath}/NIFTIs/bval DWI_Preprocessed.mif ${SavePath}/NIFTIs/DWI_Preprocessed.nii
###############################################################################



### CLEAN UP EXTRA FILES ######################################################
cd ${SavePath}/MIFs
#rm APB0.mif
#rm DWI.mif
#rm DWI_top_eddy.mif
#rm DWI_top_eddy_denoised.mif
#rm noise.mif
#rm res.mif
#rm RPE.mif
###############################################################################
