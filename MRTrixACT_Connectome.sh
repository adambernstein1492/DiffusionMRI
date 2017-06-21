#! /bin/bash
DiffDataPath=$1
T1DataPath=$2
SubjectID=$3
MRtrixPath=$4
OutPath=$5

mkdir $OutPath
cd $OutPath
mkdir temp



### REGISTER T1 TO DWI ########################################################
# Create Average B0 Image
dwiextract $DiffDataPath/data.nii -fslgrad $DiffDataPath/bvec $DiffDataPath/bval - -bzero | mrmath - mean temp/MeanB0.nii -axis 3

SUBJECTS_DIR=$T1DataPath
# Run bbregister from freesurfer
bbregister --s $SubjectID --mov temp/MeanB0.nii --reg temp/B02T1.lta --T2

# Convert Registration Matrix to MRTrix Format
lta_convert --inlta temp/B02T1.lta --outfsl temp/B02T1.mat

transformconvert temp/B02T1.mat temp/MeanB0.nii $T1DataPath/$SubjectID/mri/orig.mgz flirt_import temp/B02T1_mrtrix.txt

# Apply Inverse Transform to T1 and Parcellation
transformcalc temp/B02T1_mrtrix.txt invert temp/T12B0_mrtrix.txt

mrtransform -linear temp/T12B0_mrtrix.txt $T1DataPath/$SubjectID/mri/orig.mgz T1_Reg.mif
mrtransform -interp nearest -linear temp/T12B0_mrtrix.txt $T1DataPath/$SubjectID/mri/aparc+aseg.mgz aparc+aseg.mif
###############################################################################



### CREATE SEGMENTATION MASK ##################################################
5ttgen freesurfer aparc+aseg.mif temp/5TT.mif
###############################################################################



### CONVERT LABELS ############################################################
labelconvert aparc+aseg.mif $FREESURFER_HOME/FreeSurferColorLUT.txt $MRtrixPath/src/connectome/tables/fs_default.txt nodes.mif
###############################################################################



### DWI PROCESSING ############################################################
# Convert to mif
mrconvert $DiffDataPath/data.nii -fslgrad $DiffDataPath/bvec $DiffDataPath/bval DWI.mif

# Estimate Response Function
dwi2response msmt_5tt DWI.mif temp/5TT.mif temp/RF_WM.txt temp/RF_GM.txt temp/RF_CSF.txt -voxels temp/RF_voxels.mif

# Estimate FODs
dwi2fod msmt_csd DWI.mif temp/RF_WM.txt temp/WM_FODs.mif temp/RF_GM.txt temp/GM.mif temp/RF_CSF.txt temp/CSF.mif -mask $DiffDataPath/mask.nii

# Tractography
tckgen temp/WM_FODs.mif temp/10M.tck -act temp/5TT.mif -backtrack -crop_at_gmwmi -seed_dynamic temp/WM_FODs.mif -maxlength 250 -number 10M -cutoff 0.06

tcksift temp/10M.tck temp/WM_FODs.mif 5M_SIFT.tck -act temp/5TT.mif -term_number 5M
###############################################################################



### CONNECTOME GENERATION #####################################################
tck2connectome 1M_SIFT.tck nodes.mif connectome.csv -out_assignments tcklabels.txt

# Parse out Tracts between regions
mkdir Tracts
connectome2tck 1M_SIFT.tck tcklabels.txt Tracts/IsoTracts
###############################################################################



### CLEAN UP ##################################################################
rm -rf temp
###############################################################################
