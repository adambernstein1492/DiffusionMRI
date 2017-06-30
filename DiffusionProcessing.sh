#! /bin/bash
DataPath=$1
SavePath=$2
MRTrixPath=$3
SubjectID=$4





### CONVERT TO MRTRIX FORMAT ##################################################
mkdir $MRTrixPath
mrconvert -fslgrad $DataPath/bvec $DataPath/bval $DataPath/data.nii $MRTrixPath/DWI.mif
mrconvert mask.nii mask.mif
cd $MRTrixPath
###############################################################################



### MRTRIX FOD ESTIMATION #####################################################
# Estimate the Response Function
dwi2response tournier DWI.mif response.txt

# Estimate the FODs for each voxel within the mask
dwi2fod csd DWI.mif response.txt FOD.mif -mask MASK.mif
###############################################################################



### RUN FREESURFER ############################################################
export SUBJECT_DIR=$SavePath
recon-all -i $DataPath/T1.nii -subject $SubjectID -all
bbregiter --s $SubjectID --mov meanb0_brain.nii.gz --reg epi2T1.lta --dti

### MRTRIX ANATOMICALLY CONTRAINED TRACTOGRAPHY ###############################
# Isolate B0 images
dwiextract -bzero DWI.mif - | mrmath - mean meanb0.mif -axis 3

# Register T1 to DWI
flirt -in T1.nii -ref B0Avg.nii -omat T12DWI.mat -dof 6

###############################################################################

