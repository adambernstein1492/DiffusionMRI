import numpy
import dipy.io
import nibabel as nib

def load_diffusion_data(dwi_file, bvals_file, bvecs_file, mask_file):
    dwi = nib.load(dwi_file)
    data = dwi.get_data()
    data[numpy.isnan(data)] = 0.0
    dwi = nib.Nifti1Image(data, dwi.affine, dwi.header)
    mask = nib.load(mask_file)
    bvals,bvecs = dipy.io.read_bvals_bvecs(bvals_file, bvecs_file)

    return dwi,mask,bvals,bvecs

def progress_update(message, percent):
    print "\r" + message + "%10.1f %%" % percent,

    if(percent == 100):
        print "\n"

def factn(number, n):
    fact_n = 1.0

    if (number < n):
        fact_n = 1.0
    else:
        for i in range(number,0,-n):
            fact_n *= i

    return fact_n

def check_diffusion_input(dwi_path, bval_path, bvec_path, mask_path):
    try:
        dwi = nib.load(dwi_path)
    except:
        print "Error: Image must be a 4-D NIFTI file"
        quit()

    try:
        length = dwi.shape[3]
    except:
        print "Error: Image must be a 4-D NIFTI file"
        quit()

    try:
        bval,bvec = dipy.io.read_bvals_bvecs(bval_path, bvec_path)
    except:
        print "Error: Cannot read in bval and bvecs files"
        quit()

    if dwi.shape[3] != bval.shape[0] or dwi.shape[3] != bvec.shape[0] or bvec.shape[1] != 3:
        print "Error: bvals, bvecs and dwi dimensions are not in agreement"
        quit()

    if mask_path != "None":
        try:
            mask = nib.load(mask_path)
        except:
            print "Error: Mask must be a 3-D NIFTI file"
            quit()

        if mask.shape != dwi.shape[0:3]:
            print "Error: Mask must have same dimensions as the dwi"
            quit()
