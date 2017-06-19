import numpy
import dipy.io
import nibabel as nib

def load_diffusion_data(dwi_file, bvals_file, bvecs_file, mask_file):
    dwi = nib.load(dwi_file)
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
