import numpy as np
import dipy.io
import nibabel as nib
import os
import shutil

def load_diffusion_data(dwi_file, bvals_file, bvecs_file, mask_file):
    dwi = nib.load(dwi_file)
    data = dwi.get_data()
    data[np.isnan(data)] = 0.0
    dwi = nib.Nifti1Image(data, dwi.affine, dwi.header)
    mask = nib.load(mask_file)
    bvals,bvecs = dipy.io.read_bvals_bvecs(bvals_file, bvecs_file)

    signs = get_bvec_signs(dwi)
    for i in range(3):
        bvecs[:,i] *= signs[i]

    return dwi,mask,bvals,bvecs
    
def remove_nan(image_file):
    img = nib.load(image_file)
    data = img.get_data()
    data[np.isnan(data)] = 0.0
    img = nib.Nifti1Image(data, img.affine, img.header)
    
    nib.save(img, os.path.split(image_file)[0] + "tmp.nii")
    os.remove(image_file)
    shutil.move(os.path.split(image_file)[0] + "tmp.nii", image_file)

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

def read_direction_file(path_to_file):

    file = open(path_to_file, 'r')
    lines = file.readlines()

    directions = []
    for i in range(len(lines)):
        directions.append(lines[i].split())

    for line in range(len(directions)):
        for number in range(len(directions[line])):
            directions[line][number] = float(directions[line][number])

    return directions

def b_to_q(bvals, bvecs, big_delta, little_delta):
    # q = 1/(2*pi) * gyromagnetic_ratio * sqrt(bvals/diffuson_time)

    diffusion_time = big_delta - little_delta / 3.0

    qvals = 1 / (2 * np.pi) * np.sqrt(bvals / diffusion_time)

    # Scale vectors by q-value
    qvectors = np.zeros(bvecs.shape)
    for i in range(bvals.shape[0]):
        qvectors[i,:] = bvecs[i,:] * qvals[i]

    return qvectors

def select_largest_shell(dwi,bval,bvec,tolerance=30):
    largest_bval = np.amax(bval)

    count = 0
    for i in range(len(bval)):
        if np.abs(bval[i]-largest_bval) <= tolerance:
            count += 1

    dwi_large = np.zeros((dwi.shape[0], dwi.shape[1], dwi.shape[2], count))
    bval_large = np.zeros((count))
    bvec_large = np.zeros((count,3))

    index = 0
    for i in range(dwi.shape[3]):
        if np.abs(bval[i]-largest_bval) <= 30:
            dwi_large[:,:,:,index] = dwi[:,:,:,i]
            bval_large[index] = bval[i]
            bvec_large[index,:] = bvec[i,:]
            index += 1

    return dwi_large, bval_large, bvec_large

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
            
def get_bvec_signs(dwi):
    signs = np.zeros(3)
    affine = dwi.affine
    
    for i in range(3):
        index = 0
        max_val = 0
        for j in range(3):
            if np.abs(affine[i,j]) > np.abs(max_val):
                max_val = affine[i,j]
                index = j
                
        if max_val < 0:
            signs[index] = -1
        else:
            signs[index] = 1
            
    return signs    
