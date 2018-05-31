import os

def scan_for_files(filename, start_directory, files):
    file_and_dir_names = os.listdir(start_directory)

    if (start_directory[-1] != '/'):
        start_directory += '/'

    for i in file_and_dir_names:
        # Look for matching files in current directory
        if ((i == filename) and os.path.isfile(start_directory + i)):
            files.append(os.path.abspath(start_directory + filename))

        # Look in sub-directories
        if os.path.isdir(start_directory + i):
            files = scan_for_files(filename, os.path.abspath(start_directory + i), files)

    return files

def scan_for_nifti(start_directory, files):
    file_and_dir_names = os.listdir(start_directory)

    if (start_directory[-1] != '/'):
        start_directory += '/'

    for i in file_and_dir_names:
        # Look for matching files in current directory
        if ((i[-4:] == '.nii') and os.path.isfile(start_directory + i)):
            files.append(os.path.abspath(start_directory + i))

        # Look in sub-directories
        if os.path.isdir(start_directory + i):
            files = scan_for_nifti(os.path.abspath(start_directory + i), files)

    return files

def scan_for_diffusion_imgs(start_directory, bvals):
    file_and_dir_names = os.listdir(start_directory)

    if (start_directory[-1] != '/'):
        start_directory += '/'

    for i in file_and_dir_names:
        # Look for bval file
        if ((i[-5:] == '.bval') and os.path.isfile(start_directory + i)):
            bvals.append(os.path.abspath(start_directory + i))

        if os.path.isdir(start_directory + i):
            bvals = scan_for_diffusion_imgs((start_directory + i), bvals)

    return bvals

def scan_for_tracts(start_directory, files):
    file_and_dir_names = os.listdir(start_directory)

    if (start_directory[-1] != '/'):
        start_directory += '/'

    for i in file_and_dir_names:
        # Look for matching files in current directory
        if ((i[-4:] == '.tck') and os.path.isfile(start_directory + i)):
            files.append(os.path.abspath(start_directory + i))

        # Look in sub-directories
        if os.path.isdir(start_directory + i):
            files = scan_for_nifti(os.path.abspath(start_directory + i), files)

    return files
