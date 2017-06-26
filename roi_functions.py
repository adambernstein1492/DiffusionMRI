import numpy as np
import nibabel as nib

def calc_roi_stats(param_map_path, outpath):
    # Load in Paramter Map
    params = nib.load(param_map_path)
    params = params.get_data()

    # Load in aparc+aseg ROIs
    rois1 = nib.load(freesurfer_directory + outpath + "aparc+aseg.nii")
    rois1 = rois1.get_data()
    rois2 = nib.load(freesurfer_directory + outpath + "aparc.a2009s+aseg.nii")
    rois2 = rois2.get_data()
    rois3 = nib.load(freesurfer_directory + outpath + "wmparc.nii")
    rois3 = rois3.get_data()

    # Read in Region Key
    file_location = os.path.dirname(__file__)
    fs_color_lut = file_location + "/config/FreeSurferColorLUT.txt"

    with open(fs_color_lut) as f:
        lines = f.readlines()

    stats_aparc_aseg = [['StructureName', 'Mean', 'Median', 'Min', 'Max', 'STD', 'Number of Voxels']]
    stats_aparc_aseg_2009 = [['StructureName', 'Mean', 'Median', 'Min', 'Max', 'STD', 'Number of Voxels']]
    stats_aparc_aseg_wm = [['StructureName', 'Mean', 'Median', 'Min', 'Max', 'STD', 'Number of Voxels']]
    index_1 = 1
    index_2 = 2
    index_3 = 3
    for i in range(len(lines)):
        if lines[i][0] != '#' and lines[i][0] != "\n" and lines[i][0] != "0":
            line_values = lines[i].split()

            param_values = params[rois1 == int(line_values[0])]
            if len(param_values) != 0:
                stats_aparc_aseg.append([])
                stats_aparc_aseg[index_1].append(line_values[1])
                stats_aparc_aseg[index_1].append(np.mean(param_values))
                stats_aparc_aseg[index_1].append(np.median(param_values))
                stats_aparc_aseg[index_1].append(np.amax(param_values))
                stats_aparc_aseg[index_1].append(np.amin(param_values))
                stats_aparc_aseg[index_1].append(np.std(param_values))
                stats_aparc_aseg[index_1].append(len(param_values))
                index_1 += 1

            param_values = params[rois2 == int(line_values[0])]
            if len(param_values) != 0:
                stats_aparc_aseg_2009.append([])
                stats_aparc_aseg_2009[index_2].append(line_values[1])
                stats_aparc_aseg_2009[index_2].append(np.mean(param_values))
                stats_aparc_aseg_2009[index_2].append(np.median(param_values))
                stats_aparc_aseg_2009[index_2].append(np.amax(param_values))
                stats_aparc_aseg_2009[index_2].append(np.amin(param_values))
                stats_aparc_aseg_2009[index_2].append(np.std(param_values))
                stats_aparc_aseg_2009[index_2].append(len(param_values))
                index_2 += 1

            param_values = params[rois3 == int(line_values[0])]
            if len(param_values) != 0:
                stats_aparc_aseg_wm.append([])
                stats_aparc_aseg_wm[index_3].append(line_values[1])
                stats_aparc_aseg_wm[index_3].append(np.mean(param_values))
                stats_aparc_aseg_wm[index_3].append(np.median(param_values))
                stats_aparc_aseg_wm[index_3].append(np.amax(param_values))
                stats_aparc_aseg_wm[index_3].append(np.amin(param_values))
                stats_aparc_aseg_wm[index_3].append(np.std(param_values))
                stats_aparc_aseg_wm[index_3].append(len(param_values))
                index_3 += 1

    return stats_aprac_aseg, stats_aparc_aseg_2009, stats_aparc_aseg_wm
