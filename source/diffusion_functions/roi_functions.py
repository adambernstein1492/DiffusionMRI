import numpy as np
import nibabel as nib
import os

def calc_roi_stats(param_map_path):
    # Load in Paramter Map
    params = nib.load(param_map_path)
    params = params.get_data()

    # Load in aparc+aseg ROIs and save voxel size
    rois1 = nib.load("aparc+aseg.nii")
    voxel_size = rois1.header.get_zooms()
    voxel_volume = voxel_size[0] * voxel_size[1] * voxel_size[2]
    rois1 = rois1.get_data()
    rois2 = nib.load("aparc.a2009s+aseg.nii")
    rois2 = rois2.get_data()
    rois3 = nib.load("wmparc.nii")
    rois3 = rois3.get_data()

    # Read in Region Key
    file_location = os.path.dirname(__file__)
    fs_color_lut = file_location + "/config/FreeSurferColorLUT.txt"

    with open(fs_color_lut) as f:
        lines = f.readlines()

    stats_aparc_aseg = [['StructureName', 'Mean', 'Median', 'Min', 'Max', 'STD', 'Number of Voxels', 'Volume (mm^3)']]
    stats_aparc_aseg_2009 = [['StructureName', 'Mean', 'Median', 'Min', 'Max', 'STD', 'Number of Voxels', 'Volume (mm^3)']]
    stats_aparc_aseg_wm = [['StructureName', 'Mean', 'Median', 'Min', 'Max', 'STD', 'Number of Voxels', 'Volume (mm^3)']]
    index_1 = 1
    index_2 = 1
    index_3 = 1
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
                stats_aparc_aseg[index_1].append(len(param_values) * voxel_volume)
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
                stats_aparc_aseg_2009[index_2].append(len(param_values) * voxel_volume)
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
                stats_aparc_aseg_wm[index_3].append(len(param_values) * voxel_volume)
                index_3 += 1

    return stats_aparc_aseg, stats_aparc_aseg_2009, stats_aparc_aseg_wm

def parse_freesurfer_stats(filename, output):
    with open(filename) as f:
        lines = f.readlines()

    stats = [['StructureName', 'NumberVoxels', 'Volume', 'Mean', 'STD', 'Min', 'Max', 'Range']]
    index = 1
    for i in range(len(lines)):
        if lines[i][0] != "#" and lines[i][0] != "\n":
            line_stats = lines[i].split()
            stats.append([])
            stats[index].append(line_stats[4])
            stats[index].append(line_stats[2])
            stats[index].append(line_stats[3])
            stats[index].append(line_stats[5])
            stats[index].append(line_stats[6])
            stats[index].append(line_stats[7])
            stats[index].append(line_stats[8])
            stats[index].append(line_stats[9])
            index += 1

    return stats
    with open(output, "wb") as f:
        writer = csv.writer(f)
        writer.writerows(stats)
