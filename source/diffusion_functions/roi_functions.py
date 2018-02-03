import numpy as np
import nibabel as nib
import os

def calc_roi_stats(param_map_path, lut):
    # Load in Paramter Map
    params = nib.load(param_map_path)
    params = params.get_data()
    params[np.isnan(params)] = 0.0

    # Load in aparc+aseg ROIs and save voxel size
    rois = nib.load("parc.nii")
    voxel_size = rois.header.get_zooms()
    voxel_volume = voxel_size[0] * voxel_size[1] * voxel_size[2]
    rois = rois.get_data()

    # Read in Region Key
    with open(lut) as f:
        lines = f.readlines()

    lines_split = []
    for line in lines:
        lines_split.append(line.split())

    # Create Dictionary
    region_stats = {'region_name': [],
                    'region_mean': np.zeros(len(lines_split)),
                    'region_median': np.zeros(len(lines_split)),
                    'region_min': np.zeros(len(lines_split)),
                    'region_max': np.zeros(len(lines_split)),
                    'region_std': np.zeros(len(lines_split)),
                    'region_volume': np.zeros(len(lines_split))}

    for i in range(len(lines_split)):
        param_values = params[rois == int(lines_split[i][0])]

        region_stats['region_name'].append(lines_split[i][1])

        if len(param_values) != 0:
            region_stats['region_mean'][i] = np.mean(param_values)
            region_stats['region_median'][i] = np.median(param_values)
            region_stats['region_min'][i] = np.median(param_values)
            region_stats['region_max'][i] = np.median(param_values)
            region_stats['region_std'][i] = np.std(param_values)
            region_stats['region_volume'][i] = len(param_values) * voxel_volume

    return region_stats

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
