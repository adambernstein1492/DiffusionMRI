import numpy as np
import nibabel as nib
import os

def calc_roi_stats(param_map_path, gm_lut, wm_lut):
    # Load in Paramter Map
    params = nib.load(param_map_path)
    params = params.get_data()
    params[np.isnan(params)] = 0.0

    # Load in aparc+aseg ROIs and save voxel size
    gm_rois = nib.load("GM_parc_reg.nii")
    voxel_size = gm_rois.header.get_zooms()
    voxel_volume = voxel_size[0] * voxel_size[1] * voxel_size[2]
    gm_rois = gm_rois.get_data()

    wm_rois = nib.load("WM_parc_reg.nii")
    wm_rois = wm_rois.get_data()

    # Read in Region Key
    with open(gm_lut) as f:
        lines = f.readlines()

    gm_regions = []
    for line in lines:
        gm_regions.append(line.split())

    with open(wm_lut) as f:
        lines = f.readlines()

    wm_regions = []
    for line in lines:
        wm_regions.append(line.split())

    # Create Dictionary
    region_stats = {'region_name': [],
                    'region_mean': np.zeros(len(gm_regions) + len(wm_regions)),
                    'region_median': np.zeros(len(gm_regions) + len(wm_regions)),
                    'region_min': np.zeros(len(gm_regions) + len(wm_regions)),
                    'region_max': np.zeros(len(gm_regions) + len(wm_regions)),
                    'region_std': np.zeros(len(gm_regions) + len(wm_regions)),
                    'region_volume': np.zeros(len(gm_regions) + len(wm_regions))}

    for i in range(len(gm_regions) + len(wm_regions)):
        if i < len(gm_regions):
            param_values = params[gm_rois == int(gm_regions[i][0])]

            region_stats['region_name'].append(gm_regions[i][1])

            if len(param_values) != 0:
                region_stats['region_mean'][i] = np.mean(param_values)
                region_stats['region_median'][i] = np.median(param_values)
                region_stats['region_min'][i] = np.amin(param_values)
                region_stats['region_max'][i] = np.amax(param_values)
                region_stats['region_std'][i] = np.std(param_values)
                region_stats['region_volume'][i] = len(param_values) * voxel_volume

        if i >= len(gm_regions):
            j = i - len(gm_regions)

            param_values = params[wm_rois == int(wm_regions[j][0])]

            region_stats['region_name'].append(wm_regions[j][1])

            if len(param_values) != 0:
                region_stats['region_mean'][i] = np.mean(param_values)
                region_stats['region_median'][i] = np.median(param_values)
                region_stats['region_min'][i] = np.amin(param_values)
                region_stats['region_max'][i] = np.amax(param_values)
                region_stats['region_std'][i] = np.std(param_values)
                region_stats['region_volume'][i] = len(param_values) * voxel_volume

    return region_stats
