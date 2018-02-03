import numpy as np

def track_stats(values):
    with open(values) as f:
        parameter_values = f.readlines()

    # Separate each Track and make it a float
    data_interp = np.zeros((len(parameter_values), 101))
    index = 0
    for tract in parameter_values:
        indiv_points = tract.split()
        indiv_points = [float(i) for i in indiv_points]

        xp = np.linspace(0,100,num=len(indiv_points))
        x = np.linspace(0,100,num=101)
        # Interpolate
        data_interp[index,:] = np.interp(x, xp, indiv_points)
        index += 1

    # Calculate and return mean and std along track
    mean_along_tract = np.mean(data_interp, axis=0)
    std_along_tract = np.std(data_interp, axis=0)

    return mean_along_tract, std_along_tract
