import numpy as np
import scipy.special
import util

def eval_spherical_harmonics(directions, order):
    # Convert to spherical coordinates
    dirs_sphere = cart_to_sphere(directions)
    num_harmonics = (order + 1) * (order + 2) / 2

    B = np.zeros((directions.shape[0], num_harmonics))

    index = 0
    for L in range(0,order+1,2):
        for m in range(-L, 0):
            scale = np.sqrt((2*L + 1) / (4 * np.pi) * util.factn(L-m, 1) / util.factn(L+m, 1))
            for i in range(directions.shape[0]):
                P1, _ = scipy.special.lpmn(m, L, np.cos(dirs_sphere[i,0]))
                B[i,index] = np.sqrt(2) * np.real(scale * P1[m,L] * np.exp(1j * m * dirs_sphere[i,1]))
            index += 1

        m = 0
        scale = np.sqrt((2*L + 1) / (4 * np.pi) * util.factn(L-m, 1) / util.factn(L+m, 1))
        for i in range(directions.shape[0]):
            P1, _ = scipy.special.lpmn(m, L, np.cos(dirs_sphere[i,0]))
            B[i,index] = np.real(scale * P1[m,L] * np.exp(1j * m * dirs_sphere[i,1]))
        index += 1

        for m in range(1, L+1):
            scale = np.sqrt(2) * np.sqrt((2*L + 1) / (4 * np.pi) * util.factn(L-m, 1) / util.factn(L+m, 1))
            for i in range(directions.shape[0]):
                P1, _ = scipy.special.lpmn(m, L, np.cos(dirs_sphere[i,0]))
                B[i,index] = np.sqrt(2) * np.imag(scale * P1[m,L] * np.exp(1j * m * dirs_sphere[i,1]))
            index += 1

    return B

def fit_to_SH(signal, directions, mask, order, reg=0.006):

    L = calc_normalization_matrix(order)
    B = eval_spherical_harmonics(directions, order)

    # Used for Progress update
    count = 0.0
    percent_prev = 0.0
    num_vox = np.sum(mask)

    # Fit Coeffs for all signal values
    coeffs = np.zeros((signal.shape[0], signal.shape[1], signal.shape[2], B.shape[1]))
    for x in range(signal.shape[0]):
        for y in range(signal.shape[1]):
            for z in range(signal.shape[2]):
                 if mask[x,y,z] != 0:
                     first_term = np.matmul(np.transpose(B), B) + reg * L
                     second_term = np.matmul(np.transpose(B), signal[x,y,z,:])

                     coeffs[x,y,z,:] = np.matmul(np.linalg.inv(first_term), second_term)

                     # Update Progress
                     count += 1.0
                     percent = np.around((count / num_vox * 100), decimals = 1)
                     if(percent != percent_prev):
                         util.progress_update("Fitting Spherical Harmonics: ", percent)
                         percent_prev = percent

    return coeffs


def eval_SH_basis(coeffs, directions):
    pass

def calc_normalization_matrix(order):

    num_harmonics = (order + 1) * (order + 2) / 2

    L = np.zeros((num_harmonics, num_harmonics))

    i = 0
    for l in range(0,order+1,2):
        for m in range(-l,l+1):
            L[i,i] = l ** 2 * (l + 1) ** 2
            i += 1

    return L

def cart_to_sphere(directions):
    # Calc Theta
    theta = np.arccos(directions[:,2])

    # Calc Phi
    phi = np.zeros((directions.shape[0]))

    for i in range(directions.shape[0]):
        ratio = np.absolute(np.divide(directions[i,1], directions[i,0], out=np.zeros_like(directions[i,1]), where=directions[i,0]!=0))
        if ((directions[i,0] >= 0) and (directions[i,1] >= 0)): # 1st Quadrant
            phi[i] = np.arctan(ratio)

        if ((directions[i,0] < 0) and (directions[i,1] >= 0)):  # 2nd Quadrant
            phi[i] = np.pi - np.arctan(ratio)

        if ((directions[i,0] < 0) and (directions[i,1] < 0)):   # 3rd Quadrant
            phi[i] = np.pi + np.arctan(ratio)

        if ((directions[i,0] >= 0) and (directions[i,1] < 0)):  # 4th Quadrant
            phi[i] = (2.0 * np.pi) - np.arctan(ratio)

    sphere_coords = np.zeros((directions.shape[0], 2))
    sphere_coords[:,0] = theta
    sphere_coords[:,1] = phi

    return sphere_coords
