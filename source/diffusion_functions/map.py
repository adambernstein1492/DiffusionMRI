import numpy as np
import nibabel as nib
import scipy.special
import util
import dti


def main_map(dwi_file, bval_file, bvec_file, mask_file, little_delta, big_delta,
             out_path, order=6, b_thresh_dti=2100, calc_rtps=True, calc_ng=True,
             calc_pa=True, calc_dki=False, return_dti=False):

    # Load in Data
    dwi, mask, bvals, bvecs = util.load_diffusion_data(dwi_file, bval_file, bvec_file, mask_file)
    data = dwi.get_data()
    mask = mask.get_data()

    # Fit DTI
    if return_dti:
        eigen_values, eigen_vectors = dti.main_dti(dwi_file, bval_file, bvec_file,
                            mask_file, (out_path + "DTI_"), b_thresh_dti, True, True, True, True)
        eigen_values[eigen_values <= 0] = 1e-5
    else:
        eigen_values, eigen_vectors = dti.main_dti(dwi_file, bval_file, bvec_file,
                            mask_file, "", b_thresh_dti, False, False, False, False)
        eigen_values[eigen_values <= 0] = 1e-5

    # Determine Diffusion Time
    diffusion_time = big_delta - little_delta / 3

    # Calculate u-vectors
    uvectors = np.sqrt(2 * eigen_values * diffusion_time)

    # Convert b values and vectors to q-vectors
    qvectors = b_to_q(bvals, bvecs, diffusion_time)

    # Invert Eigenvectors (Inverse of a rotation matrix is its transpose)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                eigen_vectors[i,j,k,:,:] = eigen_vectors[i,j,k,:,:].T

    # Fit MAP
    coeffs = fit_map(data, qvectors, mask, diffusion_time, uvectors, eigen_vectors,
                     order, lam=0.2)

    num_coeffs = int(round(1.0/6 * (order/2 + 1) * (order/2 + 2) * (2*order + 3)))
    uvectors = np.sort(uvectors, axis=3)

    # Save Coefficients and UVectors
    img = nib.Nifti1Image(coeffs, dwi.affine, dwi.header)
    nib.save(img, (out_path + 'coeffs.nii'))
    img = nib.Nifti1Image(uvectors, dwi.affine, dwi.header)
    nib.save(img, (out_path + 'uvecs.nii'))

    # Calculate and Save Scalar MAPs
    if calc_pa:
        print "Calculating Propagator Anisotropy"
        pa_dti, pa_dti_theta, pa, pa_theta = calc_propagator_anisotropy(coeffs, uvectors, order, mask)

        img = nib.Nifti1Image(pa_dti, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'PA_DTI.nii'))

        img = nib.Nifti1Image(pa_dti_theta, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'PA_DTI_theta.nii'))

        img = nib.Nifti1Image(pa, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'PA.nii'))

        img = nib.Nifti1Image(pa_theta, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'PA_theta.nii'))

    if calc_ng:
        print "Calculating Non-Gaussianity"
        ng, ng_par, ng_perp = calc_non_gaussianity(coeffs, order, num_coeffs)

        img = nib.Nifti1Image(ng, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'NonGaussianity.nii'))

        img = nib.Nifti1Image(ng_par, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'AxialNonGaussianity.nii'))

        img = nib.Nifti1Image(ng_perp, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'PerpindicularNonGaussianity.nii'))

    if calc_rtps:
        print "Calculating Return origin, axis, and plane probabilities"
        rtop, rtap, rtpp = calc_return_to_probabilities(coeffs, uvectors, order, num_coeffs)

        img = nib.Nifti1Image(rtop, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'RTOP.nii'))

        img = nib.Nifti1Image(rtap, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'RTAP.nii'))

        img = nib.Nifti1Image(rtpp, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'RTPP.nii'))

    if calc_dki:
        pass

    if return_dti:
        img = nib.Nifti1Image(eigen_vectors, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'dti_eigen_vectors.nii'))
        img = nib.Nifti1Image(eigen_values, dwi.affine, dwi.header)
        nib.save(img, (out_path + 'dti_eigen_values.nii'))

def fit_map(data, qvectors, mask, diffusion_time, uvectors, eigen_vectors,
            order = 6, lam = 0.2):

    # Calculate the number of coefficients for a give order
    num_coeffs = int(round(1.0/6 * (order/2 + 1) * (order/2 + 2) * (2*order + 3)))

    # Allocate space for Outputs
    data[data <= 0] = np.finfo(float).eps
    coeffs = np.zeros((data.shape[0], data.shape[1], data.shape[2], num_coeffs))

    # Calculate "B" for estimating B0
    B = calc_b(order, num_coeffs)

    # Calculate Coefficients for Regularization Matrix
    reg_matrix_coeffs = calc_reg_matrix_coeffs(order, num_coeffs)

    count = 0.0
    percent_prev = 0.0
    num_vox = np.sum(mask)

    for i in range(data.shape[0]):
        for j in range (data.shape[1]):
            for k in range(data.shape[2]):
                if mask[i,j,k] != 0:
                    # Rotate qvectors
                    qvec = np.dot(eigen_vectors[i,j,k,:,:], qvectors.T).T

                    # Evaluate Gaussian Hermite Basis Functions
                    Q = gaussian_hermite_q_space(order, num_coeffs, qvec, uvectors[i,j,k,:])

                    # Create Regularization for least squares fit
                    R = create_regularization_matrix(order, num_coeffs, uvectors[i,j,k,:], reg_matrix_coeffs)

                    # Perform Linear Fit
                    coeffs[i,j,k,:] = (np.linalg.lstsq((np.dot(Q.T, Q) + lam * R),
                                                   (np.dot(Q.T, data[i,j,k,:])))[0])

                    # Normalize coefficients by estimated B0
                    EstimatedB0 = np.dot(coeffs[i,j,k,:], B)
                    coeffs[i,j,k,:] /= EstimatedB0

                    # Update Progress
                    count += 1.0
                    percent = np.around((count / num_vox * 100), decimals = 1)
                    if(percent != percent_prev):
                        util.progress_update("Fitting MAP: ", percent)
                        percent_prev = percent

    return coeffs

def gaussian_hermite_q_space(order, num_coeffs, qvec, uvector):
    # Allocate Space for output
    Q = np.zeros((qvec.shape[0], num_coeffs), dtype=complex)
    H = np.zeros((qvec.shape[0], order+1, 3), dtype=complex)

    # Points to evaluate hermite polynomials at
    x1 = 2 * np.pi * uvector[0] * qvec[:,0]
    x2 = 2 * np.pi * uvector[1] * qvec[:,1]
    x3 = 2 * np.pi * uvector[2] * qvec[:,2]

    # Evaluate separated basis functions at each order
    for i in range(order+1):
        coeffs = np.zeros(order+1)
        coeffs[i] = 1

        H[:,i,0] = (1j**(-i) / np.sqrt(2**i * np.math.factorial(i)) *
                    np.exp(-(x1**2)/2) * np.polynomial.hermite.hermval(x1, coeffs))
        H[:,i,1] = (1j**(-i) / np.sqrt(2**i * np.math.factorial(i)) *
                    np.exp(-(x2**2)/2) * np.polynomial.hermite.hermval(x2, coeffs))
        H[:,i,2] = (1j**(-i) / np.sqrt(2**i * np.math.factorial(i)) *
                    np.exp(-(x3**2)/2) * np.polynomial.hermite.hermval(x3, coeffs))

    # Multiply separable functions together to get 3-D function
    index = 0
    for N in range(0,order+1,2):
        for n1 in range(N+1):
            for n2 in range(N+1):
                for n3 in range(N+1):
                    if((n1+n2+n3) == N):
                        Q[:,index] = H[:,n1,0] * H[:,n2,1] * H[:,n3,2]
                        index += 1

    return np.real(Q)

def calc_b(order, num_coeffs):
    # Allocate Space
    B = np.zeros(num_coeffs)

    index = 0
    for N in range(0,order+1,2):
        for n1 in range(order+1):
            for n2 in range(order+1):
                for n3 in range(order+1):
                    if((n1+n2+n3) == N):
                        if((n1 % 2) == 0):
                            B1 = np.sqrt(util.factn(n1,1)) / util.factn(n1,2)
                            B2 = np.sqrt(util.factn(n2,1)) / util.factn(n2,2)
                            B3 = np.sqrt(util.factn(n3,1)) / util.factn(n3,2)

                            B[index] = B1 * B2 * B3
                            index += 1

    return B

def b_to_q(bvals, bvecs, diffusion_time):
    # q = 1/(2*pi) * gyromagnetic_ratio * sqrt(bvals/diffuson_time)

    qvals = 1 / (2 * np.pi) * np.sqrt(bvals / diffusion_time)

    # Scale vectors by q-value
    qvectors = np.zeros(bvecs.shape)
    for i in range(bvals.shape[0]):
        qvectors[i,:] = bvecs[i,:] * qvals[i]

    return qvectors

def create_regularization_matrix(order, num_coeffs, uvector, reg_matrix_coeffs):
    R = np.zeros((num_coeffs, num_coeffs))
    scale = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    ux = uvector[0]
    uy = uvector[1]
    uz = uvector[2]

    scale[0] = ux ** 3 / (uy * uz)
    scale[1] = 2 * ux * uy / uz
    scale[2] = uy ** 3 / (uz * ux)
    scale[3] = 2 * uy * uz / ux
    scale[4] = uz ** 3 / (ux * uy)
    scale[5] = 2 * ux * uz / uy

    for i in range(6):
        R = R + scale[i] * reg_matrix_coeffs[:,:,i]

    return R

def calc_reg_matrix_coeffs(order, num_coeffs):
    reg_matrix_coeffs = np.zeros((num_coeffs, num_coeffs, 6))

    U = calc_u_matrix(order)
    T = calc_t_matrix(order)
    S = calc_s_matrix(order)

    indexN = 0
    for N in range(0,order+1,2):
        for n1 in range(order+1):
            for n2 in range(order+1):
                for n3 in range(order+1):
                    if((n1+n2+n3) == N):
                        indexM = 0

                        for M in range(0,order+1,2):
                            for m1 in range(order+1):
                                for m2 in range(order+1):
                                    for m3 in range(order+1):
                                        if((m1+m2+m3) == M):
                                            reg_matrix_coeffs[indexN, indexM, 0] = S[n1,m1] * U[n2,m2] * U[n3,m3]
                                            reg_matrix_coeffs[indexN, indexM, 1] = T[n1,m1] * T[n2,m2] * U[n3,m3]
                                            reg_matrix_coeffs[indexN, indexM, 2] = U[n1,m1] * S[n2,m2] * U[n3,m3]
                                            reg_matrix_coeffs[indexN, indexM, 3] = U[n1,m1] * T[n2,m2] * T[n3,m3]
                                            reg_matrix_coeffs[indexN, indexM, 4] = U[n1,m1] * U[n2,m2] * S[n3,m3]
                                            reg_matrix_coeffs[indexN, indexM, 5] = T[n1,m1] * U[n2,m2] * T[n3,m3]

                                            indexM += 1
                        indexN += 1

    return reg_matrix_coeffs

def calc_u_matrix(order):
    U = np.zeros((order+1, order+1))

    for n in range(order+1):
        for m in range(order+1):
            if(m == n):
                U[n,m] = (-1) ** n / (2 * np.sqrt(np.pi))

    return U

def calc_t_matrix(order):
    T = np.zeros((order+1, order+1))
    mn_factors = [0.0, 0.0, 0.0]

    for n in range(order+1):
        for m in range(order+1):
            if(m == n):
                mn_factors[0] = (1+2*n)
            else:
                mn_factors[0] = 0

            if((m+2) == n):
                mn_factors[1] = np.sqrt((n * (n-1)))
            else:
                mn_factors[1] = 0

            if((m == (n+2))):
                mn_factors[2] = np.sqrt((m * (m-1)))
            else:
                mn_factors[2] = 0

            T[n,m] = (-1) ** (n+1) * np.pi ** (1.5) * np.sum(mn_factors)

    return T

def calc_s_matrix(order):
    S = np.zeros((order+1, order+1))
    mn_factors = [0.0, 0.0, 0.0, 0.0, 0.0]

    for n in range(order+1):
        for m in range(order+1):
            if(m == n):
                mn_factors[0] = 3 * (2*n**2 + 2*n + 1)
            else:
                mn_factors[0] = 0

            if(m == (n+2)):
                mn_factors[1] = (6+4*n) * np.sqrt((util.factn(m,1) / util.factn(n,1)))
            else:
                mn_factors[1] = 0

            if(m == (n+4)):
                mn_factors[2] = np.sqrt((util.factn(m,1) / util.factn(n,1)))
            else:
                mn_factors[2] = 0

            if((m+2) == n):
                mn_factors[3] = (6+4*m) * np.sqrt((util.factn(n,1) / util.factn(m,1)))
            else:
                mn_factors[3] = 0

            if((m+4) == n):
                mn_factors[4] = np.sqrt((util.factn(n,1) / util.factn(m,1)))
            else:
                mn_factors[4] = 0

            S[n,m] = 2*(-1)**n * np.pi**3.5 * np.sum(mn_factors)

    return S

def calc_return_to_probabilities(coeffs, uvectors, order, num_coeffs):
    # Calculate sign factors
    b, signs = calc_b_and_signs(order, num_coeffs)
    b = b[0,:] * b[1,:] * b[2,:]

    # Reshape Data for Faster Processing
    uvectors_lin = np.reshape(uvectors, (uvectors.shape[0] * uvectors.shape[1] * uvectors.shape[2], uvectors.shape[3]))
    coeffs_lin = np.reshape(coeffs, (coeffs.shape[0] * coeffs.shape[1] * coeffs.shape[2], coeffs.shape[3]))

    # Calculate Coefficients
    rtop_scale = (8 * np.pi**3 * uvectors_lin[:,0]**2 * uvectors_lin[:,1]**2 * uvectors_lin[:,2]**2) ** (-0.5)
    rtap_scale_x = (4 * np.pi**2 * uvectors_lin[:,1]**2 * uvectors_lin[:,2]**2) ** (-0.5)
    rtap_scale_y = (4 * np.pi**2 * uvectors_lin[:,0]**2 * uvectors_lin[:,2]**2) ** (-0.5)
    rtap_scale_z = (4 * np.pi**2 * uvectors_lin[:,0]**2 * uvectors_lin[:,1]**2) ** (-0.5)
    rtpp_scale_x = (2 * np.pi * uvectors_lin[:,0]**2)**(-0.5)
    rtpp_scale_y = (2 * np.pi * uvectors_lin[:,1]**2)**(-0.5)
    rtpp_scale_z = (2 * np.pi * uvectors_lin[:,2]**2)**(-0.5)

    # Estimate Probabilites
    rtop = np.zeros(coeffs_lin.shape[0])
    rtap_x = np.zeros(coeffs_lin.shape[0])
    rtap_y = np.zeros(coeffs_lin.shape[0])
    rtap_z = np.zeros(coeffs_lin.shape[0])
    rtpp_x = np.zeros(coeffs_lin.shape[0])
    rtpp_y = np.zeros(coeffs_lin.shape[0])
    rtpp_z = np.zeros(coeffs_lin.shape[0])
    rtap = np.zeros((coeffs.shape[0], coeffs.shape[1], coeffs.shape[2], 3))
    rtpp = np.zeros((coeffs.shape[0], coeffs.shape[1], coeffs.shape[2], 3))

    rtop = rtop_scale * np.dot(signs[0,:]*signs[1,:]*signs[2,:]*b, coeffs_lin.T)
    rtap_x = rtap_scale_x * np.dot(signs[1,:]*signs[2,:]*b, coeffs_lin.T)
    rtap_y = rtap_scale_y * np.dot(signs[0,:]*signs[2,:]*b, coeffs_lin.T)
    rtap_z = rtap_scale_z * np.dot(signs[0,:]*signs[1,:]*b, coeffs_lin.T)
    rtpp_x = rtpp_scale_x * np.dot(signs[0,:]*b, coeffs_lin.T)
    rtpp_y = rtpp_scale_y * np.dot(signs[1,:]*b, coeffs_lin.T)
    rtpp_z = rtpp_scale_z * np.dot(signs[2,:]*b, coeffs_lin.T)


    # Filter Data for Crazy Outliers
    rtop[rtop < 0] = 0
    rtap_x[rtap_x < 0] = 0
    rtap_y[rtap_y < 0] = 0
    rtap_z[rtap_z < 0] = 0
    rtpp_x[rtpp_x < 0] = 0
    rtpp_y[rtpp_y < 0] = 0
    rtpp_z[rtpp_z < 0] = 0

    rtop = rtop ** (1.0/3)
    rtap_x = rtap_x ** 0.5
    rtap_y = rtap_y ** 0.5
    rtap_z = rtap_z ** 0.5

    rtop = np.reshape(rtop, coeffs.shape[0:3])
    rtap[:,:,:,0] = np.reshape(rtap_x, coeffs.shape[0:3])
    rtap[:,:,:,1] = np.reshape(rtap_y, coeffs.shape[0:3])
    rtap[:,:,:,2] = np.reshape(rtap_z, coeffs.shape[0:3])
    rtpp[:,:,:,0] = np.reshape(rtpp_x, coeffs.shape[0:3])
    rtpp[:,:,:,1] = np.reshape(rtpp_y, coeffs.shape[0:3])
    rtpp[:,:,:,2] = np.reshape(rtpp_z, coeffs.shape[0:3])

    rtop = np.real(rtop)
    rtap = np.real(rtap)
    rtpp = np.real(rtpp)

    rtop[np.isnan(rtop)] = 0
    rtap[np.isnan(rtap)] = 0
    rtpp[np.isnan(rtpp)] = 0

    return rtop, rtap, rtpp


def calc_non_gaussianity(coeffs, order, num_coeffs):
    # Calculate sign factors
    b, signs = calc_b_and_signs(order, num_coeffs)

    #Reshape for Faster Processing
    coeffs_lin = np.reshape(coeffs, (coeffs.shape[0] * coeffs.shape[1] * coeffs.shape[2], num_coeffs))
    a1 = np.zeros((coeffs_lin.shape[0], order+1))
    a2 = np.zeros((coeffs_lin.shape[0], order+1))
    a3 = np.zeros((coeffs_lin.shape[0], order+1))
    a23 = np.zeros((coeffs_lin.shape[0], (order+1)**2))
    a13 = np.zeros((coeffs_lin.shape[0], (order+1)**2))
    a12 = np.zeros((coeffs_lin.shape[0], (order+1)**2))
    perp_23 = np.zeros((coeffs_lin.shape[0], num_coeffs))
    perp_13 = np.zeros((coeffs_lin.shape[0], num_coeffs))
    perp_12 = np.zeros((coeffs_lin.shape[0], num_coeffs))

    ng = np.zeros(coeffs_lin.shape[0])
    ng_par_x = np.zeros(coeffs_lin.shape[0])
    ng_par_y = np.zeros(coeffs_lin.shape[0])
    ng_par_z = np.zeros(coeffs_lin.shape[0])
    ng_perp_x = np.zeros(coeffs_lin.shape[0])
    ng_perp_y = np.zeros(coeffs_lin.shape[0])
    ng_perp_z = np.zeros(coeffs_lin.shape[0])

    ng_par = np.zeros((coeffs.shape[0], coeffs.shape[1], coeffs.shape[2], 3))
    ng_perp = np.zeros((coeffs.shape[0], coeffs.shape[1], coeffs.shape[2], 3))

    # Parallel Coefficients
    index = 0
    for N in range(0,order+1,2):
        for n1 in range(N+1):
            for n2 in range(N+1):
                for n3 in range(N+1):
                    if((n1+n2+n3) == N):
                        a1[:,n1] += signs[1,index] * signs[2,index] * b[1,index] * b[2,index] * coeffs_lin[:,index]
                        a2[:,n2] += signs[0,index] * signs[2,index] * b[0,index] * b[2,index] * coeffs_lin[:,index]
                        a3[:,n3] += signs[0,index] * signs[1,index] * b[0,index] * b[1,index] * coeffs_lin[:,index]
                        index += 1

    # Perpindicular Coefficients
    n1_index = []
    n2_index = []
    n3_index = []
    index = 0
    for N in range(0,order+1,2):
        for n1 in range(N+1):
            for n2 in range(N+1):
                for n3 in range(N+1):
                    if((n1+n2+n3) == N):
                        perp_23[:,index] = signs[0,index] * b[0,index] * coeffs_lin[:,index]
                        perp_13[:,index] = signs[1,index] * b[1,index] * coeffs_lin[:,index]
                        perp_12[:,index] = signs[2,index] * b[2,index] * coeffs_lin[:,index]
                        n1_index.append(n1);
                        n2_index.append(n2);
                        n3_index.append(n3);
                        index += 1

    index = 0
    for N in range(order+1):
        for M in range(order+1):

            for i in range(num_coeffs):
                if (n2_index[i] == N and n3_index[i] == M):
                    a23[:,index] += perp_23[:,i]

            for i in range(num_coeffs):
                if (n1_index[i] == N and n3_index[i] == M):
                    a13[:,index] += perp_13[:,i]

            for i in range(num_coeffs):
                if (n1_index[i] == N and n2_index[i] == M):
                    a12[:,index] += perp_12[:,i]

            index += 1


    # Calculate Nongaussianities
    ng = np.sqrt(1 - coeffs_lin[:,0]**2 / np.sum(coeffs_lin**2, axis=1))
    ng_par_x = np.sqrt(1 - a1[:,0]**2 / np.sum(a1**2, axis = 1))
    ng_par_y = np.sqrt(1 - a2[:,0]**2 / np.sum(a2**2, axis = 1))
    ng_par_z = np.sqrt(1 - a3[:,0]**2 / np.sum(a3**2, axis = 1))
    ng_perp_x = np.sqrt(1 - a23[:,0]**2 / np.sum(a23**2, axis = 1))
    ng_perp_y = np.sqrt(1 - a13[:,0]**2 / np.sum(a13**2, axis = 1))
    ng_perp_z = np.sqrt(1 - a12[:,0]**2 / np.sum(a12**2, axis = 1))

    # Change Back to Image Matrix
    ng = np.reshape(ng, coeffs.shape[0:3])
    ng_par[:,:,:,0] = np.reshape(ng_par_x, coeffs.shape[0:3])
    ng_par[:,:,:,1] = np.reshape(ng_par_y, coeffs.shape[0:3])
    ng_par[:,:,:,2] = np.reshape(ng_par_z, coeffs.shape[0:3])
    ng_perp[:,:,:,0] = np.reshape(ng_perp_x, coeffs.shape[0:3])
    ng_perp[:,:,:,1] = np.reshape(ng_perp_y, coeffs.shape[0:3])
    ng_perp[:,:,:,2] = np.reshape(ng_perp_z, coeffs.shape[0:3])

    # Clean up and div by zeros
    ng[np.isnan(ng)] = 0
    ng_par[np.isnan(ng_par)] = 0
    ng_perp[np.isnan(ng_perp)] = 0

    return ng, ng_par, ng_perp


def calc_b_and_signs(order, num_coeffs):
    b = np.zeros((3,num_coeffs))
    signs = np.zeros((3,num_coeffs))

    index = 0
    for N in range(0,order+1,2):
        for n1 in range(N+1):
            for n2 in range(N+1):
                for n3 in range(N+1):
                    if((n1+n2+n3) == N):
                        if((n1 % 2) == 0):
                            b[0,index] = np.sqrt(util.factn(n1,1)) / util.factn(n1,2)
                            signs[0,index] = (-1)**(n1/2)
                        if((n2 % 2) == 0):
                            b[1,index] = np.sqrt(util.factn(n2,1)) / util.factn(n2,2)
                            signs[1,index] = (-1)**(n2/2)
                        if((n3 % 2) == 0):
                            b[2,index] = np.sqrt(util.factn(n3,1)) / util.factn(n3,2)
                            signs[2,index] = (-1)**(n3/2)

                        index += 1

    return b, signs


def calc_propagator_anisotropy(coeffs, uvectors, order, mask):
    # Calculate Isotropic Coefficients
    u0 = calc_u0(uvectors)

    # Caluclate the Anisotropy of the DTI signal
    pa_dti, pa_dti_theta = calculate_pa_dti(u0, uvectors, mask)

    # Calculate the General Propagator Anisotropy
    pa = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2]))
    pa_theta = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2]))
    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            for z in range(mask.shape[2]):
                if(mask[x,y,z] != 0):
                    # Estimate Transformation Matrix
                    tmn = calc_tmn(uvectors[x,y,z,:], u0[x,y,z], order, order/2)

                    # Calculate PA
                    pa[x,y,z], pa_theta[x,y,z] = calculate_pa(coeffs[x,y,z,:], uvectors[x,y,z,:], u0[x,y,z], tmn)

    # Apply scaling function
    pa_dti = pa_scaling(pa_dti, 0.4)
    pa = pa_scaling(pa, 0.4)

    return pa_dti, pa_dti_theta, pa, pa_theta


def calc_u0(uvectors):
    # Compute factors for polynomial
    uvectors = uvectors ** 2
    u3 = -3
    u2 = -np.sum(uvectors, axis=3)
    u1 = (uvectors[:,:,:,0] * uvectors[:,:,:,1] + uvectors[:,:,:,0] * uvectors[:,:,:,2] +
          uvectors[:,:,:,1] * uvectors[:,:,:,2])
    u = 3 * uvectors[:,:,:,0] * uvectors[:,:,:,1] * uvectors[:,:,:,2]

    # Allocate Space for Output
    u0 = np.zeros((uvectors.shape[0], uvectors.shape[1], uvectors.shape[2]))
    for x in range(uvectors.shape[0]):
        for y in range(uvectors.shape[1]):
            for z in range(uvectors.shape[2]):
                poly = np.array([u3, u2[x,y,z], u1[x,y,z], u[x,y,z]])

                # Find roots of polynomial
                r = np.roots(poly)

                # Find the single real, positive root
                for i in r:
                    if(i >= 0 and np.isreal(i)):
                        u0[x,y,z] = np.sqrt(np.real(i))
                        break

    return u0


def calculate_pa_dti(u0, uvectors, mask):
    pa_dti = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2]))
    pa_dti_theta = np.zeros((mask.shape[0], mask.shape[1], mask.shape[2]))

    for x in range(u0.shape[0]):
        for y in range(u0.shape[1]):
            for z in range(u0.shape[2]):
                if(mask[x,y,z] != 0):
                    pa_dti[x,y,z] = ((8 * u0[x,y,z]**3 * uvectors[x,y,z,0] * uvectors[x,y,z,1] * uvectors[x,y,z,2]) /
                                     ((uvectors[x,y,z,0]**2 + u0[x,y,z]**2) * (uvectors[x,y,z,1]**2 + u0[x,y,z]**2) * (uvectors[x,y,z,2]**2 + u0[x,y,z]**2)))

                    pa_dti_theta[x,y,z] = np.arccos(np.sqrt(pa_dti[x,y,z]))
                    pa_dti[x,y,z] = np.real(np.sin(pa_dti_theta[x,y,z]))


    return pa_dti, pa_dti_theta


def calculate_pa(coeffs, uvectors, u0, tmn):

    b, norm_factor = change_basis(coeffs, tmn, u0)

    anorm =  1 / (8 * np.pi**(1.5) * uvectors[0] * uvectors[1] * uvectors[2])

    cosPA = np.sqrt(np.dot(b**2, norm_factor) / (np.sum(coeffs**2) * anorm))

    pa_theta = np.arccos(cosPA)
    pa = np.sin(pa_theta)

    return pa, pa_theta


def change_basis(coeffs, tmn, u0):
    norm_factor = np.zeros((tmn.shape[1]))

    for i in range(tmn.shape[1]):
        norm_factor[i] = scipy.special.gamma(2*i + 1.5) / (util.factn(2*i, 1) * 4 * np.pi**2 * u0**3)

    b = np.dot(coeffs,tmn) / norm_factor

    return b, norm_factor


def calc_tmn(uvectors, u0, order_a, order_b):
    tmn = np.zeros((int(np.round(1/6.0 * (order_a/2 + 1) * (order_a/2 + 2) * (2*order_a + 3))), 1 + order_a / 2))

    tx = np.zeros((order_a+1, order_b+1))
    ty = np.zeros((order_a+1, order_b+1))
    tz = np.zeros((order_a+1, order_b+1))

    psi_x = uvectors[0] / u0
    psi_y = uvectors[1] / u0
    psi_z = uvectors[2] / u0

    for m in range(order_a+1):
        for n in range(0,2*order_b+1,2):
            if((m+n)%2 == 0):
                coeff = (np.sqrt(util.factn(m,1) * util.factn(n,1)) / np.sqrt(2 * np.pi) *
                        np.sqrt(util.factn(n,1)) / (util.factn(n,2) * u0))

                for r in range(0,m+1,2):
                    for s in range(0,n+1,2):
                        tx[m,n/2] += (1 + psi_x**2) ** (-(m+n-r-s+1)/2.0) * psi_x**(n-s) * (-1)**((r+s)/2.0) * 2**((m+n-r-s)/2.0) * util.factn(m+n-r-s-1,2) / (util.factn(m-r,1) * util.factn(n-s,1) * util.factn(r,2) * util.factn(s,2))
                        ty[m,n/2] += (1 + psi_y**2) ** (-(m+n-r-s+1)/2.0) * psi_y**(n-s) * (-1)**((r+s)/2.0) * 2**((m+n-r-s)/2.0) * util.factn(m+n-r-s-1,2) / (util.factn(m-r,1) * util.factn(n-s,1) * util.factn(r,2) * util.factn(s,2))
                        tz[m,n/2] += (1 + psi_z**2) ** (-(m+n-r-s+1)/2.0) * psi_z**(n-s) * (-1)**((r+s)/2.0) * 2**((m+n-r-s)/2.0) * util.factn(m+n-r-s-1,2) / (util.factn(m-r,1) * util.factn(n-s,1) * util.factn(r,2) * util.factn(s,2))

                tx[m,n/2] *= coeff
                ty[m,n/2] *= coeff
                tz[m,n/2] *= coeff


    m_index = 0
    for M in range(0,order_a+1,2):
        for m1 in range(M+1):
            for m2 in range(M+1):
                for m3 in range(M+1):
                    if((m1+m2+m3) == M):
                        for N in range(order_b+1):
                            for n1 in range(N+1):
                                for n2 in range(N+1):
                                    for n3 in range(N+1):
                                        if((n1+n2+n3) == N):
                                            tmn[m_index, N] += tx[m1, n1] * ty[m2, n2] * tz[m3,n3]
                        m_index += 1

    return tmn


def pa_scaling(pa, param):
    pa = pa**(3*param) / (1 - 3*pa**param + 3*pa**(2*param))

    return pa
