import numpy as np

def wave_prep(N, M, lambda_k):
    """
    Accepts number of harmonics, N and M, and min wavelength, lambda_k,
    then returns randomized amplitude, A, and phase shift, phi
    """
    
    A = np.zeros([2*N+1,2*M+1])
    phi = np.zeros([2*N+1,2*M+1])
    for i, n in enumerate(np.arange(-N, N+1)):
        for j, m in enumerate(np.arange(-M, M+1)):
            if ((n**2 + m**2) < (N**2 + 1)):
                A[i][j] = lambda_k / 2 * (n**2 + m**2)/(N**2)
            else:
                A[i][j] = 0
            phi[i][j] = lambda_k * np.random.randn() / 4
    
    return A, phi



def axi_fourier_series(pts, phi, A, lambda_k, sample_points, N, M):
    """
    Accepts x y z points in the pts array, with other
    parameters for the fourier series perturbation process
    """
    
    # convert coordinates from inches to mils
    x = pts[:,0] * 1000
    y = pts[:,1] * 1000
    z = pts[:,2] * 1000
    
    # Set axial boundaries for linear tapering
    bounds = [50, 550, 6250, 6750]
    
    # Map cartesian coordinates into cylindrical coordinate system
    th = np.zeros(sample_points)
    r  = np.zeros(sample_points)
    for i in np.arange(sample_points):
      th[i] = np.arctan2(z[i], y[i])
      r[i]  = np.sqrt(y[i]**2 + z[i]**2)
    
    # Create N_mat and M_mat to avoid a triple do loop
    N_mat = np.zeros([2*N+1, 2*M+1])
    for i, n in enumerate(np.arange(-N,N+1)):
        N_mat[i,:] = n
    
    M_mat = np.zeros([2*N+1, 2*M+1])
    for j, m in enumerate(np.arange(-M,M+1)):
        M_mat[:,j] = m
    
    # Set limits and constants for personalization and speed
    # max_height limits perturbation amplitudes to this quantity
    # trans_point is a location upstream of internal geometry that should not be perturbed (units are mils)
    # min_radius is a location outside of internal geometry that should not be perturbed
    # inv_maxr and inv_lambda are used for speed due to cost of division
    # scale_height is a scaling factor to reduce the double sum
    
    max_height    = 1.25 * lambda_k 
    trans_point   = 5000
    min_radius    = 650
    inv_maxr      = 1 / np.max(r[:])
    inv_maxx      = 1 / bounds[1]
    inv_lambda_th = 1 / (lambda_k * N)
    inv_lambda_x  = 1 / (lambda_k * M)
    
    # Perturbation Loop
    # Perturb all sampled points in the radial direction between bounds
    # Else if statement avoids perturbing internal geometry
    # r(i) * inv_maxr scales perturbations linearly from the min to max radius
    B = np.zeros([2*N+1,2*M+1])
    double_sum = np.zeros(sample_points)
    for i in np.arange(sample_points):
        
        if  bounds[0] < x[i] and x[i] <= bounds[1]:
            B = np.cos((2 * np.pi * N_mat * th[i] * r[i] * inv_lambda_th) + (2 * np.pi * M_mat * x[i] * inv_lambda_x) + phi)
            double_sum[i] = np.sum(A*B) / (2*N+1)
            if (double_sum[i] > max_height):
                double_sum[i] = max_height
            r[i] = r[i] + (double_sum[i] + lambda_k/2) * (x[i] - bounds[0]) / (bounds[1] - bounds[0])
        
        elif bounds[1] < x[i] and x[i] <= trans_point:
            B = np.cos((2 * np.pi * N_mat * th[i] * r[i] * inv_lambda_th) + (2 * np.pi * M_mat * x[i] * inv_lambda_x) + phi)
            double_sum[i] = np.sum(A*B) / (2*N+1)
            if (double_sum[i] > max_height):
                double_sum[i] = max_height
            r[i] = r[i] + (double_sum[i] + lambda_k/2)
        
        elif trans_point < x[i] and x[i] <= bounds[2] and r[i] > min_radius:
            B = np.cos((2 * np.pi * N_mat * th[i] * r[i] * inv_lambda_th) + (2 * np.pi * M_mat * x[i] * inv_lambda_x) + phi)
            double_sum[i] = np.sum(A*B) / (2*N+1)
            if (double_sum[i] > max_height):
                double_sum[i] = max_height
            r[i] = r[i] + (double_sum[i] + lambda_k/2)
        
        elif bounds[2] < x[i] and x[i] <= bounds[3] and r[i] > min_radius:
            B = np.cos((2 * np.pi * N_mat * th[i] * r[i] * inv_lambda_th) + (2 * np.pi * M_mat * x[i] * inv_lambda_x) + phi)
            double_sum[i] = np.sum(A*B) / (2*N+1)
            if (double_sum[i] > max_height):
                double_sum[i] = max_height
            r[i] = r[i] + (double_sum[i] + lambda_k/2) * (bounds[3] - x[i]) / (bounds[3] - bounds[2])
    
    # Convert coordinates back to cartesian
    x[:] = x[:]                 / 1000 
    y[:] = r[:] * np.cos(th[:]) / 1000 
    z[:] = r[:] * np.sin(th[:]) / 1000 
    
    pts_perturbed = np.zeros([sample_points,3])
    pts_perturbed[:,0] = x[:]
    pts_perturbed[:,1] = y[:]
    pts_perturbed[:,2] = z[:]
    
    return pts_perturbed
    
    
    
def cart_fourier_series(pts, phi, A, lambda_k, sample_points, N, M):
    """
    Accepts the pts array, with other parameters for the fourier series 
    perturbation process and returns x y z for the perturbed solution
    """
    
    x = pts[:,0]
    y = pts[:,1]
    
    # Create N_mat and M_mat to avoid a triple do loop
    N_mat = np.zeros([2*N+1, 2*M+1])
    for i, n in enumerate(np.arange(-N,N+1)):
        N_mat[i,:] = n
    
    M_mat = np.zeros([2*N+1, 2*M+1])
    for j, m in enumerate(np.arange(-M,M+1)):
        M_mat[:,j] = m
    
    inv_lambda_x   = 1 / (lambda_k * N)
    inv_lambda_y   = 1 / (lambda_k * M)

    # Perturbation Loop
    B = np.zeros([2*N+1,2*M+1])
    z = np.zeros(sample_points)
    for i in np.arange(sample_points):
        B = np.cos((2 * np.pi * N_mat * x[i] * inv_lambda_x) + (2 * np.pi * M_mat * y[i] * inv_lambda_y) + phi)
        z[i] = np.sum(A*B) / (2*N+1)
    
    # Convert coordinates back to inches
    x[:] = x[:] / 1000
    y[:] = y[:] / 1000
    z[:] = z[:] / 1000

    return x, y, z