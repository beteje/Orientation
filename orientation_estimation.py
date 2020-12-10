#!/usr/bin/env python3
import numpy as np
import cupy, math, os
from scipy import ndimage as ndi
from numba import cuda
from vispy import gloo
from matplotlib import pyplot as plt
from matplotlib.colors import hsv_to_rgb
from generate_fibre_volume import FibreVolume
import plot_orientation as po

#------------------------------------------------------------------------------
# Create simulated data for each of the names corresponding to a set of angle ranges and PSNR 
def simulate_Data(volume_size, n_fibres, elvtn_rng, azth_rng, radius_lim, length_lim, gap_lim, PSNR):
    # create the volume
    volume = FibreVolume(volume_size=volume_size, n_fibres=n_fibres, PSNR=PSNR)
    volume.make_volume(elvtn_rng=elvtn_rng, azth_rng=azth_rng, radius_lim=radius_lim, 
                       length_lim=length_lim, gap=gap_lim)
    
    # plot slices of original and noisy data
    volume.plot_slices(int(min(volume_size)/2))
    
    return volume

#------------------------------------------------------------------------------
# Create a canvas for the given colours at the specified locations
def create_Canvas(clrs, X, Y, Z, n, volume_size):
    c = po.Canvas()
    ps = c.pixel_scale
    data = np.zeros(n, [('a_position', np.float32, 3),
                        ('a_color', np.float32, 4),
                        ('a_size', np.float32)])
    data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
    data['a_color'] = np.concatenate((clrs, np.ones((n,1))), axis=1)
    data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
    c.data_prog.bind(gloo.VertexBuffer(data))
    c.app.run()
    return c

#-----------------------------------------------------------------------------
# Define the colours for each of the identified fibre locations
def set_Colours(volume, X, Y, Z, n):
    clrs = np.ones((volume.volume_size[0], volume.volume_size[1], volume.volume_size[2], 3), dtype=np.float32)
    clrs_vec = np.zeros((n,3), dtype=np.float32)
    for i, (x, y, z) in enumerate(zip(X, Y, Z)):
        if volume.elevation[x, y, z] / np.pi <= 0.5:
            clrs[x, y, z] = hsv_to_rgb([volume.azimuth[x, y, z] / np.pi, 2 * volume.elevation[x, y, z] / np.pi, 1.0])
            clrs_vec[i] = clrs[x, y, z]
        else:
            clrs[x, y, z] = hsv_to_rgb([volume.azimuth[x, y, z] / np.pi, 1.0, 1.5 - volume.elevation[x, y, z] / np.pi])
            clrs_vec[i] = clrs[x, y, z]
    return clrs, clrs_vec

#-----------------------------------------------------------------------------
# Estimate the orientation information (circularity, azimuth angles, elevation angles)
def orientation_est(vol, W_est, W_circ, Thresh_circ, res, TPB, BPG):
    # Define size and sigma for Gaussian filters
    s = math.floor(W_est/2)             # Size of 3D Gaussian filters (filters are s by s by 2 voxels in size)
    sigma = (s+2)/4                     # Sigma value of the 3D Gaussians (assuming symmetrical Gaussians)
    
    # Compute gradients using filtering:
    Gx = ndi.filters.gaussian_filter(vol, sigma, order=[1,0,0])
    Gy = ndi.filters.gaussian_filter(vol, sigma, order=[0,1,0])
    Gz = ndi.filters.gaussian_filter(vol, sigma, order=[0,0,1])
    
    # Calculate magnitude and angles for gradient assuming it is a 3D vector:
    G = np.sqrt(abs(Gx)**2 + abs(Gy)**2 + abs(Gz)**2)       # Magnitude
    # Azimuth angle (in radians):
    phi = np.arctan2(Gy, Gx)
    phi[phi<0] = phi[phi<0] + 2*np.pi
    # Elevation angle (in radians):
    theta = np.arccos(Gz/G)
    
    # Determine circularity
    circularity = circularity_est(Gx, Gy, Gz, W_circ, Thresh_circ, TPB, BPG)
    
    # Determine dominant angles
    azimuth, elevation = angle_search(G, theta, phi, res, W_est, TPB, BPG)
    
    return circularity, azimuth, elevation
 
#-----------------------------------------------------------------------------
# Calculate the per voxel circularity 
def circularity_est(Gx, Gy, Gz, W, Thresh, TPB, BPG):
    # Determine size of volume:
    [M, N, P] = Gx.shape
    W2 = math.floor(W/2)
    
    # Determine whether each of the planes is circular in local window:
    Gxy = cupy.array(Gx + 1j*Gy)
    Gxz = cupy.array(Gx + 1j*Gz)
    Gyz = cupy.array(Gy + 1j*Gz)
    
    Rxy = cupy.zeros(Gxy.shape, dtype='float')
    Rxz = cupy.zeros(Gxz.shape, dtype='float')
    Ryz = cupy.zeros(Gyz.shape, dtype='float')
    
    local_average[BPG, TPB](abs(Gxy)**2, W, W2, Rxy)
    local_average[BPG, TPB](abs(Gxz)**2, W, W2, Rxz)
    local_average[BPG, TPB](abs(Gyz)**2, W, W2, Ryz)

    R1xy = cupy.zeros(Gxy.shape, dtype='complex')
    R1xz = cupy.zeros(Gxz.shape, dtype='complex')
    R1yz = cupy.zeros(Gyz.shape, dtype='complex')
    
    local_average[BPG, TPB](Gxy**2, W, W2, R1xy)
    local_average[BPG, TPB](Gxz**2, W, W2, R1xz)
    local_average[BPG, TPB](Gyz**2, W, W2, R1yz)
    
    Rxy = cupy.asnumpy(Rxy)
    Rxz = cupy.asnumpy(Rxz)
    Ryz = cupy.asnumpy(Ryz)
    R1xy = cupy.asnumpy(R1xy)
    R1xz = cupy.asnumpy(R1xz)
    R1yz = cupy.asnumpy(R1yz)

    Circularity_XY = abs(R1xy)**2/Rxy**2
    Circularity_XZ = abs(R1xz)**2/Rxz**2
    Circularity_YZ = abs(R1yz)**2/Ryz**2

    circularity = np.zeros((M,N,P))
    circularity[(Circularity_XY>Thresh) + (Circularity_XZ>Thresh) + (Circularity_YZ>Thresh)] = 1 
    
    cupy._default_memory_pool.free_all_blocks()
    
    return circularity

#-----------------------------------------------------------------------------
# Calculate the volume wide directionality measures 
def directionality_est(azimuth, elevation, circularity):
    X, Y, Z = circularity.nonzero()
    
    # convert to radians:
    az = azimuth[X, Y, Z] 
    el = elevation[X, Y, Z] 

    # Generate 3D vector (adjusting for axial data):
    Vx = np.sin(el) * np.cos(2 * az)
    Vy = np.sin(el) * np.sin(2 * az)
    Vz = np.cos(el)

    # Generate complex numbers:
    Vxy = Vx + 1j * Vy
    Vyz = Vy + 1j * Vz
    Vxz = Vx + 1j * Vz

    # Calculate R values:
    Rxy = np.mean(abs(Vxy.flatten())**2)
    Rxz = np.mean(abs(Vyz.flatten())**2)
    Ryz = np.mean(abs(Vxz.flatten())**2)

    # Calculate R1 values:
    R1xy = np.mean(Vxy.flatten()**2)
    R1xz = np.mean(Vyz.flatten()**2)
    R1yz = np.mean(Vxz.flatten()**2)

    # Determine the three circularity numbers:
    circularity_XY = abs(R1xy)**2 / Rxy**2
    circularity_XZ = abs(R1xz)**2 / Rxz**2
    circularity_YZ = abs(R1yz)**2 / Ryz**2

    # Spherical mean:
    x = np.zeros(3)
    x[0] = np.mean(Vx.flatten())
    x[1] = np.mean(Vy.flatten())
    x[2] = np.mean(Vz.flatten())

    # Mean resultant length:
    R_length = np.sqrt(sum(abs(x)**2))

    # Mean direction:
    x_mean = x / R_length
    az_mean = np.arctan2(x_mean[1], x_mean[0])    
    if az_mean<0: az_mean += 2 * np.pi
    az_mean *= 180 / np.pi / 2
    el_mean = np.arccos(x_mean[2]) * 180 / np.pi
    
    return circularity_XY, circularity_YZ, circularity_XZ, R_length, az_mean, el_mean
    

#-----------------------------------------------------------------------------
# Perform a 2D search to find the dominant angles for the azimuth and elevation
def angle_search(G, theta, phi, res, W, TPB, BPG):
    # Determine a vector of possible orientations:
    angle_vector2 = np.concatenate(([0, res/2], np.arange(res, 180, res), [180-res/2]))
    N_angle2 = np.size(angle_vector2)
    angle_vector1 = np.concatenate(([0, res/2], np.arange(res, 180, res), [180-res/2]))
    N_angle1 = np.size(angle_vector1)
    
    # Compute a 2D search
    min_values = cupy.ones(G.shape, dtype='float') * 100000
    dominant_theta = -cupy.ones(G.shape, dtype='float')
    dominant_phi = -cupy.ones(G.shape, dtype='float')
    
    # Calculate the individual images for each value of theta:
    W2 = math.floor(W/2)
    for i in range(N_angle1):
        print('Theta Angle: ' + str(angle_vector1[i]))
        theta_est = angle_vector1[i] * np.pi/180
        for j in range(N_angle2):
            # Obtain local value of phi and convert to radians:
            phi_est = angle_vector2[j] * np.pi/180
            # Calculate test:
            test_values = cupy.array(G * abs(np.cos(theta_est) * np.cos(theta) + np.sin(theta_est) * np.sin(theta) * np.cos(phi_est-phi)))
            return_values = cupy.zeros(G.shape, dtype='float')
            # Perform local average:
            local_average[BPG, TPB](test_values, W, W2, return_values)
            # find the values to replace
            replace_Values[BPG, TPB](min_values, return_values, dominant_theta, dominant_phi, theta_est, phi_est)
    
    azimuth = cupy.asnumpy(dominant_phi);
    elevation = cupy.asnumpy(dominant_theta);
    cupy._default_memory_pool.free_all_blocks()
    
    return azimuth, elevation

#-----------------------------------------------------------------------------
# Cuda function to determine local 3D averages
@cuda.jit
def local_average(vol, W, W2, output):
    i, j, k = cuda.grid(3)
    x_size, y_size, z_size = vol.shape
    
    # start of x
    if (i < W2):
        # start of y
        if (j < W2):
            # start of z
            if (k < W2):
                for x in range(-i, W2+1):
                    for y in range(-j, W2+1):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*(W2+1+j)*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-i, W2+1):
                    for y in range(-j, W2+1):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*(W2+1+j)*W)
            # end of z
            elif (k < z_size):
                for x in range(-i, W2+1):
                    for y in range(-j, W2+1):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*(W2+1+j)*(W2+z_size-k))
        # middle of y
        elif (j <= y_size - (W2+1)):
            # start of z
            if (k < W2):
                for x in range(-i, W2+1):
                    for y in range(-W2, W2+1):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*W*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-i, W2+1):
                    for y in range(-W2, W2+1):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*W*W)
            # end of z
            elif (k < z_size):
                for x in range(-i, W2+1):
                    for y in range(-W2, W2+1):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*W*(W2+z_size-k))
        # end of y
        elif (j < y_size):
            # start of z
            if (k < W2):
                for x in range(-i, W2+1):
                    for y in range(-W2, y_size-j):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*(W2+y_size-j)*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-i, W2+1):
                    for y in range(-W2, y_size-j):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*(W2+y_size-j)*W)
            # end of z
            elif (k < z_size):
                for x in range(-i, W2+1):
                    for y in range(-W2, y_size-j):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+1+i)*(W2+y_size-j)*(W2+z_size-k))
    # middle of x
    elif (i <= x_size - (W2+1)):
        # start of y
        if (j < W2):
            # start of z
            if (k < W2):
                for x in range(-W2, W2+1):
                    for y in range(-j, W2+1):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*(W2+1+j)*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-W2, W2+1):
                    for y in range(-j, W2+1):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*(W2+1+j)*W)
            # end of z
            elif (k < z_size):
                for x in range(-W2, W2+1):
                    for y in range(-j, W2+1):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*(W2+1+j)*(W2+z_size-k))
        # middle of y
        elif (j <= y_size - (W2+1)):
            # start of z
            if (k < W2):
                for x in range(-W2, W2+1):
                    for y in range(-W2, W2+1):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*W*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-W2, W2+1):
                    for y in range(-W2, W2+1):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*W*W)
            # end of z
            elif (k < z_size):
                for x in range(-W2, W2+1):
                    for y in range(-W2, W2+1):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*W*(W2+z_size-k))
        # end of y
        elif (j < y_size):
            # start of z
            if (k < W2):
                for x in range(-W2, W2+1):
                    for y in range(-W2, y_size-j):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*(W2+y_size-j)*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-W2, W2+1):
                    for y in range(-W2, y_size-j):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*(W2+y_size-j)*W)
            # end of z
            elif (k < z_size):
                for x in range(-W2, W2+1):
                    for y in range(-W2, y_size-j):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / (W*(W2+y_size-j)*(W2+z_size-k))
    # end of x
    elif (i < x_size):
        # start of y
        if (j < W2):
            # start of z
            if (k < W2):
                for x in range(-W2, x_size-i):
                    for y in range(-j, W2+1):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*(W2+1+j)*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-W2, x_size-i):
                    for y in range(-j, W2+1):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*(W2+1+j)*W)
            # end of z
            elif (k < z_size):
                for x in range(-W2, x_size-i):
                    for y in range(-j, W2+1):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*(W2+1+j)*(W2+z_size-k))
        # middle of y
        elif (j <= y_size - (W2+1)):
            # start of z
            if (k < W2):
                for x in range(-W2, x_size-i):
                    for y in range(-W2, W2+1):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*W*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-W2, x_size-i):
                    for y in range(-W2, W2+1):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*W*W)
            # end of z
            elif (k < z_size):
                for x in range(-W2, x_size-i):
                    for y in range(-W2, W2+1):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*W*(W2+z_size-k))
        # end of y
        elif (j < y_size):
            # start of z
            if (k < W2):
                for x in range(-W2, x_size-i):
                    for y in range(-W2, y_size-j):
                        for z in range(-k, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*(W2+y_size-j)*(W2+1+k))
            # middle of z
            elif (k <= z_size - (W2+1)):
                for x in range(-W2, x_size-i):
                    for y in range(-W2, y_size-j):
                        for z in range(-W2, W2+1):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*(W2+y_size-j)*W)
            # end of z
            elif (k < z_size):
                for x in range(-W2, x_size-i):
                    for y in range(-W2, y_size-j):
                        for z in range(-W2, z_size-k):
                            output[i, j, k] += vol[i+x, j+y, k+z]
                output[i, j, k] = output[i, j, k] / ((W2+x_size-i)*(W2+y_size-j)*(W2+z_size-k))

#-----------------------------------------------------------------------------
# Cuda function to find the values in the volume to replace with the updated values
@cuda.jit
def replace_Values(orig, test, theta, phi, theta_est, phi_est):
    i, j, k = cuda.grid(3)
    # Find values where the test data is less than the minimum data
    if test[i, j, k] < orig[i, j, k]:
        orig[i, j, k] = test[i, j, k]
        phi[i, j, k] = phi_est
        theta[i, j, k] = theta_est

if __name__ == '__main__':
    fld = os.getcwd() + '/Test/'

    # Initialize simulation parameters  
    # volume_size = (256, 256, 256)
    volume_size = (64, 64, 64)
    # n_fibres = 100
    n_fibres = 10
    radius_lim = (2, 10)
    length_lim = (0.3, 0.9)
    gap_lim = 1
    PSNR = [30, 20, 10]
    elvtn_rng = [[0, 180], [70, 100], [30, 35], [110, 115]]
    azth_rng = [[0, 180], [90, 120], [70, 75], [150, 155]]
    nme = ['random','range','single','double']
    
    # Initialize estimation parameters
    W_est = 11
    W_circ = 7
    Thresh = 0.9
    res = 5
    
    # Initialize Cuda parameters
    threadsperblock = (8, 8, 8)
    blockspergrid_x = math.ceil(volume_size[0] / threadsperblock[0])
    blockspergrid_y = math.ceil(volume_size[1] / threadsperblock[1])
    blockspergrid_z = math.ceil(volume_size[2] / threadsperblock[2])
    blockspergrid = (blockspergrid_x, blockspergrid_y, blockspergrid_z)
    
    # Choice of simulation
    print('Select a simulation:')
    print('[1] Random elevation & azimuth angles')
    print('[2] Range of elevation angles ' + str(elvtn_rng[1]) + ' & azimuth angles ' + str(azth_rng[1]))
    print('[3] Single elevation angle ' + str(elvtn_rng[2]) + ' & azimuth angle ' + str(azth_rng[2]))
    print('[4] Two elevation angles ' + str(elvtn_rng[2]) + ' & ' + str(elvtn_rng[3]) + ' & azimuth angles ' + str(azth_rng[2]) + ' & ' + str(azth_rng[3]))
    inputValid = False
    while not inputValid:
        inputRaw = input('Simulation: ')
        inputNo = int(inputRaw) - 1
        if inputNo > -1 and inputNo < len(nme):
            sim = inputNo
            print('Selected simulation: ' + nme[inputNo])
            inputValid = True
            break
        else:
            print('Please select a valid simulation number')
    
    # Simulate the data
    if nme[sim] == 'double':
        volume = simulate_Data(volume_size, n_fibres, elvtn_rng[sim-1:sim+1], azth_rng[sim-1:sim+1], radius_lim, length_lim, gap_lim, PSNR)
    else:
        volume = simulate_Data(volume_size, n_fibres, [elvtn_rng[sim]], [azth_rng[sim]], radius_lim, length_lim, gap_lim, PSNR)
    
    fn = 'Data_' + nme[sim] + '.hdf5'
    volume.save_volume(fld + fn)  
    
    # Plot histograms of angles
    volume.plot_histogram('Elevation', (0, 180), 'darkgreen')
    el_fn = 'Data_Elevation_' + nme[sim] + '.eps'
    plt.savefig(fld + el_fn, transparent=True, bbox_inches='tight')
    az_fn = 'Data_Azimuth_' + nme[sim] + '.eps'
    volume.plot_histogram('Azimuth', (0, 180), 'darkblue')
    plt.savefig(fld + az_fn, transparent=True, bbox_inches='tight')
    
    # Calculate orientation angles
    for i in range(len(PSNR)):
        results = FibreVolume()
        results.circularity, results.azimuth, results.elevation = orientation_est(volume.noisy_data['PSNR'+str(PSNR[i])], 
                                                                                  W_est, W_circ, Thresh, res, threadsperblock, blockspergrid)
        # Determine the volume directionality measures
        results.circularity_XY, results.circularity_YZ, results.circularity_XZ, results.R_length, results.azimuth_mean, results.elevation_mean = directionality_est(results.azimuth, 
                                                                                                                                                                    results.elevation, results.circularity)
        
        # Save the results
        outFn = 'Results_' + nme[sim] + '_SNR' + str(PSNR[i]) + '.hdf5'
        results.save_volume(fld + outFn)    

        # Plot histograms of angles
        el_fn = 'Results_Elevation_' + nme[sim] + '_SNR' + str(PSNR[i]) + '.eps'
        results.plot_histogram('Elevation', (0, 180), 'darkgreen')
        plt.savefig(fld + el_fn, transparent=True, bbox_inches='tight')
        az_fn = 'Results_Azimuth_' + nme[sim] + '_SNR' + str(PSNR[i]) + '.eps'
        results.plot_histogram('Azimuth', (0, 180), 'darkblue')
        plt.savefig(fld + az_fn, transparent=True, bbox_inches='tight')
    
    # Create a 3D canvas of the simulated data
    X, Y, Z = volume.data.nonzero()
    n = len(Z)

    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    sim_clrs, sim_clrs_vec = set_Colours(volume, X, Y, Z, n)

    # Create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_orig = create_Canvas(sim_clrs_vec, X, Y, Z, n, volume_size)
    canvas_orig.create_animation(fld + 'rand_orig.gif')
    
    # Create a 3D canvas of the 30dB PSNR results
    results = FibreVolume()
    results.load_volume(fld + 'Results_' + nme[sim] + '_SNR30.hdf5')
    X, Y, Z = results.circularity.nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    results_clrs, results_clrs_vec = set_Colours(results, X, Y, Z, n)
    
    # Create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_results = create_Canvas(results_clrs_vec, X, Y, Z, n, volume_size)
    canvas_results.create_animation(fld + nme[sim] + '_results30.gif')