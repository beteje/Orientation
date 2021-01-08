#!/usr/bin/env python3
import numpy as np
import cupy as cp
import math, os
from cupyx.scipy import ndimage as cndi
from vispy import gloo
from vispy.color import ColorArray
from matplotlib import pyplot as plt
from generate_fibre_volume import FibreVolume
import plot_orientation as po

#------------------------------------------------------------------------------
# Create simulated data for each of the names corresponding to a set of angle ranges and PSNR 
def simulate_Data(volume_size, n_fibres, elvtn_rng, azth_rng, radius_lim, length_lim, gap_lim, PSNR, fn):
    # create the volume
    volume = FibreVolume(volume_size=volume_size, n_fibres=n_fibres)
    volume.make_volume(elvtn_rng=elvtn_rng, azth_rng=azth_rng, radius_lim=radius_lim, 
                       length_lim=length_lim, gap=gap_lim)
    volume.save_volume(fn + '.hdf5')
    
    # plot slices of original noisy data
    volume.plot_slices(int(min(volume_size)/2))
    
    for i in range(len(PSNR)):
        volume.load_volume(fn + '.hdf5')
        volume.add_noise(PSNR[i])
        volume.save_volume(fn + '_PSNR' + str(PSNR[i]) + '.hdf5')

#------------------------------------------------------------------------------
# Create a canvas for the given colours at the specified locations
def create_Canvas(clrs, X, Y, Z, n, volume_size):
    c = po.Canvas()
    ps = c.pixel_scale
    data = np.zeros(n, [('a_position', np.float32, 3),
                        ('a_color', np.float32, 4),
                        ('a_size', np.float32)])
    data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
    data['a_color'] = clrs
    data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
    c.data_prog.bind(gloo.VertexBuffer(data))
    c.app.run()
    return c

#-----------------------------------------------------------------------------
# Define the colours for each of the identified fibre locations
def set_Colours(elevation, azimuth):
    hsv = np.zeros((n,3), dtype=np.float32)
    hsv[:,0] = 360*azimuth/np.pi
    hsv[:,1] = np.where(elevation/np.pi <=0.5, 2*elevation/np.pi, 1.0)
    hsv[:,2] = np.where(elevation/np.pi <=0.5, 1.0, 1.5-elevation/np.pi)
    clrs = ColorArray(color=hsv, color_space="hsv").rgba
    return clrs

#-----------------------------------------------------------------------------
# Estimate the orientation information (circularity, azimuth angles, elevation angles)
def orientation_est(vol, W_est, W_circ, Thresh_circ, res):
    # Define size and sigma for Gaussian filters
    s = math.floor(W_est/2)             # Size of 3D Gaussian filters (filters are s by s by 2 voxels in size)
    sigma = (s+2)/4                     # Sigma value of the 3D Gaussians (assuming symmetrical Gaussians)
    
    vol = cp.array(vol, dtype=np.float32)
    
    # Compute gradients using filtering:
    Gx = cndi.filters.gaussian_filter(vol, sigma, order=[1,0,0])
    Gy = cndi.filters.gaussian_filter(vol, sigma, order=[0,1,0])
    Gz = cndi.filters.gaussian_filter(vol, sigma, order=[0,0,1])
    
    # Calculate magnitude and angles for gradient assuming it is a 3D vector:
    G = cp.sqrt(abs(Gx)**2 + abs(Gy)**2 + abs(Gz)**2)     # Magnitude
    # Azimuth angle (in radians):
    phi = cp.arctan2(Gy, Gx)
    phi[phi<0] = phi[phi<0] + 2*np.pi
    # Elevation angle (in radians):
    theta = cp.arccos(Gz/G)
    
    # Determine circularity
    circularity = circularity_est(Gx, Gy, Gz, W_circ, Thresh_circ)
    
    # Determine dominant angles
    azimuth, elevation = angle_search(G, circularity, theta, phi, res, W_est)
    
    circularity = cp.asnumpy(circularity)
    azimuth = cp.asnumpy(azimuth)
    elevation = cp.asnumpy(elevation)
    
    return circularity, azimuth, elevation
 
#-----------------------------------------------------------------------------
# Calculate the per voxel circularity 
def circularity_est(Gx, Gy, Gz, W, Thresh):
    # Calculate average region:
    Weights = cp.ones((W,W,W), dtype=np.float32)
    
    # Determine whether each of the planes is circular in local window:
    Gxy = cp.array(Gx + 1j*Gy, dtype=np.complex64)
    Gxz = cp.array(Gx + 1j*Gz, dtype=np.complex64)
    Gyz = cp.array(Gy + 1j*Gz, dtype=np.complex64)
        
    Rxy = cp.zeros(Gxy.shape, dtype=np.float32)
    Rxz = cp.zeros(Gxz.shape, dtype=np.float32)
    Ryz = cp.zeros(Gyz.shape, dtype=np.float32)
    
    Rxy = cndi.convolve(abs(Gxy)**2, Weights)
    Rxz = cndi.convolve(abs(Gxz)**2, Weights)
    Ryz = cndi.convolve(abs(Gyz)**2, Weights)
    
    Gxy2 = Gxy**2
    Gxz2 = Gxz**2
    Gyz2 = Gyz**2
    
    R1xy_r = cndi.convolve(cp.real(Gxy2), Weights)
    R1xz_r = cndi.convolve(cp.real(Gxz2), Weights)
    R1yz_r = cndi.convolve(cp.real(Gyz2), Weights)
    
    R1xy_i = cndi.convolve(cp.imag(Gxy2), Weights)
    R1xz_i = cndi.convolve(cp.imag(Gxz2), Weights)
    R1yz_i = cndi.convolve(cp.imag(Gyz2), Weights)
    
    R1xy = R1xy_r + 1j*R1xy_i
    R1xz = R1xz_r + 1j*R1xz_i
    R1yz = R1yz_r + 1j*R1yz_i
    
    Circularity_XY = abs(R1xy)**2/Rxy**2
    Circularity_XZ = abs(R1xz)**2/Rxz**2
    Circularity_YZ = abs(R1yz)**2/Ryz**2

    circularity = cp.zeros(Gx.shape, dtype=np.int8)
    circularity[(Circularity_XY>Thresh) + (Circularity_XZ>Thresh) + (Circularity_YZ>Thresh)] = 1 
     
    del Weights, Gxy, Gxz, Gyz, Rxy, Rxz, Ryz, R1xy_r, R1xz_r, R1yz_r, R1xy_i, R1xz_i, R1yz_i, R1xy, R1xz, R1yz, Circularity_XY, Circularity_XZ, Circularity_YZ
    cp._default_memory_pool.free_all_blocks()
    
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
def angle_search(G, circularity, theta, phi, res, W):
    # Determine a vector of possible orientations:
    angle_vector2 = np.arange(0, 180, res)
    N_angle2 = angle_vector2.size
    angle_vector1 = np.arange(0, 180, res)
    N_angle1 = angle_vector1.size
    
    # Compute a 2D search    
    # W2 = math.floor(W/2)
    Weights = cp.ones((W,W,W), dtype=np.float32)/(W*W*W)
    circ_nonzero = circularity.nonzero()
    n = len(circ_nonzero[0])
    min_values = cp.ones(n, dtype=np.float32) * 100000
    dominant_theta = -cp.ones(n, dtype=np.uint8)
    dominant_phi = -cp.ones(n, dtype=np.uint8)
    new_values = cp.empty_like(min_values)
    
    # Calculate the individual images for each value of theta:
    for i in range(N_angle1):
        print('Theta Angle: {:d}'.format(angle_vector1[i]))
        theta_est = angle_vector1[i] * np.pi/180
        cos_theta = np.cos(theta_est)
        sin_theta = np.sin(theta_est) 
        for j in range(N_angle2):
            print('\tPhi Angle: {:d}'.format(angle_vector1[j]))
            # Obtain local value of phi and convert to radians:
            phi_est = angle_vector2[j] * np.pi/180
            # Calculate test:
            test_values = G * abs(cos_theta * cp.cos(theta) + sin_theta * cp.sin(theta) * cp.cos(phi_est-phi))
            # Perform local average:
            averaged_values = cndi.convolve(test_values, Weights)
            new_values = averaged_values[circ_nonzero]
            # find the values to replace
            dominant_phi = cp.where(new_values < min_values, phi_est, dominant_phi)
            dominant_theta = cp.where(new_values < min_values, theta_est, dominant_theta)
            min_values = cp.where(new_values < min_values, new_values, min_values)
        
    azimuth = cp.empty_like(circularity, dtype=np.float32)
    azimuth[:] = cp.nan
    azimuth[circ_nonzero] = dominant_phi
    elevation = cp.empty_like(circularity, dtype=np.float32)
    elevation[:] = cp.nan
    elevation[circ_nonzero] = dominant_theta
    
    del test_values, averaged_values, min_values, new_values, dominant_phi, dominant_theta, circ_nonzero 
    cp._default_memory_pool.free_all_blocks()
    
    return azimuth, elevation

#-----------------------------------------------------------------------------
if __name__ == '__main__':
    fld = os.getcwd() + '/Test/'

    # Initialize simulation parameters  
    volume_size = (256,256,256)
    n_fibres = 100
    radius_lim = (2, 10)
    length_lim = (0.3, 0.9)
    gap_lim = 1
    # PSNR = [30, 20, 10]
    PSNR = [30]
    elvtn_rng = [[0, 180], [70, 100], [30, 35], [110, 115]]
    azth_rng = [[0, 180], [90, 120], [70, 75], [150, 155]]
    nme = ['random','range','single','double']
    
    # Initialize estimation parameters
    W_est = 11
    W_circ = 7
    Thresh = 0.9
    res = 1
    
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
    fn = fld + 'Data_' + nme[sim]
    if nme[sim] == 'double':
        volume = simulate_Data(volume_size, n_fibres, elvtn_rng[sim-1:sim+1], azth_rng[sim-1:sim+1], radius_lim, length_lim, gap_lim, PSNR, fn)
    else:
        volume = simulate_Data(volume_size, n_fibres, [elvtn_rng[sim]], [azth_rng[sim]], radius_lim, length_lim, gap_lim, PSNR, fn)
    
    plt.close()
    
    # Plot histograms of angles
    volume = FibreVolume(volume_size=volume_size, n_fibres=n_fibres)
    volume.load_volume(fn + '.hdf5')
    volume.plot_histogram('Elevation', (0, 180), 'darkgreen')
    el_fn = 'Data_Elevation_' + nme[sim] + '.eps'
    plt.savefig(fld + el_fn, transparent=True, bbox_inches='tight')
    plt.close()
    az_fn = 'Data_Azimuth_' + nme[sim] + '.eps'
    volume.plot_histogram('Azimuth', (0, 180), 'darkblue')
    plt.savefig(fld + az_fn, transparent=True, bbox_inches='tight')
    plt.close()
    
    # # Create a 3D canvas of the simulated data
    X, Y, Z = volume.data.nonzero()
    n = len(Z)

    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    sim_clrs = set_Colours(volume.elevation[X, Y, Z], volume.azimuth[X, Y, Z])

    # Create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_orig = create_Canvas(sim_clrs, X, Y, Z, n, volume_size)
    canvas_orig.create_animation(fld + nme[sim] +'_orig.gif')
    canvas_orig.close()
    
    # Calculate orientation angles
    for i in range(len(PSNR)):
        volume.load_volume(fn + '_PSNR' + str(PSNR[i]) + '.hdf5')
        results = FibreVolume(volume_size=volume_size, n_fibres=n_fibres)
        results.circularity, results.azimuth, results.elevation = orientation_est(volume.data, W_est, W_circ, Thresh, res)
        # Determine the volume directionality measures
        results.circularity_XY, results.circularity_YZ, results.circularity_XZ, results.R_length, results.mean_azimuth, results.mean_elevation = directionality_est(results.azimuth, 
                                                                                                                                                                    results.elevation, results.circularity)
        # Save the results
        outFn = 'Results_' + nme[sim] + '_SNR' + str(PSNR[i]) + '.hdf5'
        results.save_volume(fld + outFn)    

        # Plot histograms of angles
        el_fn = 'Results_Elevation_' + nme[sim] + '_SNR' + str(PSNR[i]) + '.eps'
        results.plot_histogram('Elevation', (0, 180), 'darkgreen')
        plt.savefig(fld + el_fn, transparent=True, bbox_inches='tight')
        plt.close()
        az_fn = 'Results_Azimuth_' + nme[sim] + '_SNR' + str(PSNR[i]) + '.eps'
        results.plot_histogram('Azimuth', (0, 180), 'darkblue')
        plt.savefig(fld + az_fn, transparent=True, bbox_inches='tight')
        plt.close()
    
    # Create a 3D canvas of the 30dB PSNR results
    results = FibreVolume()
    results.load_volume(fld + 'Results_' + nme[sim] + '_SNR30.hdf5')
    X, Y, Z = results.circularity.nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    results_clrs = set_Colours(results.elevation[X, Y, Z], results.azimuth[X, Y, Z])
    
    # Create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_results = create_Canvas(results_clrs, X, Y, Z, n, volume_size)
    canvas_results.create_animation(fld + nme[sim] + '_results30.gif')
    canvas_results.close()