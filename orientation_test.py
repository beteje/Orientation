import numpy as np
import orientation_simulation as sim
from matplotlib import pyplot as plt
from matplotlib.colors import hsv_to_rgb
from vispy import gloo
import plot_orientation as po
import scipy.io as sio

fld = '/media/beth/DATADRIVE1/google-drive/Data & Results/Orientation/'

volume_size = (256, 256, 256)
n_fibres = 100
radius_lim = (2, 10)
length_lim = (0.3, 0.9)
gap_lim = 1
smooth_lvl = 0
PSNR = (30, 20, 10)
elvtn_rng = [[0, 180], [70, 100], [30, 35], [110, 115]]
azth_rng = [[0, 180], [90, 120], [70, 75], [150, 155]]
nme = ('random','range','single','double')

#------------------------------------------------------------------------------
# Create simulated data for each of the names corresponding to a set of angle ranges and PSNR 
def simulate_Data(nme, volume_size, n_fibres, elvtn_rng, azth_rng, radius_lim, length_lim, gap_lim, PSNR):
    for i in range(len(nme)): 
        if nme[i] == 'double':
            volume = sim.make_volume(volume_size, n_fibres, elvtn_rng[i-1:i+1], azth_rng[i-1:i+1], radius_lim, length_lim, gap_lim, PSNR)
        else:
            volume = sim.make_volume(volume_size, n_fibres, [elvtn_rng[i]], [azth_rng[i]], radius_lim, length_lim, gap_lim, PSNR)
            
        fn = 'testData_' + nme[i] + '.mat'
        sio.savemat(fld + fn, volume)
        
        # plot histograms of angles
        po.plot_histogram(volume['elevation'], 'Elevation', (0, 180), 'darkgreen')
        po.plot_histogram(volume['azimuth'], 'Azimuth', (0, 180), 'darkblue')
        
        # plot slices of original and noisy data
        po.plot_slices(volume, PSNR, 128)
        
#------------------------------------------------------------------------------
# Create a canvas for the given colours at the specified locations
def create_Canvas(clrs, X, Y, Z, volume_size):
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
def set_Colours(volume, X, Y, Z, n, volume_size):
    clrs = np.ones((volume_size[0], volume_size[1], volume_size[2], 3), dtype=np.float32)
    clrs_vec = np.zeros((n,3), dtype=np.float32)
    for i, (x, y, z) in enumerate(zip(X, Y, Z)):
        if volume['elevation'][x, y, z] / np.pi <= 0.5:
            clrs[x, y, z] = hsv_to_rgb([volume['azimuth'][x, y, z] / np.pi, 2 * volume['elevation'][x, y, z] / np.pi, 1.0])
            clrs_vec[i] = clrs[x, y, z]
        else:
            clrs[x, y, z] = hsv_to_rgb([volume['azimuth'][x, y, z] / np.pi, 1.0, 1.5 - volume['elevation'][x, y, z] / np.pi])
            clrs_vec[i] = clrs[x, y, z]
    return clrs, clrs_vec

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # simulate_Data(nme, volume_size, n_fibres, elvtn_rng, azth_rng, radius_lim, length_lim, gap_lim, PSNR)
    
    # plot colour wheel
    # po.plot_color_wheel((0,180), (0,180))
    # plt.savefig(fld + 'Colour_Wheel.eps', transparent=True, bbox_inches='tight')

    #-----------------------------------------------------------------------------
    # for i in range(len(PSNR)):
    #     for j in range(len(nme)):
    #         fn = 'Results_' + nme[j] + '_SNR' + str(PSNR[i]) + '.mat'
    #         el_fn = 'Elevation_' + nme[j] + '_SNR' + str(PSNR[i]) + '.eps'
    #         az_fn = 'Azimuth_' + nme[j] + '_SNR' + str(PSNR[i]) + '.eps'
    #         results = sio.loadmat(fld+fn)
    #         X, Y, Z = results['circularity'].nonzero()
    
    #         # plot histograms of angles
    #         po.plot_histogram(results['elevation'][X, Y, Z], 'Elevation', (0, 180), 'darkgreen')
    #         plt.savefig(fld + el_fn, transparent=True, bbox_inches='tight')
    #         po.plot_histogram(results['azimuth'][X, Y, Z], 'Azimuth', (0, 180), 'darkblue')
    #         plt.savefig(fld + az_fn, transparent=True, bbox_inches='tight')

    #-----------------------------------------------------------------------------
    print('Loading Random Test Data')
    volume_rand = sio.loadmat(fld + 'testData_random.mat')
    X, Y, Z = volume_rand['data'].nonzero()
    n = len(Z)

    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    orig_rand_clrs, orig_rand_clrs_vec = set_Colours(volume_rand, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(orig_rand_clrs[128,:,:,:])
    # ax2.imshow(orig_rand_clrs[:,128,:,:])
    # ax3.imshow(orig_rand_clrs[:,:,128,:])

    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_orig_rand = create_Canvas(orig_rand_clrs_vec, X, Y, Z, volume_size)
    canvas_orig_rand.create_animation('rand_orig.gif')

    #-----------------------------------------------------------------------------
    print('Loading Random SNR 30dB Results')
    results_rand_30 = sio.loadmat(fld + 'Results_random_SNR30.mat')
    X, Y, Z = results_rand_30['circularity'].nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    results_rand_30_clrs, results_rand_30_clrs_vec = set_Colours(results_rand_30, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(results_rand_30_clrs[128,:,:,:])
    # ax2.imshow(results_rand_30_clrs[:,128,:,:])
    # ax3.imshow(results_rand_30_clrs[:,:,128,:])
    
    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_results_rand_30 = create_Canvas(results_rand_30_clrs_vec, X, Y, Z, volume_size)
    canvas_results_rand_30.create_animation('rand_results30.gif')

    #-----------------------------------------------------------------------------
    print('Loading Range Test Data')
    volume_range = sio.loadmat(fld + 'testData_range.mat')
    X, Y, Z = volume_range['data'].nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    orig_range_clrs, orig_range_clrs_vec = set_Colours(volume_range, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(orig_range_clrs[128,:,:,:])
    # ax2.imshow(orig_range_clrs[:,128,:,:])
    # ax3.imshow(orig_range_clrs[:,:,128,:])
    
    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_orig_range = create_Canvas(orig_range_clrs_vec, X, Y, Z, volume_size)
    canvas_orig_range.create_animation('range_orig.gif')

    # -----------------------------------------------------------------------------
    print('Loading Range SNR 30 dB Results')
    results_range_30 = sio.loadmat(fld + 'Results_range_SNR30.mat')
    X, Y, Z = results_range_30['circularity'].nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    results_range_30_clrs, results_range_30_clrs_vec = set_Colours(results_range_30, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(results_range_30_clrs[128,:,:,:])
    # ax2.imshow(results_range_30_clrs[:,128,:,:])
    # ax3.imshow(results_range_30_clrs[:,:,128,:])
    
    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_results_range_30 = create_Canvas(results_range_30_clrs_vec, X, Y, Z, volume_size)
    canvas_results_range_30.create_animation('range_results30.gif')

    #-----------------------------------------------------------------------------
    print('Loading Double Test Data')
    volume_double = sio.loadmat(fld + 'testData_double.mat')
    X, Y, Z = volume_double['data'].nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    orig_double_clrs, orig_double_clrs_vec = set_Colours(volume_double, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(orig_double_clrs[128,:,:,:])
    # ax2.imshow(orig_double_clrs[:,128,:,:])
    # ax3.imshow(orig_double_clrs[:,:,128,:])
    
    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_orig_double = create_Canvas(orig_double_clrs_vec, X, Y, Z, volume_size)
    canvas_orig_double.create_animation('double_orig.gif')

    # -----------------------------------------------------------------------------
    print('Loading Double SNR 30 dB Results')
    results_double_30 = sio.loadmat(fld + 'Results_double_SNR30.mat')
    X, Y, Z = results_double_30['circularity'].nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    results_double_30_clrs, results_double_30_clrs_vec = set_Colours(results_double_30, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')   
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(results_double_30_clrs[128,:,:,:])
    # ax2.imshow(results_double_30_clrs[:,128,:,:])
    # ax3.imshow(results_double_30_clrs[:,:,128,:])
    
    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_results_double_30 = create_Canvas(results_double_30_clrs_vec, X, Y, Z, volume_size)
    canvas_results_double_30.create_animation('double_results30.gif')

    #-----------------------------------------------------------------------------
    print('Loading Single Test Data')
    volume_single = sio.loadmat(fld + 'testData_single.mat')
    X, Y, Z = volume_single['data'].nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    orig_single_clrs, orig_single_clrs_vec = set_Colours(volume_single, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')  
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(orig_single_clrs[128,:,:,:])
    # ax2.imshow(orig_single_clrs[:,128,:,:])
    # ax3.imshow(orig_single_clrs[:,:,128,:])
    
    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_orig_single = create_Canvas(orig_single_clrs_vec, X, Y, Z, volume_size)
    canvas_orig_single.create_animation('single_orig.gif')

    #-----------------------------------------------------------------------------
    print('Loading Single SNR 30 dB Results')
    results_single_30 = sio.loadmat(fld + 'Results_single_SNR30.mat')
    X, Y, Z = results_single_30['circularity'].nonzero()
    n = len(Z)
    
    # Calculate the colours of the fibres based on their angles
    print('Calculating Colours')
    results_single_30_clrs, results_single_30_clrs_vec = set_Colours(results_single_30, X, Y, Z, n, volume_size)
    
    # print('Plotting Slices')  
    # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
    # ax1.imshow(results_single_30_clrs[128,:,:,:])
    # ax2.imshow(results_single_30_clrs[:,128,:,:])
    # ax3.imshow(results_single_30_clrs[:,:,128,:])
    
    # create a 3d canvas to display the fibre orientation
    print('Generating Canvas')
    canvas_results_single_30 = create_Canvas(results_single_30_clrs_vec, X, Y, Z, volume_size)
    canvas_results_single_30.create_animation('single_results30.gif')
