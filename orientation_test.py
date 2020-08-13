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
intersect = False
gap_lim = 1
smooth_lvl = 0
PSNR = (30, 20, 10)
elvtn_rng = [[0, 180], [70, 100], [30, 35], [110, 115]]
azth_rng = [[0, 180], [90, 120], [70, 75], [150, 155]]
nme = ('random','range','single','double')

# for i in range(len(nme)): 
#     if nme[i] == 'double':
#         volume = sim.make_volume(volume_size, n_fibres, elvtn_rng[i-1:i+1], azth_rng[i-1:i+1], radius_lim, length_lim, gap_lim, intersect, PSNR, smooth_lvl)
#     else:
#         volume = sim.make_volume(volume_size, n_fibres, [elvtn_rng[i]], [azth_rng[i]], radius_lim, length_lim, gap_lim, intersect, PSNR, smooth_lvl)
        
#     fn = 'testData_' + nme[i] + '.mat'
#     sio.savemat(fn, volume)
    
#     # plot histograms of angles
#     po.plot_histogram(volume['elvtn'], 'Elevation', (0, 180), 'darkgreen')
#     po.plot_histogram(volume['azth'], 'Azimuth', (0, 180), 'darkblue')
    
#     # plot slices of original and noisy data
#     po.plot_slices(volume, PSNR, 128)

#------------------------------------------------------------------------------
# plot colour wheel
# po.plot_color_wheel((0,180), (0,180))
plt.savefig(fld+'Colour_Wheel.eps', transparent=True, bbox_inches='tight')

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
#         plt.savefig(fld+el_fn, transparent=True, bbox_inches='tight')
#         po.plot_histogram(results['azimuth'][X, Y, Z], 'Azimuth', (0, 180), 'darkblue')
#         plt.savefig(fld+az_fn, transparent=True, bbox_inches='tight')

#-----------------------------------------------------------------------------
# volume_rand = sio.loadmat(fld + 'testData_random.mat')
# X, Y, Z = volume_rand['data'].nonzero()
# n = len(Z)
# orig_rand_clrs = np.ones((256,256,256,3), dtype=np.float32)
# orig_rand_clrs_vec = np.zeros((n,3), dtype=np.float32)

# for i, (x, y, z) in enumerate(zip(X, Y, Z)):
#     if volume_rand['elvtn'][x, y, z] / np.pi <= 0.5:
#         orig_rand_clrs[x, y, z] = hsv_to_rgb([volume_rand['azth'][x, y, z] / np.pi, 2 * volume_rand['elvtn'][x, y, z] / np.pi, 1.0])
#         orig_rand_clrs_vec[i] = orig_rand_clrs[x, y, z]
#     else:
#         orig_rand_clrs[x, y, z] = hsv_to_rgb([volume_rand['azth'][x, y, z] / np.pi, 1.0, 1.5 - volume_rand['elvtn'][x, y, z] / np.pi])
#         orig_rand_clrs_vec[i] = orig_rand_clrs[x, y, z]
    
# # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# # ax1.imshow(orig_rand_clrs[128,:,:,:])
# # ax2.imshow(orig_rand_clrs[:,128,:,:])
# # ax3.imshow(orig_rand_clrs[:,:,128,:])

# # create a 3d canvas to display the fibre orientation
# canvas_orig_rand = po.Canvas()
# ps = canvas_orig_rand.pixel_scale
# orig_rand_data = np.zeros(n, [('a_position', np.float32, 3),
#                     ('a_color', np.float32, 4),
#                     ('a_size', np.float32)])
# orig_rand_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
# orig_rand_data['a_color'] = np.concatenate((orig_rand_clrs_vec, np.ones((n,1))), axis=1)
# orig_rand_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
# canvas_orig_rand.data_prog.bind(gloo.VertexBuffer(orig_rand_data))
# canvas_orig_rand.app.run()

# -----------------------------------------------------------------------------
# results_rand_30 = sio.loadmat(fld + 'Results_random_SNR30.mat')
# X, Y, Z = results_rand_30['circularity'].nonzero()
# n = len(Z)
# results_rand_30_clrs = np.ones((256,256,256,3), dtype=np.float32)
# results_rand_30_clrs_vec = np.zeros((n,3), dtype=np.float32)

# for i, (x, y, z) in enumerate(zip(X, Y, Z)):
#     if results_rand_30['elevation'][x, y, z] / np.pi <= 0.5:
#         results_rand_30_clrs[x, y, z] = hsv_to_rgb([results_rand_30['azimuth'][x, y, z] / np.pi, 2 * results_rand_30['elevation'][x, y, z] / np.pi, 1.0])
#         results_rand_30_clrs_vec[i] = results_rand_30_clrs[x, y, z]
#     else:
#         results_rand_30_clrs[x, y, z] = hsv_to_rgb([results_rand_30['azimuth'][x, y, z] / np.pi, 1.0, 1.5 - results_rand_30['elevation'][x, y, z] / np.pi])
#         results_rand_30_clrs_vec[i] = results_rand_30_clrs[x, y, z]
    
# # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# # ax1.imshow(results_rand_30_clrs[128,:,:,:])
# # ax2.imshow(results_rand_30_clrs[:,128,:,:])
# # ax3.imshow(results_rand_30_clrs[:,:,128,:])

# # create a 3d canvas to display the fibre orientation
# canvas_results_rand_30 = po.Canvas()
# ps = canvas_results_rand_30.pixel_scale
# results_rand_30_data = np.zeros(n, [('a_position', np.float32, 3),
#                     ('a_color', np.float32, 4),
#                     ('a_size', np.float32)])
# results_rand_30_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
# results_rand_30_data['a_color'] = np.concatenate((results_rand_30_clrs_vec, np.ones((n,1))), axis=1)
# results_rand_30_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
# canvas_results_rand_30.data_prog.bind(gloo.VertexBuffer(results_rand_30_data))
# canvas_results_rand_30.app.run()

#-----------------------------------------------------------------------------
# volume_range = sio.loadmat(fld + 'testData_range.mat')
# X, Y, Z = volume_range['data'].nonzero()
# n = len(Z)
# orig_range_clrs = np.ones((256,256,256,3), dtype=np.float32)
# orig_range_clrs_vec = np.zeros((n,3), dtype=np.float32)

# for i, (x, y, z) in enumerate(zip(X, Y, Z)):
#     if volume_range['elvtn'][x, y, z] / np.pi <= 0.5:
#         orig_range_clrs[x, y, z] = hsv_to_rgb([volume_range['azth'][x, y, z] / np.pi, 2 * volume_range['elvtn'][x, y, z] / np.pi, 1.0])
#         orig_range_clrs_vec[i] = orig_range_clrs[x, y, z]
#     else:
#         orig_range_clrs[x, y, z] = hsv_to_rgb([volume_range['azth'][x, y, z] / np.pi, 1.0, 1.5 - volume_range['elvtn'][x, y, z] / np.pi])
#         orig_range_clrs_vec[i] = orig_range_clrs[x, y, z]
    
# # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# # ax1.imshow(orig_range_clrs[128,:,:,:])
# # ax2.imshow(orig_range_clrs[:,128,:,:])
# # ax3.imshow(orig_range_clrs[:,:,128,:])

# # create a 3d canvas to display the fibre orientation
# canvas_orig_range = po.Canvas()
# ps = canvas_orig_range.pixel_scale
# orig_range_data = np.zeros(n, [('a_position', np.float32, 3),
#                     ('a_color', np.float32, 4),
#                     ('a_size', np.float32)])
# orig_range_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
# orig_range_data['a_color'] = np.concatenate((orig_range_clrs_vec, np.ones((n,1))), axis=1)
# orig_range_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
# canvas_orig_range.data_prog.bind(gloo.VertexBuffer(orig_range_data))
# canvas_orig_range.app.run()

# -----------------------------------------------------------------------------
# results_range_30 = sio.loadmat(fld + 'Results_range_SNR30.mat')
# X, Y, Z = results_range_30['circularity'].nonzero()
# n = len(Z)
# results_range_30_clrs = np.ones((256,256,256,3), dtype=np.float32)
# results_range_30_clrs_vec = np.zeros((n,3), dtype=np.float32)

# for i, (x, y, z) in enumerate(zip(X, Y, Z)):
#     if results_range_30['elevation'][x, y, z] / np.pi <= 0.5:
#         results_range_30_clrs[x, y, z] = hsv_to_rgb([results_range_30['azimuth'][x, y, z] / np.pi, 2 * results_range_30['elevation'][x, y, z] / np.pi, 1.0])
#         results_range_30_clrs_vec[i] = results_range_30_clrs[x, y, z]
#     else:
#         results_range_30_clrs[x, y, z] = hsv_to_rgb([results_range_30['azimuth'][x, y, z] / np.pi, 1.0, 1.5 - results_range_30['elevation'][x, y, z] / np.pi])
#         results_range_30_clrs_vec[i] = results_range_30_clrs[x, y, z]
    
# # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# # ax1.imshow(results_range_30_clrs[128,:,:,:])
# # ax2.imshow(results_range_30_clrs[:,128,:,:])
# # ax3.imshow(results_range_30_clrs[:,:,128,:])

# # create a 3d canvas to display the fibre orientation
# canvas_results_range_30 = po.Canvas()
# ps = canvas_results_range_30.pixel_scale
# results_range_30_data = np.zeros(n, [('a_position', np.float32, 3),
#                     ('a_color', np.float32, 4),
#                     ('a_size', np.float32)])
# results_range_30_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
# results_range_30_data['a_color'] = np.concatenate((results_range_30_clrs_vec, np.ones((n,1))), axis=1)
# results_range_30_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
# canvas_results_range_30.data_prog.bind(gloo.VertexBuffer(results_range_30_data))
# canvas_results_range_30.app.run()

#-----------------------------------------------------------------------------
# volume_double = sio.loadmat(fld + 'testData_double.mat')
# X, Y, Z = volume_double['data'].nonzero()
# n = len(Z)
# orig_double_clrs = np.ones((256,256,256,3), dtype=np.float32)
# orig_double_clrs_vec = np.zeros((n,3), dtype=np.float32)

# for i, (x, y, z) in enumerate(zip(X, Y, Z)):
#     if volume_double['elvtn'][x, y, z] / np.pi <= 0.5:
#         orig_double_clrs[x, y, z] = hsv_to_rgb([volume_double['azth'][x, y, z] / np.pi, 2 * volume_double['elvtn'][x, y, z] / np.pi, 1.0])
#         orig_double_clrs_vec[i] = orig_double_clrs[x, y, z]
#     else:
#         orig_double_clrs[x, y, z] = hsv_to_rgb([volume_double['azth'][x, y, z] / np.pi, 1.0, 1.5 - volume_double['elvtn'][x, y, z] / np.pi])
#         orig_double_clrs_vec[i] = orig_double_clrs[x, y, z]
    
# # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# # ax1.imshow(orig_double_clrs[128,:,:,:])
# # ax2.imshow(orig_double_clrs[:,128,:,:])
# # ax3.imshow(orig_double_clrs[:,:,128,:])

# # create a 3d canvas to display the fibre orientation
# canvas_orig_double = po.Canvas()
# ps = canvas_orig_double.pixel_scale
# orig_double_data = np.zeros(n, [('a_position', np.float32, 3),
#                     ('a_color', np.float32, 4),
#                     ('a_size', np.float32)])
# orig_double_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
# orig_double_data['a_color'] = np.concatenate((orig_double_clrs_vec, np.ones((n,1))), axis=1)
# orig_double_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
# canvas_orig_double.data_prog.bind(gloo.VertexBuffer(orig_double_data))
# canvas_orig_double.app.run()

# -----------------------------------------------------------------------------
# results_double_30 = sio.loadmat(fld + 'Results_double_SNR30.mat')
# X, Y, Z = results_double_30['circularity'].nonzero()
# n = len(Z)
# results_double_30_clrs = np.ones((256,256,256,3), dtype=np.float32)
# results_double_30_clrs_vec = np.zeros((n,3), dtype=np.float32)

# for i, (x, y, z) in enumerate(zip(X, Y, Z)):
#     if results_double_30['elevation'][x, y, z] / np.pi <= 0.5:
#         results_double_30_clrs[x, y, z] = hsv_to_rgb([results_double_30['azimuth'][x, y, z] / np.pi, 2 * results_double_30['elevation'][x, y, z] / np.pi, 1.0])
#         results_double_30_clrs_vec[i] = results_double_30_clrs[x, y, z]
#     else:
#         results_double_30_clrs[x, y, z] = hsv_to_rgb([results_double_30['azimuth'][x, y, z] / np.pi, 1.0, 1.5 - results_double_30['elevation'][x, y, z] / np.pi])
#         results_double_30_clrs_vec[i] = results_double_30_clrs[x, y, z]
    
# # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# # ax1.imshow(results_double_30_clrs[128,:,:,:])
# # ax2.imshow(results_double_30_clrs[:,128,:,:])
# # ax3.imshow(results_double_30_clrs[:,:,128,:])

# # create a 3d canvas to display the fibre orientation
# canvas_results_double_30 = po.Canvas()
# ps = canvas_results_double_30.pixel_scale
# results_double_30_data = np.zeros(n, [('a_position', np.float32, 3),
#                     ('a_color', np.float32, 4),
#                     ('a_size', np.float32)])
# results_double_30_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
# results_double_30_data['a_color'] = np.concatenate((results_double_30_clrs_vec, np.ones((n,1))), axis=1)
# results_double_30_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
# canvas_results_double_30.data_prog.bind(gloo.VertexBuffer(results_double_30_data)) 158, 68, 495
# canvas_results_double_30.app.run()

#-----------------------------------------------------------------------------
# volume_single = sio.loadmat(fld + 'testData_single.mat')
# X, Y, Z = volume_single['data'].nonzero()
# n = len(Z)
# orig_single_clrs = np.ones((256,256,256,3), dtype=np.float32)
# orig_single_clrs_vec = np.zeros((n,3), dtype=np.float32)

# for i, (x, y, z) in enumerate(zip(X, Y, Z)):
#     if volume_single['elvtn'][x, y, z] / np.pi <= 0.5:
#         orig_single_clrs[x, y, z] = hsv_to_rgb([volume_single['azth'][x, y, z] / np.pi, 2 * volume_single['elvtn'][x, y, z] / np.pi, 1.0])
#         orig_single_clrs_vec[i] = orig_single_clrs[x, y, z]
#     else:
#         orig_single_clrs[x, y, z] = hsv_to_rgb([volume_single['azth'][x, y, z] / np.pi, 1.0, 1.5 - volume_single['elvtn'][x, y, z] / np.pi])
#         orig_single_clrs_vec[i] = orig_single_clrs[x, y, z]
    
# # fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# # ax1.imshow(orig_single_clrs[128,:,:,:])
# # ax2.imshow(orig_single_clrs[:,128,:,:])
# # ax3.imshow(orig_single_clrs[:,:,128,:])

# # create a 3d canvas to display the fibre orientation
# canvas_orig_single = po.Canvas()
# ps = canvas_orig_single.pixel_scale
# orig_single_data = np.zeros(n, [('a_position', np.float32, 3),
#                     ('a_color', np.float32, 4),
#                     ('a_size', np.float32)])
# orig_single_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
# orig_single_data['a_color'] = np.concatenate((orig_single_clrs_vec, np.ones((n,1))), axis=1)
# orig_single_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
# canvas_orig_single.data_prog.bind(gloo.VertexBuffer(orig_single_data))
# canvas_orig_single.app.run()

#-----------------------------------------------------------------------------
results_single_30 = sio.loadmat(fld + 'Results_single_SNR30.mat')
X, Y, Z = results_single_30['circularity'].nonzero()
n = len(Z)
results_single_30_clrs = np.ones((256,256,256,3), dtype=np.float32)
results_single_30_clrs_vec = np.zeros((n,3), dtype=np.float32)

for i, (x, y, z) in enumerate(zip(X, Y, Z)):
    if results_single_30['elevation'][x, y, z] / np.pi <= 0.5:
        results_single_30_clrs[x, y, z] = hsv_to_rgb([results_single_30['azimuth'][x, y, z] / np.pi, 2 * results_single_30['elevation'][x, y, z] / np.pi, 1.0])
        results_single_30_clrs_vec[i] = results_single_30_clrs[x, y, z]
    else:
        results_single_30_clrs[x, y, z] = hsv_to_rgb([results_single_30['azimuth'][x, y, z] / np.pi, 1.0, 1.5 - results_single_30['elevation'][x, y, z] / np.pi])
        results_single_30_clrs_vec[i] = results_single_30_clrs[x, y, z]
    
# fig, [ax1, ax2, ax3] = plt.subplots(1, 3)
# ax1.imshow(results_single_30_clrs[128,:,:,:])
# ax2.imshow(results_single_30_clrs[:,128,:,:])
# ax3.imshow(results_single_30_clrs[:,:,128,:])

# create a 3d canvas to display the fibre orientation
canvas_results_single_30 = po.Canvas()
ps = canvas_results_single_30.pixel_scale
results_single_30_data = np.zeros(n, [('a_position', np.float32, 3),
                    ('a_color', np.float32, 4),
                    ('a_size', np.float32)])
results_single_30_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
results_single_30_data['a_color'] = np.concatenate((results_single_30_clrs_vec, np.ones((n,1))), axis=1)
results_single_30_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
canvas_results_single_30.data_prog.bind(gloo.VertexBuffer(results_single_30_data))
canvas_results_single_30.app.run()