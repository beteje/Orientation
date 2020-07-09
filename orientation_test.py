import numpy as np
import orientation_simulation as sim
from matplotlib.colors import hsv_to_rgb
from vispy import gloo
import plot_orientation as po

# %gui qt

elvtn_rng = (10, 179)
azth_rng = (0, 179)
volume_size = (128, 128, 128)
n_fibres = 5
radius_lim = (3, 10)
length_lim = (0.4, 0.9)
gap_lim = 1
intersect = False
PSNR = 10
smooth_lvl = 0
volume = sim.make_volume(volume_size, n_fibres, elvtn_rng, azth_rng, radius_lim, length_lim, gap_lim, intersect, PSNR, smooth_lvl)

# plot histograms of angles
po.plot_histogram(volume['elvtn'], 'Elevation', (0, 180), 'darkgreen')
po.plot_histogram(volume['azth'], 'Azimuth', (0, 180), 'darkblue')

# plot slices of original and noisy data
po.plot_slices(volume['data'], volume['noisy_data'], 50)

# plot colour wheel
po.plot_color_wheel((0,180), (0,180))

# create a 3d canvas to display the fibres
X, Y, Z = volume['data'].nonzero()
canvas_volume = po.Canvas()
n = len(Z)
ps = canvas_volume.pixel_scale
volume_data = np.zeros(n, [('a_position', np.float32, 3),
                    ('a_color', np.float32, 4),
                    ('a_size', np.float32)])
volume_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
volume_data['a_color'] = np.concatenate((np.zeros((n, 3)), np.ones((n,1))), axis=1)
volume_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
canvas_volume.data_prog.bind(gloo.VertexBuffer(volume_data))
canvas_volume.app.run()

# generate the rgb colors for each of the points
clrs = np.zeros((n,3), dtype=np.float32)
    
for i, (x, y, z) in enumerate(zip(X, Y, Z)):
    if volume['elvtn'][x, y, z] / np.pi <= 0.5:
        clrs[i] = hsv_to_rgb([volume['azth'][x, y, z] / np.pi, 2 * volume['elvtn'][x, y, z] / np.pi, 1.0])
    else:
        clrs[i] = hsv_to_rgb([volume['azth'][x, y, z] / np.pi, 1.0, 1.5 - volume['elvtn'][x, y, z] / np.pi])

# create a 3d canvas to display the fibre orientation
canvas_orientation = po.Canvas()
ps = canvas_orientation.pixel_scale
orientation_data = np.zeros(n, [('a_position', np.float32, 3),
                    ('a_color', np.float32, 4),
                    ('a_size', np.float32)])
orientation_data['a_position'] = np.array([2*X/volume_size[0]-1, 2*Y/volume_size[1]-1, 2*Z/volume_size[2]-1]).T
orientation_data['a_color'] = np.concatenate((clrs, np.ones((n,1))), axis=1)
orientation_data['a_size'] = np.ones(n) * ps * 2 * np.sqrt(2)
canvas_orientation.data_prog.bind(gloo.VertexBuffer(orientation_data))
canvas_orientation.app.run()

