import numpy as np
from skimage.draw import circle
from scipy import ndimage as ndi
from random import choice
import multiprocessing
from itertools import repeat
import fill_voids

def make_volume(volume_size, n_fibres, elvtn_rng, azth_rng, radius_lim, length_lim, gap, PSNR):
    # set limits for the data generation
    max_fails = 100
    # median_rad = 3
    max_len_loss = 0.5
    
    # initialize data
    volume = np.zeros(volume_size, dtype=np.uint8)
    elvtn_data = np.zeros_like(volume, dtype=np.float32)
    azth_data = np.zeros_like(volume, dtype=np.float32)
    diameter = np.zeros_like(volume, dtype=np.float32)
    
    n_generated = 0
    n_fails = 0
    while n_generated < n_fibres and n_fails < max_fails:
        # create fibres
        n_cores = multiprocessing.cpu_count()
        n = max(n_fibres-n_generated, n_cores)
        params = [volume_size, length_lim, radius_lim, azth_rng, elvtn_rng, gap, max_len_loss]
        with multiprocessing.Pool(processes=n_cores) as pool:
            fibres = pool.starmap(generate_fibre, repeat(params, n))

        for f in fibres:
            # if intersection is not allowed check if any fibres already exist in the gap region
            if np.any(volume[f['gap_fibre'][:, 0], f['gap_fibre'][:, 1], f['gap_fibre'][:, 2]]):
                n_fails += 1
                print("The number of fails is {}".format(n_fails))
                if n_fails == max_fails:
                    print("Maximum number of fails exceeded. Generated {} fibres".format(n_generated))

                continue
            
            if n_generated < n_fibres:
                # add the fibre to the volume
                volume[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = 1
                elvtn_data[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = f['elvtn']
                azth_data[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = f['azth']
                diameter[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = f['radius'] * 2
                n_generated += 1
                n_fails = 0
                print("The number of generated fibres is {}".format(n_generated))
    
    out = {'data': volume,
           'elevation': elvtn_data,
           'azimuth': azth_data,
           'diameter': diameter}
    
    for i in range(len(PSNR)):
        # add noise to the data
        nme = 'noisy_data_' + str(PSNR[i])
        out[nme] = add_noise(volume, PSNR[i])

    return out

def generate_fibre(volume_size, length_lim, radius_lim, azth_rng, elvtn_rng, gap, max_len_loss):
    np.random.seed()
    length = min(volume_size)
    gap_fibre_len = 0
    
    while float(gap_fibre_len) / length < max_len_loss:
        # set the length of the fibre
        length = min(volume_size)
        length = np.floor(length * np.random.uniform(length_lim[0], length_lim[1])).astype(np.int32)
    
        # set the location of the fibre
        offset = [np.random.uniform(olim[0], olim[1]) for olim in zip(-np.array(volume_size)/2, np.array(volume_size)/2)]
        offset = np.round(offset).astype(np.int32)
        # set the angles of the fibre
        ch = choice(range(len(azth_rng)))
        azth = np.random.uniform(np.deg2rad(azth_rng[ch][0]), np.deg2rad(azth_rng[ch][1]))
        elvtn = np.random.uniform(np.deg2rad(elvtn_rng[ch][0]), np.deg2rad(elvtn_rng[ch][1]))
        # set the radius of the fibre
        radius = np.random.uniform(radius_lim[0], radius_lim[1])
          
        # calculate the rotation matrices    
        rot_x = np.array([[1., 0., 0],
                       [0., np.cos(elvtn), -np.sin(elvtn)],
                       [0., np.sin(elvtn), np.cos(elvtn)]])
    
        rot_z = np.array([[np.cos(azth+np.pi/2), -np.sin(azth+np.pi/2), 0],
                       [np.sin(azth+np.pi/2), np.cos(azth+np.pi/2), 0],
                       [0., 0., 1.]])
    
        # generate the orientation vector
        orient_vec = np.array([0, 0, 1])
        orient_vec = np.dot(rot_x, orient_vec)
        orient_vec = np.dot(rot_z, orient_vec)
    
        # calculate the steps along the orientation vector
        n_steps = np.round(length)
        half_steps = int(np.ceil(n_steps / 2.))
        steps = range(half_steps - int(n_steps), half_steps)
    
        # draw a circle perpendicular to the orientation vector
        X, Y = circle(0, 0, radius)
        Z = np.repeat(0, len(Y))
        circle_pts = np.array([X, Y, Z])
        circle_pts = np.dot(rot_x, circle_pts)
        circle_pts = np.dot(rot_z, circle_pts)
        
        # draw the equivalent circle for the gap radius
        X, Y = circle(0, 0, radius + gap)
        Z = np.repeat(0, len(Y))
        gap_circle_pts = np.array([X, Y, Z])
        gap_circle_pts = np.dot(rot_x, gap_circle_pts)
        gap_circle_pts = np.dot(rot_z, gap_circle_pts)
        
        # propogate both circles along the orientation vector
        step_shifts = np.array([step * orient_vec for step in steps])
        center_shift = np.array([np.array(volume_size) * 0.5 + offset])
        
        slices = np.round(np.array([circle_pts.T + (step_shift + center_shift)
                                            for step_shift in step_shifts]))
        gap_slices = np.round(np.array([gap_circle_pts.T + (step_shift + center_shift)
                                            for step_shift in step_shifts]))
        
        # filter all the points which are outside the boundary
        pt_filter = lambda pt: np.all(np.greater_equal(pt, (0, 0, 0))) and \
                               np.all(np.less(np.array(pt), volume_size))
        
        # for each of the slices add those points within the volume to the fibre
        fibre = None
        for slc in slices:
            slice_mask = [pt_filter(pt) for pt in slc]
            slice_pts = slc[slice_mask].astype(np.int32)
            if len(slice_pts) > 0:
                fibre = slice_pts if fibre is None else \
                                                    np.concatenate((fibre, slice_pts))
        
        # fill in any holes
        volume = np.zeros(volume_size, dtype=np.uint8)
        volume[fibre[:, 0], fibre[:, 1], fibre[:, 2]] = 1
        volume = fill_voids.fill(volume)
        fibre = np.empty((np.count_nonzero(volume), 3), dtype=np.uint32)
        fibre[:, 0], fibre[:, 1], fibre[:, 2] = volume.nonzero()
        
        # repeat for the gap fibre
        n_slices = 0
        gap_fibre = None
        for slc in gap_slices:
            slice_mask = [pt_filter(pt) for pt in slc]
            slice_pts = slc[slice_mask].astype(np.int32)
            if len(slice_pts) > 0:
                n_slices += 1
                gap_fibre = slice_pts if gap_fibre is None else \
                                                    np.concatenate((gap_fibre, slice_pts))
        
        # calculate the length of the gap fibre                                                
        gap_fibre_len = np.round(n_slices).astype(np.int32)
    
    return {'fibre': fibre, 'gap_fibre': gap_fibre, 'elvtn': elvtn, 'azth': azth, 'radius': radius}

def add_noise(data, PSNR):    
    data = data.astype(np.float32)
    # generate Gaussian noise
    noise = np.random.randn(*data.shape)
    
    # calculate the required standard deviation of the noise
    sigma = np.sqrt(10 ** (-PSNR / 10) * (np.max(np.abs(data)) **2) / np.mean(noise**2))
    
    # add the noise to the data
    noisy_data = data + sigma * noise
    
    return noisy_data


if __name__ == '__main__':
    elvtn_rng = (30, 30)
    azth_rng = (45, 45)
    volume_size = (128, 128, 128)
    n_fibres = 4
    radius_lim = (2, 4)
    length_lim = (0.2, 0.8)
    gap = 3 
    intersect = False
    PSNR = 30
    smooth_lvl = 0
    out = make_volume(volume_size, n_fibres, elvtn_rng, azth_rng, radius_lim, length_lim, gap, intersect, PSNR)