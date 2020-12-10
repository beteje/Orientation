#!/usr/bin/env python3
import numpy as np
from skimage.draw import circle
from random import choice
from itertools import repeat
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import multiprocessing, fill_voids, h5py

class FibreVolume:
    def __init__(self, volume_size=(128, 128, 128), n_fibres=10, PSNR=[30]):
        # Size of the volume
        self.volume_size = volume_size
        # Number of fibres in the volume
        self.n_fibres = n_fibres
        # PSNR value or values for the noisy data
        self.PSNR = PSNR
        # initialize data
        self.data = np.zeros(volume_size, dtype=np.uint8)
        self.elevation = np.zeros_like(self.data, dtype=np.float32)
        self.azimuth = np.zeros_like(self.data, dtype=np.float32)
        self.diameter = np.zeros_like(self.data, dtype=np.float32)
        self.circularity = []
        self.circularity_XY = []
        self.circularity_XZ = []
        self.circularity_YZ = []
        self.mean_elevation = []
        self.mean_azimuth = []
        self.R_length = []
        self.noisy_data = {}
    
    # -------------------------------------------------------------------------           
    def make_volume(self, elvtn_rng=[[30, 30]], azth_rng=[[45, 45]],  
                 radius_lim=(2, 4), length_lim=(0.2, 0.8), gap=3):
        # set limits for the data generation
        max_fails = 100
        # median_rad = 3
        max_len_loss = 0.5
        
        n_generated = 0
        n_fails = 0
        while n_generated < self.n_fibres and n_fails < max_fails:
            # create fibres
            n_cores = multiprocessing.cpu_count()
            n = max(self.n_fibres-n_generated, n_cores)
            params = [self.volume_size, length_lim, radius_lim, azth_rng, elvtn_rng, gap, max_len_loss]
            with multiprocessing.Pool(processes=n_cores) as pool:
                fibres = pool.starmap(self.generate_fibre, repeat(params, n))
    
            for f in fibres:
                # if intersection is not allowed check if any fibres already exist in the gap region
                if np.any(self.data[f['gap_fibre'][:, 0], f['gap_fibre'][:, 1], f['gap_fibre'][:, 2]]):
                    n_fails += 1
                    print("The number of fails is {}".format(n_fails))
                    if n_fails == max_fails:
                        print("Maximum number of fails exceeded. Generated {} fibres".format(n_generated))
    
                    continue
                
                if n_generated < self.n_fibres:
                    # add the fibre to the volume
                    self.data[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = 1
                    self.elevation[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = f['elvtn']
                    self.azimuth[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = f['azth']
                    self.diameter[f['fibre'][:, 0], f['fibre'][:, 1], f['fibre'][:, 2]] = f['radius'] * 2
                    n_generated += 1
                    n_fails = 0
                    print("The number of generated fibres is {}".format(n_generated))
        
        # add noise to the data
        self.add_noise()
    
    # -------------------------------------------------------------------------
    def add_noise(self):    
        data = self.data.astype(np.float32)
        
        for i in range(len(self.PSNR)):
            # generate Gaussian noise
            noise = np.random.randn(*data.shape)
        
            # calculate the required standard deviation of the noise
            sigma = np.sqrt(10** (-self.PSNR[i] / 10) * (np.max(np.abs(data))**2) / np.mean(noise**2))
        
            # add the noise to the data
            nme = 'PSNR' + str(self.PSNR[i])
            self.noisy_data.update({nme : data + sigma * noise})
    
    # -------------------------------------------------------------------------
    def generate_fibre(self, volume_size, length_lim, radius_lim, azth_rng, elvtn_rng, gap, max_len_loss):
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
        
    # -------------------------------------------------------------------------
    def plot_histogram(self, xlabel, xlims, barcolor):
        # identify the data to be plotted
        if len(self.circularity) != 0:
            X, Y, Z = self.circularity.nonzero()
        else:
            X, Y, Z = self.data.nonzero()
            
        if xlabel.lower() == 'elevation':
            data = self.elevation[X, Y, Z]
        elif xlabel.lower() == 'azimuth':
            data = self.azimuth[X, Y, Z]
        else:
            print('Please select either the Elevation or Azimuth data to plot')
            return
        
        # convert data to be plotted to an array in degrees
        data = np.rad2deg(data.flatten())

        # set options for the histogram
        weights = np.ones_like(data)/float(len(data))
        num_bins = int((xlims[1] - xlims[0])/5)
        
        # plot the histogram
        fig = plt.figure(figsize=(14,8))
        ax = fig.add_subplot(111)
        n, bins, patches = ax.hist(data, num_bins, xlims, color=barcolor, rwidth=0.8, weights=weights)
        
        # set the x-axis ticks, limits & labels
        ax.set_xlim(xlims)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.tick_params(axis='x', labelsize=20, which='major', direction='out', length=8, width=2)
        ax.tick_params(axis='x', which='minor', direction='out', length=4, width=2)
        ax.set_xlabel(xlabel + ' (deg.)', labelpad=2, fontsize=30, color='black')
        
        # set the y-axis ticks & labels
        ax.set_ylim((0, 1))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.tick_params(axis='y', labelsize=20, which='major', direction='out', length=8, width=2)
        ax.tick_params(axis='y', which='minor', direction='out', length=4, width=2)
        vals = ax.get_yticks()
        ax.set_yticklabels(['{:3.0f}'.format(x*100) for x in vals])
        ax.set_ylabel('Frequency of Occurrence (%)', labelpad=2, fontsize=30, color='black')
        
        # set general figure options
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)
        
        plt.tight_layout()
    
    # -----------------------------------------------------------------------------
    def plot_slices(self, slice_no):
        # Plot slices of the original data from each dimension
        fig, axs = plt.subplots(3, len(self.PSNR)+1)
        axs[0, 0].imshow(self.data[:, :, slice_no], cmap=plt.cm.get_cmap('Greys'), origin='lower')
        axs[0, 0].set(xlabel='x axis', ylabel='y axis')
        
        axs[1, 0].imshow(self.data[:,slice_no,:], cmap=plt.cm.get_cmap('Greys'), origin='lower')
        axs[1, 0].set(xlabel='x axis', ylabel='z axis')
    
        axs[2, 0].imshow(self.data[slice_no,:,:], cmap=plt.cm.get_cmap('Greys'), origin='lower')
        axs[2, 0].set(xlabel='x axis', ylabel='z axis')
    
        # Repeat for each of the different noise levels
        for i in range(len(self.PSNR)):
            nme = 'PSNR' + str(self.PSNR[i])
            tmp = self.noisy_data[nme][:, :, slice_no]
            mn = np.amin(tmp)
            mx = np.amax(tmp)
            tmp = (tmp - mn)/(mx-mn)
            axs[0, i+1].imshow(tmp, cmap=plt.cm.get_cmap('Greys'),  origin='lower')
            axs[0, i+1].set(xlabel='x axis', ylabel='y axis')
    
            tmp = self.noisy_data[nme][:,slice_no,:]
            mn = np.amin(tmp)
            mx = np.amax(tmp)
            tmp = (tmp - mn)/(mx-mn)
            im = axs[1, i+1].imshow(tmp, cmap=plt.cm.get_cmap('Greys'),  origin='lower')
            axs[1, i+1].set(xlabel='x axis', ylabel='z axis')
        
            tmp = self.noisy_data[nme][slice_no,:,:]
            mn = np.amin(tmp)
            mx = np.amax(tmp)
            tmp = (tmp - mn)/(mx-mn)
            im = axs[2, i+1].imshow(tmp, cmap=plt.cm.get_cmap('Greys'),  origin='lower')
            axs[2, i+1].set(xlabel='x axis', ylabel='z axis')
    
        fig.colorbar(im, ax=axs[:,:])
    
        for ax in axs.flat:
            ax.label_outer()
        
        plt.show()
        
    # -----------------------------------------------------------------------------
    def save_volume(self, fname):     
        # Save the data related to the current fibre volume to an hdf5 file
        with h5py.File(fname, "w") as f:
            f.create_dataset('volume_size', data=self.volume_size)
            f.create_dataset('n_fibres', data=self.n_fibres)
            f.create_dataset('PSNR', data=self.PSNR)
            f.create_dataset('data', data=self.data)
            f.create_dataset('diameter', data=self.diameter)
            f.create_dataset('azimuth', data=self.azimuth)
            f.create_dataset('elevation', data=self.elevation)
            f.create_dataset('circularity', data=self.circularity)
            f.create_dataset('circularity_XY', data=self.circularity_XY)
            f.create_dataset('circularity_XZ', data=self.circularity_XZ)
            f.create_dataset('circularity_YZ', data=self.circularity_YZ)
            f.create_dataset('mean_elevation', data=self.mean_elevation)
            f.create_dataset('mean_azimuth', data=self.mean_azimuth)
            f.create_dataset('R_length', data=self.R_length)
            if self.noisy_data:
                for i in range(len(self.PSNR)):
                    nme = 'PSNR' + str(self.PSNR[i])
                    f.create_dataset('noisy_data/'+nme, data=self.noisy_data[nme])
            
    # -----------------------------------------------------------------------------
    def load_volume(self, fname):
        # Load the data from the hdf5 file into the current fibre volume
        f = h5py.File(fname, "r") 
        self.volume_size = f['volume_size'][()]
        self.n_fibres = f['n_fibres'][()]
        self.PSNR = f['PSNR'][()]
        self.data = f['data'][()]
        self.diameter = f['diameter'][()]
        self.azimuth = f['azimuth'][()]
        self.elevation = f['elevation'][()]
        self.circularity = f['circularity'][()]
        self.circularity_XY = f['circularity_XY'][()]
        self.circularity_XZ = f['circularity_XZ'][()]
        self.circularity_YZ = f['circularity_YZ'][()]
        self.mean_elevation = f['mean_elevation'][()]
        self.mean_azimuth = f['mean_azimuth'][()]
        self.R_length = f['R_length'][()]
        if 'noisy_data' in f.keys():
            for i in range(len(self.PSNR)):
                nme = 'PSNR' + str(self.PSNR[i])
                self.noisy_data[nme] = f['noisy_data/'+nme][()]
           
if __name__ == '__main__':
    out = FibreVolume()
    out.make_volume()
    # plot histograms of angles
    out.plot_histogram('Elevation', (0, 180), 'darkgreen')
    out.plot_histogram('Azimuth', (0, 180), 'darkblue')
    
    # plot slices of original and noisy data
    out.plot_slices(64)