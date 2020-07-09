import numpy as np
import math
import itertools
from skimage import io, filters, morphology, feature
from scipy import ndimage as ndi
import vigra
from multiprocessing import Pool
import matplotlib as mpl
import matplotlib.pyplot as plt

def prepare_data(data):
    if data.ndim == 2:
        # threshold data using Otsu's method
        th_val = filters.threshold_otsu(data)
        data_seg = (data > th_val).astype(np.uint8)
        
        # create a 1 pixel wide skeleton
        skel = morphology.skeletonize(data_seg)
    elif data.ndim == 3:
        # threshold data using Otsu's method
        data_seg = np.zeros(data.shape)
        for i in range(data.shape[0]):
            if np.any(data[i] != data[i][0]):
                data_seg[i] = (data[i] > filters.threshold_otsu(data[i])).astype(data_seg.dtype)
                
            
        # create a 1 pixel wide skeleton
        skel = morphology.skeletonize_3d(data_seg)
    
    return data_seg, skel

def calc_orientation(skel, window_radius):
    # set the amount to pad by
    padding = window_radius*2 + 1
    # pad the skeleton
    padded_skel = np.pad(skel, pad_width=(padding,), mode='constant', constant_values=0) 
    
    # determine coordinates of non-zero values in the skeleton
    ycords, xcords = np.nonzero(padded_skel)
    # initialize the orientation map
    orientation_map = np.zeros_like(padded_skel, dtype=np.float32)
    for i, (yy, xx) in enumerate(zip(ycords, xcords)):
        # select patch surrounding the current point
        patch_skel = padded_skel[(yy-window_radius):(yy+window_radius+1),
                                 (xx-window_radius):(xx+window_radius+1)]
        
        # estimate structure tensor elements
        Axx, Axy, Ayy = feature.structure_tensor(patch_skel, sigma=0.1)
        tensor_vals = np.array([[np.mean(Axx), np.mean(Axy)], [np.mean(Axy), np.mean(Ayy)]])
        
        # perform eigen-analysis
        w, v = np.linalg.eig(tensor_vals)
        
        # calculate orientation
        orientation = math.atan2(*v[:,np.argmax(w)])
        orientation_map[yy, xx] = orientation if orientation >= 0.0 else (np.pi - np.abs(orientation))
    
    # remove padded values
    orientation_map = orientation_map[padding:-padding, padding:-padding]

    return orientation_map

def calc_orientation_3d(skel, data, window_size, sigma):
    # find the points corresponding to non-zero elements of the skeleton
    Z, Y, X = skel.nonzero()
    pts = zip(Z, Y, X)
    
    # set the half window size
    ws2 = np.uint32(window_size/2)
    
    patches = []
    for pt in pts:
        # extract patch surounding the current point
        lim0 = pt - ws2
        lim1 = pt + ws2
    
        if any(np.array(lim0) < 0) or any(np.array(lim1) > (data.shape[0] - 1)):
            patches.append(None)
        else:
            z0, y0, x0 = lim0
            z1, y1, x1 = lim1
            patches.append(data[z0:z1, y0:y1, x0:x1])
    
    # use parallel pool to estimate the structure tensor of each patch
    args = zip(patches, itertools.repeat(sigma))
    proc_pool = Pool(processes=12)
    results = np.array(proc_pool.map(unpack_estimate_tensor, args))
    proc_pool.close()
    proc_pool.join()
    proc_pool.terminate()
    
    # extract the individual orientation components from the results
    lat_arr, azth_arr = results.T
    lat = np.zeros_like(skel, dtype=np.float32)
    azth = np.zeros_like(skel, dtype=np.float32)
    lat[Z, Y, X] = lat_arr
    azth[Z, Y, X] = azth_arr
         
    return lat, azth

def unpack_estimate_tensor(args):
    # for use with parallel pool to call estimate_tensor
    return estimate_tensor(*args)

def estimate_tensor(patch, sigma):
    # If the current point has no patch return empty orientation
    if patch is None:
        return (0, 0)
    
    # estimate the components of the structure tensor
    img = vigra.filters.structureTensor(patch, 1, 1, sigma_d=sigma)
    Axx = img[:, :, :, 0]
    Axy = img[:, :, :, 1]
    Axz = img[:, :, :, 2]
    Ayy = img[:, :, :, 3]
    Ayz = img[:, :, :, 4]
    Azz = img[:, :, :, 5]
    tensor_vals = np.array([[np.mean(Azz), np.mean(Ayz), np.mean(Axz)],
                              [np.mean(Ayz), np.mean(Ayy), np.mean(Axy)],
                              [np.mean(Axz), np.mean(Axy), np.mean(Axx)]])
    # tensor_vals = tensor_vals[::-1, ::-1]
    
    # perform eigen-analysis
    eps = 1e-8
    w, v = np.linalg.eig(tensor_vals)  
    mv = v[:, np.argmin(w)]
    mv[np.abs(mv) < eps] = 0
    
    # extract the orientation angles
    G = np.sqrt(mv[2]**2 + mv[1]**2)
    lat = np.arcsin(np.around(G, decimals=3))
    azth = np.arctan(mv[1] / mv[2]) if mv[2] else np.pi/2.
    
    return (lat, azth)

def plot_orientation_map(orient_map, skel):
    # produce a circular grayscale dilation of the orientation &  bindary dilation of skeleton
    orient_map = ndi.grey_dilation(orient_map, structure=morphology.disk(1)).astype(np.float32)
    skel = ndi.binary_dilation(skel, structure=morphology.disk(1)).astype(np.float32)
    # generate a mask of the orientation map
    masked_orient_map = np.ma.masked_where(skel == 0, orient_map)
    
    # create a color map object and normalize colours between 0 & pi
    cmap_obj = mpl.cm.get_cmap('hsv', 2056)
    cmap_obj.set_bad(color='black')
    norm = mpl.colors.Normalize(0.0, np.pi)
    
    # plot the masked orientation map using the color map
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_axis_off()
    ax.imshow(masked_orient_map, cmap=cmap_obj, interpolation=None)
    
    # create a color wheel for the different angles
    display_axes = fig.add_axes([0.770, -0.05, 0.2, 0.2], projection='polar')
    cb = mpl.colorbar.ColorbarBase(display_axes, cmap=cmap_obj, norm=norm, orientation='horizontal')
    display_axes.text(0.05, 0.3,'0', color='white', fontsize=20, weight='bold', horizontalalignment='center',
                      verticalalignment='center', transform=display_axes.transAxes)
    display_axes.text(0.85, 0.3, '180',color='white', fontsize=20, weight='bold', horizontalalignment='center',
                      verticalalignment='center', transform=display_axes.transAxes)
    cb.outline.set_visible(False)
    display_axes.set_axis_off()
    
    plt.show()

def plot_orientation_map_3d(azth_data, lat_data):
    # find the limits of the region where data exists
    r = np.any(azth_data, axis=(1, 2))
    c = np.any(azth_data, axis=(0, 2))
    z = np.any(azth_data, axis=(0, 1))
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    
    # adjust the limits of the azimuth & latitude
    azth, lat = azth_data[rmin:rmax, cmin:cmax, zmin:zmax], \
                np.abs(lat_data[rmin:rmax, cmin:cmax, zmin:zmax])
     
    # create a skeleton from the azimuth
    skel = azth.copy().astype(np.float32)
    skel[skel.nonzero()] = 1.
    
    # produce a spherical grey dilation of the azimuth & latitude & binary dilation of skeleton
    azth = ndi.grey_dilation(azth, structure=morphology.ball(1))
    lat = ndi.grey_dilation(lat, structure=morphology.ball(1))
    skel = ndi.grey_dilation(skel, structure=morphology.ball(1))
    
    # generate the rgb colors for each of the points in the skeleton
    Z, Y, X = skel.nonzero()
    vol_orient = np.zeros(skel.shape + (3,), dtype=np.float32)
    for z, y, x in zip(Z, Y, X):
        vol_orient[z, y, x] = mpl.colors.hsv_to_rgb([azth[z, y, x]/np.pi, lat[z, y, x]/np.pi, 1.0])
  
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.voxels(skel, facecolors = vol_orient, edgecolors = np.clip(2*vol_orient - 0.5, 0, 1), linewidth=0.5)
    plt.show

if __name__ == '__main__':
    # load image
    img = io.imread('/media/DATADRIVE1/Code/quanfima/data/polymer_slice.tif')
    # img = np.random.random((64,64,64))
    
    # segment data & create skeleton
    img_seg, skeleton = prepare_data(img)

    if img.ndim == 2:
        # calculate the orientation map
        orientation = calc_orientation(skeleton, 12)

        # plot the orientation map
        plot_orientation_map(orientation, skeleton)

    elif img.ndim == 3:
        # calculate the orientation map
        lat, azth = calc_orientation_3d(skeleton, img, 22, 0.025)
    
        plot_orientation_map_3d(azth, lat)
    
    
    