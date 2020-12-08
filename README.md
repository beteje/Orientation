# Orientation
Estimation of the orientation of fibres in a volume

We estimate 3D cellular alignment based on image gradients and directional statistics.

The intensity gradients of the volumetric image are used to construct a 3D vector field and the local dominant orientations of this vector field then determined.

- _orientation_estimation_ performs the simulation and estimation of the orientation angles (requires cuda libraries)
- _generate_fibre_volume_ provides the FibreVolume class which produces the volume, positions the fibres within the volume and adds noise. Also provides additional functions to plot slices of the volume and histograms of the angles. The positioning of the fibres is a simple version based on the [Quanfima package](https://github.com/rshkarin/quanfima/tree/master/quanfima) which provides a wider range of functionality than the one here
- _plot_orientation_ Produces 3D visualizations of the volumes using the [VisPy library](http://vispy.org/)
