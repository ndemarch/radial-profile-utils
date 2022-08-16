# radial-profile-utils
A collection of code to create radial profiles and measure various quantities. Requires [Numpy](https://numpy.org/install/), [Astropy](https://docs.astropy.org/en/stable/install.html), [Matplotlib](https://matplotlib.org/stable/users/installing/index.html#:~:text=If%20you%20are%20using%20the,sudo%20yum%20install%20python3%2Dmatplotlib), [Scipy](https://scipy.org/install/), [Astroquery],(https://astroquery.readthedocs.io/en/latest/) [Photutils](https://photutils.readthedocs.io/en/stable/install.html#)

This repository contains the following functions:

1. galaxy_cutout() --> To create thumbnail from our original galaxy image
2. MaskSources() --> To mask objects/edges in our image
3. PlateauFinder() --> To find where our profile or CDF becomes constant (i.e. total flux cutoff).
                   --> This function can also return r50 and r90 estimates.
4. sky_background() --> To measure the contribution of the sky in our image using a circular annulus
