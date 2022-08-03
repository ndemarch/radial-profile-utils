# radial-profile-utils
A collection of code to create radial profiles and measure various quantities. 

This repository contains the following functions:

1. galaxy_cutout() --> to create thumbnail from our original galaxy image
2. MaskSources() --> to mask objects/edges in our image
3. PlateauFinder() --> find where our profile or CDF becomes constant (i.e. total flux cutoff)
4. sky_background() --> measure the contribution of the sky in our image using a circular annulus
