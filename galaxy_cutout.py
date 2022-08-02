'''
To create smaller thumnails for our galaxy images
Created November 2021: Nick DeMarchi
'''

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle, Skycoord
from astroquery.vizier import Vizier
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel

def galaxy_cutout(data, header, ra, dec, name, size=10.0, centre_and_radius=False):
    """
    Two-dimensional cutout for galaxy image.

    Parameters:
    -----------
    data : (N X M array)
          The galaxy image.
    header : (FITS file header)
            The fits file header for our image.
    ra : (float)
         The right ascension of our target galaxy.
    dec : (float)
          The declination of our target galaxy.
    name : (str)
           The galaxy name
    size : (float)
           The factor of R25 (optical size of our galaxy) we want the cutout to be. Default is 10.
    centre_and_radius : (bool)
                        If we want our function to return the galaxy centre and R25 radius.
                        Default is False.
    Returns:
    --------
         cutout : (N X M array)
                  The new cutout image
         centre : (2D array, optional)
                  The galaxy centre in pixels
         r : (float)
             The galaxy radius R25 in arc-seconds

    """
    # remove bad pixels from the data

    is_good = np.isfinite(data)

    data[~is_good] = 0

    # sky coordinate for galaxy centre in RA and DEC

    galaxy_position = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')

    # load WCS

    w = WCS(header)

    # load an optical radius, R25

    leda_query = Vizier.query_region(name, radius='0d0m20s', catalog='HyperLEDA')

    r = Angle((10 ** leda_query[0]['logD25'][0] * 0.1 / 60 / 2) * u.deg).arcsec

    cutout = Cutout2D(data, galaxy_position, size=size * r * u.arcsec, wcs=w)

    xp, yp = skycoord_to_pixel(galaxy_position, wcs=cutout.wcs)

    centre = (xp, yp)

    if not centre_and_radius:

        return cutout

    else:

        return cutout, centre, r
