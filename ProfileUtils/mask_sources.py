#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
To mask sources in 2-dimensional photometric galaxy images.
Created: April 2022 - Nick DeMarchi
'''

import warnings
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.wcs import WCS
from astroquery.vizier import Vizier


class MaskSources(object):
    """
    Class for creating a mask for a galaxy image.
    
    Parameters
    ----------
    data : (N X M array)
           The galaxy image.
           
    header : (astropy.io.fits.header.Header)
             The FITS file header for our image.
             
    name : (string)
           The galaxy name.
           
    path_to_table : (string)
                    This should be a path to a table that contains the columns 'RA', 'DEC', 'pa', 'inclination'
                    and 'Galaxy'.
                    NOTE: This was created for my purpose to use the VERTICO tables from Toby Brown 2021 paper. 
                    The __init__ can be edited to load these values itself rather than point to a directory 
                    where the table is stored.
                    
    ps : (float)
         The pixel scale four our image.
         
    include_galaxy : (bool)
                     If True we also mask our galaxy. Default is False.
                     
    large_incl_corr : (bool)
                      Large inclination correction. If True we correct for inclination based on parameter 'incl_limit'.
                      The default is True.
                      
    incl_limit : (float)
                 The cos(inclination) limit. Default is cos(inclination) = 0.2.
    Returns
    -------
            The boolean mask for our image; under ProfileMask.disk_mask
    """

    def __init__(self, data, header, name, path_to_table, ps=None, include_galaxy=False, large_incl_corr=True, incl_limit=0.2):
        
        if data.ndim != 2:
            raise ValueError('data must be 2 dimensional image')
        if ps is None:
            raise ValueError('pixel-scale must be specified')

        self.data = data
        self.header = header
        self.name = name
        self.path_to_table = path_to_table
        self.ps = ps
        self.include_galaxy = include_galaxy
        # load WCS
        self.w = WCS(header)
        self.large_incl_corr = large_incl_corr
        self.incl_limit = incl_limit
        # run our nested function
        self.disk_mask = self.calc_mask()
    
    
    def calc_mask(self):
        '''
        
        Main function that will calculate the boolean mask. Nested is a function that produces a radius map
        to select sources to mask. 

        Returns
        -------
        (bool)
            Boolean mask of sources to mask

        '''
        
        def radius_map(shape, ra, dec, pa, incl, w, large_incl_corr=True,incl_limit=0.2):
            '''
            
            A function to create a radius map. We can select specific sources and mask them
            to a certain size. We choose 1.5 times the optical size R25 for the mask.
            
            '''

            # All inputs assumed as Angle
            if large_incl_corr and (np.isnan(pa.rad + incl.rad)):
                pa = Angle(0 * u.rad)
                incl = Angle(0 * u.rad)
                # Not written to the header
                msg = 'PA or INCL is NaN in radius calculation, setting both to zero'
                warnings.warn(msg, UserWarning)
                # Warning ends
            cos_pa, sin_pa = np.cos(pa.rad), np.sin(pa.rad)
            cos_incl = np.cos(incl.rad)
            if large_incl_corr and (cos_incl < incl_limit):
                cos_incl = incl_limit
            xcm, ycm = ra.rad, dec.rad
            coordinates = np.zeros(list(shape) + [2])
            # Original coordinate is (y, x)
            # :1 --> x, RA --> the one needed to be divided by cos(incl)
            # :0 --> y, Dec
            coordinates[:, :, 0], coordinates[:, :, 1] = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
            # Now, value inside coordinates is (x, y)
            # :0 --> x, RA --> the one needed to be divided by cos(incl)
            # :1 --> y, Dec
            for i in range(shape[0]):
                coordinates[i] = Angle(w.wcs_pix2world(coordinates[i], 1) * u.deg).rad
            coordinates[:, :, 0] = 0.5 * (coordinates[:, :, 0] - xcm) * (np.cos(coordinates[:, :, 1]) + np.cos(ycm))
            coordinates[:, :, 1] -= ycm
            # Now, coordinates is (dx, dy) in the original coordinate
            # cos_pa*dy-sin_pa*dx is new y
            # cos_pa*dx+sin_pa*dy is new x
            radius = np.sqrt((cos_pa * coordinates[:, :, 1] + sin_pa * coordinates[:, :, 0]) ** 2 + (
                    (cos_pa * coordinates[:, :, 0] - sin_pa * coordinates[:, :, 1]) / cos_incl) ** 2)
            radius = Angle(radius * u.rad).arcsec
            return radius
        
        # initialize boolean mask
        disk_mask = np.ones(self.data.shape, dtype='int64') < 0
        # first step is to mask any bad pixels or edges
        initial_mask = ~np.isfinite(self.data)
        disk_mask |= initial_mask
        # next mask all other known galaxies in sample
        load_table = fits.open(self.path_to_table)[1].data
        if self.include_galaxy:
            table = load_table
        else:
            table = load_table[~np.isin(load_table['Galaxy'], self.name)]
        pos = table['RA'] * u.deg, table['DEC'] * u.deg, table['pa'] * u.deg, table['inclination'] * u.deg
        ra, dec, pa, incl = Angle(pos[0]), Angle(pos[1]), Angle(pos[2]), Angle(pos[3])
        table_length = len(table['Galaxy'])
        for i in range(table_length):
            query = Vizier.query_region(table['Galaxy'][i], radius='0d0m20s', catalog='HyperLEDA')
            r25_extra = Angle((10 ** query[0]['logD25'][0] * 0.1 / 60 / 2) * u.deg).arcsec
            radius_map_extra = radius_map(self.data.shape, ra[i], dec[i], pa[i], incl[i], self.w)
            disk_mask |= radius_map_extra < (1.5 * r25_extra)
        # now load our galaxy ra,dec
        n_load_table = fits.open(self.path_to_table)[1].data
        n_table = n_load_table[np.isin(n_load_table['Galaxy'], self.name)]
        pos = n_table['RA'] * u.deg, n_table['DEC'] * u.deg
        n_ra, n_dec = Angle(pos[0]), Angle(pos[1])
        # Mask out other objects in the field of view
        c0 = SkyCoord(n_ra, n_dec, frame='icrs', equinox="J2000")
        # choose field of view to be size of input data
        leda_query = Vizier.query_region(c0,
                                         width=self.ps * self.data.shape[0] * u.arcsec,
                                         height=self.ps * self.data.shape[1] * u.arcsec,
                                         catalog='HyperLEDA')
        n_objects = len(leda_query[0]['RAJ2000'])
        print('number of masked objects is: ', n_objects)
        for j in range(n_objects):
            # need to make sure my galaxy exists in field of view and is not masked
            if len(leda_query[0]['ANames'][j].split()) > 0:
                n_names = len(leda_query[0]['ANames'][j].split())
                for k in range(n_names):
                    if leda_query[0]['ANames'][j].split()[k] == self.name:
                        print(leda_query[0]['ANames'][j].split()[k])
                        print('found you: ', self.name)
                        continue
            else:
                gal_ra_parsed = leda_query[0]['RAJ2000'][j].split(' ')
                ra_j = Angle(gal_ra_parsed[0] + 'h' + gal_ra_parsed[1] + 'm' + gal_ra_parsed[2] + 's')
                gal_dec_parsed = leda_query[0]['DEJ2000'][j].split(' ')
                dec_j = Angle(gal_dec_parsed[0] + 'd' + gal_dec_parsed[1] + 'm' + gal_dec_parsed[2] + 's')
                # 'logR25' # Axis ratio in log scale
                incl_j = Angle((np.arccos(1. / (10 ** leda_query[0]['logR25'][j])) * 180. / np.pi) * u.deg)
                # 'PA' # Position angle
                pa_j = Angle(leda_query[0]['PA'][j] * u.deg)
                # 'logD25' # Apparent diameter
                r25_j = Angle((10 ** leda_query[0]['logD25'][j] * 0.1 / 60 / 2) * u.deg).arcsec
                # create radius map
                radius_map_j = radius_map(self.data.shape, ra_j, dec_j, pa_j, incl_j, self.w)
                # add elliptical mask for specified extra source
                disk_mask |= radius_map_j < (1.5 * r25_j)
        # the final return
        return (~disk_mask)
    
