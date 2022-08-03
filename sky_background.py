#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 14:53:13 2022

@author: ndemo1220
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from photutils import CircularAnnulus


def sky_background(data, centre, r_in, r_out, mask = None, add_plot = False):
    '''
    Measures the background in a sky annulus. Units are whatever units the map provided are in.
    NOTE: It is preferable to have units of FLUX.

    Parameters
    ----------
    data : (N X M array)
           Our galaxy image data.
           
    centre : (float)
             The galaxy centre in pixels (x,y).
             
    r_in : (float)
           The inner radius value in pixels.
           
    r_out : (float)
            The outer radius value in pixels.
            
    mask : (bool), optional
           A boolean mask for sources/edges not to be considered. The default is None.
           
    add_plot : (bool), optional
               If True, the code will return a plot of your galaxy with the background annulus. 
               The default is False.

    Raises
    ------
    ValueError
        If the inner radius is larger than the outer radius.

    Returns
    -------
    intensity : (float)
                The iaverage value of the background measured in our annulus
                
    area : (float)
           The area of our annulus in pixels 

    '''
    if not r_out > r_in:
        raise ValueError("'r_out' must be greater than 'r_in'.")
        
    
        
    if mask is not None and np.dtype(mask[0][0]) != 'bool':
        print("[WARNING]: mask is not True/False boolean. Converting to True/False array.")
        mask = np.array(mask, dtype = bool)
        
    
    ap = CircularAnnulus(positions=centre, r_in=r_in, r_out=r_out)
            
    sum_flux = ap.do_photometry(data, method = 'exact')[0][0]
    
    area = ap.area_overlap(data, method='exact', mask = mask) # just incase our ellipses are larger than the 2D image
    
    # now intensity
    
    intensity = sum_flux / area
    
    if add_plot == True:
        
        #initialize plot
        
        plt.figure(figsize = (9,8))
        
        norm = colors.LogNorm(vmax = data.max(), vmin = data[data>0].min())
        plt.imshow(data, cmap='afmhot_r', origin='lower',norm=norm)
        plt.xlabel('$x$ / pixel',fontsize = 15)
        plt.ylabel('$y$ / pixel',fontsize=15)
        cbar = plt.colorbar()
        cbar.set_label('Intensity', fontsize = 15)
        ap.plot(color = 'b', lw = 8)
    
    return intensity, area