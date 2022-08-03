# -*- coding: utf-8 -*-
''' 
To find total flux, especially for noisy data. Also can calculate r50, and r90. 
Created May 2022: Nick DeMarchi
'''
import numpy as np
from scipy.interpolate import interp1d

class PlateauFinder(object):
    
    def __init__(self, data, plateau_length = None, ratio = 1.0, interp_length = int(1e5), method = 'first'):

        '''
        Takes our x and y variables for a function that has maximum(s)/minimum(s) and 
        finds where we become constant. We assume (x,y) is a simple function and can be
        interpolated. The process is to compute the first derivative on our interpolated data
        and find regions that approach zero (i.e. constant).
        
        Parameters
        ----------
        data : (2D array)
                A 2D array of our x and y values. x and y must be same length.
        
        plateau_length : (int)
                         The amount of consecutive x-values that must exist below a calculated threshold
                         in order for the iteration to stop. If None, the default is to chose the average 
                         spacing between the original data points and find the amount of interpolated data
                         that exist in this spacing.
                         
        ratio : (float)
                If plateau_length is None, this sets whether we want the avergae spacing between a set of 
                original data points or more. The default is 1.0 which is the average spacing between two points.
                
       interp_length : (int)
                       Unless the user specifies a function (and edits the code) we obtain a smoothed function
                       for our data by interpolating (this was used on noisy data with <20 points). This sets the 
                       number of data points we wish to have in our interpolated function. The default is 10,000.
                       
       method : (string)
                'first' --> takes the first point for which our algorithm detects the 'plateau'. Returns a single
                            cutoff value
                'all'  --> considers all points which exist below the calculated threshold parameter.
                
       Returns
       -------
       
       The interpolated x and y values (PlateauFinder.interp_x, PlateauFinder.interp_y) along with the cutoff value(s)
       for where we reach our plateau (PlateauFinder.cutoff).
       
       NOTE: This was used for radial profiles for nearby galaxies and the inputs were flux and radius or intensity 
             and radius which have a relatively smooth CDF.
            
        
        '''
        
        # separate x and y
        if data.shape[0] == data.ndim:
            x = data[0]
            y = data[1]
        else:
            print("[WARNING]: make sure that 'data' input is a 2 X N array. Taking transpose of input array")
            data = data.T
            x = data[0]
            y = data[1]
        
        if isinstance(x, np.ndarray) == False:
            x = np.asarray(x,dtype = type(x))
        if isinstance(y, np.ndarray) == False:
            y = np.asarray(y,dtype = type(y))
        
        if len(x) != len(y):
            raise ValueError('Array lengths do not match')
        
        ndim_x = len(x.shape)
        ndim_y = len(y.shape)
        if (ndim_x != 1) or (ndim_y != 1):
            raise ValueError('Input for data must be a 2 X N array')
        
        if plateau_length is None:
            print('setting plateau_length to be '+str(ratio)+' times the spacing between two data points along x-axis')
            # average spacing between two points
            spacing = np.sum(x[1:] - x[:-1]) / (len(x)-1)
            # total length of data
            length = x[-1] - x[0]
            # percentage length of spacing between two points
            percent = spacing / length
            plateau_length = int(ratio * percent * interp_length)
            
        self.x = x
        self.y = y
        self.plateau_length = plateau_length
        self.interp_length = interp_length
        self.method = method
        # run the algorithm
        self.run()
        
    def run(self):
        
        '''
        run our algorithm to find the plateau locations. This function does not return anything
        itself but runs the algorithm and attached several atti=ributes to our class instance.
        
        Attributes
        ----------
        
        interp_x : (array)
                   The interpolated x-values
                   
        interp_y : (array)
                   The interpolated y-values
                   
        derivative : (2D array)
                     The x and y values for our derivative. This has size interp_length - 2.          
                   
        cutoff_index : (integer or array of integers)
                       The indexed value in which the algorithm detects the plateau. This can then
                       be used to index the x and y values at which we reach the plateau. 
        
        
        '''
        
        def f_prime(f,var,step):
            '''
            compute the first derivative using the central difference method.
            f'(x)' = (f(x+h) - f(x-h)) / (2*h)

            Parameters
            ----------
            f : (function)
                The function that takes 'var' as a variable
            var : (array)
                  Thevariable passed to our function 'f'
            step : (float)
                   The step size for our derivative

            Returns
            -------
            (array)
                The first derivative of our function 'f'

            '''
            return (f(var+step) - f(var-step))/(2*step)
        
        # interpolate
        func = interp1d(self.x, self.y, kind = 'cubic')
        r = np.linspace(self.x[0], self.x[-1], self.interp_length)
        F = func(r)
    
        # take derivative
        diff = f_prime(func,r[1:-1],1/len(r[1:-1]))
        
        self.derivative = np.array((r[1:-1],diff))
        
        # now set initial alpha by finding minimum. We can arrange our derivative values to do this
        #first keep only positive values for where our cumulative y goes to zero
        diff_pos = diff[diff>0]
        dx = np.sort(diff_pos, kind = 'stable')
        assert dx[0] == diff_pos.min(), "Minimum is not correct"
        
        # now iterate about this minimum dx[0] until we reach length desired for 'plateau'
        index = 0
        threshold = dx[index]
        cutoff_values = np.where((dx <= threshold) == True)[0]
        while len(cutoff_values) < self.plateau_length:
            index += 1
            threshold = dx[index]
            cutoff_values = np.where((diff_pos < threshold) == True)[0]
        print('Did we find our plateau ?', len(cutoff_values) >= self.plateau_length )
        if len(cutoff_values) >= self.plateau_length:
            if self.method == 'first':
                cutoff = int(cutoff_values[0])
                self.threshold = threshold
            elif self.method == 'all':
                cutoff = int(cutoff_values)
                self.threshold = threshold
# =============================================================================
#                 elif method == 'mid':
#                     ind = int(self.plateau_length / 2)
#                     cutoff = int(cutoff_values[ind])
#                 elif method == 'end':
#                     cutoff = int(cutoff_values[-1])
# =============================================================================
            else:
                self.threshold = threshold
                raise ValueError("Method must be 'first' or 'all'")
        else:
            print('[WARNING]: Check plateau_length value. Code was never able to reach this size')
            cutoff = int(cutoff_values[int(len(threshold/2))])
        print('our threshold value is: ', threshold)
    
        self.interp_x = r
        self.interp_y = F
        self.cutoff = cutoff
        
    def return_totals(self):
        '''
        This function can be called to return some statistics related to the CDF of our data.

        Raises
        ------
        ValueError
            If our algorithm failed to procude a cutoff value, a ValueError will be raised here.

        Returns
        -------
        x50 : (float)
              The x-value which contains 50 percent of our CDF.
        x90 : (float)
              The x-value which contains 90 percent of our CDF.
        x_max : (float)
                The x-value our algorithm cutoff at.
        total_y : (float)
                  The total y-value. This is our 'maximum'

        '''
        
        
    
        if self.cutoff is None:
            raise ValueError('The index which we cut off at must be provided')
            
        if self.cutoff not in range(len(self.interp_x)):
            raise ValueError('cutoff value must be in the range 0 <= cutoff <= len(x)')
        
        total_y = self.interp_y[self.cutoff]
        x_max = self.interp_x[self.cutoff]
        
        condition = total_y > 2*self.interp_y[0]
        if not condition:
            print('[WARNING]: Initial y-value is larger than half the total y-value. Not able to measure meaningful x_50')
            x50 = self.x[0]
            x90 = self.interp_x[np.where(self.interp_y <= total_y*0.9)[0][-1]]
        else:
            # Find R50 and R90
            index_50 = np.where((self.interp_y <= total_y / 2) & (self.interp_x < x_max))[0][-1]
            index_90 = np.where((self.interp_y <= total_y * 0.9) & (self.interp_x < x_max))[0][-1]
        
            x50 = self.interp_x[index_50]
            x90 = self.interp_x[index_90]
            
        return x50, x90, x_max, total_y
