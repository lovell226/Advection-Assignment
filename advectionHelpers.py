# Helper functions for testing advection schemes.
# Authors: Hilary Weller <h.weller@reading.ac.uk>
#          & Becks Lovell <r.lovell@pgr.reading.ac.uk>

import numpy as np

def cosBell(x, a, b):
    '''
    Returns a cosine bell shape between parameters a and b
    Author: Hilary Weller <h.weller@reading.ac.uk>
    Parameters:
        x (1darray): The x-locations for calculating the cosine bell.
        a (float): The start of the non-zero bell
        b (float): The end of the non-zero bell
    Returns:
        phi : 1darray. A cosine bell shape between x=a and x=b
    '''
    return np.where((x > a) & (x < b), 
                    0.5*(1 - np.cos(2*np.pi*(x - a)/(b-a))), 0.)


def l2error(phi, phiA):
    '''
    Returns the l2 error norm of phi in comparison to exact solution phiA.
    Author: Hilary Weller <h.weller@reading.ac.uk>
    Parameters:
        phi (1darray): The numerical solution
        phiA (1darray): The analytic solution
    Returns
        l2 (float): The l2 error norm
    '''
    phiError = phi - phiA
    l2 = np.sqrt(np.sum(phiError**2)/np.sum(phiA**2))
    return l2


def amp_ratio(phi, phiA):
    '''
    Returns the max amp ratio of phi in comparison to exact solution phiA.
    Parameters:
        phi (1darray): The numerical solution
        phiA (1darray): The analytic solution
    Returns
        amp (float): The amplitude ratio
    '''
    amp = np.max(phi)/np.max(phiA)
    return amp


def squareWave(x, a, b):
    '''
    Returns a square 'wave' between parameters a and b
    Parameters:
        x (1darray): The x-locations for calculating the shape.
        a (float): The start of the non-zero section
        b (float): The end of the non-zero section
    Returns:
        phi : 1darray. A square wave between x = a and x = b
    '''
    #ie 0 at all points except between a and b where = 1
    return np.where((x > a) & (x < b), 1., 0.)

