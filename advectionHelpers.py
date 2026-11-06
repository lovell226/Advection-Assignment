# Helper functions for testing advection schemes.
# Authors: Hilary Weller <h.weller@reading.ac.uk>
#          Becks Lovell <r.lovell@pgr.reading.ac.uk>

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

def l_inf_error(phi, phiA):
    '''
    Returns the l infinity error norm of phi in comparison to exact solution phiA.
    Parameters:
        phi (1darray): The numerical solution
        phiA (1darray): The analytic solution
    Returns
        l_inf (float): The l infinity error norm
    '''
    phiError = phi - phiA
    l_inf = np.max(np.abs(phiError))
    return l_inf

def Mass_diff(phi, phiA, dx):
    '''
    Returns the change in mass of phi in comparison to exact solution phiA.
    Parameters:
        phi (1darray): The numerical solution
        phiA (1darray): The analytic solution /or initial condition
    Returns
        Mass (float): The relative change in mass
    '''
    #mass0 = np.sum(np.abs(phi))*dx
    #massA = np.sum(np.abs(phiA))*dx
    mass0 = np.sum(phi)*dx
    massA = np.sum(phiA)*dx
    Mass = np.abs(mass0-massA)/massA
    return Mass

def amp_ratio(phi, phiA):
    '''
    Returns the amp ratio of phi in comparison to exact solution phiA.
    Parameters:
        phi (1darray): The numerical solution
        phiA (1darray): The analytic solution
    Returns
        amp (float): The lamplitude ratio
    '''
    amp = np.max(phi)/np.max(phiA)
    return amp

def peak_shift(x, phi, phiA):
    '''
    Distance between numerical and exact peak locations
    '''

    peak_num = x[np.argmax(phi)]
    peak_exact = x[np.argmax(phiA)]

    return peak_num - peak_exact

#def error_all_calc(phi, phiA):
#    errors = []


def squareWave(x, a, b):
    '''
    Returns a uniform distrbution type shape between parameters a and b
    Parameters:
        x (1darray): The x-locations for calculating the shape.
        a (float): The start of the non-zero jump
        b (float): The end of the non-zero jump
    Returns:
        phi : 1darray. A unifrom distribution type shape between x=a and x=b
    '''
    #ie 0 at all points except between a and b when = 1
    return np.where((x > a) & (x < b), 
                    1., 0.)

