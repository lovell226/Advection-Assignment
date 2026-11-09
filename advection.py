###############################################################################
# Date Revised: 17 Sept 2025
# Author: Hilary Weller <h.weller@reading.ac.uk>
# Description: Basic tests of advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os

def compareSchemeSnapshots(init, schemes, outFile,
                           nx = 40, nt = 10, u = 1, endTime = 0.4,
                           xLim = [0,1], yLim = [None,None],a = 0, b = 0.5):
    '''
    Author: Hilary Weller <h.weller@reading.ac.uk>
    Compares snapshots of results of advection schemes using periodic boundaries
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The minimum (xLim[0]) and maximim (xLim[1]) x-values of
                      the domain. xLim[1] is not part of the domain when periodic
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If [None,None] is specified then default limits
                      are used.
    Returns:
        None
    '''
    # changed init to call a and b as well so that I can change this for courant tests
    # Derived parameters
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    dt = endTime/nt
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    print('Courant number =', c)
    
    # Set the initial conditions
    phi0 = init(x,a,b)
    
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']

    cols = {'FTBS': 'blue',
                'CTCS': 'red',
                'FTCS': 'green',
                'LinearSemiLagrangian': 'cyan',
                'CubicSemiLagrangian': 'magenta'}
    
    # Find and plot the analytic solution
    phiA = init((x - u*endTime)%(xLim[1] - xLim[0]), a, b)
    plt.figure()
    plt.plot(x, phiA, 'k', label = "Analytic", linewidth=2)
    
    # Calculate the results of all of the schemes and plot
    for scheme in schemes:
        # The solution, phi, is stored as a list of arrays at different time
        # steps, phiLevels. Initially just the initial conditions. The latest
        # time level array is always phiLevels[-1]
        phiLevels = [phi0.copy()]
        for n in range(nt):
            scheme(phiLevels, c)
        plt.plot(x, phiLevels[-1], color=cols.get(scheme.__name__, 'black'), label = scheme.__name__)
    
    plt.xlim(xLim)
    plt.ylim(yLim)
    plt.xlabel('x')
    plt.ylabel(r'$\phi$')
    plt.legend()
    plt.savefig("plots/"+outFile)
