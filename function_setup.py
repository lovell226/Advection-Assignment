###############################################################################
# Date Revised: 17 Sept 2025
# Author: Hilary Weller <h.weller@reading.ac.uk>
# Description: Basic tests of advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, l2error, amp_ratio

# First make sure that the plots folder exists
if not os.path.exists("plots"):
    os.makedirs("plots")

def evolveScheme(init, schemes,
                           nx = 40, nt = 10, u = 1, endTime = 0.3,
                           xLim = [0,1],a = 0, b = 0.5):
    '''
    Modified versrion of ComparesSchemeSnapshots to run the numerical schemes and store final solutions
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
        a (float): The start of the non-zero section in ICs
        b (float): The end of the non-zero section in ICs
    Returns:
        Output (array): x, initial conditions, analytic solution and final state for all schemes
    '''
    Output = {}
    # Derived parameters
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx

    Output['x'] = x
    
    # Set the initial conditions
    phi0 = init(x,a,b)
    Output['Initial'] = phi0.copy()
    
    # Find the analytic solution
    phiA = init((x - u*endTime)%(xLim[1] - xLim[0]), a, b)
    Output['Analytic'] = phiA.copy()
    
    # Calculate the results of all of the schemes
    for scheme in schemes:
        # The solution, phi, is stored as a list of arrays at different time
        # steps, phiLevels. Initially just the initial conditions. The latest
        # time level array is always phiLevels[-1]
        phiLevels = [phi0.copy()]
        for n in range(nt):
            scheme(phiLevels, c)
        # add final to stored list
        Output[scheme.__name__] = phiLevels[-1].copy()

    return Output


def selectCourants(nx = 40, u = 1, endTime = 0.1, xLim = [0,1]):
    Output = {}
    dx = (xLim[1] - xLim[0])/nx
    desiredCourants = np.linspace(0.05, 5, 50)
    nts = np.round(u*endTime/(desiredCourants*dx)).astype(int)   
    # calculate courant numbers to ensure all nt divide endTime exactly
    courants = u*(endTime/nts)/dx 
    # to remove courant numbers outside our desired range
    mask = (courants >= 0.05) & (courants <= 5)
    courants = courants[mask]
    nts = nts[mask]

    #remove any duplicates as a result of rounding
    nts, indices = np.unique(nts, return_index=True)
    courants = courants[indices]
    
    Output['nts'] = nts
    Output['courants'] = courants
    return(Output)


def MetricStore(schemes):
    # To store all errors for our schemes
    metrics = {}
    for scheme in schemes:
        metrics[scheme.__name__] = {
            "L2": [],
            "Amp": []
        }
    return metrics

def parameterSweep(
        parameterValues,
        schemes,
        nx=40,
        nt=10,
        u=1,
        endTime=0.3,
        xLim=[0,1],
        varyWidth=False,
        varyTime=False,
        fixedCourant=0.8,
        init = cosBell):

    metrics = MetricStore(schemes)

    for value in parameterValues:
        # for width experiments:
        if varyWidth:
            centre = 0.3
            a = centre - value/2
            b = centre + value/2
            currentNt = nt
            currentEndTime = endTime
        # for longer time testing with fixed courant:
        elif varyTime:
            currentEndTime = value
            dx = (xLim[1]-xLim[0])/nx
            # keep Courant number fixed
            currentNt = int(np.round(currentEndTime*u/(fixedCourant*dx)))
            a = 0.1
            b = 0.5
        else:
            currentNt = value
            currentEndTime = endTime
            a = 0.1
            b = 0.5
        # the run the schemes for new parameters
        A = evolveScheme(
            init,
            schemes,
            nx=nx,
            nt=currentNt,
            u=u,
            endTime=currentEndTime,
            xLim=xLim,
            a=a,
            b=b
        )

        dx = (xLim[1]-xLim[0])/nx

        for scheme in schemes:
            # add the calculated errors to the stored metric list
            name = scheme.__name__
            
            metrics[name]["L2"].append(l2error(A[name],A["Analytic"]))

            metrics[name]["Amp"].append(amp_ratio(A[name],A["Analytic"],))

    return metrics


def plotMetric(
        x,
        outFile,
        metricStorage,
        metricName,
        xlabel,
        ylabel,
        title,
        yLim = [None, None],
        xLim = [None, None],
        logy=False 
        ):
    # seperate function to plot all metric plots
    # so colours are consistent across report
    cols = {'FTBS': 'blue',
            'CTCS': 'red',
            'FTCS': 'green',
            'LinearSemiLagrangian': 'cyan',
            'CubicSemiLagrangian': 'magenta'}
    
    plt.figure()

    for schemeName in metricStorage:
        
        plt.plot(
            x,
            metricStorage[schemeName][metricName],
            color=cols.get(schemeName, "black"),
            label=schemeName
        )

    if logy:
        plt.yscale("log")
    plt.xlim(xLim)
    plt.ylim(yLim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    plt.legend()

    plt.savefig("plots/"+outFile)
