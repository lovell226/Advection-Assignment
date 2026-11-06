###############################################################################
# Date Revised: 17 Sept 2025
# Author: Hilary Weller <h.weller@reading.ac.uk>
# Description: Basic tests of advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, squareWave, l2error, amp_ratio, Mass_diff
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2

# First make sure that the plots folder exists
if not os.path.exists("plots"):
    os.makedirs("plots")

def evolveScheme(init, schemes,
                           nx = 40, nt = 10, u = 1, endTime = 0.3,
                           xLim = [0,1],a = 0, b = 0.5):
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
    Output = {}
    #changed init to call a and b as well so that I can change this for corrant number tests
    # Derived parameters
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx

    Output['x'] = x
    #dt = endTime/nt
    print('Courant number =', c)
    
    # Set the initial conditions
    phi0 = init(x,a,b)
    Output['Initial'] = phi0.copy()
    # Colours for each line
    #cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
    #        'pink', 'orange', 'gray', 'orange']
    
    # Find and plot the analytic solution
    phiA = init((x - u*endTime)%(xLim[1] - xLim[0]), a, b)
    Output['Analytic'] = phiA.copy()
    #plt.figure()
    #plt.plot(x, phiA, 'k', label = "Analytic", linewidth=2)
    
    # Calculate the results of all of the schemes and plot
    for scheme in schemes:
        # The solution, phi, is stored as a list of arrays at different time
        # steps, phiLevels. Initially just the initial conditions. The latest
        # time level array is always phiLevels[-1]
        phiLevels = [phi0.copy()]
        for n in range(nt):
            scheme(phiLevels, c)
        #plt.plot(x, phiLevels[-1], cols[i], label = schemes[i].__name__)
        #want to add final to stored list
        Output[scheme.__name__] = phiLevels[-1].copy()

    return Output
    #plt.xlim(xLim)
    #plt.ylim(yLim)
    #plt.xlim(xLim)
    #plt.ylim(-0.2, 2)
    #plt.xlabel('x')
    #plt.ylabel(r'$\phi$')
    #plt.legend()
    #plt.savefig("plots/"+outFile)

# List of functions to do the advection
schemes = [FTBS,CTCS,SemiLagrangian,SemiLagrangian2]


#Need to keep this one for section 3
#test = (evolveScheme(cosBell, schemes,
#                       nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1], a = 0.1, b = 0.5))

#plt.figure()
#plt.plot(test['x'],test['Analytic'], label='Analytic')
#plt.plot(test['x'],test['FTBS'], label='FTBS')
#plt.plot(test['x'],test['FTCS'], label='FTCS')
#plt.plot(test['x'],test['CTCS'], label='CTCS')
#plt.show()

def selectCourants(nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1]):
    #courants = np.linspace(cMin, cMax, nC)
    Output = {}
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    #dt = endTime/nt
    #c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    #print('Courant number =', c)
        
    #nt=np.zeros(len(courants))
    desiredCourants = np.linspace(0.05, 5, 50)
    nts = np.round(
    u*endTime/(desiredCourants*dx)
    ).astype(int)   
    #nts = np.arange(5, 100)

    #courants = nx*u*endTime/(nts*(xLim[1]-xLim[0]))
    courants = u*(endTime/nts)/dx 
    mask = (courants >= 0.05) & (courants <= 5)

    courants = courants[mask]
    nts = nts[mask]
    
    for i in range(len(courants)):
        nts[i]=np.ceil((u*endTime)/(courants[i]*dx))
        courants[i]=u*(endTime/nts[i])/dx
    #remove any duplicates as a result of rounding
    nts, indices = np.unique(nts, return_index=True)
    courants = courants[indices]
    
    #print(nts)
    #print(courants)
    Output['nts'] = nts
    Output['courants'] = courants
    return(Output)

Courant_test = selectCourants(nx = 40, nt = 10, u = 1, endTime = 2.0, xLim = [0,1])

def initialiseMetricStorage(schemes):

    metrics = {}

    for scheme in schemes:

        metrics[scheme.__name__] = {
            "L2": [],
            "Mass": [],
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
        RealAnalytical = 0,
        init = cosBell):

    metrics = initialiseMetricStorage(schemes)

    for value in parameterValues:

        if varyWidth:

            centre = 0.3

            a = centre - value/2
            b = centre + value/2

            currentNt = nt
            currentEndTime = endTime
        elif varyTime:
            currentEndTime = value

            dx = (xLim[1]-xLim[0])/nx

            # keep Courant number fixed
            currentNt = int(
                np.round(
                    currentEndTime*u/(fixedCourant*dx)
                )
            )

            a = 0.1
            b = 0.5
        else:

            currentNt = value
            currentEndTime = endTime

            a = 0.1
            b = 0.5

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

            name = scheme.__name__
            if RealAnalytical == 0:
                metrics[name]["L2"].append(
                    l2error(
                        A[name],
                        A["Analytic"]
                    )
                )

                metrics[name]["Mass"].append(
                    Mass_diff(
                        A[name],
                        A["Analytic"],
                        dx
                    )
                )

                metrics[name]["Amp"].append(
                    amp_ratio(
                        A[name],
                        A["Analytic"],
                    )
                )
            #else:
            #    #want to run against a much higher resolution analytical solution (much smaller grid spacing)
            #    nx = nx*4
            #    B = evolveScheme(
            #                cosBell,
            #                schemes = [],
            #                nx=nx,
            #                nt=currentNt,
            #                u=u,
            #                endTime=endTime,
            #                xLim=xLim,
            #                a=a,
            #                b=b
            #            )
                #phiA = cosBell((x - u*endTime)%(xLim[1] - xLim[0]), a, b)
            #    metrics[name]["L2"].append(
            #        l2error(
            #            A[name],
            #            B["Analytic"]
            #        )
            #    )

            #    metrics[name]["Mass"].append(
            #        Mass_diff(
            #            A[name],
            #            B["Analytic"],
            #            dx
            #        )
            #    )

            #    metrics[name]["Amp"].append(
            #        amp_ratio(
            #            A[name],
            #            B["Analytic"],
            #        )
            #    )                

    return metrics


def plotMetric(
        x,
        outFile,
        metricStorage,
        metricName,
        xlabel,
        ylabel,
        title,
        ymax = [None, None],
        logy=False 
        ):

    plt.figure()

    for schemeName in metricStorage:

        plt.plot(
            x,
            metricStorage[schemeName][metricName],
            label=schemeName
        )

    if logy:
        plt.yscale("log")
    else:
        plt.ylim(ymax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    #plt.grid()

    plt.legend()

    plt.savefig("plots/"+outFile)

Courant_test = selectCourants(
    nx=40,
    endTime=2.0
)

courantMetrics = parameterSweep(
    Courant_test["nts"],
    schemes,
    endTime=2.0
)

plotMetric(
    Courant_test["courants"],
    "L2_Courant.pdf",
    courantMetrics,
    "L2",
    "Courant Number",
    "L2 Error",
    "L2 Error vs Courant Number",
    logy=True
)

plotMetric(
    Courant_test["courants"],
    "Mass_Conservation_Courant.pdf",
    courantMetrics,
    "Mass",
    "Courant Number",
    "Mass Difference",
    "Mass Difference vs Courant Number",
    logy=True
    #ymax= [0,0.04] 
)

plotMetric(
    Courant_test["courants"],
    "Amp_Courant.pdf",
    courantMetrics,
    "Amp",
    "Courant Number",
    "Amplitude Ratio",
    "Amplitude Ratio vs Courant Number",
    ymax = [0,2]
)

widthCells = np.array([20,15,10,8,6,4,3,2,1.5,1])

widths = widthCells*(1/40)

widthMetrics = parameterSweep(
    widths,
    schemes,
    varyWidth=True,
    nx = 40,
    nt = 15,
    endTime = 0.3,
)

plotMetric(
    widthCells,
    "L2_Width.pdf",
    widthMetrics,
    "L2",
    "Bell Width",
    "L2 Error",
    "L2 Error vs Width",
    logy=True
)

plotMetric(
    widthCells,
    "Amp_Width.pdf",
    widthMetrics,
    "Amp",
    "Bell Width",
    "Amplitude Ratio",
    "Amplitude Ratio vs Width"
)

#for a time test
times = [0.5, 1, 2, 5, 10, 20]
timeMetrics = parameterSweep(
    times,
    schemes,
    varyTime=True
)
plotMetric(
    times,
    "L2_Time.pdf",
    timeMetrics,
    "L2",
    "Time",
    "L2 Error",
    "L2 Error vs Time (C=0.8)",
    logy=True
)

plotMetric(
    times,
    "Mass_Time.pdf",
    timeMetrics,
    "Mass",
    "Time",
    "Mass Error",
    "Mass Error vs Time (C=0.8)",
    logy=True
)
plotMetric(
    times,
    "Amp_Time.pdf",
    timeMetrics,
    "Amp",
    "Time",
    "Amplitude Ratio",
    "Amplitude Ratio vs Time (C=0.8)"
)

times2 = [0.5, 1, 2, 5, 10, 20]
timeMetrics2 = parameterSweep(
    times,
    schemes,
    varyTime=True,
    init=squareWave
)
plotMetric(
    times2,
    "L2_Time_square.pdf",
    timeMetrics2,
    "L2",
    "Time",
    "L2 Error",
    "L2 Error vs Time (C=0.8)",
    logy=True
)