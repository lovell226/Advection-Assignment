###############################################################################
# Date Revised: 19 Sept 2025
# Author: Hilary Weller <h.weller@reading.ac.uk>
# Description: Create a graph of errors growing over time for advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, l2error
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2

# First make sure that the plots folder exists
if not os.path.exists("plots"):
    os.makedirs("plots")

def growingErrors(init, schemes, outFile,
                  nx = 40, nt = 10, u = 1, endTime = 0.4, 
                  xLim = [0,1], yLim = None, a = 0, b = 0.5):
    '''
    Author: Hilary Weller <h.weller@reading.ac.uk>
    Creates a graph of error growth over time for various advection schemes
    Parameters:
        nx (int): The number of points in space
        init (callable): Function to set the initial conditions.
        schemes (list of callable): List of functions to evaluate the schemes.
        nt (int): The number of time steps.
        u (float): The wind speed
        endTime (float) : The final time
        outFile (string): The name of the output file.
        xLim (tuple): The maximum and minimun spacial dimensions
        yLim (tuple): The minimum (yLim[0]) and maximum (yLim[1]) y values of
                      the graph. If None is specified then default limits
                      are used.
    Returns:
        None
    '''
    # Derived parameters
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    dt = endTime/nt
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    print('Courant number =', c)
    
    # Set the initial conditions for all schemes
    phi0 = init(x, a, b)
    
    # Store an array of all times
    times = np.linspace(0, endTime, nt+1)
    
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']
    plt.figure()
    # Loop over all the schemes and plot
    for i in range(len(schemes)):
        # The l2 error for each time step
        l2errors = np.zeros(nt+1)
        # The solution for this scheme
        phiLevels = [phi0.copy()]
        # Loop through all time steps
        for n in range(nt):
            # Update phi
            schemes[i](phiLevels, c)
            # Calculate the analytic solution
            phiA = init((x - u*times[n+1])%(xLim[1] - xLim[0]), a, b)
            # Calculate the l2 error
            l2errors[n+1] = l2error(phiLevels[-1], phiA)

        # Plot the error for this scheme
        plt.semilogy(times, l2errors, cols[i], label = schemes[i].__name__)
    
    plt.xlim([0, endTime])
    plt.ylim(yLim)
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel(r'$\ell_2$ error')
    plt.savefig("plots/"+outFile)

schemes = [FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2]

growingErrors(cosBell, schemes, "basicSchemesCosErrors.pdf",
              nx = 40, nt = 40, u = 1, endTime = 0.4,xLim = [0,1], a = 0, b = 0.5)
#Ive changed the above to match the advection params used - check through notes to make sure they aren't different


def Errors_vs_Courant(init, schemes, outFile,
                  nx = 40, u = 1, endTime = 0.4, 
                  xLim = [0,1], yLim = None, a = 0, b = 0.5, cMin=0.1, cMax=10, nC=50):
    
    # Derived parameters
    courants = np.linspace(cMin, cMax, nC)
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    #dt = endTime/nt
    #c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    #print('Courant number =', c)
    
    nt=np.zeros(len(courants))

    for i in range(len(courants)):
        nt[i]=np.ceil((u*endTime)/(courants[i]*dx))
        courants[i]=u*(endTime/nt[i])/dx
    print(nt)

    # Set the initial conditions for all schemes
    phi0 = init(x, a, b)
    
    # Store an array of all times
    #times = np.linspace(0, endTime, nt+1)
    
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']
    plt.figure()
    # Loop over all the schemes and plot
    for i in range(len(schemes)):
        # The l2 error for each courant number
        l2errors = np.zeros(nC)
        # The solution for this scheme
        #phiLevels = [phi0.copy()]
        # Loop through all time steps
        for j in range(len(nt)):
            phiLevels = [phi0.copy()]
            #for each j in nt run for nt and store final error
            print(len(nt))
            print(j)
            times = np.linspace(0, endTime, int(nt[j])+1)
            for n in range(int(nt[j])):
                # Update phi
                schemes[i](phiLevels, courants[j])
                # Calculate the analytic solution
            phiA = init((x - u*times[n+1])%(xLim[1] - xLim[0]), a, b)
                # Calculate the l2 error
            l2errors[j] = l2error(phiLevels[-1], phiA)#not sure this is right??

        # Plot the error for this scheme
        print(np.shape(courants))
        print(np.shape(l2errors))
        plt.semilogy(courants, l2errors, cols[i], label = schemes[i].__name__)
        print(courants)
    plt.xlim([0, max(courants)])
    plt.ylim(-10,100)
    plt.legend()
    plt.xlabel('Courant Number')
    plt.ylabel(r'$\ell_2$ error')
    plt.savefig("plots/"+outFile)
    

schemes = [FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2]

def Profiles_vs_Courant(init, schemes,
                        nx=40, u=1, endTime=0.4,
                        xLim=[0,1], a=0, b=0.5,
                        cVals=[0.5,1,2,5]):

    dx = (xLim[1]-xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)

    phi0 = init(x,a,b)

    for scheme in schemes:

        plt.figure(figsize=(8,6))

        phiA = init(
            (x-u*endTime)%(xLim[1]-xLim[0]),
            a,b
        )

        plt.plot(x, phiA, 'k', lw=3,
                 label='Analytic')

        for c in cVals:

            nt = int(np.ceil(
                (u*endTime)/(c*dx)
            ))

            cActual = u*(endTime/nt)/dx

            phiLevels = [phi0.copy()]

            for n in range(nt):
                scheme(phiLevels, cActual)

            plt.plot(
                x,
                phiLevels[-1],
                label=f'C={cActual:.2f}'
            )

        plt.title(scheme.__name__)
        plt.xlabel('x')
        plt.ylabel(r'$\phi$')
        plt.legend()
        plt.show()

def Error_vs_Courant(init, schemes,
                     nx=40, u=1, endTime=0.4,
                     xLim=[0,1], a=0, b=0.5,
                     cMin=0.1, cMax=5, nC=50):

    dx = (xLim[1]-xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)

    phi0 = init(x,a,b)

    requestedC = np.linspace(cMin,cMax,nC)

    plt.figure(figsize=(8,6))

    for scheme in schemes:

        errors = np.zeros(nC)
        actualC = np.zeros(nC)

        for j,c in enumerate(requestedC):

            nt = int(np.ceil(
                (u*endTime)/(c*dx)
            ))

            cActual = u*(endTime/nt)/dx
            actualC[j] = cActual

            phiLevels = [phi0.copy()]

            for n in range(nt):
                scheme(phiLevels, cActual)

            phiA = init(
                (x-u*endTime)%(xLim[1]-xLim[0]),
                a,b
            )

            errors[j] = l2error(
                phiLevels[-1],
                phiA
            )

        plt.semilogy(
            actualC,
            errors,
            label=scheme.__name__
        )

    plt.xlabel('Courant Number')
    plt.ylabel(r'$L_2$ Error')
    plt.legend()
    plt.grid()
    plt.show()

def Growth_vs_Courant(init, schemes,
                      nx=40, u=1, endTime=0.4,
                      xLim=[0,1], a=0, b=0.5,
                      cMin=0.1, cMax=5, nC=50):

    dx = (xLim[1]-xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)

    phi0 = init(x,a,b)

    amp0 = np.max(np.abs(phi0))

    requestedC = np.linspace(cMin,cMax,nC)

    plt.figure(figsize=(8,6))

    for scheme in schemes:

        growth = np.zeros(nC)
        actualC = np.zeros(nC)

        for j,c in enumerate(requestedC):

            nt = int(np.ceil(
                (u*endTime)/(c*dx)
            ))

            cActual = u*(endTime/nt)/dx
            actualC[j] = cActual

            phiLevels = [phi0.copy()]

            for n in range(nt):
                scheme(phiLevels, cActual)

            growth[j] = (
                np.max(np.abs(phiLevels[-1]))
                / amp0
            )

        plt.plot(
            actualC,
            growth,
            label=scheme.__name__
        )

    plt.axhline(
        1,
        color='k',
        linestyle='--'
    )

    plt.xlabel('Courant Number')
    plt.ylabel('Amplitude Ratio')
    plt.legend()
    plt.grid()
    plt.show()

schemes = [
    FTBS,
    FTCS,
    CTCS,
    SemiLagrangian,
    SemiLagrangian2
]

Profiles_vs_Courant(
    cosBell,
    schemes,
    cVals=[0.5,1,2,5]
)

Error_vs_Courant(
    cosBell,
    schemes,
    cMin=0.1,
    cMax=5,
    nC=50
)

Growth_vs_Courant(
    cosBell,
    schemes,
    cMin=0.1,
    cMax=5,
    nC=50
)