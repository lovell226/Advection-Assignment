###############################################################################
# Date Revised: 17 Sept 2025
# Author: Hilary Weller <h.weller@reading.ac.uk>
# Description: Basic tests of advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, l2error, amp_ratio, Mass_diff
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2

# First make sure that the plots folder exists
if not os.path.exists("plots"):
    os.makedirs("plots")

def evolveScheme(init, schemes,
                           nx = 40, nt = 10, u = 1, endTime = 0.4,
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
    Output['x'] = x
    dt = endTime/nt
    c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    print('Courant number =', c)
    
    # Set the initial conditions
    phi0 = init(x,a,b)
    Output['Initial'] = phi0.copy()
    # Colours for each line
    cols = ['blue', 'red', 'green', 'cyan', 'magenta', 'brown', 'purple',
            'pink', 'orange', 'gray', 'orange']
    
    # Find and plot the analytic solution
    phiA = init((x - u*endTime)%(xLim[1] - xLim[0]), a, b)
    Output['Analytic'] = phiA.copy()
    #plt.figure()
    #plt.plot(x, phiA, 'k', label = "Analytic", linewidth=2)
    
    # Calculate the results of all of the schemes and plot
    for i in range(len(schemes)):
        # The solution, phi, is stored as a list of arrays at different time
        # steps, phiLevels. Initially just the initial conditions. The latest
        # time level array is always phiLevels[-1]
        phiLevels = [phi0.copy()]
        for n in range(nt):
            schemes[i](phiLevels, c)
        #plt.plot(x, phiLevels[-1], cols[i], label = schemes[i].__name__)
        #want to add final to stored list
        Output[schemes[i].__name__] = phiLevels[-1].copy()

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
schemes = [FTBS,FTCS,CTCS,SemiLagrangian,SemiLagrangian2]


#Need to keep this one for section 3
test = (evolveScheme(cosBell, schemes,
                       nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1], a = 0.1, b = 0.5))

#plt.figure()
#plt.plot(test['x'],test['Analytic'], label='Analytic')
#plt.plot(test['x'],test['FTBS'], label='FTBS')
#plt.plot(test['x'],test['FTCS'], label='FTCS')
#plt.plot(test['x'],test['CTCS'], label='CTCS')
#plt.show()

def selectCourants(nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1], cMin=0.1, cMax=5, nC=10):
    #courants = np.linspace(cMin, cMax, nC)
    Output = {}
    dx = (xLim[1] - xLim[0])/nx
    x = np.arange(xLim[0], xLim[1], dx)
    #dt = endTime/nt
    #c = nx*u*endTime/(nt*(xLim[1] - xLim[0])) # Courant number, dt*u/dx
    #print('Courant number =', c)
        
    #nt=np.zeros(len(courants))

    nts = np.arange(5, 50)

    courants = nx*u*endTime/(nts*(xLim[1]-xLim[0])) 
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

Courant_test = selectCourants(nx = 40, nt = 10, u = 1, endTime = 0.3, xLim = [0,1])

def calc_errors_times(init, schemes,
                       nx = 40, u = 1, endTime = 0.3, xLim = [0,1], a = 0.1, b = 0.5):
#want to run through schemes and 
    l2errors = {}
    mass_diffs = {}
    amp_ratios ={}
    for scheme in schemes:
        l2errors[scheme.__name__] = []
        mass_diffs[scheme.__name__] = []
        amp_ratios[scheme.__name__] = []
    for i in Courant_test['nts']:
        #xLim = [0,1]
        #nx = 40
        A=evolveScheme(init, schemes,
                       nx = nx, nt = i, u = u, endTime = endTime, xLim = xLim, a = a, b = b)
        for scheme in schemes:
            l2errors[scheme.__name__].append(l2error(A['Analytic'],A[scheme.__name__]))
            mass_diffs[scheme.__name__].append(Mass_diff(A['Analytic'],A[scheme.__name__], dx = (xLim[1] - xLim[0])/nx))
            amp_ratios[scheme.__name__].append(amp_ratio(A['Analytic'],A[scheme.__name__]))
#print(np.shape(l2errors))
    return l2errors, mass_diffs, amp_ratios

l2errors, mass_diffs, amp_ratios = calc_errors_times(cosBell, schemes, nx = 40, u =1, endTime =0.3,
                                                     xlim=[1,0], a=0.1, b=0.5)

plt.figure()

for scheme in schemes:

    plt.plot(
        Courant_test['courants'],
        l2errors[scheme.__name__],
        label=scheme.__name__
    )


plt.xlabel('Courant Number')
plt.ylabel('L2 Error')
plt.title('L2 Error vs Courant Number')
plt.legend()
plt.grid()
plt.show()

plt.figure()

for scheme in schemes:

    plt.plot(
        Courant_test['courants'],
        mass_diffs[scheme.__name__],
        label=scheme.__name__
    )


plt.xlabel('Courant Number')
plt.ylabel('Mass Difference')
plt.title('Mass Difference vs Courant Number')
plt.legend()
plt.grid()
plt.show()


plt.figure()

for scheme in schemes:

    plt.plot(
        Courant_test['courants'],
        amp_ratios[scheme.__name__],
        label=scheme.__name__
    )


plt.xlabel('Courant Number')
plt.ylabel('Amp ratio')
plt.title('Amp Ratio vs Courant Number')
plt.legend()
plt.grid()
plt.show()


courants = Courant_test['courants']
#sort so courant numbers are acending instead
order = np.argsort(Courant_test['courants'])

courants = np.array(Courant_test['courants'])[order]
data = np.array([
    np.array(l2errors[scheme.__name__])[order]
    #np.array(amp_ratios[scheme.__name__])[order]
    for scheme in schemes
])
plt.figure(figsize=(8,4))

plt.imshow(
    data,
    interpolation='none',
    aspect='auto'
)

plt.colorbar(label='L2 Error')

step = 5

plt.xticks(
    np.arange(0, len(courants), step),
    np.round(courants[::step], 2),
    rotation=45
)

plt.yticks(
    np.arange(len(schemes)),
    [scheme.__name__ for scheme in schemes]
)

plt.xlabel('Courant Number')
plt.ylabel('Scheme')

plt.show()

widthCells = np.array([
    20,
    15,
    10,
    8,
    6,
    4,
    3,
    2
])

dx = (xLim[1] - xLim[0])/nx
widths = widthCells * dx

width_l2errors = {}
width_amp_ratios = {}

for scheme in schemes:
    width_l2errors[scheme.__name__] = []
    width_amp_ratios[scheme.__name__] = []

for width in widths:

    centre = 0.3

    a = centre - width/2
    b = centre + width/2

    A = evolveScheme(
        cosBell,
        schemes,
        nx=40,
        nt=10,
        u=1,
        endTime=0.3,
        xLim=[0,1],
        a=a,
        b=b
    )

    

    for scheme in schemes:

        width_l2errors[scheme.__name__].append(
            l2error(
                A['Analytic'],
                A[scheme.__name__]
            )
        )

        width_amp_ratios[scheme.__name__].append(
            amp_ratio(
                A['Analytic'],
                A[scheme.__name__]
            )
        )

plt.figure()

for scheme in schemes:

    plt.plot(
        widths,
        width_l2errors[scheme.__name__],
        label=scheme.__name__
    )

plt.xlabel("Bell Width")
plt.ylabel("L2 Error")
plt.yscale('log')

plt.legend()
plt.grid()

plt.show()

plt.figure()

for scheme in schemes:

    plt.plot(
        widths,
        width_amp_ratios[scheme.__name__],
        label=scheme.__name__
    )

plt.xlabel("Bell Width")
plt.ylabel("Amplitude Ratio")

plt.legend()
plt.grid()

plt.show()