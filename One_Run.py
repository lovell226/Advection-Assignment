import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, l2error, l_inf_error, Mass_diff, amp_ratio, peak_shift
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2
#from advection import compareSchemeSnapshots
#from growingErrors import growingErrors


def diagnostics(x, dx, phi0, phi, phiA):
    '''
    Collect all diagnostics into a single dictionary.

    Parameters
    ----------
    x : ndarray
        Grid locations

    dx : float
        Grid spacing

    phi0 : ndarray
        Initial condition

    phi : ndarray
        Numerical solution

    phiA : ndarray
        Exact solution

    Returns
    -------
    dict
    '''

    metrics = {}

    metrics['L2'] = l2error(phi, phiA)

    metrics['Linf'] = l_inf_error(phi, phiA)

    metrics['massError'] = Mass_diff(phi0, phi, dx)

    metrics['ampRatio'] = amp_ratio(phi0, phi)

    metrics['peakShift'] = peak_shift(x, phi, phiA)

    metrics['maxVal'] = np.max(phi)

    metrics['minVal'] = np.min(phi)

    metrics['overshoot'] = np.max(phi) - np.max(phi0)

    metrics['undershoot'] = np.min(phi) - np.min(phi0)

    metrics['widthInitial'] = field_width(x, phi0)

    metrics['widthFinal'] = field_width(x, phi)

    metrics['widthGrowth'] = (
        metrics['widthFinal']
        - metrics['widthInitial']
    )

    return metrics

def field_width(x, phi):

    total = np.sum(phi)

    if np.abs(total) < 1e-12:
        return np.nan

    mean = np.sum(x*phi)/total

    sigma = np.sqrt(
        np.sum((x-mean)**2 * phi)/total
    )

    return sigma

def runExperiment(
        scheme,
        phiInit,
        c,
        nt,
        x,
        u=1.0):
    '''
    Run a complete advection experiment.

    Parameters
    ----------
    scheme : function
        FTBS, FTCS, CTCS, SemiLagrangian etc.

    phiInit : ndarray
        Initial condition

    c : float
        Courant number

    nt : int
        Number of timesteps

    x : ndarray
        Grid points

    u : float
        Advection velocity

    Returns
    -------
    dict
    '''

    dx = x[1] - x[0]

    dt = c*dx/u

    phiLevels = [phiInit.copy()]

    for n in range(nt):
        scheme(phiLevels, c)

    phiFinal = phiLevels[-1]

    # exact solution

    t = 1.0

    L = x[-1] - x[0] + dx

    xDeparture = (x - u*t) % L

    phiExact = np.interp(
        xDeparture,
        x,
        phiInit,
        period=L
    )

    metrics = diagnostics(
        x,
        dx,
        phiInit,
        phiFinal,
        phiExact
    )

    return {
        "phi0": phiInit,
        "phiFinal": phiFinal,
        "phiExact": phiExact,
        "metrics": metrics,
        "courant": c,
        "nt": nt
    }
def courantSweep(
        scheme,
        phiInit,
        courants,
        nt,
        x):

    results = []

    for c in courants:

        results.append(
            runExperiment(
                scheme,
                phiInit,
                c,
                nt,
                x
            )
        )

    return results

courants = np.linspace(0.1,5,40)
xLim=[0,1]
dx = (xLim[1] - xLim[0])/40
x = np.arange(xLim[0], xLim[1], dx)
results = courantSweep(
    SemiLagrangian2,
    cosBell(x,0.2,0.4),
    courants,
    100,
    x
)

L2 = [r['metrics']['L2'] for r in results]

plt.plot(courants, L2)
plt.xlabel('Courant Number')
plt.ylabel('L2 Error')
plt.show()

def gaussian(x, sigma, xc=0.5):

    return np.exp(
        -(x-xc)**2/(2*sigma**2)
    )

def widthSweep(
        scheme,
        sigmas,
        c,
        nt,
        x):

    results = []

    for sigma in sigmas:

        phi0 = gaussian(x,sigma)

        results.append(
            runExperiment(
                scheme,
                phi0,
                c,
                nt,
                x
            )
        )

    return results

sigmas = [0.2,0.1,0.05,0.025,0.0125]

results = widthSweep(
    SemiLagrangian2,
    sigmas,
    0.8,
    100,
    x
)

errors = [r['metrics']['L2'] for r in results]

plt.semilogx(sigmas, errors)
plt.xlabel('width')
plt.ylabel('L2 Error')
plt.show()