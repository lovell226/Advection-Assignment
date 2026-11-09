###############################################################################
# Date Revised: 08 Nov 2026
# Author: Becks Lovell <r.lovell@pgr.reading.ac.uk>
# Description: All tests of Advection schemes
###############################################################################
import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, squareWave, l2error, amp_ratio
from schemesBasic import FTBS, CTCS, FTCS, LinearSemiLagrangian, CubicSemiLagrangian
from advection import compareSchemeSnapshots
from function_setup import evolveScheme, selectCourants, MetricStore, parameterSweep, plotMetric

# First make sure that the plots folder exists
if not os.path.exists("plots"):
    os.makedirs("plots")

#plots for section 3

schemesbasic = [FTBS, CTCS, FTCS]

compareSchemeSnapshots(cosBell, schemesbasic, "basicSchemesCosBell.pdf",
                       nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1], a = 0.1, b = 0.5)

compareSchemeSnapshots(squareWave, schemesbasic, "basicSchemesSquare.pdf",
                       nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1], a = 0.1, b = 0.5)

#plots for section 5.1

schemes = [FTBS, CTCS, LinearSemiLagrangian, CubicSemiLagrangian]

compareSchemeSnapshots(cosBell, schemes, "NewSchemesCos.pdf",
                       nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1], a = 0.1, b = 0.5)

compareSchemeSnapshots(squareWave, schemes, "NewSchemesSquare.pdf",
                       nx = 40, nt = 10, u = 1, endTime = 0.1, xLim = [0,1], a = 0.1, b = 0.5)


# plots for section 5.2 - large courant number testing

Courant_test = selectCourants(nx=40,endTime=2.0)

courantMetrics = parameterSweep(Courant_test["nts"],schemes,endTime=2.0, init=cosBell)

plotMetric(
    Courant_test["courants"],
    "L2_Courant.pdf",
    courantMetrics,
    "L2",
    "Courant Number",
    r'$\ell_2$ error',
    "L2 Error vs Courant Number - CosBell",
    logy=True,  
    yLim=[0,10]
)


# for sqaure wave:
courantMetrics2 = parameterSweep(Courant_test["nts"],schemes,endTime=2.0, init=squareWave)

plotMetric(
    Courant_test["courants"],
    "Square_L2_Courant.pdf",
    courantMetrics2,
    "L2",
    "Courant Number",
    r'$\ell_2$ error',
    "L2 Error vs Courant Number - Square Wave",
    logy=True,
    yLim=[0,10]
)


# example plots for section 5.2 - courant number of 2.0
compareSchemeSnapshots(cosBell, schemes, "cos_courant_large.pdf",
                       nx = 40, nt = 50, u = 1, endTime = 2.0, xLim = [0,1], a = 0.1, b = 0.5, yLim=[-0.1,1.1])

compareSchemeSnapshots(squareWave, schemes, "square_courant_large.pdf",
                       nx = 40, nt = 50, u = 1, endTime = 2.0, xLim = [0,1], a = 0.1, b = 0.5, yLim=[-0.1,1.1])

# plots for section 5.3 - spatial scale testing

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
    "Number of Spatial Cells",
    r'$\ell_2$ error',
    "L2 Error vs Width - CosBell",
    logy=True
)

plotMetric(
    widthCells,
    "Amp_Width.pdf",
    widthMetrics,
    "Amp",
    "Number of Spatial Cells",
    "Amplitude Ratio",
    "Amplitude Ratio vs Width - CosBell"
)

widthMetrics2 = parameterSweep(
    widths,
    schemes,
    varyWidth=True,
    nx = 40,
    nt = 15,
    endTime = 0.3,
    init=squareWave
)

plotMetric(
    widthCells,
    "Square_L2_Width.pdf",
    widthMetrics2,
    "L2",
    "Number of Spatial Cells",
    r'$\ell_2$ error',
    "L2 Error vs Width - Square Wave",
    logy=True
)

plotMetric(
    widthCells,
    "Square_Amp_Width.pdf",
    widthMetrics2,
    "Amp",
    "Number of Spatial Cells",
    "Amplitude Ratio",
    "Amplitude Ratio vs Width - Square Wave"
)

# snapshots for section 5.3
compareSchemeSnapshots(cosBell, schemes, "cos_width_large.pdf",
                       nx = 40, nt = 15, u = 1, endTime = 0.3, xLim = [0,1], a = 0.05, b = 0.55, yLim=[-0.1,1.1])
compareSchemeSnapshots(cosBell, schemes, "cos_width_med.pdf",
                       nx = 40, nt = 15, u = 1, endTime = 0.3, xLim = [0,1], a = 0.2, b = 0.4, yLim=[-0.1,1.1])
compareSchemeSnapshots(cosBell, schemes, "cos_width_small.pdf",
                       nx = 40, nt = 15, u = 1, endTime = 0.3, xLim = [0,1], a = 0.275, b = 0.325, yLim=[-0.1,1.1])

compareSchemeSnapshots(squareWave, schemes, "square_width_large.pdf",
                       nx = 40, nt = 15, u = 1, endTime = 0.3, xLim = [0,1], a = 0.05, b = 0.55, yLim=[-0.1,1.1])
compareSchemeSnapshots(squareWave, schemes, "square_width_med.pdf",
                       nx = 40, nt = 15, u = 1, endTime = 0.3, xLim = [0,1], a = 0.2, b = 0.4, yLim=[-0.1,1.1])
compareSchemeSnapshots(squareWave, schemes, "square_width_small.pdf",
                       nx = 40, nt = 15, u = 1, endTime = 0.3, xLim = [0,1], a = 0.275, b = 0.325, yLim=[-0.1,1.1])
